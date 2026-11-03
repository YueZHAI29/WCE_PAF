
# estimate PAF with the selected model with one knot and right-constrained

library(tidyverse)
library(conflicted)
library(table1)
library(splines)
library(survival)
library(WCE)
library(doParallel)
library(foreach)
library(ggplot2)

conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("filter", "dplyr")

load("/e3n_database.RData") 

thm_k_td <- readRDS("~/longdata.rds")
thm_k_td$duree_regle <- splines::ns(thm_k_td$agedeb - thm_k_td$age_regle, df = 3)

complete_data <- thm_k_td
ev <- sort(unique(complete_data[complete_data$outcome == 1, "tfin_thm", drop = TRUE]))
Id <- unique(complete_data$ident)

# extract small_data containing all individuals at all event time
small_data <- complete_data %>%
  filter(tfin_thm %in% ev)  %>%
  arrange(ident, tfin_thm)

###############################################################################################
##   Step 1. Compute (and save) the WCE matrix (time-varying WCE covariates) (~20 minutes)   ##
###############################################################################################
generate_WCE_matrix <- function(complete_data,
                                exposure_var = "treatment_lag1year",
                                cutoff = 120,
                                nknots = 1,
                                knots = 60,
                                spline_order = 4,
                                num_cores = 11) {
  # Build augmented knot sequence and basis
  knots_augm <- c((-spline_order + 1):0, knots, cutoff:(cutoff + spline_order - 1))
  Bbasis <- splineDesign(knots = knots_augm, x = 1:cutoff, ord = spline_order)
  # Define event times and subject IDs
  ev <- sort(unique(complete_data[complete_data$outcome == 1, "tfin_thm", drop = TRUE]))
  Id <- unique(complete_data$ident)
  
  # Parallel backend setup
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Parallel WCE matrix computation
  WCEvar_matrix <- foreach(i = 1:length(Id), .combine = rbind, .packages = "WCE") %dopar% {
    WCE:::wcecalc(
      ev,
      complete_data[[exposure_var]][complete_data$ident == Id[i]],
      complete_data$tfin_thm[complete_data$ident == Id[i]],
      Bbasis,
      cutoff
    )
  }
  stopCluster(cl)
  # Clean output
  WCEvar_matrix <- as.data.frame(WCEvar_matrix)
  WCEvar_matrix <- WCEvar_matrix[!is.na(WCEvar_matrix[, 1]), ]
  colnames(WCEvar_matrix) <- paste0("WCEvar", 1:ncol(WCEvar_matrix))
  
  return(WCEvar_matrix)
}

###  1 knot, cutoff = 120 months
WCEvar_matrix_lag_1knot <- generate_WCE_matrix(
  complete_data = complete_data,
  exposure_var = "treatment_lag1year",
  nknots = 1,
  cutoff = 120,
  knots = 60,
  spline_order = 4,
  num_cores = 20
)

## function to fit wce models
fit_wce_model <- function(complete_data = complete_data, #modified by Ana 
                          cutoff = 120,
                          n.knots = 2,
                          knots = c(12, 60),
                          spline_order = 4,
                          time.start = "mdeb_thm",
                          time.stop = "mfin_thm",
                          Constrained = "right",
                          covariates = c("imcq", "kcsein_famq", "mbsq", "duree_regle", "menotypfin", "etude"),
                          weighted = TRUE,
                          id_var = "ident",
                          WCE_var_matrix) { # Added by Ana
  
  # extract small_data for event time only
  #data <- complete_data #removed by Ana as never used - we work with complete_data in the function
  ev <- sort(unique(complete_data[complete_data$outcome == 1, "tfin_thm", drop = TRUE]))
  Id <- unique(complete_data$ident)
  smalldat <- complete_data[complete_data$tfin_thm %in% ev, ]
  
  stopifnot(nrow(WCE_var_matrix) == nrow(smalldat)) # Check added by Ana - the error before was caused by merging the complete data with the WCE matrix estimated at event times only
  smalldat <- cbind(smalldat, WCE_var_matrix) # Added by Ana within the function as the design matrix only exist for event times
  
  # create spline basis
  knots_augm <- c((-spline_order + 1):0, knots, cutoff:(cutoff + spline_order - 1))
  basis_matrix <- splines::splineDesign(knots = knots_augm, x = 1:cutoff, ord = spline_order)
  
  # string of TD WCE variables according to constrains
  # cubic: drop the first or last two, quadratic: drop the first or last one
  if (Constrained == "left") {
    string_of_var <- paste0("WCEvar", (spline_order - 1):ncol(basis_matrix))
  } else if (Constrained == "right") {
    string_of_var <- paste0("WCEvar", 1:(n.knots + 2))
  } else {
    string_of_var <- paste0("WCEvar", 1:ncol(basis_matrix))
  }
  # build formula with covariates
  formula <- as.formula(paste0(
    "Surv(", time.start, ", ", time.stop, ", outcome) ~ ",
    paste(c(string_of_var, covariates), collapse = " + ")
  ))
  # Fit IP weighted or not weighted models
  if (weighted == TRUE) {
    weights <- smalldat[["ipw"]]
    cox_wce <- coxph(
      formula = formula, data = smalldat, cluster = smalldat[[id_var]],
      weights = weights, control = coxph.control(timefix = FALSE)
    )
  } else if (weighted == FALSE) {
    cox_wce <- coxph(
      formula = formula, data = smalldat, cluster = smalldat[[id_var]],
      control = coxph.control(timefix = FALSE)
    )
  }
  
  coefs <- cox_wce$coefficients
  basis_col <- gsub("WCEvar", "",string_of_var) # Added by Ana - extract the column index of the targeted basis
  basis_col <- as.numeric(basis_col) # Added by Ana - extract the column index of the targeted basis
  coef_expo <- basis_matrix[, basis_col] %*% coefs[c(string_of_var)] # Ana : correcting the left constrained estimation
  
  return(list(
    cox_model = cox_wce,
    coef_expo = as.numeric(coef_expo),
    spline_basis = basis_matrix,
    string_of_var = string_of_var
  ))
}

small_data_lag_1knot <- cbind(small_data, WCEvar_matrix_lag_1knot)

bic_model <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "right",
  n.knots = 1, knots = 60, 
  WCE_var_matrix = WCEvar_matrix_lag_1knot
)

# get survival estimates under observed exposure
null_ipw <- coxph(Surv(mdeb_thm, mfin_thm, outcome) ~ 1, weights = smalldat$ipw, data = smalldat)
surv_curve <- survfit(null_ipw)

# fup_length <- 120

calculate_paf <- function(fup_length, 
                          complete_data, 
                          string_of_var,
                          bic_model, 
                          surv_curve) {
  
  predict_mat <- complete_data %>%
    group_by(ident) %>%
    dplyr::filter(tdeb_thm == 0) # take the baseline row for everyone --> Q Ana : assuming time fix covariates? 
  
  # Added by Ana : Set exposure to 0 to get the counter factual --> if the entire population was not exposed
  for (i in seq_along(string_of_var)) {
    predict_mat[[string_of_var[i]]] <- 0
  }
  predict_mat$mfin_thm <- predict_mat$mdeb_thm + fup_length    # add fup time to time zero
  
  predict_mat$predicted_surv_0 <- predict(bic_model, newdata = predict_mat, type = "survival") # predict survival in the absence of exposure
  Proba_predict_never_exposed <- 1 - mean(predict_mat$predicted_surv_0) #Average survival if the entire population was not exposed
 
  ### Get observed survival
  surv_begin <- summary(surv_curve, times = predict_mat$mdeb_thm)$surv
  surv_end <- summary(surv_curve, times = predict_mat$mfin_thm)$surv
  predict_mat$proba_observed <- surv_end / surv_begin
  Proba_predict_observed <- 1 - mean(predict_mat$proba_observed)
  
  # PAF:
  paf <- (Proba_predict_observed - Proba_predict_never_exposed) / Proba_predict_observed
  
  return(c(
    fup_length = fup_length,
    Proba_predict_observed = Proba_predict_observed,
    Proba_predict_never_exposed = Proba_predict_never_exposed,
    PAF = paf
  ))
}

# calculate PAF at different time points
fup_seq <- seq(1, 180, by = 5)

paf_results <- sapply(fup_seq, function(fup) {
  calculate_paf(fup, 
                complete_data, 
                bic_model$string_of_var,
                bic_model$cox_model, 
                surv_curve)
})

paf_df <- as.data.frame(t(paf_results))

# Set the output file
pdf("/PAF_plots_bic.pdf", width = 10, height = 5)   # adjust size 

# Make two plots side-by-side
par(mfrow = c(1, 2))  # 1 row, 2 columns

# --- Plot 1: PAF over follow-up ---
plot(
  paf_df$fup_length, paf_df$PAF,
  type = "l", col = "blue3", lwd = 2,
  xlab = "Follow-up length (months)",
  ylab = "Estimated PAF",
  main = "Estimated PAF over follow-up time"
)

# --- Plot 2: Observed vs Never Exposed cumulative risk ---
matplot(
  paf_df$fup_length,
  paf_df[, c("Proba_predict_observed", "Proba_predict_never_exposed")],
  type = "l", lty = 1, col = c("red4", "green4"), lwd = 2,
  xlab = "Follow-up length (months)",
  ylab = "Cumulative risk",
  main = "Observed vs Counterfactual (Never Exposed)"
)
legend("topleft",
       legend = c("Observed", "Never Exposed"),
       col = c("red4", "green4"), lty = 1, lwd = 2)

# Close the PDF device — this actually writes the file
dev.off()



