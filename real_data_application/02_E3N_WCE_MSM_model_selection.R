# Read-data application for WCE-MSM
# E3N cohort, estimating PAF of breast cancer due to HRT
# Part one Model Selection

rm(list = ls()) 

library(tidyverse)
library(conflicted)
library(table1)
library(splines)
library(survival)
library(WCE)
library(doParallel)
library(foreach)
library(ggplot2)
load("~/e3n_database.RData") 

#################################################
#####    1. Construct propensity score    #######  
#################################################
### current age
thm_k_td$age <- thm_k_td$mdeb_thm/12
thm_k_td$before_2002 <- ifelse(thm_k_td$id_q <=7, 1, 0) # binary indicator for calendar time

### lag HRT for different length
thm_k_td <- thm_k_td %>%
  arrange(ident, mdeb_thm) %>%
  group_by(ident) %>%
  mutate(treatment_lag1year = dplyr::lag(thm_use, n = 12))
thm_k_td$treatment_lag1year[is.na(thm_k_td$treatment_lag1year)] <-0

### missing covariables - impute with median/mean
thm_k_td$imcq[is.na(thm_k_td$imcq)] <- median(thm_k_td$imcq, na.rm = TRUE)
thm_k_td$etude[is.na(thm_k_td$etude)] <- median(as.numeric(thm_k_td$etude), na.rm = TRUE)
thm_k_td$age_regle[is.na(thm_k_td$age_regle)] <- median(thm_k_td$age_regle, na.rm = TRUE)

### Model for treatment and mammography
# Check descriptives and to select covariates for the inverse probability weight model
# table1( ~ as.factor(mamq) + as.factor(treatment_lag3months) + as.factor(treatment_lag2years)
#         + pillq  +  age + imcq + kcsein_famq + mbsq +
#           agedeb  + as.factor(menotypfin) + as.factor(generation)| as.factor(thm_use), 
#         data = thm_k_td )

###  Estimate stabilized weights
### denominator(predicted treatment probability given confounders and previous treatment)
model_treatment <- glm(
  thm_use ~ mamq + treatment_lag3months +
    splines::ns(age, df = 3) + imcq + as.factor(kcsein_famq) + as.factor(mbsq) +
    splines::ns(agedeb, df = 3) + as.factor(menotypfin) + as.factor(generation)
    + as.factor(before_2002),
  family = binomial, data = thm_k_td
)
# summary(model_treatment)$coefficients
proba_hrt_predicted <- predict(model_treatment, newdata = thm_k_td, type = "response")
# summary(proba_hrt_predicted)
thm_k_td$weight <- ifelse(thm_k_td$thm_use == 1, proba_hrt_predicted, 1 - proba_hrt_predicted)
thm_k_td$weight.cum <- ave(thm_k_td$weight, thm_k_td$ident, FUN = cumprod)

### numerator
model_treatment_num <- glm(
  thm_use ~ treatment_lag3months +
    splines::ns(age, df = 3) + imcq + as.factor(kcsein_famq) + as.factor(mbsq) +
    splines::ns(agedeb, df = 3) + as.factor(menotypfin) + as.factor(generation)
    + as.factor(before_2002),
  family = binomial, data = thm_k_td
)
proba_hrt_predicted_num <- predict(model_treatment_num, newdata = thm_k_td, type = "response")
thm_k_td$weight_num <- ifelse(thm_k_td$thm_use == 1, proba_hrt_predicted_num, 1 - proba_hrt_predicted_num)
thm_k_td$weight.cum.num <- ave(thm_k_td$weight_num, thm_k_td$ident, FUN = cumprod)

### stablized weights, truncated at 1st and 99th percentile:
thm_k_td$ipw <- thm_k_td$weight.cum.num / thm_k_td$weight.cum

truncated_ipw <- ifelse(thm_k_td$ipw < quantile(thm_k_td$ipw, 0.01),
  quantile(thm_k_td$ipw, 0.01), thm_k_td$ipw
)
truncated_ipw <- ifelse(thm_k_td$ipw > quantile(thm_k_td$ipw, 0.99),
  quantile(thm_k_td$ipw, 0.99), thm_k_td$ipw
)
thm_k_td$ipw <- truncated_ipw
rm(model_treatment)
rm(model_treatment_num)
rm(proba_hrt_predicted)
rm(proba_hrt_predicted_num)

thm_k_td$kcsein_famq <- as.factor(thm_k_td$kcsein_famq)
thm_k_td$mbsq <- as.factor(thm_k_td$mbsq)
thm_k_td$menotypfin <- as.factor(thm_k_td$menotypfin)
thm_k_td$etude <- as.factor(thm_k_td$etude)
saveRDS(thm_k_td, file = "~/R/abc_meht/WP1/data/20250718_longdata.rds") # database with weights and transformed variates

############################################
#####     2. IP weighted Cox model     #####
############################################
thm_k_td$duree_regle <- splines::ns(thm_k_td$agedeb - thm_k_td$age_regle, df = 3)

### IP weighted cox model
## 1) with current treatment
cox_ipw <- coxph(
  Surv(mdeb_thm, mfin_thm, outcome) ~ thm_use + imcq + kcsein_famq + mbsq +
    splines::ns(agedeb - age_regle, df = 3) + menotypfin + etude,
  weights = thm_k_td$ipw,
  data = thm_k_td
)
print(summary(cox_ipw))

## 2) with treatment lagged one year
cox_ipw_lag1year <- coxph(
  Surv(mdeb_thm, mfin_thm, outcome) ~ treatment_lag1year + imcq + kcsein_famq + mbsq +
    splines::ns(agedeb - age_regle, df = 3) + menotypfin + etude,
  weights = thm_k_td$ipw,
  data = thm_k_td
)
print(summary(cox_ipw_lag1year))

# time since menopause as time scale
# cox_fup <- coxph(Surv(tdeb_thm, tfin_thm, outcome) ~ thm_use + imcq + kcsein_famq + mbsq +
#                    splines::ns(age, df = 4) + menotypfin + etude,
#                  weights = thm_k_td$ipw,  data = thm_k_td)
# summary(cox_fup)

#####           WCE package         #####
###  check dataframe structure for wce package
thm_k_td <- as.data.frame(thm_k_td)
checkWCE(data = thm_k_td, id = "ident", event = "outcome", start = "tdeb_thm", stop = "tfin_thm", expos = "thm_use")

################################################
###    With modified code from WCE package   ###
# necessary variables for wce models
complete_data <- thm_k_td
ev <- sort(unique(complete_data[complete_data$outcome == 1, "tfin_thm", drop = TRUE]))
Id <- unique(complete_data$ident)
###############################################################################################
##   Step 1. Compute (and save) the WCE matrix (time-varying WCE covariates)    ##
###############################################################################################
generate_WCE_matrix <- function(complete_data,
                                exposure_var = "treatment_lag1year",
                                cutoff = 120,
                                nknots = 2,
                                knots = c(12, 60),
                                spline_order = 4,
                                num_cores = 9,
                                output_path = "~/R/abc_meht/WP1/data/WCE/") {
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

  # Save if output_path is provided
  saveRDS(WCEvar_matrix,
    file = paste0(
      output_path, Sys.Date(), "_", nknots, "_knots_",
      paste(knots, collapse = "_"), "_", exposure_var, ".rds"
    )
  )
  return(WCEvar_matrix)
}

### Model 1 (database 1): 2 knots, cutoff = 120 months
WCEvar_matrix_lag_2knots <- generate_WCE_matrix(
  complete_data = complete_data,
  exposure_var = "treatment_lag1year",
  nknots = 2,
  cutoff = 120,
  knots = c(12, 60),
  spline_order = 4,
  num_cores = 13,
  output_path = "~/R/abc_meht/WP1/data/WCE/"
)

### Model 2 (database 2): 1 knot, cutoff = 120 months
WCEvar_matrix_lag_1knot <- generate_WCE_matrix(
  complete_data = complete_data,
  exposure_var = "treatment_lag1year",
  nknots = 1,
  cutoff = 120,
  knots = 60,
  spline_order = 4,
  num_cores = 9,
  output_path = "~/R/abc_meht/WP1/data/WCE/"
)

### Model 3 (database 3): 2 knots at 40 and 80, cutoff = 120 months
WCEvar_matrix_lag_2knots_even <- generate_WCE_matrix(
  complete_data = complete_data,
  exposure_var = "treatment_lag1year",
  nknots = 2,
  cutoff = 120,
  knots = c(40, 80),
  spline_order = 4,
  num_cores = 13,
  output_path = "~/R/abc_meht/WP1/data/WCE/"
)

### Model 4 (database 4): 3 knots at 30, 60 and 90, cutoff = 120 months
WCEvar_matrix_lag_3knots_even <- generate_WCE_matrix(
  complete_data = complete_data,
  exposure_var = "treatment_lag1year",
  nknots = 3,
  cutoff = 120,
  knots = c(30, 60, 90),
  spline_order = 4,
  num_cores = 13,
  output_path = "~/R/abc_meht/WP1/data/WCE/"
)

gc()

################################################
#####       Fit different WCE models       #####
################################################
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
                          WCE_var_matrix) {
  # extract small_data for event time only
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
  
  # Remove large elements of the cow_wce output list
  cox_wce[c("residuals", "y", "linear.predictors", "terms")] <- NULL
  
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

#####################################################
### Database 1: knots = c(12,60), cutoff = 120:   ###

###  COX MSM WCE database 1 2knots, right/left/no constrain:
## higher variance for ipw models, so it's better to constrain it , especially when it's biological plausible
wce_ipw_lag_2knots_right <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE,
  Constrained = "right",
  cutoff = 120,
  n.knots = 2,
  knots = c(12, 60),
  spline_order = 4,
  time.start = "mdeb_thm",
  time.stop = "mfin_thm",
  covariates = c("imcq", "kcsein_famq", "mbsq", "duree_regle", "menotypfin", "etude"),
  id_var = "ident", 
  WCE_var_matrix = WCEvar_matrix_lag_2knots # parameter added by Ana 
)

wce_ipw_lag_2knots_left <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE,
  Constrained = "left",
  cutoff = 120,
  n.knots = 2,
  knots = c(12, 60), 
  WCE_var_matrix = WCEvar_matrix_lag_2knots # parameter added by Ana 
)

wce_ipw_lag_2knots <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE,
  Constrained = "no",
  cutoff = 120,
  n.knots = 2,
  knots = c(12, 60), 
  WCE_var_matrix = WCEvar_matrix_lag_2knots # parameter added by Ana 
)
####    Database 2 , one knot
wce_ipw_lag_1knot_right <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, 
  Constrained = "right",
  n.knots = 1, 
  knots = 60, 
  WCE_var_matrix = WCEvar_matrix_lag_1knot # parameter added by Ana 
)
wce_ipw_lag_1knot_left <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "left",
  n.knots = 1, knots = 60, 
  WCE_var_matrix = WCEvar_matrix_lag_1knot # parameter added by Ana 
)
wce_ipw_lag_1knot <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "no",
  n.knots = 1, knots = 60, 
  WCE_var_matrix = WCEvar_matrix_lag_1knot # parameter added by Ana 
)

####    Database 3 , evenly spaced two knots
wce_ipw_lag_2knots_even_right <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "right",
  n.knots = 2, knots = c(40, 80), 
  WCE_var_matrix = WCEvar_matrix_lag_2knots_even # parameter added by Ana 
  
)
wce_ipw_lag_2knots_even_left <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "left",
  n.knots = 2, knots = c(40, 80), 
  WCE_var_matrix = WCEvar_matrix_lag_2knots_even # parameter added by Ana 
  
)
wce_ipw_lag_2knots_even <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "no",
  n.knots = 2, knots = c(40, 80), 
  WCE_var_matrix = WCEvar_matrix_lag_2knots_even # parameter added by Ana 
  
)

####    Database 4 , evenly spaced THREE knots
wce_ipw_lag_3knots_even_right <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "right",
  n.knots = 3,  knots = c(30, 60, 90),
  WCE_var_matrix = WCEvar_matrix_lag_3knots_even 
)

wce_ipw_lag_3knots_even_left <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "left",
  n.knots = 3,  knots = c(30, 60, 90),
  WCE_var_matrix = WCEvar_matrix_lag_3knots_even 
)

wce_ipw_lag_3knots_even <- fit_wce_model(
  complete_data = complete_data,
  weighted = TRUE, Constrained = "no",
  n.knots = 3,  knots = c(30, 60, 90),
  WCE_var_matrix = WCEvar_matrix_lag_3knots_even 
)

model_names <- c(
  "wce_ipw_lag_2knots_right", "wce_ipw_lag_2knots_left", "wce_ipw_lag_2knots",
  "wce_ipw_lag_1knot_right", "wce_ipw_lag_1knot_left", "wce_ipw_lag_1knot",
  "wce_ipw_lag_2knots_even_right", "wce_ipw_lag_2knots_even_left", "wce_ipw_lag_2knots_even",
  "wce_ipw_lag_3knots_even_right", "wce_ipw_lag_3knots_even_left", "wce_ipw_lag_3knots_even"
)

wce_models <- mget(model_names)

# wce_ipw_2knots_right$cox_model %>% BIC()
# plot(exp(wce_ipw_2knots_right$coef_expo), xlab = "month since exposure(lagged by 1 year)", main = "2 knots, right constrained, weighted")
# exp(sum(wce_ipw_2knots_right$coef_expo))

##############################################
# Loop over model names and extract AIC and BIC
results_mat <- matrix(NA,
  nrow = length(wce_models), ncol = 3,
  dimnames = list(model_names, c("AIC", "BIC", "HR"))
)

for (m in model_names) {
  fit <- wce_models[[m]]
  results_mat[m, "AIC"] <- AIC(fit$cox_model)
  results_mat[m, "BIC"] <- BIC(fit$cox_model)
  results_mat[m, "HR"] <- exp(sum(fit$coef_expo)) # area under the curve / always vs. never
}
print(results_mat)

# Save if output_path is provided
saveRDS(results_mat, 
        file = paste0("~/", Sys.Date(), "_WCE_results_mat.rds"))


saveRDS(list(wce_models, results_mat), 
        file = paste0("~/", Sys.Date(), "_WCE_models.rds"))


#############################################
######       03_add_null_models       #######
#############################################
required_pkgs <- c(
  "tidyverse",
  "conflicted",
  "table1",
  "splines",
  "survival",
  "ggplot2",
  "slider"
)

conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")

# ---- Install any missing packages ----
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

thm_k_td <- readRDS("~/longdata.rds")

############################################
#####     1. IP weighted Cox models    #####
############################################
thm_k_td$duree_regle <- splines::ns(thm_k_td$agedeb - thm_k_td$age_regle, df = 3)

# add cumulative exposure
thm_k_td <- thm_k_td %>%
  arrange(ident, mdeb_thm) %>%   # ensure correct order
  group_by(ident) %>%
  mutate(
    cum_expo_120m = slide_dbl(
      .x = treatment_lag1year,              # the exposure variable
      .f = ~ sum(.x, na.rm = TRUE),
      .before = 120,              # includes current month + 100 previous
      .complete = FALSE           # allow partial windows (missing months treated as zero)
    )
  ) %>%
  ungroup()

## Exposure as treatment lagged one year("current" only)
cox_ipw_lag1year <- coxph(
  Surv(mdeb_thm, mfin_thm, outcome) ~ treatment_lag1year + imcq + kcsein_famq + mbsq +
    duree_regle + menotypfin + etude,
  weights = thm_k_td$ipw,
  data = thm_k_td
)

##  "Null" model
cox_ipw_null <- coxph(
  Surv(mdeb_thm, mfin_thm, outcome) ~ imcq + kcsein_famq + mbsq +
    duree_regle + menotypfin + etude,
  weights = thm_k_td$ipw,
  data = thm_k_td
)

# unweighted cumulative exposure
cox_ipw_cum <- coxph(
  Surv(mdeb_thm, mfin_thm, outcome) ~ cum_expo_120m + imcq + kcsein_famq + mbsq +
    duree_regle + menotypfin + etude,
  weights = thm_k_td$ipw,
  data = thm_k_td
)

################
null_models <- data.frame(
  model = c("lag1year", "cum_expo", "null"),
  BIC = c( BIC(cox_ipw_lag1year), BIC(cox_ipw_cum), BIC(cox_ipw_null)  ),
  AIC = c(AIC(cox_ipw_lag1year), AIC(cox_ipw_cum), AIC(cox_ipw_null)   ),
  coef = c(cox_ipw_lag1year$coef["treatment_lag1year"], cox_ipw_cum$coef["cum_expo_120m"], NA_real_  )
)

saveRDS(null_models, "null_models_BIC.rds")
