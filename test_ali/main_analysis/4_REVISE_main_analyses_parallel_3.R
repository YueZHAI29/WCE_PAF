# ============================================================
# master file for main analyses, cluster - PARALLEL VERSION (standalone)
# All function definitions (from 1_function_to_generate_data.R and the
# FIXED/lite 2_function_to_fit_models.R) are pasted directly below,
#
# Structure mirrors parallel_simulation.r / Copy_parallel_simulation2.r:
#   makeCluster -> clusterEvalQ (load pkgs) -> clusterExport (export fns)
#   -> foreach %dopar% (.errorhandling="pass") -> stopCluster
# ============================================================

# Load necessary libraries
# [REV] Keep numerical libraries at one thread per parallel worker.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(tidyverse)
library(survival)
library(MASS)
library(stats)
library(stringr)
library(matrixStats)
library(splines)
library(zoo)

# Parallel processing libraries
library(foreach)
library(doParallel)
library(parallel)


# ============================================================
# ---- BEGIN: contents of 1_function_to_generate_data.R ----
# ============================================================
# Function to simulate one dataset based on Marginal Structural Cox with WCE modeling
# Data generation mechanism according to Young et al.2010 and Xiao et al. 2014
# effect assumption and weight assumption made on marginal effect
# with binary time-varying confounder affected by previous exposure
# with missing confounder 
# Output: dat.complet, dat.obs, dat.locf
# No other time-fixed confounder was included 

library(tidyverse)

library(MASS)
library(stats)
library(stringr)
library(matrixStats)
library(splines2)
library(zoo)
library(dplyr)

sim_MSM_WCE <- function(weight_hypo = "decreasing", true_cutoff = 50, n=1000, n.visit = 100, 
                        LtoA = 1, AtoL = 0.5, logitpL = 0, LtoL = 1, lambda0 = 0.0005, HRA = 4,  
                        T0toL = 0.5, logitpA = -3, AtoA = 4, n.obs = 5,censor_max = NULL,
                        seed = NA){
  
  if (is.na(seed)) {
    set.seed(sample(10000, 1))
  } else if (is.numeric(seed)) {
    set.seed(seed)
  } else {
    stop("Invalid seed: must be numeric or NA.")
  }
  
  delta0 = -2.5 # intercept of missingness logitP
  deltaT = -0.06 # dependency of missingness on T.obs
  
  #############################################
  ##          Data generation                ##
  ##    Output: dat.obs and dat.complet      ##
  #############################################
  expit=function(x){exp(x)/(1+exp(x))}
  
  # Measurement schedule
  J <-	n.visit  # number of measurement in the "complete" dataset
  Step <- 1
  tj <- rep(seq(0,J-1,Step),n) # vector of visit times for all n individuals (0,...,9)
  Jt <- length(seq(0,J-1,Step)) # number of planned visit time (=10)
  
  # generate survival time T0 (counterfactual for A=0), for each individual
  # from a exponential distribution
  V = runif(n,0,1) 
  T0 <- -log(V)/lambda0
  T0_threshold <- quantile(T0, 0.3)
  
  L = matrix(nrow=n,ncol=n.visit)
  A = matrix(nrow=n,ncol=n.visit)
  
  # A binary
  L[,1]=rbinom(n,1,expit(logitpL 
                         # + T0toL * I(T0 < T0_threshold)
  ) )
  A[,1]=rbinom(n,1,expit(logitpA + LtoA*L[,1]) )# A0|L0: logit(P) = -1 + ..*L0
  for(k in 2:n.visit){
    L[,k] = rbinom(n,1, expit(logitpL + LtoL*L[, k-1] + AtoL * A[, k-1] + 
                                T0toL * I(T0 < T0_threshold)) )
    A[,k]= rbinom(n,1,expit(logitpA + LtoA*L[,k] + AtoA * A[,k-1] ) )# A1|A0,L0,L1,T>=k
  }
  #colMeans(L);  colMeans(A); matrixStats::colSds(L)
  #summary(rowSums(A))
  #hist(rowSums(A))
  
  #####################################################################
  ###       Calculate WCE  using a certain weight assumption        ###
  
  .generate_true_weight <- function(true_cutoff, weight_hypo) {
    x <- 0:true_cutoff
    
    if (weight_hypo == "decreasing") {
      # exponentially decreasing
      y <- exp(-0.1 * x)
      weights <- y / sum(y)
    } else if (weight_hypo == "middle_peak") {
      # as described in Sylvestre et al. 2009 paper (supplementary)
      mean <- 25
      sd <- 0.1
      x_for_density <- seq(mean - 4 * sd, mean + 4 * sd, length.out = true_cutoff + 1)
      y <- dnorm(x_for_density, mean = mean, sd = sd)
      weights <- y / sum(y)
    } else if (weight_hypo == "constant") {
      y <- rep(1 / length(x), length(x))
      weights <- y
    } else if (weight_hypo == "only_current") {
      y <- c(1, rep(0, length(x) - 1))
      weights <- y
    } else {
      stop("Invalid weight_hypo value")
    }
    return(weights)
  }
  weights <- .generate_true_weight(true_cutoff, weight_hypo )
  # plot(weights ~ x)
  
  # generate matrix of exposure and lagged values (long format)
  id =c(1:n)
  Amat = as.data.frame(cbind(id, A))
  Amat = Amat %>% tidyr::pivot_longer(!id, names_to = "visit", values_to = "A")
  Amat$visit <- rep(1:n.visit, n)
  for (i in 1:true_cutoff) {
    Amat = Amat %>%
      group_by(id) %>%
      dplyr::mutate("Alag{i}":= dplyr::lag(A, i))
    Amat[is.na(Amat[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
  }
  ## weighted exposure at each time
  #  = weights for each lagged time exposure *exposure value
  wceA = as.matrix(Amat[, 3:ncol(Amat)]) %*% as.vector(weights)
  ## weighted cumulative exposure = cumsum of the above at each time 
  wceA = cbind(as.data.frame(Amat[, 1:2]), wceA)
  colnames(wceA)[3] <- "wce"
  # pivot back to wide format for survival time generation
  wceA = wceA %>% tidyr::pivot_wider(values_from = "wce", names_from = "visit")
  wceA <- as.matrix(wceA[,-1])
  rm(Amat)
  
  # Generation of survival time
  T.obs <- rep(NA, n)
  # function to calculate integral of  exp(log(HRA)*wceA[,t])dt (Xiao et al.2014)
  v=1    # visit No.1
  exp_sum <- exp(log(HRA) * wceA[, 1])
  exp_sum_lag1 <- 0
  new.t = (T0 - exp_sum_lag1 ) * exp(-log(HRA) * wceA[, 1])
  T.obs = ifelse(T0 <= exp_sum , v-1+new.t, T.obs) 
  
  for(v in 2:n.visit){
    # check if T0 < integral of  exp(log(HRA)*wceA[,t])dt
    #integral:
    exp_sum <- 0
    for (icol in 1:v) {
      exp_sum <- exp_sum + exp(log(HRA) * wceA[, icol])
    }
    exp_sum_lag1 <- exp_sum - exp(log(HRA) * wceA[, v])
    new.t = (T0 - exp_sum_lag1 ) * exp(-log(HRA) * wceA[, v])
    T.obs = ifelse(is.na(T.obs) & T0 <= exp_sum , v-1+new.t, T.obs)
  } 
  D.obs=ifelse(is.na(T.obs),0,1)
  T.obs=ifelse(is.na(T.obs),n.visit,T.obs)
  
  # ----- ADD UNIFORM CENSORING TO TARGET ~10% CENSORING RATE -----
  # Simulate censoring time from Uniform(0, max_time)
  if (!is.null(censor_max)) {
    max_time <- n.visit * censor_max
    censoring_time <- runif(n, min = 0, max = max_time)
    
    # Apply censoring
    T.obs.censored <- pmin(T.obs, censoring_time)
    D.obs.censored <- as.integer(T.obs <= censoring_time & T.obs < 100)  # 1=event observed, 0=censored
    
    # Replace original T.obs and D.obs
    T.obs <- T.obs.censored
    D.obs <- D.obs.censored
  }
  # Check actual censoring rate (optional)
  # cat("Censoring rate: ", mean(D.obs == 0), "\n")
  
  # summary(T.obs)
  # summary(D.obs); table(D.obs)
  
  # Create data frame and reshape into 'long' format (1 row for each visit)
  colnames(A)=paste0("A.",0:(n.visit-1))
  colnames(L)=paste0("L.",0:(n.visit-1))
  dat=data.frame(id=1:n,T.obs,D.obs,A,L)
  
  dat.long=reshape(data = dat, varying=c(paste0("A.",0:(n.visit-1)), paste0("L.",0:(n.visit-1)) ), 
                   direction="long",idvar="id")
  dat.long=dat.long[order(dat.long$id,dat.long$time),]
  
  #generate start and stop times for each row
  dat.long$time.stop=dat.long$time+1
  dat.long = dat.long[dat.long$time < dat.long$T.obs, ]
  dat.long$time.stop = ifelse(dat.long$time.stop > dat.long$T.obs,
                              dat.long$T.obs, dat.long$time.stop)
  dat.long$event = ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)
  
  #visit number
  dat.long$visit = ave(rep(1, nrow(dat.long)), dat.long$id, FUN=cumsum)
  
  #generate lagged A values
  for (i in 1:n.visit) {
    dat.long=dat.long %>%
      group_by(id) %>%
      dplyr::mutate("Alag{i}":= dplyr::lag(A, i))
    dat.long[is.na(dat.long[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
  }
  
  #generate lagged L values # in case past values of L are needed (especially for gformula)
  for (i in 1:10) {
    dat.long=dat.long %>%
      group_by(id) %>%
      dplyr::mutate("Llag{i}" := dplyr::lag(L, i)) 
    dat.long[is.na(dat.long[,paste0("Llag",i)] ) , paste0("Llag",i)] <- 0
  }
  # baseline L
  dat.long=dat.long %>%
    group_by(id) %>%
    dplyr::mutate(L.baseline = dplyr::first(L))
  dat.complet <- as.data.frame(dat.long)
  
  #######  NOT NEEDED  #################
  #### Generate missing value
  dat <- dat.long
  # MISSING AT-RISK MEASUREMENTS
  deltaVis_1 <- deltaVis_2 <- deltaVis_3 <- 0.0877 # coefficient of time
  deltaYother1_1 <- deltaYother1_2 <- deltaYother1_3 <- 0.005 # depend on other variables
  deltaYother2_1 <- deltaYother2_2 <- deltaYother2_3 <- 0.005
  delta0_1 <- delta0_2 <- delta0_3 <- delta0   # intercept of missingness logitP
  deltaZ1_1<-deltaZ1_2<-deltaZ1_3<-(-0.1113)
  deltaZ2_1<-deltaZ2_2<-deltaZ2_3<- 0.0255
  deltaD_1<-deltaD_2<-deltaD_3 <- 0.1609
  deltaT_1<-deltaT_2<-deltaT_3<- deltaT
  
  #############################################
  #   GENERATE MISSING AT RISK MEASUREMENTS  ##
  # Marker 1 = L
  dat$Proba_1 <- (dat$time!=0)*expit(delta0_1 + deltaVis_1 * dat$time 
                                     + deltaD_1*dat$event + deltaT_1*dat$T.obs )
  
  dat$mis_1 <- rbinom(nrow(dat),1,dat$Proba_1)
  
  # only n.obs (five) values are observed
  visit.observed <- seq(1,n.visit, n.visit/n.obs) # observed visit time
  
  dat[!(dat$visit %in% visit.observed), "L"] <- NA
  dat[dat$mis_1 == 1,"L"] <- NA
  
  #generate lagged L values #
  for (i in 1:(n.visit-1)) {
    dat.long=dat.long %>%
      group_by(id) %>%
      dplyr::mutate("Llag{i}" := dplyr::lag(L, i)) 
    dat.long[is.na(dat.long[,paste0("Llag",i)] ) , paste0("Llag",i)] <- 0
  }
  dat.obs <- as.data.frame(dat)
  
  ##############################
  ##    locf
  dat <- dat.obs
  for(ind in 1:n){
    dat[dat$id == ind, "L"] <- na.locf(dat[dat$id == ind, "L"],na.rm = FALSE)
  }
  dat.locf <- as.data.frame(dat)
  
  return(list(dat.complet = dat.complet, dat.locf = dat.locf, dat.obs = dat.obs))
}

# dat.list <- sim_MSM_WCE(weight_hypo = "decreasing", true_cutoff = 50, n=1000, n.visit = 100,
#                         LtoA = 1, AtoL = 0.5, logitpL = 0, LtoL = 1, lambda0 = 0.0005, HRA = 4,
#                         T0toL = 0.5, logitpA = -3, AtoA = 4, n.obs = 5,
#                         seed = 1)
#nrow(dat.list$dat.complet)
#dat.complet <- dat.list$dat.complet


# ============================================================
# ---- END: contents of 1_function_to_generate_data.R ----
# ---- BEGIN: contents of 2_function_to_fit_models.R (FIXED / lite) ----
# ============================================================
# Function to fit all models and return model objects on an input dataset
# input a simulated by sim_MSM_WCE_cond() function

# @import function sim_MSM_WCE()

library(tidyverse)
library(survival)

library(MASS)
library(stats)
library(stringr)
library(matrixStats)
library(splines)

compute_ipw <- function(data) {
  # denominator model: Pr(Ak | Lk, A(k-1))
  wt.mod <- glm(A ~ L + Alag1, family = "binomial", data = data)
  pred.wt <- predict(wt.mod, type = "response")
  data$wt <- ifelse(data$A == 1, pred.wt, 1 - pred.wt)
  data$wt.cum <- ave(data$wt, data$id, FUN = cumprod)
  
  # numerator model: Pr(Ak | A(k-1))
  wt.mod.num <- glm(A ~ Alag1, family = "binomial", data = data)
  pred.wt.num <- predict(wt.mod.num, type = "response")
  data$wt.num <- ifelse(data$A == 1, pred.wt.num, 1 - pred.wt.num)
  data$wt.cum.num <- ave(data$wt.num, data$id, FUN = cumprod)
  
  # stabilized weights
  data$ipw.s <- data$wt.cum.num / data$wt.cum
  
  # display summary and 99th percentile
  # summary(data$ipw.s)
  # quantile(data$ipw.s, 0.99)
  
  return(data)
}

# an internal function to fit weighted or unweighted WCE models 
.ipw_wce <- function(nknots = 1, cutoff = 50, weight_var = "ipw.s", covariates = NULL,
                     Constrained = 'no', spline_order = 4, n.visit = 100, n = NULL, 
                     data = NULL, n_full_expo = FALSE){
  dat.complet = data
  time_since_expo <- seq(0,cutoff,1)
  knots_quantile <- seq(1/(nknots + 1), nknots/(nknots + 1), length.out = nknots)
  inner_knots = round(quantile(time_since_expo, knots_quantile), 0)    
  bsknots_augm <- c((-spline_order+1):0, inner_knots, cutoff: (cutoff + spline_order -1))
  basis_matrix <- splines::splineDesign(x = time_since_expo, knots = bsknots_augm,  ord = spline_order)
  
  # calculated time varying wce variable
  for (j in 1:ncol(basis_matrix)){
    for (i in 0:cutoff){
      if (i==0){
        wcej <- dat.complet[, "A"] * basis_matrix[i+1,j]
      }
      else if (i>=1){
        expo_temp <- dat.complet[,paste0("Alag",i)] * basis_matrix[i+1,j]
        wcej <- wcej + expo_temp
      }
      dat.complet[, paste0("WCEvar",j)] <- wcej
    } # end of i loop for each basis function
  } # end of j loop for all basis
  
  # cubic: drop the first or last two, quadratic: drop the first or last one
  if (Constrained == 'left') {
    string_of_var <- paste0("WCEvar", (spline_order-1):ncol(basis_matrix))
  } else if (Constrained == 'right') {
    string_of_var <- paste0("WCEvar", 1:(nknots + 2))
  } else {
    string_of_var <- paste0("WCEvar", 1: ncol(basis_matrix))
  }
  
  formula <- as.formula(paste0("Surv(time, time.stop, event) ~ ", 
                               paste(string_of_var, collapse = " + ")) )
  
  if (is.null(weight_var) == FALSE){
    weights <- dat.complet[[weight_var]]
    cox_wce <- coxph(formula = formula, data = dat.complet, cluster = id, 
                     weights = weights, control = coxph.control(timefix = FALSE))
  }  else if(is.null(weight_var) == TRUE & is.null(covariates) != TRUE){
    formula_adj <- as.formula(paste0("Surv(time, time.stop, event) ~ ", 
                                     paste(string_of_var, collapse = " + "), "+ " , covariates) )
    cox_wce <- coxph(formula = formula_adj, data = dat.complet, cluster = id, 
                     control = coxph.control(timefix = FALSE))
  } else if (is.null(weight_var) == TRUE & is.null(covariates) == TRUE){
    cox_wce <- coxph(formula = formula, data = dat.complet, cluster = id, 
                     control = coxph.control(timefix = FALSE))
  }
  
  ###    PREDICTED P_Y(A=0) (at the end of follow-up)  ########
  predmat <- data.frame(id = unique(dat.complet[, c("id")]))
  predmat$A <- 0
  for (i in 1:3) {
    predmat=predmat %>%
      group_by(id) %>%
      mutate("Alag{i}":= dplyr::lag(A, i))
    predmat[is.na(predmat[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
  }
  predmat$event <- 0
  predmat$time <- 0
  predmat$time.stop <- n.visit
  
  # set time-varying variables all equals 0
  for (i in 1:ncol(basis_matrix)) {
    predmat=predmat %>%
      group_by(id) %>%
      mutate("WCEvar{i}":= 0)
  }
  # predmat$predict0_ipw <- predict(cox_wce, predmat, type = "survival")
  # Pr_A0_ipw <- 1- mean(predmat$predict0_ipw)
  
  survcurv <- survfit(cox_wce, newdata = predmat[1,], se.fit=TRUE, conf.type = "log" )
  survcurv <- summary(survcurv, times = 100)
  Pr_A0_ipw <- data.frame(
    Pr_A0_ipw  = 1-survcurv$surv,
    lower  = 1-survcurv$upper,
    upper  = 1-survcurv$lower
  )
  
  ### Predict Pr_Y(A=1)  -- SKIPPED in the "lite" version.
  # PAF only needs Pr_Y_A0 (see file 5, PAF_est_mat <- (P_observed - P0_predicted)/P_observed).
  # This block built an n*n.visit row prediction frame and re-did a per-id
  # dplyr::lag() loop of length `cutoff` -- profiling showed this dominates
  # .ipw_wce()'s runtime (~84% of total fit_models() time). Set n_full_expo = TRUE
  # in the function call to restore the original behaviour if Pr_Y_A1 is needed
  # again (e.g. for the main manuscript simulation, outside the bootstrap).
  if (isTRUE(n_full_expo)) {
    predpop1 <- data.frame(matrix(nrow = n*n.visit, ncol = 4))
    # id 
    predpop1[,1] <- rep(1:n, each=n.visit)
    # time
    predpop1[,2] <- rep(0:(n.visit-1), n)
    # timestop
    predpop1[,3] <- rep(1:n.visit, n)
    # event
    predpop1[,4] <- rep(0, nrow(predpop1))
    colnames(predpop1) = c("id", "time","time.stop", "event")
    
    predpop1$A <- 1
    for (i in 1:cutoff) {
      predpop1=predpop1 %>%
        group_by(id) %>%
        mutate("Alag{i}":= dplyr::lag(A, i))
      predpop1[is.na(predpop1[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
    }  
    # calculate time varying wce variable
    for (j in 1:ncol(basis_matrix)){
      for (i in 0:cutoff){
        if (i==0){
          wcej <- predpop1[, "A"] * basis_matrix[i+1,j]
        }
        else if (i>=1){
          expo_temp <- predpop1[,paste0("Alag",i)] * basis_matrix[i+1,j]
          wcej <- wcej + expo_temp
        }
        predpop1[, paste0("WCEvar",j)] <- wcej
      } # end of i loop for each basis function
    } # end of j loop for all basis
    cumhaz <- predict(cox_wce, predpop1, type = "expected")
    predpop1$predict1 <- exp(-cumhaz)
    predpop1$surv1 <- ave(predpop1$predict1, predpop1$id,FUN=cumprod)
    predpop1$surv1cum <- ave(predpop1$surv1, predpop1$id, FUN = min)
    Pr_A1_ipw <- 1- mean(predpop1$surv1cum) 
  } else {
    Pr_A1_ipw <- NA_real_
  }
  
  # reconstruct estimated weights
  numbers_wce <- as.numeric(gsub("\\D", "", string_of_var))
  basis_matrix_used <- basis_matrix[, numbers_wce]
  
  coef_wce_var <- summary(cox_wce)$coef[1:length(string_of_var),1]
  weight_predict <- basis_matrix_used %*% coef_wce_var
  
  # Delta method for confidence interval estimation
  varcov <- vcov(cox_wce)
  varcov <- as.matrix(varcov[1:length(string_of_var), 1:length(string_of_var)])
  sd_weight <- vector(length = cutoff+1)
  for (time in 1:(cutoff+1)) {
    sd_weight[time] = sqrt(t(basis_matrix_used[time, ]) %*% varcov %*% basis_matrix_used[time,])
  }
  
  ### Make a plot of the predicted effect(HR)
  weight_lower  = weight_predict - 1.96 * sd_weight
  weight_upper  = weight_predict + 1.96 * sd_weight
  weight_mat_plot <- data.frame(time_since_expo, weight_predict, weight_lower, weight_upper, sd_weight)
  weightplot <- ggplot(weight_mat_plot, aes(x = time_since_expo)) + 
    geom_line(aes(y = exp(weight_predict)), color = "black", linewidth = 1.2)+ # Point estimate (main line)
    geom_ribbon(aes(ymin = exp(weight_lower), ymax =exp(weight_upper)), fill = "lightgray", alpha = 0.5)+ # Shade area between confidence interval (optional)
    geom_line(aes(y = exp(weight_lower)), color = "skyblue4", linetype = "dashed")+ # Upper and lower bound lines (dashed for distinction)
    geom_line(aes(y = exp(weight_upper)), color = "skyblue4", linetype = "dashed")+ # Axis labels and title
    geom_line(aes(y = 1), color = "Orange", linetype = "dashed")+ # Line in between the 95% CI 
    labs(x = "Time since exposure", y = "Hazard ratio", title = "Estimate hazard ratio by time since exposure")
  
  cox_null <- coxph(Surv(time, time.stop, event) ~ 1, data = dat.complet, 
                    cluster = id, control = coxph.control(timefix = FALSE))
  KM_estimator <- summary(survfit(cox_null), times = n.visit)$surv
  
  return(list(
    coefficients = summary(cox_wce)$coefficients, 
    Pr_Y_A1 = Pr_A1_ipw,
    Pr_Y_A0 = Pr_A0_ipw,
    BIC = BIC(cox_wce),
    AIC = AIC(cox_wce),
    coef_weight_predict = weight_predict,
    coef_sd = sd_weight,    
    loglik = cox_wce$loglik,
    vcov = cox_wce$var,
    score_LRT = summary(cox_wce)$logtest, 
    score_logrank_test = summary(cox_wce)$sctest,
    score_waldtest = summary(cox_wce)$waldtest,
    nevent = cox_wce$nevent,
    proba_observed = KM_estimator,
    concordance = cox_wce$concordance
  ))
}

.predict_P0_ipw <- function(cox_ipw = cox_ipw ){
  predmat <- data.frame(id = unique(dat.complet[, c("id")]))
  predmat$A <- 0
  for (i in 1:3) {
    predmat=predmat %>%
      group_by(id) %>%
      mutate("Alag{i}":= dplyr::lag(A, i))
    predmat[is.na(predmat[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
  }
  predmat$event <- 0
  predmat$time <- 0
  # get the maximum follow-up time from cox model
  n.visit <- max(cox_ipw$y[, "stop"])
  predmat$time.stop <- n.visit
  
  # time-varying variables all equals 0
  for (i in 1:length(cox_ipw$coefficients)) {
    predmat=predmat %>%
      group_by(id) %>%
      mutate("WCEvar{i}":= 0)
  }
  predmat$Asum_half <- 0
  predmat$Asum_whole <- 0
  predmat$predict0_ipw <- predict(cox_ipw, predmat, type = "survival")
  Pr_A0_ipw <- 1- mean(predmat$predict0_ipw)
  # predict(cox_ipw, predmat[1,], type = "survival", se.fit=TRUE)
  
  survcurv <- survfit(cox_ipw, newdata = predmat[1,], se.fit=TRUE, conf.type = "log" )
  survcurv <- summary(survcurv, times = 100)
  predicted_P0 <- data.frame(
    Pr_A0_ipw  = 1-survcurv$surv,
    lower  = 1-survcurv$upper,
    upper  = 1-survcurv$lower
  )
  return(predicted_P0)
}


fit_models <- function(dat.complet = dat.complet, n=1000, n.visit = 100){
  # Calculate weight
  compute_ipw <- function(data) {
    # denominator model: Pr(Ak | Lk, A(k-1))
    wt.mod <- glm(A ~ L + Alag1, family = "binomial", data = data)
    pred.wt <- predict(wt.mod, type = "response")
    data$wt <- ifelse(data$A == 1, pred.wt, 1 - pred.wt)
    data$wt.cum <- ave(data$wt, data$id, FUN = cumprod)
    
    # numerator model: Pr(Ak | A(k-1))
    wt.mod.num <- glm(A ~ Alag1, family = "binomial", data = data)
    pred.wt.num <- predict(wt.mod.num, type = "response")
    data$wt.num <- ifelse(data$A == 1, pred.wt.num, 1 - pred.wt.num)
    data$wt.cum.num <- ave(data$wt.num, data$id, FUN = cumprod)
    
    # stabilized weights
    data$ipw.s <- data$wt.cum.num / data$wt.cum
    
    # display summary and 99th percentile
    # summary(data$ipw.s)
    # quantile(data$ipw.s, 0.99)
    
    return(data)
  }
  dat.complet <- compute_ipw(dat.complet)
  
  # calculate cumulated exposure for half and whole of the fup period
  dat.complet <- dat.complet %>%
    mutate(Asum_half = rowSums(across(Alag1:Alag50)))
  dat.complet$Asum_half <- dat.complet$Asum_half + dat.complet$A
  
  dat.complet <- dat.complet %>%
    mutate(Asum_whole = rowSums(across(Alag1:Alag100)))
  dat.complet$Asum_whole <- dat.complet$Asum_whole + dat.complet$A
  
  
  #######################################################
  ###    PREDICTED P_Y(A=0) (at the end of follow-up) for cox ipw model  ########
  .predict_P0_ipw <- function(cox_ipw = cox_ipw ){
    predmat <- data.frame(id = unique(dat.complet[, c("id")]))
    predmat$A <- 0
    for (i in 1:3) {
      predmat=predmat %>%
        group_by(id) %>%
        mutate("Alag{i}":= dplyr::lag(A, i))
      predmat[is.na(predmat[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
    }
    predmat$event <- 0
    predmat$time <- 0
    # get the maximum follow-up time from cox model
    n.visit <- max(cox_ipw$y[, "stop"])
    predmat$time.stop <- n.visit
    
    # time-varying variables all equals 0
    for (i in 1:length(cox_ipw$coefficients)) {
      predmat=predmat %>%
        group_by(id) %>%
        mutate("WCEvar{i}":= 0)
    }
    predmat$Asum_half <- 0
    predmat$Asum_whole <- 0
    predmat$predict0_ipw <- predict(cox_ipw, predmat, type = "survival")
    Pr_A0_ipw <- 1- mean(predmat$predict0_ipw)
    # predict(cox_ipw, predmat[1,], type = "survival", se.fit=TRUE)
    
    survcurv <- survfit(cox_ipw, newdata = predmat[1,], se.fit=TRUE, conf.type = "log" )
    survcurv <- summary(survcurv, times = 100)
    predicted_P0 <- data.frame(
      Pr_A0_ipw  = 1-survcurv$surv,
      lower  = 1-survcurv$upper,
      upper  = 1-survcurv$lower
    )
    return(predicted_P0)
  }
  
  ###    Predicted P_Y(A=1) (at the end of follow-up)
  ###    For everyone always treated senario
  .predict_P1_ipw <- function(cox_ipw = cox_ipw, cutoff = NULL, nknots = NULL,
                              Constrained = 'no', spline_order = 4, n.visit = 100,
                              n_full_expo = FALSE){
    if (!isTRUE(n_full_expo)) return(NA_real_)
    predpop1 <- data.frame(matrix(nrow = n*n.visit, ncol = 4))
    # id 
    predpop1[,1] <- rep(1:n, each=n.visit)
    # time
    predpop1[,2] <- rep(0:(n.visit-1), n)
    # timestop
    predpop1[,3] <- rep(1:n.visit, n)
    # event
    predpop1[,4] <- rep(0, nrow(predpop1))
    colnames(predpop1) = c("id", "time","time.stop", "event")
    
    predpop1$A <- 1
    if (is.null(cutoff) == FALSE & is.null(nknots) == FALSE){
      for (i in 1:cutoff) {
        predpop1=predpop1 %>%
          group_by(id) %>%
          mutate("Alag{i}":= dplyr::lag(A, i))
        predpop1[is.na(predpop1[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
      }  
      # get the spline basis in the original model
      time_since_expo <- seq(0,cutoff,1)
      knots_quantile <- seq(1/(nknots + 1), nknots/(nknots + 1), length.out = nknots)
      inner_knots = round(quantile(time_since_expo, knots_quantile), 0)    
      bsknots_augm <- c((-spline_order+1):0, inner_knots, cutoff: (cutoff + spline_order -1))
      basis_matrix <- splines::splineDesign(x = time_since_expo, knots = bsknots_augm,  ord = spline_order)
      # calculated time varying wce variable
      for (j in 1:ncol(basis_matrix)){
        for (i in 0:cutoff){
          if (i==0){
            wcej <- predpop1[, "A"] * basis_matrix[i+1,j]
          }
          else if (i>=1){
            expo_temp <- predpop1[,paste0("Alag",i)] * basis_matrix[i+1,j]
            wcej <- wcej + expo_temp
          }
          predpop1[, paste0("WCEvar",j)] <- wcej
        } # end of i loop for each basis function
      } # end of j loop for all basis
    }
    else if (is.null(cutoff) == FALSE){
      for (i in 1:n.visit) {
        predpop1=predpop1 %>%
          group_by(id) %>%
          mutate("Alag{i}":= dplyr::lag(A, i))
        predpop1[is.na(predpop1[,paste0("Alag",i)] ) , paste0("Alag",i)] <- 0
      }  
      predpop1 <- predpop1 %>%
        mutate(Asum_half = rowSums(across(Alag1:Alag50)))
      predpop1$Asum_half <- predpop1$Asum_half + predpop1$A
      
      predpop1 <- predpop1 %>%
        mutate(Asum_whole = rowSums(across(Alag1:Alag100)))
      predpop1$Asum_whole <- predpop1$Asum_whole + predpop1$A
    }
    cumhaz <- predict(cox_ipw, predpop1, type = "expected")
    predpop1$predict1 <- exp(-cumhaz)
    predpop1$surv1 <- ave(predpop1$predict1, predpop1$id,FUN=cumprod)
    predpop1$surv1cum <- ave(predpop1$surv1, predpop1$id, FUN = min)
    Pr_A1_ipw <- 1- mean(predpop1$surv1cum) 
    return(Pr_A1_ipw)
  } 
  
  # cox ipw null
  cox_null <- coxph(Surv(time, time.stop, event) ~ 1, data = dat.complet, 
                    cluster = id, control = coxph.control(timefix = FALSE))
  surv_null <- survfit(cox_null)
  KM_estimator <- summary(surv_null, times = n.visit)$surv
  KM_lower <- summary(survfit(cox_null), times = n.visit)$lower
  KM_upper <-summary(survfit(cox_null), times = n.visit)$upper
  nevent <- KM_estimator*n
  
  cox_ipw_null <- coxph(Surv(time, time.stop, event) ~ 1, data = dat.complet, 
                        cluster = id, weights = ipw.s, control = coxph.control(timefix = FALSE))
  cumhz_end_fup <- basehaz(cox_ipw_null)[nrow(basehaz(cox_ipw_null)),"hazard"]
  null_coefficients <- matrix(0, nrow = 1, ncol = 6)
  colnames(null_coefficients)<- c("coef", "exp(coef)", "se(coef)",  "robust se" ,"z","Pr(>|z|)")
  cox_ipw_null <- list(
    coefficients = null_coefficients,
    BIC = BIC(cox_ipw_null),
    AIC = AIC(cox_ipw_null),
    nevent = cox_ipw_null$nevent,
    proba_observed = KM_estimator,
    Pr_Y_A1 = 1-exp(-cumhz_end_fup),
    Pr_Y_A0 = 1-exp(-cumhz_end_fup),
    KM_lower = KM_lower,
    KM_upper = KM_upper
  )
  # cox unadjusted, current
  cox_unadj_current <- coxph(Surv(time, time.stop, event) ~ A, data = dat.complet, 
                             cluster = id, control = coxph.control(timefix = FALSE))
  cox_unadj_current <- list(
    coefficients = summary(cox_unadj_current)$coefficients, 
    Pr_Y_A1 = .predict_P1_ipw(cox_unadj_current),
    Pr_Y_A0 = .predict_P0_ipw(cox_unadj_current),
    BIC = BIC(cox_unadj_current),
    AIC = AIC(cox_unadj_current)
  )
  # cox unadjusted, cum half
  cox_unadj_cum_half <- coxph(Surv(time, time.stop, event) ~ Asum_half, data = dat.complet, 
                              cluster = id, control = coxph.control(timefix = FALSE))
  cox_unadj_cum_half <- list(
    coefficients = summary(cox_unadj_cum_half)$coefficients, 
    Pr_Y_A1 = .predict_P1_ipw(cox_unadj_cum_half, cutoff = 50),
    Pr_Y_A0 = .predict_P0_ipw(cox_unadj_cum_half),
    BIC = BIC(cox_unadj_cum_half),
    AIC = AIC(cox_unadj_cum_half)
  )
  ##################
  ## The same as BIC()
  BIC_volinsky <- function(fit) {
    ll <- fit$loglik[2]
    k  <- length(coef(fit))
    d  <- fit$nevent
    -2 * ll + log(d) * k
  }
  # BIC_volinsky(cox_unadj_current)  # debug leftover: cox_unadj_current was
  # already overwritten into a list at line ~414, so fit$nevent is NULL here
  # and this line kills fit_models() with "non-numeric argument to log()".
  ###########################
  
  # cox unadjusted, wce
  cox_unadj_wce_1knots <- .ipw_wce(nknots = 1, cutoff = 50, weight_var = NULL, covariates = NULL, data = dat.complet, n = n)
  
  # cox adjusted, wce  (difficult to decide how to get the predicted probabilities for now)
  #  cox_adj_wce_1knots <- .ipw_wce(nknots = 1, cutoff = true_cutoff, weight_var = NULL, covariates = "L", dat.complet = dat.complet)
  #
  # cox_adj_current <- coxph(Surv(time, time.stop, event) ~ A + L, data = dat.complet, 
  #                          cluster = id, control = coxph.control(timefix = FALSE))
  #   
  # cox_L_effect <- coxph(Surv(time, time.stop, event) ~ L, data = dat.complet, 
  #                       cluster = id, weights = ipw.s, control = coxph.control(timefix = FALSE))
  
  # cox ipw current
  cox_ipw_current <- coxph(Surv(time, time.stop, event) ~ A, data = dat.complet, 
                           cluster = id, weights = ipw.s, control = coxph.control(timefix = FALSE))
  cox_ipw_current <- list(
    coefficients = summary(cox_ipw_current)$coefficients, 
    Pr_Y_A1 = .predict_P1_ipw(cox_ipw_current),
    Pr_Y_A0 = .predict_P0_ipw(cox_ipw_current),
    BIC = BIC(cox_ipw_current),
    AIC = AIC(cox_ipw_current),
    loglik = cox_ipw_current$loglik,
    vcov = cox_ipw_current$var,
    score_LRT = summary(cox_ipw_current)$logtest, 
    score_logrank_test = summary(cox_ipw_current)$sctest,
    score_waldtest = summary(cox_ipw_current)$waldtest,
    nevent = cox_ipw_current$nevent,
    proba_observed = KM_estimator,
    concordance = cox_ipw_current$concordance
  )
  # cox ipw cum
  cox_ipw_cum_half <- coxph(Surv(time, time.stop, event) ~ Asum_half, data = dat.complet, 
                            cluster = id, weights = ipw.s, control = coxph.control(timefix = FALSE))
  cox_ipw_cum_half <- list(
    coefficients = summary(cox_ipw_cum_half)$coefficients, 
    BIC = BIC(cox_ipw_cum_half),
    AIC = AIC(cox_ipw_cum_half),
    Pr_Y_A1 = .predict_P1_ipw(cox_ipw_cum_half, cutoff = 50),
    Pr_Y_A0 = .predict_P0_ipw(cox_ipw_cum_half),
    loglik = cox_ipw_cum_half$loglik,
    vcov = cox_ipw_cum_half$var,
    score_LRT = summary(cox_ipw_cum_half)$logtest, 
    score_logrank_test = summary(cox_ipw_cum_half)$sctest,
    score_waldtest = summary(cox_ipw_cum_half)$waldtest,
    proba_observed = KM_estimator,
    concordance = cox_ipw_cum_half$concordance
  )
  cox_ipw_cum_whole <- coxph(Surv(time, time.stop, event) ~ Asum_whole, data = dat.complet, 
                             cluster = id, weights = ipw.s, control = coxph.control(timefix = FALSE))
  cox_ipw_cum_whole <- list(
    coefficients = summary(cox_ipw_cum_whole)$coefficients, 
    BIC = BIC(cox_ipw_cum_whole),
    AIC = AIC(cox_ipw_cum_whole),
    Pr_Y_A1 = .predict_P1_ipw(cox_ipw_cum_whole, cutoff = 100),
    Pr_Y_A0 = .predict_P0_ipw(cox_ipw_cum_whole),
    loglik = cox_ipw_cum_whole$loglik,
    vcov = cox_ipw_cum_whole$var,
    score_LRT = summary(cox_ipw_cum_whole)$logtest, 
    score_logrank_test = summary(cox_ipw_cum_whole)$sctest,
    score_waldtest = summary(cox_ipw_cum_whole)$waldtest,
    proba_observed = KM_estimator,
    concordance = cox_ipw_cum_whole$concordance
  )
  
  # Cox ipw WCE half, unconstrained
  cox_ipw_wce_1knots_half <- .ipw_wce(nknots = 1, cutoff = 50, weight_var = "ipw.s", data = dat.complet, n = n)
  cox_ipw_wce_2knots_half <- .ipw_wce(nknots = 2, cutoff = 50, weight_var = "ipw.s", data = dat.complet, n = n)
  cox_ipw_wce_3knots_half <- .ipw_wce(nknots = 3, cutoff = 50, weight_var = "ipw.s", data = dat.complet, n = n)
  
  # Cox ipw WCE half right constrained
  cox_ipw_wce_1knots_half_right <- .ipw_wce(nknots = 1, cutoff = 50, weight_var = "ipw.s",
                                            Constrained = 'right', data = dat.complet, n = n)
  cox_ipw_wce_2knots_half_right <- .ipw_wce(nknots = 2, cutoff = 50, weight_var = "ipw.s",
                                            Constrained = 'right', data = dat.complet, n = n)
  cox_ipw_wce_3knots_half_right <- .ipw_wce(nknots = 3, cutoff = 50, weight_var = "ipw.s",
                                            Constrained = 'right', data = dat.complet, n = n)
  
  # Cox ipw WCE half left constrained
  cox_ipw_wce_1knots_half_left <- .ipw_wce(nknots = 1, cutoff = 50, weight_var = "ipw.s",
                                           Constrained = 'left', data = dat.complet, n = n)
  cox_ipw_wce_2knots_half_left <- .ipw_wce(nknots = 2, cutoff = 50, weight_var = "ipw.s",
                                           Constrained = 'left', data = dat.complet, n = n)
  cox_ipw_wce_3knots_half_left <- .ipw_wce(nknots = 3, cutoff = 50, weight_var = "ipw.s",
                                           Constrained = 'left', data = dat.complet, n = n)
  # Cox ipw WCE whole, unconstrained
  cox_ipw_wce_1knots_whole <- .ipw_wce(nknots = 1, cutoff = n.visit, weight_var = "ipw.s", data = dat.complet, n = n)
  cox_ipw_wce_2knots_whole <- .ipw_wce(nknots = 2, cutoff = n.visit, weight_var = "ipw.s", data = dat.complet, n = n)
  cox_ipw_wce_3knots_whole <- .ipw_wce(nknots = 3, cutoff = n.visit, weight_var = "ipw.s", data = dat.complet, n = n)
  
  # 17-19 WCE whole, right constrained
  cox_ipw_wce_1knots_whole_right <- .ipw_wce(nknots = 1, cutoff = n.visit, weight_var = "ipw.s",
                                             Constrained = 'right', data = dat.complet, n = n)
  cox_ipw_wce_2knots_whole_right <- .ipw_wce(nknots = 2, cutoff = n.visit, weight_var = "ipw.s",
                                             Constrained = 'right', data = dat.complet, n = n)
  cox_ipw_wce_3knots_whole_right <- .ipw_wce(nknots = 3, cutoff = n.visit, weight_var = "ipw.s",
                                             Constrained = 'right', data = dat.complet, n = n)
  # Cox ipw WCE whole left constrained
  cox_ipw_wce_1knots_whole_left <- .ipw_wce(nknots = 1, cutoff = n.visit, weight_var = "ipw.s",
                                            Constrained = 'left', data = dat.complet, n = n)
  cox_ipw_wce_2knots_whole_left <- .ipw_wce(nknots = 2, cutoff = n.visit, weight_var = "ipw.s",
                                            Constrained = 'left', data = dat.complet, n = n)
  cox_ipw_wce_3knots_whole_left <- .ipw_wce(nknots = 3, cutoff = n.visit, weight_var = "ipw.s",
                                            Constrained = 'left', data = dat.complet, n = n)
  
  list_model <- list(
    cox_ipw_null = cox_ipw_null,
    cox_unadj_current = cox_unadj_current, 
    cox_unadj_cum_half = cox_unadj_cum_half, 
    cox_unadj_wce_1knots = cox_unadj_wce_1knots, 
    cox_ipw_current = cox_ipw_current,
    cox_ipw_cum_half = cox_ipw_cum_half,
    cox_ipw_cum_whole = cox_ipw_cum_whole,
    cox_ipw_wce_1knots_half = cox_ipw_wce_1knots_half,
    cox_ipw_wce_2knots_half = cox_ipw_wce_2knots_half,
    cox_ipw_wce_3knots_half = cox_ipw_wce_3knots_half,
    cox_ipw_wce_1knots_half_right = cox_ipw_wce_1knots_half_right,
    cox_ipw_wce_2knots_half_right = cox_ipw_wce_2knots_half_right,
    cox_ipw_wce_3knots_half_right = cox_ipw_wce_3knots_half_right,
    cox_ipw_wce_1knots_half_left = cox_ipw_wce_1knots_half_left,
    cox_ipw_wce_2knots_half_left = cox_ipw_wce_2knots_half_left,
    cox_ipw_wce_3knots_half_left = cox_ipw_wce_3knots_half_left,
    cox_ipw_wce_1knots_whole = cox_ipw_wce_1knots_whole,
    cox_ipw_wce_2knots_whole = cox_ipw_wce_2knots_whole,
    cox_ipw_wce_3knots_whole = cox_ipw_wce_3knots_whole,
    cox_ipw_wce_1knots_whole_right = cox_ipw_wce_1knots_whole_right,
    cox_ipw_wce_2knots_whole_right = cox_ipw_wce_2knots_whole_right,
    cox_ipw_wce_3knots_whole_right = cox_ipw_wce_3knots_whole_right,
    cox_ipw_wce_1knots_whole_left = cox_ipw_wce_1knots_whole_left,
    cox_ipw_wce_2knots_whole_left = cox_ipw_wce_2knots_whole_left,
    cox_ipw_wce_3knots_whole_left = cox_ipw_wce_3knots_whole_left
  )
  
  return(list_model)
}

# fit <- fit_models(dat.complet = dat.complet, n=1000, n.visit = 100)

# ============================================================
# ---- END: contents of 2_function_to_fit_models.R ----
# ---- BEGIN: parrelel version of 3_simulate_and_stock_model_outputs - censor ----
# ============================================================

# Helper function to extract only essential results from model fits
extract_essential_results <- function(model_list) {
  essential_results <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract only essential elements, not full model objects
    essential_results[[model_name]] <- list(
      coefficients = model$coefficients,
      Pr_Y_A1 = model$Pr_Y_A1,
      Pr_Y_A0 = model$Pr_Y_A0,
      BIC = model$BIC,
      AIC = model$AIC,
      proba_observed = model$proba_observed
    )
    
    # Include additional elements if they exist
    if (!is.null(model$coef_weight_predict)) {
      essential_results[[model_name]]$coef_weight_predict <- model$coef_weight_predict
    }
    if (!is.null(model$coef_sd)) {
      essential_results[[model_name]]$coef_sd <- model$coef_sd
    }
    if (!is.null(model$loglik)) {
      essential_results[[model_name]]$loglik <- model$loglik
    }
    if (!is.null(model$concordance)) {
      essential_results[[model_name]]$concordance <- model$concordance
    }
  }
  
  return(essential_results)
}

#  Save each repetition before returning it to the master process.
# A temporary file in the same directory prevents partial files being reused.
atomic_save_rds <- function(object, path) {
  tmp <- tempfile(
    pattern = paste0(basename(path), ".tmp_"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = FALSE)
  
  if (!file.rename(tmp, path)) {
    stop("Could not finalize result file: ", path)
  }
  
  invisible(path)
}

read_valid_result <- function(path, config, sim_id) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(x) || !identical(x$status, "success") ||
      !identical(x$config, config) || !identical(x$sim_id, as.integer(sim_id)) ||
      !is.list(x$model_results) || length(x$model_results) != 25L) return(NULL)
  x
}

run_main_repetition <- function(sim_id, config, replicate_dir) {
  sim_id <- as.integer(sim_id)
  seed <- as.integer(sim_id * config$seed_multiplier)
  path <- file.path(replicate_dir, sprintf("main_%s_sim%03d_results.rds",
                                           config$weight_hypo, sim_id))
  do.call(RNGkind, as.list(config$rng_kind))
  started <- Sys.time()
  warning_messages <- character()
  gc()
  result <- tryCatch(withCallingHandlers({
    # [REV] Explicit seed bypasses the random seed drawn inside sim_MSM_WCE.
    args <- config$data_args
    args$seed <- seed
    dat.list <- do.call(sim_MSM_WCE, args)
    dat.complet <- dat.list$dat.complet
    rm(dat.list)
    gc()
    fitted <- fit_models(dat.complet = dat.complet,
                         n = args$n, n.visit = args$n.visit)
    rm(dat.complet)
    model_results <- extract_essential_results(fitted)
    rm(fitted)
    gc()
    if (length(model_results) != 25L) stop("Expected 25 fitted model summaries.")
    list(status = "success", model_results = model_results, error_message = NA_character_)
  }, warning = function(w) {
    warning_messages <<- c(warning_messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    list(status = "failed", model_results = NULL, error_message = conditionMessage(e))
  })
  result$sim_id <- sim_id
  result$seed <- seed
  result$config <- config
  result$warnings <- unique(warning_messages)
  result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  # Small, uncompressed checkpoint: no full Cox objects or simulated datasets.
  atomic_save_rds(result, path)
  list(sim_id = sim_id, status = result$status, elapsed_seconds = result$elapsed_seconds,
       error_message = result$error_message)
}

# [REV] Parallel driver with fixed worker count, checkpoints and resume.
sim_and_stock_model_outputs <- function(n_sim = 500, weight_hypo = "decreasing", true_cutoff = 50, n = 1000,
                                        n.visit = 100, LtoA = 1, AtoL = 0.5, logitpL = 0, LtoL = 1,
                                        lambda0 = 0.005, HRA = 4, T0toL = 1,
                                        logitpA = -3, AtoA = 4, n.obs = 5, n_cores = 30,
                                        censor_max = NULL, output_root = "main_analysis_output",
                                        scenario_id = 1L, chunk_size = n_cores,
                                        seed_multiplier = 12345L) {
  for (value in list(n_sim, n_cores, chunk_size, scenario_id, seed_multiplier)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value))
      stop("Repetition, worker, chunk, scenario and seed settings must be positive integers.")
  }
  if (n_sim * seed_multiplier > .Machine$integer.max) stop("Seed exceeds R integer range.")
  sim_ids <- seq_len(as.integer(n_sim))
  hr_tag <- gsub(".", "", format(HRA, scientific = FALSE, trim = TRUE), fixed = TRUE)
  censor_tag <- if (is.null(censor_max)) "no_censor" else "censor"
  run_dir <- file.path(output_root, censor_tag,
                       sprintf("scenario%d_hr%s_n%d_T0L%s_lowrisk", scenario_id, hr_tag, n,
                               format(T0toL, scientific = FALSE, trim = TRUE)))
  replicate_dir <- file.path(run_dir, "replicates")
  log_dir <- file.path(run_dir, "logs")
  dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  replicate_dir <- normalizePath(replicate_dir, mustWork = TRUE)
  data_args <- list(weight_hypo = weight_hypo, true_cutoff = true_cutoff, n = n,
                    n.visit = n.visit, LtoA = LtoA, AtoL = AtoL, logitpL = logitpL, LtoL = LtoL,
                    lambda0 = lambda0, HRA = HRA, T0toL = T0toL, logitpA = logitpA, AtoA = AtoA,
                    n.obs = n.obs, censor_max = censor_max)
  # Code and data settings must match before old checkpoints can be reused.
  # n_sim and worker count are excluded so the run can be extended/reallocated.
  config <- list(format_version = 1L, scenario_id = scenario_id, weight_hypo = weight_hypo,
                 data_args = data_args, seed_multiplier = seed_multiplier,
                 rng_kind = RNGkind(), R_version = as.character(getRversion()),
                 package_versions = vapply(c("survival", "dplyr", "tidyr", "MASS", "splines", "zoo"),
                                           function(pkg) as.character(packageVersion(pkg)), character(1)),
                 code = lapply(list(sim_MSM_WCE, fit_models, .ipw_wce, extract_essential_results),
                               function(fn) paste(deparse(fn), collapse = "\n")))
  config_path <- file.path(run_dir, "run_config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), config))
      stop("Existing output uses different code/settings. Choose a new output_root: ", run_dir)
  } else atomic_save_rds(config, config_path)
  files <- file.path(replicate_dir, sprintf("main_%s_sim%03d_results.rds", weight_hypo, sim_ids))
  ready <- vapply(seq_along(files), function(j)
    !is.null(read_valid_result(files[j], config, sim_ids[j])), logical(1))
  pending <- sim_ids[!ready]
  cat(sprintf("Scenario %d: %d requested, %d reusable, %d pending; %d workers requested.\n",
              scenario_id, n_sim, sum(ready), length(pending), n_cores))
  if (length(pending)) {
    # [REV] One PSOCK cluster per scenario. Workers are local to the allocated node.
    cl <- parallel::makeCluster(min(as.integer(n_cores), length(pending)), type = "PSOCK",
                               outfile = file.path(log_dir, sprintf("workers_R%d.log", n_sim)))
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, {
      library(tidyverse); library(survival); library(MASS)
      library(stats); library(stringr); library(matrixStats); library(splines); library(zoo)
      NULL
    })
    parallel::clusterExport(cl, c("sim_MSM_WCE", "fit_models", ".ipw_wce", ".predict_P0_ipw",
                                  "compute_ipw", "extract_essential_results", "atomic_save_rds",
                                  "run_main_repetition"), envir = environment())
    chunks <- split(pending, ceiling(seq_along(pending) / as.integer(chunk_size)))
    for (chunk_i in seq_along(chunks)) {
      batch <- chunks[[chunk_i]]
      cat(sprintf("Chunk %d/%d: sim %d to %d\n", chunk_i, length(chunks), min(batch), max(batch)))
      chunk_status <- foreach(i_sim = batch, .errorhandling = "pass",
                              .options.snow = list(preschedule = FALSE)) %dopar% {
        run_main_repetition(i_sim, config, replicate_dir)
      }
      for (j in seq_along(chunk_status)) {
        x <- chunk_status[[j]]
        if (inherits(x, "error")) cat("Worker/checkpoint error for sim", batch[j], ":", conditionMessage(x), "\n")
        else cat("sim", x$sim_id, x$status, "elapsed", round(x$elapsed_seconds, 1), "seconds\n")
      }
      rm(chunk_status)
      gc()
    }
  }
  # [REV] Re-read checkpoints in sim_id order; never silently drop failures.
  results <- lapply(seq_along(files), function(j) read_valid_result(files[j], config, sim_ids[j]))
  status_rows <- lapply(seq_along(files), function(j) {
    x <- tryCatch(readRDS(files[j]), error = function(e) NULL)
    valid <- !is.null(results[[j]])
    data.frame(sim_id = sim_ids[j], seed = sim_ids[j] * seed_multiplier,
               status = if (valid) "success" else if (is.list(x) && identical(x$status, "failed")) "failed" else "missing_or_invalid",
               elapsed_seconds = if (is.list(x) && !is.null(x$elapsed_seconds)) x$elapsed_seconds else NA_real_,
               error_message = if (is.list(x) && !is.null(x$error_message)) x$error_message else NA_character_,
               warnings = if (is.list(x)) paste(x$warnings, collapse = " | ") else "",
               stringsAsFactors = FALSE)
  })
  status_table <- do.call(rbind, status_rows)
  write.csv(status_table, file.path(run_dir, sprintf("run_status_R%d.csv", n_sim)), row.names = FALSE)
  good <- vapply(results, function(x) !is.null(x), logical(1))
  if (!all(good)) stop("No combined RDS written: ", sum(!good),
                        " repetitions failed/missing. Inspect run_status and rerun this script.")
  # Preserve the original downstream interface: list_scenario$list_model_output.
  combined <- list(list_model_output = lapply(results, function(x) x$model_results))
  combined_name <- sprintf("weight%d_half_hr%s_n_%d_nsim_%d_T0L%s_lowrisk_CI.rds",
                           scenario_id, hr_tag, n, n_sim,
                           format(T0toL, scientific = FALSE, trim = TRUE))
  saveRDS(combined, file.path(run_dir, combined_name), compress = TRUE)
  cat("All", n_sim, "repetitions successful. Combined RDS:", file.path(run_dir, combined_name), "\n")
  invisible(list(output_file = file.path(run_dir, combined_name), n_success = sum(good)))
}

# ============================================================
# Run from the same parent working directory used for bootstrap output.
# ============================================================

N_CORES <- 40
N_SIM <- 500
N_SUBJ <- 10000
HRA <- 1.5

SCENARIO_IDS <- 1:4             
OUTPUT_ROOT <- "main_analysis_output"
CHUNK_SIZE <- N_CORES
SEED_MULTIPLIER <- 12345

scenario_names <- c("decreasing", "middle_peak", "constant", "only_current")
for (scenario_id in SCENARIO_IDS) {
  if (!(scenario_id %in% seq_along(scenario_names))) stop("Invalid scenario ID.")
  start_time <- Sys.time()
  tryCatch({
    sim_and_stock_model_outputs(
      n_sim = N_SIM, n = N_SUBJ, weight_hypo = scenario_names[scenario_id],
      lambda0 = 0.0005, true_cutoff = 50, HRA = HRA, T0toL = 1,
      n_cores = N_CORES, censor_max = 5, output_root = OUTPUT_ROOT,
      scenario_id = as.integer(scenario_id), chunk_size = CHUNK_SIZE,
      seed_multiplier = SEED_MULTIPLIER
    )
  }, error = function(e) {
    message("Scenario ", scenario_id, " incomplete: ", conditionMessage(e))
    stop(e)  # nonzero Rscript exit; completed checkpoints remain available
  })
  cat("Scenario", scenario_id, "total time:",
      round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "minutes\n")
}
