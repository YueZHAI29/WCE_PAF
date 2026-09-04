# ============================================================
# master bootstrap file for cluster - PARALLEL VERSION (SELF-CONTAINED)
# All function definitions (from 1_function_to_generate_data.R and the
# FIXED/lite 2_function_to_fit_models.R) are pasted directly below,
# no source() calls -- so this single file is all your colleague needs
# to launch the job.
#
# Structure mirrors parallel_simulation.r / Copy_parallel_simulation2.r:
#   makeCluster -> clusterEvalQ (load pkgs) -> clusterExport (export fns)
#   -> foreach %dopar% (.errorhandling="pass") -> stopCluster
# with per-replicate gc() hygiene as in the n>=50000 chunked version.
# ============================================================

# Load necessary libraries
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
# ---- BEGIN: bootstrap driver functions and parallel execution ----
# ============================================================


# ============================================================
# ---- resampler: fixes the id-collision bug (fresh sequential ids per draw) ----
# ============================================================
resample_subjects <- function(data, ids) {
  idx_by_id <- split(seq_len(nrow(data)), data$id)
  drawn     <- idx_by_id[as.character(ids)]
  out       <- data[unlist(drawn, use.names = FALSE), ]
  out$id    <- rep(seq_along(drawn), lengths(drawn))
  out
}

# ============================================================
# ---- fit the non-null candidate set, record both BIC and AIC selections ----
# ============================================================
fit_and_summarize <- function(data, n, n.visit = 100) {
  fit <- fit_models(dat.complet = data, n = n, n.visit = n.visit)
  all_ipw_names <- grep("ipw", names(fit), value = TRUE)

  # cox_ipw_null is retained in fit because its proba_observed and nevent
  # components are needed below. It is excluded ONLY from model selection.
  # All constant/non-null and WCE candidates remain eligible.
  ipw_names <- setdiff(all_ipw_names, "cox_ipw_null")
  if ("cox_ipw_null" %in% ipw_names || length(ipw_names) == 0L) {
    stop("Failed to construct the non-null model-selection candidate set.")
  }

  bic_vals <- sapply(fit[ipw_names], function(m) m[["BIC"]])
  aic_vals <- sapply(fit[ipw_names], function(m) m[["AIC"]])
  names(bic_vals) <- names(aic_vals) <- ipw_names

  best_bic_name <- names(which.min(bic_vals))
  best_aic_name <- names(which.min(aic_vals))

  if (any(c(best_bic_name, best_aic_name) == "cox_ipw_null")) {
    stop("Internal error: the null model entered BIC/AIC model selection.")
  }

  extract_P0 <- function(model_name) {
    tmp <- fit[[model_name]][["Pr_Y_A0"]]
    if ("Pr_A0_ipw" %in% names(tmp)) {
      c(P0 = as.numeric(tmp["Pr_A0_ipw"]),
        P0_lower = as.numeric(tmp["lower"]),
        P0_upper = as.numeric(tmp["upper"]))
    } else {
      c(P0 = as.numeric(tmp), P0_lower = NA_real_, P0_upper = NA_real_)
    }
  }

  P_observed <- 1 - fit[["cox_ipw_null"]][["proba_observed"]]
  n_events   <- fit[["cox_ipw_null"]][["nevent"]]

  paf_from <- function(model_name) {
    p0 <- extract_P0(model_name)
    PAF       <- (P_observed - p0["P0"])       / P_observed
    PAF_lower <- (P_observed - p0["P0_upper"]) / P_observed
    PAF_upper <- (P_observed - p0["P0_lower"]) / P_observed
    c(p0, PAF = unname(PAF), PAF_lower = unname(PAF_lower), PAF_upper = unname(PAF_upper))
  }

  extract_diag <- function(model_name) {
    coefs <- fit[[model_name]][["coefficients"]]
    c(max_abs_coef = if (!is.null(coefs) && nrow(coefs) > 0) max(abs(coefs[, "coef"])) else NA_real_,
      max_se       = if (!is.null(coefs) && nrow(coefs) > 0) max(coefs[, "se(coef)"]) else NA_real_)
  }

  bic_res <- paf_from(best_bic_name);  bic_diag <- extract_diag(best_bic_name)
  aic_res <- paf_from(best_aic_name);  aic_diag <- extract_diag(best_aic_name)

  data.frame(
    model_BIC         = best_bic_name,
    BIC_value         = min(bic_vals),
    P0_BIC            = bic_res["P0"],
    P0_lower_BIC      = bic_res["P0_lower"],
    P0_upper_BIC      = bic_res["P0_upper"],
    PAF_BIC           = bic_res["PAF"],
    PAF_lower_BIC     = bic_res["PAF_lower"],
    PAF_upper_BIC     = bic_res["PAF_upper"],
    max_abs_coef_BIC  = bic_diag["max_abs_coef"],
    max_se_BIC        = bic_diag["max_se"],
    flag_unstable_BIC = isTRUE(bic_diag["max_abs_coef"] > 5 | bic_diag["max_se"] > 1),
    flag_implausible_PAF_BIC = isTRUE(abs(bic_res["PAF"]) > 1),
    model_AIC         = best_aic_name,
    AIC_value         = min(aic_vals),
    P0_AIC            = aic_res["P0"],
    P0_lower_AIC      = aic_res["P0_lower"],
    P0_upper_AIC      = aic_res["P0_upper"],
    PAF_AIC           = aic_res["PAF"],
    PAF_lower_AIC     = aic_res["PAF_lower"],
    PAF_upper_AIC     = aic_res["PAF_upper"],
    max_abs_coef_AIC  = aic_diag["max_abs_coef"],
    max_se_AIC        = aic_diag["max_se"],
    flag_unstable_AIC = isTRUE(aic_diag["max_abs_coef"] > 5 | aic_diag["max_se"] > 1),
    flag_implausible_PAF_AIC = isTRUE(abs(aic_res["PAF"]) > 1),
    P_observed        = P_observed,
    n_events          = n_events,
    same_model_selected = (best_bic_name == best_aic_name),
    n_selection_candidates = length(ipw_names),
    null_model_excluded = TRUE,
    row.names = NULL
  )
}

# ============================================================
# bootstrap_PAF_parallel: parallel inner bootstrap loop for ONE dataset
# ============================================================
bootstrap_PAF_parallel <- function(d, scenario, sim_id, B, n.visit = 100, n_cores = NULL,
                                   out_dir = "output_n1000_bootstrap") {

  if (is.null(n_cores)) n_cores <- min(parallel::detectCores() - 1, B)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  subject_ids <- unique(d$id)
  n_id        <- length(subject_ids)

  # point estimate on the original (unresampled) dataset
  point_res <- fit_and_summarize(d, n = n_id, n.visit = n.visit)
  point_res$scenario <- scenario
  point_res$sim_id   <- sim_id
  point_res$boot_id  <- 0L
  point_res$type     <- "point_estimate"
  cat(sprintf(paste0(
    "[%s | sim %d] point estimates | ",
    "BIC PAF: %.4f (%s) | AIC PAF: %.4f (%s)\n"
  ), scenario, sim_id,
  point_res$PAF_BIC, point_res$model_BIC,
  point_res$PAF_AIC, point_res$model_AIC))

  set.seed(sim_id * 1000L)
  boot_id_list <- lapply(1:B, function(b) sample(subject_ids, size = n_id, replace = TRUE))

  cat(sprintf("[%s | sim %d] setting up parallel processing with %d cores...\n", scenario, sim_id, n_cores))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    library(tidyverse); library(survival); library(MASS); library(stats)
    library(stringr); library(matrixStats); library(splines); library(zoo)
  })
  # Only export FUNCTION NAMES (see note): .ipw_wce is a top-level helper
  # that fit_models() reaches via lexical scoping, so it must be exported
  # explicitly -- everything else fit_models() needs is nested inside it
  # and travels with it automatically. Local loop variables (d, n_id,
  # boot_id_list, n.visit, scenario, sim_id) are picked up automatically
  # by foreach's own variable auto-capture -- do not clusterExport them.
  clusterExport(cl, c("sim_MSM_WCE", "fit_models", ".ipw_wce",
                      "resample_subjects", "fit_and_summarize"))

  prog <- file.path(out_dir, sprintf("bootstrap_%s_sim%03d_B%d_progress.log", scenario, sim_id, B))
  cat(paste0(
    "scenario\tsim_id\tboot_id\ttime_stamp\t",
    "PAF_BIC\tmodel_BIC\tPAF_AIC\tmodel_AIC\telapsed_min\n"
  ), file = prog)

  t_start <- Sys.time()
  cat(sprintf("[%s | sim %d] starting parallel bootstrap (B = %d)...\n", scenario, sim_id, B))

  boot_rows <- foreach(b = 1:B,
                       .packages = c("tidyverse", "survival", "MASS", "stats",
                                    "stringr", "matrixStats", "splines", "zoo"),
                       .errorhandling = "pass") %dopar% {

    gc()

    out <- tryCatch({
      boot_data <- resample_subjects(d, boot_id_list[[b]])
      r <- fit_and_summarize(boot_data, n = n_id, n.visit = n.visit)
      rm(boot_data); gc()
      r$scenario <- scenario
      r$sim_id   <- sim_id
      r$boot_id  <- b
      r$type     <- "bootstrap"
      list(row = r, error = NA_character_)
    }, error = function(e) list(row = NULL, error = conditionMessage(e)))

    paf_bic_str   <- if (!is.null(out$row)) round(out$row$PAF_BIC, 4) else NA
    model_bic_str <- if (!is.null(out$row)) out$row$model_BIC else paste0("FAILED: ", out$error)
    paf_aic_str   <- if (!is.null(out$row)) round(out$row$PAF_AIC, 4) else NA
    model_aic_str <- if (!is.null(out$row)) out$row$model_AIC else paste0("FAILED: ", out$error)
    cat(sprintf("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%.1f\n",
                scenario, sim_id, b, format(Sys.time(), "%H:%M:%S"),
                paf_bic_str, model_bic_str, paf_aic_str, model_aic_str,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
        file = prog, append = TRUE)
    out
  }

  stopCluster(cl)
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("[%s | sim %d] parallel bootstrap completed in %.1f minutes (%d cores)\n",
              scenario, sim_id, elapsed, n_cores))

  ok       <- sapply(boot_rows, function(x) is.null(x$error) || is.na(x$error))
  n_failed <- sum(!ok)
  if (n_failed > 0) {
    cat(sprintf("WARNING: %d of %d replicates failed\n", n_failed, B))
    print(head(sapply(boot_rows[!ok], function(x) x$error), 5))
  }

  boot_results <- do.call(rbind, lapply(boot_rows[ok], function(x) x$row))
  all_results  <- rbind(point_res, boot_results)

  out_rds <- file.path(out_dir, sprintf("bootstrap_%s_sim%03d_B%d_results.rds", scenario, sim_id, B))
  out_csv <- file.path(out_dir, sprintf("bootstrap_%s_sim%03d_B%d_results.csv", scenario, sim_id, B))
  saveRDS(list(scenario = scenario, sim_id = sim_id, all_results = all_results, point = point_res,
               elapsed_min = elapsed, n_failed = n_failed, n_cores = n_cores,
               boot_rows_raw = boot_rows),
          file = out_rds, compress = TRUE)
  write.csv(all_results, out_csv, row.names = FALSE)
  cat("Saved:", out_rds, "and .csv\n")

  paf_bic <- boot_results$PAF_BIC
  paf_aic <- boot_results$PAF_AIC
  cat(sprintf(paste0(
    "  -> %s sim %d | BIC: SE=%.4f, mean=%.4f, 95%%CI [%.4f, %.4f] | ",
    "AIC: SE=%.4f, mean=%.4f, 95%%CI [%.4f, %.4f] | ",
    "events=%.0f | failed=%d/%d\n\n"
  ), scenario, sim_id,
  sd(paf_bic), mean(paf_bic), quantile(paf_bic, .025), quantile(paf_bic, .975),
  sd(paf_aic), mean(paf_aic), quantile(paf_aic, .025), quantile(paf_aic, .975),
  mean(boot_results$n_events), n_failed, B))

  list(all_results = all_results, point = point_res, elapsed_min = elapsed, n_failed = n_failed)
}

# ============================================================
# ---- CONFIG: change these to scale up ----
# ============================================================
SCENARIOS      <- c("middle_peak")  # scenario 2
OUTER_SIM_IDS  <- 1:200             # STEP 1: set to 1:5 to calibrate timing, THEN switch to 1:200
INNER_B        <- 200               # B=200 is ample for the SE ratio; slightly conservative for coverage
N_CORES        <- 32
LAMBDA0        <- 0.0005   # KEPT AT THE MANUSCRIPT VALUE -- so the true PAF is unchanged (~0.45).
                            # Event count is raised by increasing n instead (see N_SUBJ), which
                            # leaves the true PAF invariant since it depends on lambda0, HRA and
                            # the weight function but NOT on sample size.
HRA            <- 4
N_SUBJ         <- 2000     # ~178 events expected; small enough to still address the "smaller n" request
TRUE_PAF       <- 0.402761 # exact manuscript value for weight_hypo="middle_peak"
                            #   P0_true = 0.048771
                            #   true_P_observed = 0.08166
                            #   true_PAF = (0.08166 - 0.048771) / 0.08166 = 0.402761
                            # Depends only on weight_hypo, lambda0, HRA -- NOT on n, so this value
                            # is unchanged from the manuscript's n=10000/50000 scenarios.
OUT_DIR        <- "output_n2000_scenario2_no_null_BIC_AIC"

cat("Starting n=2000 scenario-2 bootstrap study with non-null BIC/AIC selection...\n")
cat(sprintf("Plan: %d scenarios x %d datasets x %d inner bootstraps = %d total fits (+ %d point estimates)\n",
            length(SCENARIOS), length(OUTER_SIM_IDS), INNER_B,
            length(SCENARIOS) * length(OUTER_SIM_IDS) * INNER_B,
            length(SCENARIOS) * length(OUTER_SIM_IDS)))

start_time <- Sys.time()
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (scenario in SCENARIOS) {
  for (sim_id in OUTER_SIM_IDS) {

    check_rds <- file.path(OUT_DIR, sprintf("bootstrap_%s_sim%03d_B%d_results.rds", scenario, sim_id, INNER_B))
    if (file.exists(check_rds)) {
      cat(sprintf(">> %s sim_id=%d already done -- skipping\n", scenario, sim_id))
      next
    }

    cat(sprintf("\n========== scenario = %s | sim_id = %d ==========\n", scenario, sim_id))
    set.seed(sim_id)
    d <- sim_MSM_WCE(weight_hypo = scenario, true_cutoff = 50, n = N_SUBJ,
                     n.visit = 100, lambda0 = LAMBDA0, HRA = HRA, T0toL = 1,
                     seed = sim_id)$dat.complet
    cat("Event rate:", mean(tapply(d$event, d$id, max)),
        "| n_events:", sum(tapply(d$event, d$id, max)), "\n")

    bootstrap_PAF_parallel(d, scenario = scenario, sim_id = sim_id, B = INNER_B,
                           n.visit = 100, n_cores = N_CORES, out_dir = OUT_DIR)
  }
}

total_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\nTotal execution time: %.1f minutes (%.1f hours) for %d scenarios x %d datasets x %d inner\n",
            total_min, total_min / 60, length(SCENARIOS), length(OUTER_SIM_IDS), INNER_B))

# ============================================================
# ---- Manuscript-ready BIC/AIC summaries ----
# Both criteria use the SAME non-null candidate set. The null fit is used
# only to obtain the observed risk and event count; it cannot be selected.
# Each criterion receives its own point estimate, bootstrap SE, percentile
# interval, coverage indicator, diagnostics, and model-selection frequencies.
# ============================================================
expected_files <- unlist(lapply(SCENARIOS, function(scen) {
  file.path(OUT_DIR, sprintf("bootstrap_%s_sim%03d_B%d_results.rds",
                             scen, OUTER_SIM_IDS, INNER_B))
}))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0L) {
  stop(sprintf("Cannot summarize: %d expected result file(s) are missing. First missing file: %s",
               length(missing_files), missing_files[1]))
}

all_files <- expected_files
all_data  <- lapply(all_files, readRDS)

percentile_ci <- function(x) {
  if (sum(is.finite(x)) < 2L) return(c(lower = NA_real_, upper = NA_real_))
  unname(quantile(x[is.finite(x)], c(0.025, 0.975), na.rm = TRUE,
                  names = FALSE, type = 7))
}

# ---- per-dataset table: one row per outer simulation ----
per_dataset <- do.call(rbind, lapply(all_data, function(x) {
  boot_rows <- x$all_results$type == "bootstrap"
  bp_bic <- x$all_results$PAF_BIC[boot_rows]
  bp_aic <- x$all_results$PAF_AIC[boot_rows]
  ci_bic <- percentile_ci(bp_bic)
  ci_aic <- percentile_ci(bp_aic)

  data.frame(
    scenario       = x$scenario,
    sim_id         = x$sim_id,
    true_PAF       = TRUE_PAF,
    null_model_excluded = TRUE,

    point_model_BIC = as.character(x$point$model_BIC),
    point_PAF_BIC   = x$point$PAF_BIC,
    boot_mean_BIC   = mean(bp_bic, na.rm = TRUE),
    boot_median_BIC = median(bp_bic, na.rm = TRUE),
    boot_SE_BIC     = sd(bp_bic, na.rm = TRUE),
    ci_lower_BIC    = ci_bic[1],
    ci_upper_BIC    = ci_bic[2],
    ci_width_BIC    = ci_bic[2] - ci_bic[1],
    covered_BIC     = if (is.na(TRUE_PAF)) NA else
                        ci_bic[1] <= TRUE_PAF && TRUE_PAF <= ci_bic[2],

    point_model_AIC = as.character(x$point$model_AIC),
    point_PAF_AIC   = x$point$PAF_AIC,
    boot_mean_AIC   = mean(bp_aic, na.rm = TRUE),
    boot_median_AIC = median(bp_aic, na.rm = TRUE),
    boot_SE_AIC     = sd(bp_aic, na.rm = TRUE),
    ci_lower_AIC    = ci_aic[1],
    ci_upper_AIC    = ci_aic[2],
    ci_width_AIC    = ci_aic[2] - ci_aic[1],
    covered_AIC     = if (is.na(TRUE_PAF)) NA else
                        ci_aic[1] <= TRUE_PAF && TRUE_PAF <= ci_aic[2],

    BIC_AIC_same_point_model =
      as.character(x$point$model_BIC) == as.character(x$point$model_AIC),
    n_selection_candidates = x$point$n_selection_candidates,
    n_events       = x$point$n_events,
    n_failed       = x$n_failed,
    elapsed_min    = x$elapsed_min,
    row.names = NULL
  )
}))
per_dataset <- per_dataset[order(per_dataset$scenario, per_dataset$sim_id), ]

if (any(per_dataset$point_model_BIC == "cox_ipw_null" |
        per_dataset$point_model_AIC == "cox_ipw_null")) {
  stop("Validation failed: null model found among point selections.")
}

per_dataset_file <- file.path(OUT_DIR, "per_dataset_summary_BIC_AIC.csv")
write.csv(per_dataset, per_dataset_file, row.names = FALSE)

modal_model <- function(x) names(sort(table(x), decreasing = TRUE))[1]
modal_pct <- function(x) 100 * max(table(x)) / length(x)

# ---- scenario-level summary: BIC and AIC metrics in separate columns ----
summary_rows <- lapply(unique(per_dataset$scenario), function(scen) {
  pd        <- per_dataset[per_dataset$scenario == scen, ]
  scen_data <- Filter(function(x) x$scenario == scen, all_data)
  boot_all  <- do.call(rbind, lapply(scen_data, function(x) {
    x$all_results[x$all_results$type == "bootstrap", ]
  }))

  if (any(boot_all$model_BIC == "cox_ipw_null" |
          boot_all$model_AIC == "cox_ipw_null")) {
    stop(sprintf("Validation failed: null model found in bootstrap selections for %s.", scen))
  }

  n_out <- nrow(pd)
  emp_se_bic <- sd(pd$point_PAF_BIC)
  emp_se_aic <- sd(pd$point_PAF_AIC)
  cov_bic <- if (is.na(TRUE_PAF)) NA_real_ else 100 * mean(pd$covered_BIC)
  cov_aic <- if (is.na(TRUE_PAF)) NA_real_ else 100 * mean(pd$covered_AIC)
  cov_mcse_bic <- if (is.na(cov_bic)) NA_real_ else
    100 * sqrt((cov_bic / 100) * (1 - cov_bic / 100) / n_out)
  cov_mcse_aic <- if (is.na(cov_aic)) NA_real_ else
    100 * sqrt((cov_aic / 100) * (1 - cov_aic / 100) / n_out)

  data.frame(
    scenario                 = scen,
    n_outer_datasets         = n_out,
    n_boot_per_dataset       = INNER_B,
    n_subjects               = N_SUBJ,
    mean_events              = mean(pd$n_events),
    true_PAF                 = TRUE_PAF,
    null_model_excluded      = TRUE,
    n_selection_candidates   = unique(pd$n_selection_candidates)[1],
    n_failed_total           = sum(pd$n_failed),

    mean_point_PAF_BIC       = mean(pd$point_PAF_BIC),
    bias_BIC                 = if (is.na(TRUE_PAF)) NA_real_ else
                                 mean(pd$point_PAF_BIC) - TRUE_PAF,
    bias_MCSE_BIC            = emp_se_bic / sqrt(n_out),
    empirical_SE_BIC         = emp_se_bic,
    mean_bootstrap_SE_BIC    = mean(pd$boot_SE_BIC),
    SE_ratio_BIC             = mean(pd$boot_SE_BIC) / emp_se_bic,
    sd_bootstrap_SE_BIC      = sd(pd$boot_SE_BIC),
    coverage_pct_BIC         = cov_bic,
    coverage_MCSE_BIC        = cov_mcse_bic,
    mean_CI_width_BIC        = mean(pd$ci_width_BIC),
    pct_unstable_BIC         = 100 * mean(boot_all$flag_unstable_BIC, na.rm = TRUE),
    pct_implausible_PAF_BIC  = 100 * mean(boot_all$flag_implausible_PAF_BIC, na.rm = TRUE),
    n_unique_boot_models_BIC = length(unique(boot_all$model_BIC)),
    modal_boot_model_BIC     = modal_model(boot_all$model_BIC),
    modal_boot_model_pct_BIC = modal_pct(boot_all$model_BIC),
    n_unique_point_models_BIC = length(unique(pd$point_model_BIC)),
    modal_point_model_BIC     = modal_model(pd$point_model_BIC),
    modal_point_model_pct_BIC = modal_pct(pd$point_model_BIC),

    mean_point_PAF_AIC       = mean(pd$point_PAF_AIC),
    bias_AIC                 = if (is.na(TRUE_PAF)) NA_real_ else
                                 mean(pd$point_PAF_AIC) - TRUE_PAF,
    bias_MCSE_AIC            = emp_se_aic / sqrt(n_out),
    empirical_SE_AIC         = emp_se_aic,
    mean_bootstrap_SE_AIC    = mean(pd$boot_SE_AIC),
    SE_ratio_AIC             = mean(pd$boot_SE_AIC) / emp_se_aic,
    sd_bootstrap_SE_AIC      = sd(pd$boot_SE_AIC),
    coverage_pct_AIC         = cov_aic,
    coverage_MCSE_AIC        = cov_mcse_aic,
    mean_CI_width_AIC        = mean(pd$ci_width_AIC),
    pct_unstable_AIC         = 100 * mean(boot_all$flag_unstable_AIC, na.rm = TRUE),
    pct_implausible_PAF_AIC  = 100 * mean(boot_all$flag_implausible_PAF_AIC, na.rm = TRUE),
    n_unique_boot_models_AIC = length(unique(boot_all$model_AIC)),
    modal_boot_model_AIC     = modal_model(boot_all$model_AIC),
    modal_boot_model_pct_AIC = modal_pct(boot_all$model_AIC),
    n_unique_point_models_AIC = length(unique(pd$point_model_AIC)),
    modal_point_model_AIC     = modal_model(pd$point_model_AIC),
    modal_point_model_pct_AIC = modal_pct(pd$point_model_AIC),

    pct_BIC_AIC_agree_boot   = 100 * mean(boot_all$same_model_selected, na.rm = TRUE),
    pct_BIC_AIC_agree_point  = 100 * mean(pd$BIC_AIC_same_point_model),
    row.names = NULL
  )
})
summary_table <- do.call(rbind, summary_rows)

# ---- long-form model-selection frequency table ----
frequency_rows <- lapply(unique(per_dataset$scenario), function(scen) {
  pd        <- per_dataset[per_dataset$scenario == scen, ]
  scen_data <- Filter(function(x) x$scenario == scen, all_data)
  boot_all  <- do.call(rbind, lapply(scen_data, function(x) {
    x$all_results[x$all_results$type == "bootstrap", ]
  }))

  make_freq <- function(models, criterion, selection_level) {
    counts <- sort(table(models), decreasing = TRUE)
    data.frame(
      scenario = scen,
      criterion = criterion,
      selection_level = selection_level,
      model = names(counts),
      count = as.integer(counts),
      percent = 100 * as.integer(counts) / sum(counts),
      row.names = NULL
    )
  }

  rbind(
    make_freq(boot_all$model_BIC, "BIC", "bootstrap"),
    make_freq(boot_all$model_AIC, "AIC", "bootstrap"),
    make_freq(pd$point_model_BIC, "BIC", "point_estimate"),
    make_freq(pd$point_model_AIC, "AIC", "point_estimate")
  )
})
model_frequency <- do.call(rbind, frequency_rows)

cat("\n================ BIC/AIC MANUSCRIPT SUMMARY ================\n")
print(summary_table, digits = 4)
cat("\nModel-selection frequencies:\n")
print(model_frequency, row.names = FALSE)

summary_file <- file.path(OUT_DIR, "manuscript_summary_table_BIC_AIC.csv")
frequency_file <- file.path(OUT_DIR, "model_selection_frequency_BIC_AIC.csv")
rds_summary_file <- file.path(OUT_DIR, "manuscript_summary_BIC_AIC.rds")
write.csv(summary_table, summary_file, row.names = FALSE)
write.csv(model_frequency, frequency_file, row.names = FALSE)
saveRDS(list(summary_table = summary_table,
             per_dataset = per_dataset,
             model_frequency = model_frequency,
             all_data = all_data),
        rds_summary_file)
cat("\nSaved:\n",
    per_dataset_file, "\n",
    summary_file, "\n",
    frequency_file, "\n",
    rds_summary_file, "\n")
