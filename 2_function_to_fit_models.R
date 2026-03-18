# Function to fit all models and return model objects on an input dataset
# input: the dataset(list) simulated by sim_MSM_WCE_cond() function

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
                     data = NULL){
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
  
  ### Predict Pr_Y(A=1)
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
                             Constrained = 'no', spline_order = 4, n.visit = 100){
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
  ## BIC() gives the same output as volinsky's  adapted BIC
  BIC_volinsky <- function(fit) {
    ll <- fit$loglik[2]
    k  <- length(coef(fit))
    d  <- fit$nevent
    -2 * ll + log(d) * k
  }
  BIC_volinsky(cox_unadj_current)
  ###########################
  
  # cox unadjusted, wce
  cox_unadj_wce_1knots <- .ipw_wce(nknots = 1, cutoff = 50, weight_var = NULL, covariates = NULL, data = dat.complet)
  
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

