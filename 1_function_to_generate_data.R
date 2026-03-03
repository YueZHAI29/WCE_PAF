# Function to simulate one dataset based on Marginal Structural Cox with WCE modeling
# Data generation mechanism according to Young et al.2010 and Xiao et al. 2014
# effect assumption and weight assumption made on marginal effect
# with binary time-varying confounder affected by previous exposure
# with optional missing confounder 
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

