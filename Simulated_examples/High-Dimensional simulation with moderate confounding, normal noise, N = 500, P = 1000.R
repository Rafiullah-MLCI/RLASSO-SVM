######################################
# High-Dimensional RLASSO_SVM simulation, moderate confounding, normal noise #
######################################

rm(list = ls())


# Libraries----------------------------------------

library(tidyverse)
library(BART)
library(grf)
library(rlearner)
library(future)
library(bartCause)
library(MASS)
library(AER)
library(randomForest)
library(e1071)
library(ggplot2)
library(reshape2)
library(dplyr)
######## Functions----------------------------------------
bias <- function(x,y) mean(x-y)
PEHE <- function(x,y) sqrt(mean((x-y)^2))
MSE <- function(x, y) mean((x - y)^2)


MC_se <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)

# Correlated covariates
get_features <- function(N, P) {
  
  # Generate correlated uniforms from a Gaussian Copula
  mysigma = matrix(1, P, P)
  
  for (i in 1:P) {
    for (j in 1:P) {
      mysigma[i, j] = 0.5^abs(i - j) + ifelse(i == j, 0, 0.1)
    }
  }
  
  mycop = MASS::mvrnorm(N, rep(0, P), Sigma = mysigma)
  unif = pnorm(mycop)
  
  
  # Transform in continuous and binary covariates-----------------------------
  x = matrix(NA, N, P)
  x[, 1:400] = qnorm(unif[, 1:400])
  x[, 401:P] = qbinom(unif[, 401:P], 1, 0.3)
  
  return(x)
}


### INFO---------------------------------------------------------------
N =500
P = 1000  # We tried P = 1000 to see how the problem degenerate
B = 500
 
# Store Metrics
Bias_XLasso = c(NA); Bias_ULasso = c(NA);  Bias_RLASSO = c(NA);  Bias_SLasso = c(NA); Bias_TLasso = c(NA); Bias_CF = c(NA)
PEHE_XLasso = c(NA); PEHE_ULasso = c(NA);  PEHE_RLASSO = c(NA); PEHE_SLasso = c(NA); PEHE_TLasso = c(NA);  PEHE_CF = c(NA)
PEHE_SBART = c(NA); Bias_SBART = c(NA); 
PEHE_RLASSO_SVM = c(NA); Bias_RLASSO_SVM = c(NA) 

for (b in 1:B) {
  
  
  cat("\n\n\n\n*** Iteration ", b, "\n\n\n")
  
  set.seed(b*50)
  # Simulate correlated covariates ---------------------------------
  x <- get_features(N, P)
  
  # Generate quantities, with strong targeted selection
  mu_base = 2 + 0.5*sin(pi*x[, 1]) - 0.25*x[, 2]^2 + 0.75*x[, 3]*x[, 9]  
  mu = 5*mu_base
  ITE = 1 + 2*abs(x[, 4]) + 1*x[, 10]
  
  Pscore = 0.9*plogis(1.2 - mu_base)
  w <- rbinom(N, 1, Pscore)
  y = mu + ITE*w + rnorm(N, 0, sd(ITE)/2)
  
  y_train = y
  z_train = w
  
  train_augmX = cbind(x)
  
  Train_ITE = ITE
  
  ###### MODELS ESTIMATION  ------------------------------------------------
  
  # For thresholding propensity scores
  trim = function(x, min, max) {
    x[x>max] = max
    x[x<min] = min
    return(x)
  }
  
  sanitize_x = function(x){
    # make sure x is a numeric matrix with named columns (for caret)
    if (!is.matrix(x) | !is.numeric(x) | any(is.na(x))) {
      stop("x must be a numeric matrix with no missing values")
    }
    colnames(x) = stringr::str_c("covariate_", 1:ncol(x))
    return(x)
  }
  
  sanitize_input = function(x,w,y) {
    x = sanitize_x(x)
    
    if (!is.numeric(w)) {
      stop("the input w should be a numeric vector")
    }
    if (is.numeric(w) & all(w %in% c(0,1))) {
      w = w==1
    }
    
    # make sure y is a numeric vector
    if (!is.numeric(y)) {
      stop("y should be a numeric vector")
    }
    
    # make sure the dimensions align
    if (length(y)!=nrow(x) | length(w)!=nrow(x)) {
      stop("nrow(x), length(w), and length(y) should all be equal")
    }
    
    return(list(x=x,
                w=w,
                y=y))
  }
  ########### RLASSO-SVM------------------------------------
  rlasso_svm <- function(x, w, y,
                         alpha = 1,
                         k_folds = NULL,
                         foldid = NULL,
                         lambda_tau = NULL,
                         lambda_choice = c("lambda.min","lambda.1se"),
                         rs = FALSE,
                         p_hat = NULL,
                         m_hat = NULL,
                         penalty_factor = NULL,
                         kernel = "radial",        
                         type = "nu-regression",   
                         cost = 1,                 
                         gamma = 0.1,              
                         epsilon = 0.1)  {
    
    # Sanitize input
    input <- sanitize_input(x, w, y)
    x <- input$x; w <- input$w; y <- input$y
    
    # Standardize X
    standardization <- caret::preProcess(x, method = c("center", "scale"))
    x_scl <- predict(standardization, x)
    x_scl <- x_scl[, !is.na(colSums(x_scl)), drop = FALSE]
    
    lambda_choice <- match.arg(lambda_choice)
    nobs <- nrow(x_scl)
    pobs <- ncol(x_scl)
    
    # Cross-validation folds
    if (is.null(foldid) || length(foldid) != length(w)) {
      if (is.null(k_folds)) {
        k_folds <- floor(max(3, min(10, length(w)/4)))
      }
      foldid <- sample(rep(seq(k_folds), length = length(w)))
    }
    
    # Penalty factors
    if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
      penalty_factor_nuisance <- rep(1, pobs)
      if (rs) {
        penalty_factor_tau <- c(0, rep(1, 2 * pobs))
      } else {
        penalty_factor_tau <- c(0, rep(1, pobs))
      }
    } else {
      penalty_factor_nuisance <- penalty_factor
      if (rs) {
        penalty_factor_tau <- c(0, penalty_factor, penalty_factor)
      } else {
        penalty_factor_tau <- c(0, penalty_factor)
      }
    }
    
    # --- Outcome model with SVM ---
    if (is.null(m_hat)) {
      y_fit <- e1071::svm(y ~ ., data = cbind(y, x),
                          type = type,
                          kernel = kernel,
                          cost = cost,
                          gamma = gamma,
                          epsilon = epsilon)
      m_hat <- predict(y_fit, x)
    } else { y_fit <- NULL }
    
    # --- Propensity model with SVM ---
    if (is.null(p_hat)) {
      w_fit <- e1071::svm(w ~ ., data = cbind(w, x),
                          type = type,
                          kernel = kernel,
                          cost = cost,
                          gamma = gamma)
      p_hat <- predict(w_fit, x)
    } else { w_fit <- NULL }
    
    # Residualization
    y_tilde <- y - m_hat
    
    if (rs) {
      x_scl_tilde <- cbind(as.numeric(w - p_hat) * cbind(1, x_scl), x_scl)
      x_scl_pred  <- cbind(1, x_scl, x_scl * 0)
    } else {
      x_scl_tilde <- cbind(as.numeric(w - p_hat) * cbind(1, x_scl))
      x_scl_pred  <- cbind(1, x_scl)
    }
    
    # Tau model (treatment effect)
    tau_fit <- glmnet::cv.glmnet(x_scl_tilde,
                                 y_tilde,
                                 foldid = foldid,
                                 alpha = alpha,
                                 lambda = lambda_tau,
                                 penalty.factor = penalty_factor_tau,
                                 standardize = FALSE)
    
    tau_beta <- as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    tau_hat  <- x_scl_pred %*% tau_beta
    
    ret <- list(tau_fit = tau_fit,
                tau_beta = tau_beta,
                w_fit = w_fit,
                y_fit = y_fit,
                p_hat = p_hat,
                m_hat = m_hat,
                tau_hat = tau_hat,
                rs = rs,
                standardization = standardization)
    
    class(ret) <- "rlasso_svm"
    return(ret)
  }
  
  
  # Prediction function
  predict.rlasso_svm <- function(object, newx = NULL, ...) {
    if (!is.null(newx)) {
      newx <- sanitize_x(newx)
      newx_scl <- predict(object$standardization, newx)
      newx_scl <- newx_scl[, !is.na(colSums(newx_scl)), drop = FALSE]
      
      if (object$rs) {
        newx_scl_pred <- cbind(1, newx_scl, newx_scl * 0)
      } else {
        newx_scl_pred <- cbind(1, newx_scl)
      }
      tau_hat <- newx_scl_pred %*% object$tau_beta
    } else {
      tau_hat <- object$tau_hat
    }
    return(tau_hat)
  }
  
  
  ######################### R-LASSO-SVM
  Rlasso_svm <- rlasso_svm(x = x, w = w, y = y, type="nu-regression", kernel="polynomial", cost= 500, gamma = 0.1, epsilon = 0.1)
  rlasso_svm_pred <- predict(Rlasso_svm, newx = x)
  ##### R-LASSO---------------------------------------------------
  
  Rlasso <- rlasso(x = x, w = w, y = y)
  rlasso_pred <- predict(Rlasso, x = x)
  
  #######S-Lasso--------------------------------------------------
  SLasso <- slasso(x = x, w = w, y = y)
  SLasso_pred <- predict(SLasso, x = x)
  
  ######################### T LASSO------------------------------
  TLasso <- tlasso(x = x, w = w, y = y)
  TLasso_pred <- predict(TLasso, x = x)
  
  ######################### X LASSO------------------------------
  XLasso <- xlasso(x = x, w = w, y = y)
  XLasso_pred <- predict(XLasso, x)
  
  ######################### U LASSO------------------------------
  ULasso <- ulasso(x = x, w = w, y = y)
  ULasso_pred <- predict(ULasso, x)
  
  ######################### Causal forest------------------------
  CF <- causal_forest(X = x, Y = y, W = w)
  
  #########S-BART--------------------------------------------------------
  myBART <- wbart(x.train = cbind(train_augmX, z_train), y.train = y_train, 
                  nskip = 1000, ndpost = 2000, printevery = 3000)
  XZ0_train <- cbind(train_augmX, z_train)
  XZ0_train[,"z_train"] <- ifelse(XZ0_train[,"z_train"] == 1, 0, 1)
  Y0_train_SBART <- predict(myBART, newdata = XZ0_train)
  All_obs_SBART <- cbind(Y1_SBART = myBART$yhat.train.mean, train_augmX, z_train)
  All_count_SBART <- cbind(Y0_SBART = colMeans(Y0_train_SBART), XZ0_train)
  All_Trt_SBART <- All_obs_SBART
  All_Trt_SBART[which(All_Trt_SBART[,"z_train"] == 0),] <- All_count_SBART[which(All_count_SBART[,"z_train"] == 1),]
  All_Ctrl_SBART <- All_count_SBART
  All_Ctrl_SBART[which(All_Ctrl_SBART[,"z_train"] == 1),] <- All_obs_SBART[which(All_obs_SBART[,"z_train"] == 0),]
  
  
  
  
  ###### Compute metrics-----------------------------------------------------
  
  Tau_RLASSO_SVM <- rlasso_svm_pred
  Tau_RLASSO <- rlasso_pred
  Tau_SLasso <- SLasso_pred
  Tau_TLasso <- TLasso_pred
  Tau_XLasso <- XLasso_pred
  Tau_ULasso <- ULasso_pred
  Tau_CF <- CF$predictions
  
  mean(Tau_RLASSO)
  mean(Tau_RLASSO_SVM)
  mean(Tau_SLasso)
  mean(Tau_TLasso)
  mean(Tau_XLasso)
  mean(Tau_ULasso)
  
  
  sd(Tau_RLASSO)
  sd(Tau_RLASSO_SVM)
  sd(Tau_SLasso)
  sd(Tau_TLasso)
  sd(Tau_XLasso)
  sd(Tau_ULasso)
  sd(Tau_CF)
  
  ######## Store ITE estimates-------------------------------------------
  Bias_RLASSO[b] <- bias(Tau_RLASSO, ITE)
  Bias_RLASSO_SVM[b] <- bias(Tau_RLASSO_SVM, ITE)
  Bias_SLasso[b] <- bias(Tau_SLasso, ITE)
  Bias_TLasso[b] <- bias(Tau_TLasso, ITE)
  Bias_XLasso[b] <- bias(Tau_XLasso, ITE)
  Bias_ULasso[b] <- bias(Tau_ULasso, ITE)
  Bias_CF[b] <- bias(Tau_CF, ITE)
  Bias_SBART[b] <- bias(ITE, All_Trt_SBART[,"Y1_SBART"] - All_Ctrl_SBART[,"Y0_SBART"])
  
  PEHE_RLASSO[b] <- PEHE(Tau_RLASSO, ITE)
  PEHE_RLASSO_SVM[b] <- PEHE(Tau_RLASSO_SVM, ITE)
  PEHE_SLasso[b] <- PEHE(Tau_SLasso, ITE)
  PEHE_TLasso[b] <- PEHE(Tau_TLasso, ITE)
  PEHE_XLasso[b] <- PEHE(Tau_XLasso, ITE)
  PEHE_ULasso[b] <- PEHE(Tau_ULasso, ITE)
  PEHE_CF[b] <- PEHE(Tau_CF, ITE)
  PEHE_SBART[b] <- PEHE(ITE, All_Trt_SBART[,"Y1_SBART"] - All_Ctrl_SBART[,"Y0_SBART"])
  
}

######## Final Metrics-----------------------------------------
Bias_Final <- data.frame(
  
  RLASSO = c(mean(Bias_RLASSO), MC_se(Bias_RLASSO, B)),
  RLASSO_SVM = c(mean(Bias_RLASSO_SVM), MC_se(Bias_RLASSO_SVM, B)),
  SLASSO = c(mean(Bias_SLasso), MC_se(Bias_SLasso, B)),
  TLASSO = c(mean(Bias_TLasso), MC_se(Bias_TLasso, B)),
  XLASSO = c(mean(Bias_XLasso), MC_se(Bias_XLasso, B)),
  ULASSO = c(mean(Bias_ULasso), MC_se(Bias_ULasso, B)),
  CF = c(mean(Bias_CF), MC_se(Bias_CF, B)),
  SBART = c(mean(Bias_SBART, na.rm = TRUE), MC_se(Bias_SBART, B))
  
)
rownames(Bias_Final) <- c("Bias", "SE")

PEHE_Final <- data.frame(
  
  RLASSO = c(mean(PEHE_RLASSO), MC_se(PEHE_RLASSO, B)),
  RLASSO_SVM = c(mean(PEHE_RLASSO_SVM), MC_se(PEHE_RLASSO_SVM, B)),
  SLASSO = c(mean(PEHE_SLasso), MC_se(PEHE_SLasso, B)),
  TLASSO = c(mean(PEHE_TLasso), MC_se(PEHE_TLasso, B)),
  XLASSO = c(mean(PEHE_XLasso), MC_se(PEHE_XLasso, B)),
  ULASSO = c(mean(PEHE_ULasso), MC_se(PEHE_ULasso, B)),
  CF = c(mean(PEHE_CF), MC_se(PEHE_CF, B)),
  SBART = c(mean(PEHE_SBART, na.rm = TRUE), MC_se(PEHE_SBART, B))
)
rownames(PEHE_Final) <- c("PEHE", "SE")
######## All iterations results------------------------------------
PEHE_Single = data.frame(
  
  RLASSO = PEHE_RLASSO,
  RLASSO_SVM = PEHE_RLASSO_SVM,
  SLASSO = PEHE_SLasso,
  TLASSO = PEHE_TLasso, 
  XLASSO = PEHE_XLasso,
  ULASSO = PEHE_ULasso,
  CRF = PEHE_CF
)

Bias_Single = data.frame(
  RLASSO = Bias_RLASSO,
  RLASSO_SVM = Bias_RLASSO_SVM,
  SLASSO = Bias_SLasso, 
  TLASSO = Bias_TLasso,
  XLASSO = Bias_XLasso, 
  ULASSO = Bias_ULasso,
  CRF = Bias_CF
)

# Print Results
Bias_Final
PEHE_Final

###### Box plot, line chart of Bias and PEHE------------------------------------

###### Step 1: Combine all vectors into data frames
PEHE_Single <- data.frame(
  
  RLASSO = PEHE_RLASSO,
  RLASSO_SVM = PEHE_RLASSO_SVM,
  SLASSO = PEHE_SLasso,
  TLASSO = PEHE_TLasso, 
  XLASSO = PEHE_XLasso,
  ULASSO = PEHE_ULasso,
  CF = PEHE_CF, 
  SBART = PEHE_SBART
  
  
)

Bias_Single <- data.frame(
  
  RLASSO = Bias_RLASSO, 
  RLASSO_SVM = Bias_RLASSO_SVM,
  SLASSO = Bias_SLasso, 
  TLASSO = Bias_TLasso,
  XLASSO = Bias_XLasso, 
  ULASSO = Bias_ULasso,
  CF = Bias_CF,
  SBART = Bias_SBART
)

# Step 2: Add iteration column
PEHE_Single$Iteration <- 1:nrow(PEHE_Single)
Bias_Single$Iteration <- 1:nrow(Bias_Single)

# Step 3: Melt to long format
PEHE_long <- melt(PEHE_Single, id.vars = "Iteration", variable.name = "Model", value.name = "PEHE")
Bias_long <- melt(Bias_Single, id.vars = "Iteration", variable.name = "Model", value.name = "Bias")

# Step 4: Summary stats for line and bar plots
PEHE_summary <- PEHE_long %>% group_by(Model) %>% summarise(mean_PEHE = mean(PEHE))
Bias_summary <- Bias_long %>% group_by(Model) %>% summarise(mean_Bias = mean(Bias))

# === LINE CHARTS (mean across models) ===
ggplot(PEHE_summary, aes(x = reorder(Model, mean_PEHE), y = mean_PEHE, group = 1)) +
  geom_line() + geom_point() +
  labs(title = "Sample Size N = 500, Number of Covariates P = 1000", x = "Model", y = "Mean PEHE") +
  theme_minimal() + theme(plot.title = element_text(face = "italic", size = 10))+ coord_flip()

ggplot(Bias_summary, aes(x = reorder(Model, mean_Bias), y = mean_Bias, group = 1)) +
  geom_line() + geom_point() +
  labs(title = "Sample Size N = 500, Number of Covariates P = 1000", x = "Model", y = "Mean Bias") +
  theme_minimal() +theme(plot.title = element_text(face = "italic", size = 10))+ coord_flip()

# === BOXPLOTS ===
# Boxplot for PEHE
ggplot(PEHE_long, aes(x = Model, y = PEHE, fill = Model)) +
  geom_boxplot() +
  labs(title = "Sample Size N = 500, Number of Covariates P = 1000", x = "Model", y = "PEHE") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 10), axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot for Bias
ggplot(Bias_long, aes(x = Model, y = Bias, fill = Model)) +
  geom_boxplot() +
  labs(title = "Sample Size N = 500, Number of Covariates P = 1000", x = "Model", y = "Bias") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 10), axis.text.x = element_text(angle = 45, hjust = 1))





# 95% Confidence Interval
ci_lower_bcf <- mean(PEHE) - 1.96 *PEHE_se
ci_upper_bcf <- mean(PEHE) + 1.96 *PEHE_se

ci_lower_bcf <- -2.3712 - 1.96 *0.02
ci_upper_bcf <- -2.3712 + 1.96 *0.02

c(ci_lower_bcf, ci_upper_bcf)


