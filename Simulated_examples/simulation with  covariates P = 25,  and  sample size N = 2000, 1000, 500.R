###########################################################################
# simulation with  covariates P = 25 and  sample size N = 2000, 1000, 500 #
###########################################################################

rm(list = ls())


#### Libraries------------------------------------------------------------------

library(tidyverse)
library(BART)
library(rlearner)
library(future)
library(bartCause)
library(stan4bart)
library(MASS)
library(AER)
library(e1071)
library(dplyr)
library(ggplot2)
library(reshape2)        

####### Functions---------------------------------------------------------------
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
  
  
  # Transform in continuous and binary covariates
  x = matrix(NA, N, P)
  x[, 1:10] = qnorm(unif[, 1:10])
  x[, 11:P] = qbinom(unif[, 11:P], 1, 0.3)
  
  return(x)
}


### INFO------------------------------------------------------------------------
N =2000  # We tried N = 500, N =1000, and N = 2000
P = 25  # We tried P = 25 to see how the problem degenerate
B = 1000

##### Store Metrics-------------------------------------------------------------
Bias_XLasso = c(NA); Bias_ULasso = c(NA);  Bias_RLASSO = c(NA);  Bias_SLasso = c(NA); Bias_TLasso = c(NA);   Bias_KNN = c(NA)
PEHE_XLasso = c(NA); PEHE_ULasso = c(NA);  PEHE_RLASSO = c(NA) ; PEHE_SLasso = c(NA); PEHE_TLasso = c(NA);   PEHE_KNN = c(NA)
MSE_XLasso = c(NA); MSE_ULasso = c(NA); 
PEHE_SOLS = c(NA); Bias_SOLS = c(NA); MSE_SOLS = c(NA)
PEHE_TOLS = c(NA); Bias_TOLS = c(NA); MSE_TOLS = c(NA)
PEHE_RLASSO_SVM <- Bias_RLASSO_SVM <- MSE_RLASSO_SVM <- rep(NA, B)
MSE_RLASSO = c(NA);  MSE_SLasso = c(NA); MSE_TLasso = c(NA); MSE_KNN = c(NA); MSE_SOLS = c(NA); MSE_TOLS = c(NA);

for (b in 1:B) {
  
  
  cat("\n\n\n\n*** Iteration ", b, "\n\n\n")
  
  set.seed(b*50)
  
  # Simulate correlated covariates ---------------------------------
  x <- get_features(N, P)
  
  # Generate quantities, with strong targeted selection
  mu_base = 2 + 0.5*sin(pi*x[, 1]) - 0.25*x[, 2]^2 + 0.75*x[, 3]*x[, 9]  
  mu = 5*mu_base
  # plot(X[, 1],  mu)    # Uncomment to check quantities
  
  ITE = 1 + 2*abs(x[, 4]) + 1*x[, 10]
  
  Pscore = 0.9*plogis(1.2 - mu_base)
  w <- rbinom(N, 1, Pscore)
  
  y = mu + ITE*w + rnorm(N, 0, sd(ITE)/2)
  
  y_train = y
  z_train = w
  
  train_augmX = cbind(x)
  
  Train_ITE = ITE
  
  
  data1 = data.frame(w, x)
  ###### MODELS ESTIMATION  ------------------------------------------------
  # For thresholding propensity scores
  trim = function(x, min, max) {
    x[x>max] = max
    x[x<min] = min
    return(x)
  }
  ####### utils.R
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
  
  #######RLASSO-SVM function----------------------------------------------------
  
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
                         kernel = NULL,        
                         type = NULL,   
                         cost = NULL,                 
                         gamma = NULL,              
                         epsilon = NULL)  {
    
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
  Rlasso_svm <- rlasso_svm(x = x, w = w, y = y,  type="eps-regression", kernel="radial", cost= 1000, gamma = 0.1, epsilon = 0.1)
  rlasso_svm_pred <- predict(Rlasso_svm, newx = x)
  ##### R-LASSO
  
  Rlasso <- rlasso(x = x, w = w, y = y)
  rlasso_pred <- predict(Rlasso, x = x)
  
  #######S-Lasso
  SLasso <- slasso(x = x, w = w, y = y)
  SLasso_pred <- predict(SLasso, x = x)
  
  ######################### T LASSO
  TLasso <- tlasso(x = x, w = w, y = y)
  TLasso_pred <- predict(TLasso, x = x)
  
  ######################### X LASSO
  XLasso <- xlasso(x = x, w = w, y = y)
  XLasso_pred <- predict(XLasso, x)
  
  ######################### U LASSO
  ULasso <- ulasso(x = x, w = w, y = y)
  ULasso_pred <- predict(ULasso, x)
  
  
  ######################### S-Linear Regression
  SOLS <- lm(cbind.data.frame(y_train, z_train, train_augmX))
  Y0_train_SOLS <- predict(SOLS, newdata = cbind.data.frame(train_augmX, z_train = 0), se.fit = T)
  Y1_train_SOLS <- predict(SOLS, newdata = cbind.data.frame(train_augmX, z_train = 1), se.fit = T)
  PEHE_SOLS[b] <- PEHE(Train_ITE, Y1_train_SOLS$fit - Y0_train_SOLS$fit)
  
  ######################### T-Linear Regression
  TOLS1 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 1, ])
  TOLS0 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 0, ])
  Y0_train_TOLS <- predict(TOLS0, newdata = cbind.data.frame(train_augmX, z_train == 0), se.fit = T)
  Y1_train_TOLS <- predict(TOLS1, newdata = cbind.data.frame(train_augmX, z_train == 1), se.fit = T)
  PEHE_TOLS[b] <- PEHE(Train_ITE, Y1_train_TOLS$fit - Y0_train_TOLS$fit)
  
  # Compute metrics
  
  Tau_RLASSO_SVM <- rlasso_svm_pred
  Tau_RLASSO <- rlasso_pred
  Tau_SLasso <- SLasso_pred
  Tau_TLasso <- TLasso_pred
  Tau_XLasso <- XLasso_pred
  Tau_ULasso <- ULasso_pred
  
  
  
  
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
  
  
  
  Bias_RLASSO[b] <- bias(Tau_RLASSO, ITE)
  Bias_RLASSO_SVM[b] <- bias(Tau_RLASSO_SVM, ITE)
  Bias_SLasso[b] <- bias(Tau_SLasso, ITE)
  Bias_TLasso[b] <- bias(Tau_TLasso, ITE)
  Bias_XLasso[b] <- bias(Tau_XLasso, ITE)
  Bias_ULasso[b] <- bias(Tau_ULasso, ITE)
  Bias_SOLS[b] <- bias(Y1_train_SOLS$fit - Y0_train_SOLS$fit, ITE)
  Bias_TOLS[b] <- bias(Y1_train_TOLS$fit - Y0_train_TOLS$fit, ITE)
  
  
  MSE_RLASSO[b] <- MSE(Tau_RLASSO, ITE)
  MSE_RLASSO_SVM[b] <- MSE(Tau_RLASSO_SVM, ITE)
  MSE_SLasso[b] <- MSE(Tau_SLasso, ITE)
  MSE_TLasso[b] <- MSE(Tau_TLasso, ITE)
  MSE_XLasso[b] <- MSE(Tau_XLasso, ITE)
  MSE_ULasso[b] <- MSE(Tau_ULasso, ITE)
  MSE_SOLS[b] <- MSE(Y1_train_SOLS$fit - Y0_train_SOLS$fit, ITE)
  MSE_TOLS[b] <- MSE(Y1_train_TOLS$fit - Y0_train_TOLS$fit, ITE)
  
  
  PEHE_RLASSO[b] <- PEHE(Tau_RLASSO, ITE)
  PEHE_RLASSO_SVM[b] <- PEHE(Tau_RLASSO_SVM, ITE)
  PEHE_SLasso[b] <- PEHE(Tau_SLasso, ITE)
  PEHE_TLasso[b] <- PEHE(Tau_TLasso, ITE)
  PEHE_XLasso[b] <- PEHE(Tau_XLasso, ITE)
  PEHE_ULasso[b] <- PEHE(Tau_ULasso, ITE)
  PEHE_SOLS[b] <- PEHE(Y1_train_SOLS$fit - Y0_train_SOLS$fit, ITE)
  PEHE_TOLS[b] <- PEHE(Y1_train_TOLS$fit - Y0_train_TOLS$fit, ITE)
}

# Final Metrics
Bias_Final <- data.frame(
  
  RLASSO = c(mean(Bias_RLASSO), MC_se(Bias_RLASSO, B)),
  RLASSO_SVM = c(mean(Bias_RLASSO_SVM), MC_se(Bias_RLASSO_SVM, B)),
  SLASSO = c(mean(Bias_SLasso), MC_se(Bias_SLasso, B)),
  TLASSO = c(mean(Bias_TLasso), MC_se(Bias_TLasso, B)),
  XLASSO = c(mean(Bias_XLasso), MC_se(Bias_XLasso, B)),
  ULASSO = c(mean(Bias_ULasso), MC_se(Bias_ULasso, B)),
  SOLS = c(mean(Bias_SOLS), MC_se(Bias_SOLS, B)),
  TOLS = c(mean(Bias_TOLS), MC_se(Bias_TOLS, B))
)
rownames(Bias_Final) <- c("Bias", "SE")

MSE_Final <- data.frame(
  
  RLASSO = c(mean(MSE_RLASSO), MC_se(MSE_RLASSO, B)),
  RLASSO_SVM = c(mean(MSE_RLASSO_SVM), MC_se(MSE_RLASSO_SVM, B)),
  SLASSO = c(mean(MSE_SLasso), MC_se(MSE_SLasso, B)),
  TLASSO = c(mean(MSE_TLasso), MC_se(MSE_TLasso, B)),
  XLASSO = c(mean(MSE_XLasso), MC_se(MSE_XLasso, B)),
  ULASSO = c(mean(MSE_ULasso), MC_se(MSE_ULasso, B)),
  SOLS = c(mean(MSE_SOLS), MC_se(MSE_SOLS, B)),
  TOLS = c(mean(MSE_TOLS), MC_se(MSE_TOLS, B))
)
rownames(MSE_Final) <- c("MSE", "SE")

PEHE_Final <- data.frame(
  
  RLASSO = c(mean(PEHE_RLASSO), MC_se(PEHE_RLASSO, B)),
  RLASSO_SVM = c(mean(PEHE_RLASSO_SVM), MC_se(PEHE_RLASSO_SVM, B)),
  SLASSO = c(mean(PEHE_SLasso), MC_se(PEHE_SLasso, B)),
  TLASSO = c(mean(PEHE_TLasso), MC_se(PEHE_TLasso, B)),
  XLASSO = c(mean(PEHE_XLasso), MC_se(PEHE_XLasso, B)),
  ULASSO = c(mean(PEHE_ULasso), MC_se(PEHE_ULasso, B)),
  SOLS = c(mean(PEHE_SOLS), MC_se(PEHE_SOLS, B)),
  TOLS = c(mean(PEHE_TOLS), MC_se(PEHE_TOLS, B))
)
rownames(PEHE_Final) <- c("PEHE", "SE")
# All iterations results
PEHE_Single = data.frame(
  
  RLASSO = PEHE_RLASSO,
  RLASSO_SVM = PEHE_RLASSO_SVM,
  SLASSO = PEHE_SLasso,
  TLASSO = PEHE_TLasso, 
  XLASSO = PEHE_XLasso,
  ULASSO = PEHE_ULasso,
  SOLS = PEHE_SOLS, 
  TOLS = PEHE_TOLS
)

Bias_Single = data.frame(
  RLASSO = Bias_RLASSO,
  RLASSO_SVM = Bias_RLASSO_SVM,
  SLASSO = Bias_SLasso, 
  TLASSO = Bias_TLasso,
  XLASSO = Bias_XLasso, 
  ULASSO = Bias_ULasso,
  SOLS = Bias_SOLS, 
  TOLS = Bias_TOLS
)


MSE_Single = data.frame(
  
  RLASSO = MSE_RLASSO, 
  RLASSO_SVM = MSE_RLASSO_SVM,
  SLASSO = MSE_SLasso, 
  TLASSO = MSE_TLasso,
  XLASSO = MSE_XLasso, 
  ULASSO = MSE_ULasso,
  SOLS = MSE_SOLS, 
  TOLS = MSE_TOLS
)


# Print Results
Bias_Final
MSE_Final
PEHE_Final



###### Box plot, line chart of Bias, PEHE, and MSE------------------------------

# Step 1: Combine all vectors into data frames
PEHE_Single <- data.frame(
  
  RLASSO = PEHE_RLASSO,
  RLASSO_SVM = PEHE_RLASSO_SVM,
  SLASSO = PEHE_SLasso,
  TLASSO = PEHE_TLasso, 
  XLASSO = PEHE_XLasso,
  ULASSO = PEHE_ULasso,
  SOLS = PEHE_SOLS, 
  TOLS = PEHE_TOLS
)

Bias_Single <- data.frame(
  
  RLASSO = Bias_RLASSO, 
  RLASSO_SVM = Bias_RLASSO_SVM,
  SLASSO = Bias_SLasso, 
  TLASSO = Bias_TLasso,
  XLASSO = Bias_XLasso, 
  ULASSO = Bias_ULasso,
  SOLS = Bias_SOLS, 
  TOLS = Bias_TOLS
)

MSE_Single <- data.frame(
  RLASSO = MSE_RLASSO, 
  RLASSO_SVM = MSE_RLASSO_SVM,
  SLASSO = MSE_SLasso, 
  TLASSO = MSE_TLasso,
  XLASSO = MSE_XLasso, 
  ULASSO = MSE_ULasso,
  SOLS = MSE_SOLS, 
  TOLS = MSE_TOLS
)

# Step 2: Add iteration column
PEHE_Single$Iteration <- 1:nrow(PEHE_Single)
Bias_Single$Iteration <- 1:nrow(Bias_Single)
MSE_Single$Iteration <- 1:nrow(MSE_Single)

# Step 3: Melt to long format
PEHE_long <- melt(PEHE_Single, id.vars = "Iteration", variable.name = "Model", value.name = "PEHE")
Bias_long <- melt(Bias_Single, id.vars = "Iteration", variable.name = "Model", value.name = "Bias")
MSE_long <- melt(MSE_Single, id.vars = "Iteration", variable.name = "Model", value.name = "MSE")

# Step 4: Summary stats for line and bar plots
PEHE_summary <- PEHE_long %>% group_by(Model) %>% summarise(mean_PEHE = mean(PEHE))
Bias_summary <- Bias_long %>% group_by(Model) %>% summarise(mean_Bias = mean(Bias))
MSE_summary <- MSE_long %>% group_by(Model) %>% summarise(mean_MSE = mean(MSE))



# === LINE CHARTS (mean across models) ===
ggplot(PEHE_summary, aes(x = reorder(Model, mean_PEHE), y = mean_PEHE, group = 1)) +
  geom_line() + geom_point() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "Mean PEHE") +
  theme_minimal() + theme(plot.title = element_text(face = "italic", size = 9))+ coord_flip()

ggplot(Bias_summary, aes(x = reorder(Model, mean_Bias), y = mean_Bias, group = 1)) +
  geom_line() + geom_point() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "Mean Bias") +
  theme_minimal() +theme(plot.title = element_text(face = "italic", size = 9))+ coord_flip()

ggplot(MSE_summary, aes(x = reorder(Model, mean_MSE), y = mean_MSE, group = 1)) +
  geom_line() + geom_point() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "Mean MSE") +
  theme_minimal() + theme(plot.title = element_text(face = "italic", size = 9))+ coord_flip()

# === BOXPLOTS ===
# Boxplot for PEHE
ggplot(PEHE_long, aes(x = Model, y = PEHE, fill = Model)) +
  geom_boxplot() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "PEHE") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 9), axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot for Bias
ggplot(Bias_long, aes(x = Model, y = Bias, fill = Model)) +
  geom_boxplot() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "Bias") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 9), axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot for MSE
ggplot(MSE_long, aes(x = Model, y = MSE, fill = Model)) +
  geom_boxplot() +
  labs(title = "Sample Size N = 2000, Number of Covariates P = 10", x = "Model", y = "MSE") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 9), axis.text.x = element_text(angle = 45, hjust = 1))


