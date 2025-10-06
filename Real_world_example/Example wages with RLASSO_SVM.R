################################################################################
####### CARD data, the effect of education on wages, Example in the Paper ######
################################################################################

rm(list = ls())

# Libraries
library(tidyverse)
library(BART)
library(grf)
library(rlearner)
library(bartcs)
library(future)
library(bartCause)
library(stan4bart)
library(MASS)
library(AER)
library(e1071)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(reshape2)
library(dplyr)

####### Functions---------------------------------------------------------------
BScore <- function(x,y) mean((x-y)^2)
bias <- function(x, y) mean(x - y)
PEHE <- function(x, y) sqrt(mean((x - y)^2))
MC_se <- function(x, B) qt(0.975, B - 1) * sd(x) / sqrt(B)

###### Read the dataset---------------------------------------------------------

data<- read.csv(file.choose(),header=TRUE)

# Fill missing values 
data = na.omit(data)

######### Generate synthetic ITE------------------------------------------------
set.seed(123)
ITE <- 0.5 * data$educ +
  0.3 * sqrt(data$exper) +
  0.2 * log1p(data$fatheduc + data$motheduc) +
  0.1 * (data$iq / 10) -
  0.05 * data$age +
  rnorm(nrow(data), mean = 0, sd = 0.5)

######### Save the dataset with true ITE----------------------------------------
data = cbind(data, ITE)



######## Convert to binary based on a threshold (e.g., value > 15)
z <- data$educ
y = data$lwage

x = as.matrix(data[, -c(4, 8, 26, 33, 35)])
P = ncol(x)
### INFO
N <- 1612
P <- 30  # We tried P = 30 to see how the problem degenerates
B <- 20

######### Store Metrics---------------------------------------------------------
Bias_XLasso <- Bias_ULasso <- Bias_RLASSO<- Bias_SLasso <- Bias_TLasso <- Bias_CF <-  rep(NA, B)
PEHE_XLasso <- PEHE_ULasso <- PEHE_RLASSO<- PEHE_SLasso <- PEHE_TLasso <- PEHE_CF <- rep(NA, B)
PEHE_SOLS <- Bias_SOLS <- rep(NA, B)
PEHE_TOLS <- Bias_TOLS <- rep(NA, B)
PEHE_SBART <- Bias_SBART <- rep(NA, B)
PEHE_RLASSO_SVM <- Bias_RLASSO_SVM <- rep(NA, B)


Split_Mu <- matrix(NA, B, P + 1)
Split_Tau <- matrix(NA, B, P)

for (b in 1:B) {
  cat("\n\n\n\n*** Iteration ", b, "\n\n\n")
  set.seed(b * 50)
  
  w <- ifelse(z >= 14, 1, 0)
  y_train <- y
  z_train <- w
  train_augmX <- cbind(x)
  Train_ITE <- ITE
  data1 <- data.frame(w, x)
  PS_nn <- nnet(x = x, y = w, size = 20, maxit = 2000, 
                decay = 0.01, trace=FALSE, abstol = 1.0e-8) 
  PS_est = PS_nn$fitted.values
  
  ###### MODELS ESTIMATION  ----------------------------------------------------
  #### utils
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
  
  ############################## 
  ###### RLASSO-SVM Function ###------------------------------------------------------
  ############################## 
  
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
                         epsilon = NULL) {
    
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
  Rlasso_svm <- rlasso_svm(x = x, w = w, y = y, type="eps-regression", kernel="sigmoid", cost=1500, gamma=0.01, epsilon=0.1)
  rlasso_svm_pred <- predict(Rlasso_svm, newx = x)
  ######################## R-LASSO----------------------------------------------
  
  Rlasso <- rlasso(x = x, w = w, y = y)
  rlasso_pred <- predict(Rlasso, x = x)
  
  ######################## S-Lasso----------------------------------------------
  SLasso <- slasso(x = x, w = w, y = y)
  SLasso_pred <- predict(SLasso, x = x)
  
  ######################### T-LASSO---------------------------------------------
  TLasso <- tlasso(x = x, w = w, y = y)
  TLasso_pred <- predict(TLasso, x = x)
  
  ######################### X-LASSO---------------------------------------------
  XLasso <- xlasso(x = x, w = w, y = y)
  XLasso_pred <- predict(XLasso, x)
  
  ######################### U-LASSO
  ULasso <- ulasso(x = x, w = w, y = y)
  ULasso_pred <- predict(ULasso, x)
  
  ######################### Causal forest---------------------------------------
  CF <- causal_forest(X = x, Y = y, W = w)
  
  ######################### S-Linear Regression---------------------------------
  SOLS <- lm(cbind.data.frame(y_train, z_train, train_augmX))
  Y0_train_SOLS <- predict(SOLS, newdata = cbind.data.frame(train_augmX, z_train = 0), se.fit = T)
  Y1_train_SOLS <- predict(SOLS, newdata = cbind.data.frame(train_augmX, z_train = 1), se.fit = T)
  PEHE_SOLS[b] <- PEHE(Train_ITE, Y1_train_SOLS$fit - Y0_train_SOLS$fit)
  
  ######################### T-Linear Regression---------------------------------
  TOLS1 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 1, ])
  TOLS0 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 0, ])
  Y0_train_TOLS <- predict(TOLS0, newdata = cbind.data.frame(train_augmX, z_train == 0), se.fit = T)
  Y1_train_TOLS <- predict(TOLS1, newdata = cbind.data.frame(train_augmX, z_train == 1), se.fit = T)
  PEHE_TOLS[b] <- PEHE(Train_ITE, Y1_train_TOLS$fit - Y0_train_TOLS$fit)
  
  #########Compute metrics------------------------------------------------------
  
  Tau_RLASSO_SVM <- rlasso_svm_pred
  Tau_RLASSO <- rlasso_pred
  Tau_SLasso <- SLasso_pred
  Tau_TLasso <- TLasso_pred
  Tau_XLasso <- XLasso_pred
  Tau_ULasso <- ULasso_pred
  Tau_CF <- CF$predictions
  ##### Mean--------------------------------------------------------------------
  mean(Tau_RLASSO)
  mean(Tau_RLASSO_SVM)
  mean(Tau_SLasso)
  mean(Tau_TLasso)
  mean(Tau_XLasso)
  mean(Tau_ULasso)

  Bias_RLASSO[b] <- bias(Tau_RLASSO, ITE)
  Bias_RLASSO_SVM[b] <- bias(Tau_RLASSO_SVM, ITE)
  Bias_SLasso[b] <- bias(Tau_SLasso, ITE)
  Bias_TLasso[b] <- bias(Tau_TLasso, ITE)
  Bias_XLasso[b] <- bias(Tau_XLasso, ITE)
  Bias_ULasso[b] <- bias(Tau_ULasso, ITE)
  Bias_CF[b] <- bias(Tau_CF, ITE)
  Bias_SOLS[b] <- bias(Y1_train_SOLS$fit - Y0_train_SOLS$fit, ITE)
  Bias_TOLS[b] <- bias(Y1_train_TOLS$fit - Y0_train_TOLS$fit, ITE)
  
  
  PEHE_RLASSO[b] <- PEHE(Tau_RLASSO, ITE)
  PEHE_RLASSO_SVM[b] <- PEHE(Tau_RLASSO_SVM, ITE)
  PEHE_SLasso[b] <- PEHE(Tau_SLasso, ITE)
  PEHE_TLasso[b] <- PEHE(Tau_TLasso, ITE)
  PEHE_XLasso[b] <- PEHE(Tau_XLasso, ITE)
  PEHE_ULasso[b] <- PEHE(Tau_ULasso, ITE)
  PEHE_CF[b] <- PEHE(Tau_CF, ITE)
  PEHE_SOLS[b] <- PEHE(Y1_train_SOLS$fit - Y0_train_SOLS$fit, ITE)
  PEHE_TOLS[b] <- PEHE(Y1_train_TOLS$fit - Y0_train_TOLS$fit, ITE)
}

############# Final Metrics-----------------------------------------------------
Bias_Final <- data.frame(
  
  RLASSO = c(mean(Bias_RLASSO), MC_se(Bias_RLASSO, B)),
  RLASSO_SVM = c(mean(Bias_RLASSO_SVM), MC_se(Bias_RLASSO_SVM, B)),
  SLASSO = c(mean(Bias_SLasso), MC_se(Bias_SLasso, B)),
  TLASSO = c(mean(Bias_TLasso), MC_se(Bias_TLasso, B)),
  XLASSO = c(mean(Bias_XLasso), MC_se(Bias_XLasso, B)),
  ULASSO = c(mean(Bias_ULasso), MC_se(Bias_ULasso, B)),
  CF = c(mean(Bias_CF), MC_se(Bias_CF, B)),
  SOLS = c(mean(Bias_SOLS), MC_se(Bias_SOLS, B)),
  TOLS = c(mean(Bias_TOLS), MC_se(Bias_TOLS, B))
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
  SOLS = c(mean(PEHE_SOLS), MC_se(PEHE_SOLS, B)),
  TOLS = c(mean(PEHE_TOLS), MC_se(PEHE_TOLS, B))
)
rownames(PEHE_Final) <- c("PEHE", "SE")

######### All iterations results------------------------------------------------
PEHE_Single = data.frame(
  
  RLASSO = PEHE_RLASSO,
  RLASSO_SVM = PEHE_RLASSO_SVM,
  SLASSO = PEHE_SLasso,
  TLASSO = PEHE_TLasso, 
  XLASSO = PEHE_XLasso,
  ULASSO = PEHE_ULasso,
  CRF = PEHE_CF, 
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
  CRF = Bias_CF,
  SOLS = Bias_SOLS, 
  TOLS = Bias_TOLS
)


# Print Results
Bias_Final
PEHE_Final

####### Subgroup analysis-------------------------------------------------------
mytree <- rpart(
  Tau_RLASSO_SVM ~ ., 
  data =  cbind.data.frame(Tau_RLASSO_SVM, x), 
  control = list(maxdepth = 3)
)

# pdf("YOUR DIRECTORY", 
#     width = 12, height = 8)
rpart.plot(mytree, type = 2, extra = 101, clip.right.labs = FALSE, 
           box.palette = "Blues", # color scheme
           branch.lty = 1, # dotted branch lines
           shadow.col = "gray",
           branch.lwd = 1,
           tweak = 1.1,
           branch = 1, under = TRUE,  yesno = 2)

# Assuming you have a trained RLASSO_SVM model from your simulation
# Extract coefficients from one iteration (e.g., first iteration)
set.seed(42)  # For reproducibility
example_iteration <- 1  # Change this to see different iterations

# Get feature names (adapt to match your actual data)
feature_names <- colnames(x)  # From your x matrix
if(is.null(feature_names)) feature_names <- paste0("X", 1:ncol(x))

# Create coefficient data
coef_data <- data.frame(
  feature = factor(feature_names, levels = rev(feature_names)),
  coefficient = Rlasso_svm$tau_beta[-1]  # Exclude intercept
)

# Create SHAP-like plot
ggplot(coef_data, aes(x = coefficient, y = feature, fill = coefficient)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_gradient2(low = "#377eb8", high = "#e41a1c", mid = "white",
                       midpoint = 0, name = "Impact") +
  scale_x_continuous(breaks = seq(-0.2, 0.2, 0.1),
                     limits = c(-0.25, 0.25),
                     expand = c(0, 0)) +
  labs(x = "RLASSO-SVM Coefficient Impact", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "right",
    plot.margin = unit(c(10, 20, 10, 10), "pt")
  )


ggplot(data.frame(Tau = Tau_RLASSO_SVM), aes(x = Tau)) +
  geom_density(fill = "#377eb8", alpha = 0.6) +
  labs(x = "Treatment Effect (RLASSO-SVM)", y = "Density", 
       title = "Distribution of Heterogeneous Treatment Effects") +
  theme_minimal()

################ E-values-------------------------------------------------------
# Example: SE of outcome mean
SE_RLASSO_SVM <- sd(Tau_RLASSO_SVM) / sqrt(length(Tau_RLASSO_SVM))
SE_SLASSO <- sd(Tau_SLasso) / sqrt(length(Tau_SLasso))
SE_TLASSO <- sd(Tau_TLasso) / sqrt(length(Tau_TLasso))
SE_ULASSO <- sd(Tau_ULasso) / sqrt(length(Tau_ULasso))
SE_XLASSO <- sd(Tau_XLasso) / sqrt(length(Tau_XLasso))
SE_RLASSO <- sd(Tau_RLASSO) / sqrt(length(Tau_RLASSO))
SE_CF <- sd(Tau_CF) / sqrt(length(Tau_CF))

cat("SE of Mean Y:", round(SE_Y, 4), "\n")

# Simple base R solution without dplyr
calculate_e_value <- function(estimate) {
  abs_estimate <- abs(estimate)
  if (abs_estimate < 1) {
    abs_estimate <- 1 / abs_estimate
  }
  e_value <- abs_estimate + sqrt(abs_estimate * (abs_estimate - 1))
  return(e_value)
}

# Your data
data <- data.frame(
  method = c("RLASSO-SVM", "RLASSO", "SLASSO", "TLASSO", "ULASSO", "XLASSO", "CF"),
  ATE = c(0.4255, -0.5659, 0.0311, 0.0474, -0.1248, -0.0517, 0.0739),
  SE = c(0.008, 0.121, 0.002, 0.003, 0.001, 0.002, 0.002)
)

# Compute E-values robustly (handles RR < 1)
compute_evalues <- function(ATE, SE, sd_y = 1, p0 = 0.25) {
  # Standardized mean difference
  d <- ATE / sd_y
  
  # Convert SMD to log(OR)
  logOR <- d * (pi / sqrt(3))
  OR <- exp(logOR)
  
  # Convert OR to RR given baseline risk p0
  RR <- OR / (1 - p0 + p0 * OR)
  
  # Adjust formula for RR < 1 (protective effect)
  RR_for_eval <- ifelse(RR < 1, 1 / RR, RR)
  E_point <- RR_for_eval + sqrt(RR_for_eval * (RR_for_eval - 1))
  
  # CI bounds for ATE
  CI_lower <- ATE - 1.96 * SE
  CI_upper <- ATE + 1.96 * SE
  
  # Standardized mean differences for CI bounds
  d_low <- CI_lower / sd_y
  d_up  <- CI_upper / sd_y
  
  # Conservative CI bound (closer to 1)
  logOR_low <- d_low * (pi / sqrt(3))
  logOR_up  <- d_up  * (pi / sqrt(3))
  
  OR_low <- exp(logOR_low)
  OR_up  <- exp(logOR_up)
  
  RR_low <- OR_low / (1 - p0 + p0 * OR_low)
  RR_up  <- OR_up  / (1 - p0 + p0 * OR_up)
  
  # Pick CI closer to null (RR=1)
  RR_ci <- ifelse(abs(RR_low - 1) < abs(RR_up - 1), RR_low, RR_up)
  
  # Adjust for RR < 1
  RR_ci_for_eval <- ifelse(RR_ci < 1, 1 / RR_ci, RR_ci)
  E_ci <- RR_ci_for_eval + sqrt(RR_ci_for_eval * (RR_ci_for_eval - 1))
  
  return(list(E_point = E_point, E_ci = E_ci,
              CI_lower = CI_lower, CI_upper = CI_upper))
}

# Recompute results
results <- do.call(rbind, lapply(1:nrow(data), function(i) {
  res <- compute_evalues(data$ATE[i], data$SE[i], sd_y = 1, p0 = 0.25)
  data.frame(
    Method = data$method[i],
    ATE_95_CI = sprintf("%.3f (%.3f, %.3f)",
                        data$ATE[i], res$CI_lower, res$CI_upper),
    E_value_Point = round(res$E_point, 2),
    E_value_CI = round(res$E_ci, 2)
  )
}))

cat("E-VALUE ANALYSIS (Continuous Outcome, p0 = 0.25)\n")
print(results, row.names = FALSE)



