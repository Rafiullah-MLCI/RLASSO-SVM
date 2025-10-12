#' Residualized LASSO with Support Vector Machines (RLASSO-SVM)
#'
#' Estimate heterogeneous treatment effects by combining residualized LASSO with 
#' Support Vector Machines (SVM) for nuisance estimation.
#'
#' @param x Covariate matrix (n x p).
#' @param w Treatment vector (binary: 0/1).
#' @param y Outcome vector (numeric).
#' @param alpha Elastic net mixing parameter (default = 1 for LASSO).
#' @param k_folds Number of CV folds (default chosen automatically).
#' @param foldid Optional fold assignments for CV.
#' @param lambda_tau Optional lambda sequence for tau estimation.
#' @param lambda_choice Which lambda to use ("lambda.min" or "lambda.1se").
#' @param rs Residualization flag (default = FALSE).
#' @param p_hat Optional estimated propensity scores.
#' @param m_hat Optional estimated outcome regression.
#' @param penalty_factor Penalty factors for LASSO (default = equal weights).
#' @param kernel SVM kernel type ("linear", "polynomial", "radial")Default = "radial".
#' @param type SVM type (default = "nu-regression").
#' @param cost Cost parameter C for SVM (default = 1).
#' @param gamma Kernel parameter (default = 0.1).
#' @param epsilon Epsilon parameter for regression margin (default = 0.1).
#'
#' @return A fitted `rlasso_svm` object containing:
#' \itemize{
#'   \item tau_fit: glmnet object for treatment effect estimation
#'   \item tau_beta: estimated coefficients
#'   \item w_fit, y_fit: fitted SVM models
#'   \item p_hat, m_hat: estimated nuisance functions
#'   \item tau_hat: estimated CATEs
#'   \item standardization: preprocessing object
#' }
#' @export
#' 
#'##### R Utiles.
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


# RLASSO-SVM function
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
                       kernel = "radial",       # defaults  kernel 
                       type = "nu-regression",# Default type suitable for regression
                       cost = 1,
                       gamma = 0.1,
                       epsilon = 0.1) {
  
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


#### Example data
set.seed(123)
n <- 500
p <- 20

x = matrix(rnorm(n*p), n, p) # Covariates
w = rbinom(n, 1, 0.5) # Treatment
y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n) # Outcome
# Fit the model
Rlasso_svm <- rlasso_svm(x = x, w = w, y = y)

# Predict treatment effects
rlasso_svm_pred <- predict(Rlasso_svm, newx = x)
print(rlasso_svm_pred)