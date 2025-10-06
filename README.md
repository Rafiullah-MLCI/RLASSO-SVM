# RLASSO-SVM
An R function to implement Residualized LASSO and Support Vector Machines (SVM) for treatment effect estimation.
## Details
**Function**: RLASSO-SVM  
**Type**: R Function  
**Title**: Residualized LASSO with Support Vector Machines for Heterogeneous Treatment Effects Estimation  
**Version**: 1.0  
**Date**: 2025-09-29  
**Author**: Rafullah
**Maintainer**: Rafullah <rafiuom111@gmail.com>  
**Description**: The RLASSO-SVM function implements a hybrid model combining Residualized Least Absolute Shrinkage and Selection Operator (LASSO) with Support Vector Machines (SVM) for estimating heterogeneous treatment effects estimation. It uses LASSO (via glmnet) for feature selection and SVM (via e1071) for estimating nuisance parameters (m_hat, p_hat), with preprocessing handled by caret and stringr. The function supports cross-validation, standardization, and an optional robust specification (rs) for high-dimensional datasets.  
**License**: GPL-3  
**Depends**: R (>= 3.5.0)  
**Suggested Packages**: caret, stringr, e1071, glmnet  
**Encoding**: UTF-8  

## Usage
```R
# Example data
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
