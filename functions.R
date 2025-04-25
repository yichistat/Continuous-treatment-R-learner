#This script showcases the performance of the continuous-treatment R-learner for the simple simulation in Section 2 of the paper: "Towards R-learner with Continuous Treatments"
#by Yichi Zhang, Dehan Kong, and Shu Yang
#Available at: https://arxiv.org/pdf/2208.00872

library(mvtnorm)
library(parallel)
library(Matching)
library(FNN)
library(ggplot2)
library(splines2)
library(base)
library(dplyr)
library(np)
library(magrittr)
library(expm)
library(SuperLearner)
library(orthopolynom)
library(MASS)
library(polynom)

#construct the B-spline function
Bspline_construct <- function(X , df, range, degree, knots) {
  M <- bSpline(
    x = X,
    knots = knots,
    degree = degree,
    Boundary.knots = range,
    intercept = T
  )
  return(M)
}


#determine knot locations for B-spline basis construction using empirical quantiles

build_knot_X <- function(X, df, degree) {
  num_knot <- df - degree - 1
  knots <- quantile(X, probs = seq(0, 1, 1 / (num_knot + 1)))
  knots <- knots[2:(length(knots) - 1)]
  return(knots)
}

#the row-wise Kronecker product of two matrices

rowwise_kron <- function(M1, M2) {
  M <- matrix(0, nrow <- nrow(M1), ncol <- ncol(M1) * ncol(M2))
  for (i in 1:nrow(M1)) {
    M[i, ] <- kronecker(M1[i, ], M2[i, ])
  }
  return(M)
}

#Basis matrix construction
#in our code V represents the covariates

Basis_M_Build <- function(V,
                          Tt,
                          df_V,
                          df_T,
                          Vrange_list,
                          Trange_list,
                          V_knot_list,
                          T_knot_list) {
  n <- length(Tt)
  dim_V <- ncol(V)
  dim_all <- dim_V + 1
  basis_V <- lapply(1:ncol(V), function(i) {
    Bspline_construct(
      X = V[, i],
      df = df_V[i],
      range = Vrange_list[[i]],
      degree = degree,
      knots = V_knot_list[[i]]
    )
  })
  basis_T_matrix <- Bspline_construct(X = Tt[, 1], 
                                      df = df_T[1], 
                                      range = Trange_list[[1]], 
                                      degree = degree,
                                      knots = T_knot_list[[1]])
  basis_M <- matrix(1, nrow = n, ncol = 1)
  for (i in 1:dim_V) {
    basis_M <- rowwise_kron(basis_M, basis_V[[i]])
  }
  basis_V_matrix <- basis_M
  basis_M <- rowwise_kron(basis_M, basis_T_matrix)
  return(
    list(
      basis_M = basis_M,
      basis_V = basis_V,
      basis_T = basis_T_matrix
    )
  )
}

#obtain the B-spline coefficients for our proposed estimator in Section 3
#Gamma_diff is a matrix recording \Psi(X_i,T_i) - \hat{\Gamma}(X_i) based on the propensity score estimation
#m_e is a vector recording m(X_i) based on the outcome model estimation
#df_V and df_T are the numbers of basis on the direction of covariate and treatment respectively
#Vrange_list and Xrange_list are lists containing the ranges of all covariates and the treatment
#V_knot_list and X_knot_list are the knots for all directions of the b-spline functions

solve_coefficient = function(X,
                 Tt,
                 Y,
                 df_V,
                 df_T,
                 Gamma_diff,
                 m_e,
                 Vrange_list,
                 Trange_list,
                 V_knot_list,
                 T_knot_list,
                 rho) {
  #pre-setting
  V <- X %>% as.matrix()
  dim_V <- ncol(V)
  n <- length(Y)
  dim_all = dim_V + 1
  Kn = prod(df_V) * df_T
  #build basis
  Basis_M_Build_res <-
    Basis_M_Build(V %>% as.matrix(),
                  Tt,
                  df_V,
                  df_T,
                  Vrange_list,
                  Trange_list,
                  V_knot_list,
                  T_knot_list)
  basis_M <- Basis_M_Build_res$basis_M
  basis_T_matrix <- Basis_M_Build_res$basis_T_matrix
  basis_V <- Basis_M_Build_res$basis_V
  basis_V_matrix <- Basis_M_Build_res$basis_V_matrix
  #solve
  G <- t(Gamma_diff) %*% Gamma_diff / n + rho * t(basis_M) %*% basis_M / n
  inv_G <- base::solve(G)
  MpartII <- Gamma_diff %>% t
  phi <- inv_G %*% ((MpartII %*% (Y - m_e)) / n)
  return(list(phi = phi))
}

#based on the given covariates X (V_pred) and treatment T (Tt_pred)
#and the estimated B-spline coefficients (para)
#estimate CATE by our proposed estimator in Section 3

est_tau <- function(V_pred,
                    Tt_pred,
                    para,
                    df_V,
                    df_T,
                    Vrange_list,
                    Trange_list,
                    V_knot_list,
                    T_knot_list,
                    rho) {
  basis_M_pred <-
    Basis_M_Build(
      V = V_pred,
      Tt = Tt_pred,
      df_V,
      df_T,
      Vrange_list,
      Trange_list,
      V_knot_list,
      T_knot_list
    )$basis_M
  
  basis_M_pred_0 <-
    Basis_M_Build(
      V = V_pred,
      Tt = rep(0, length(Tt_pred)) %>% as.matrix,
      df_V,
      df_T,
      Vrange_list,
      Trange_list,
      V_knot_list,
      T_knot_list
    )$basis_M
  tau = (1+rho)*(basis_M_pred %*% para - basis_M_pred_0 %*% para)
  return(list(tau = tau))
}
