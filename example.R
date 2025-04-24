library(SuperLearner)
library(pracma)
library(MASS)
library(dplyr)
library(earth)
library(plotly)
source("functions.R")
numWorkers <- 16
ranseed <- (1:500) + 2222
sim <- length(ranseed)
mc = 1

#plot the true and the estimated CATE function
plot_fun = function(Ttt,V1,tau){
  plot_ly(x=V1,y=Ttt,z=tau, showscale=FALSE, type="surface",
          colors='RdBu',cmin=-1,cmax=2,
          contours = list(
            x = list(show = TRUE,color = "black", start = 0, end = 1, size = 0.05,width = 1),
            y = list(show = TRUE,color = "black", start = 0, end = 1, size = 0.05,width = 1)
          )
  ) %>% layout(
    font = list(family = "Raleway"),
    scene = list(
      xaxis = list(
        title = list(text = "X",font = list(size = 25)),
        nticks = 4,
        range=c(0,1),
        tickfont = list(size = 18)),
      yaxis = list(
        title = list(text = "T",font = list(size = 25)),
        nticks = 4,
        range=c(0,1),
        tickfont = list(size = 18)
      ),
      zaxis = list(
        title = "",
        nticks = 4,
        range=c(-0.8,2.5),
        tickfont = list(size = 18)
      ),
      camera=list(
        eye = list(x=0.7*1.5, y=-1.2*1.5, z=0.3*1.5)
      ),
      aspectratio = list(x = 0.9, y = 0.9, z = 1)
    )
  ) %>% config(mathjax = 'cdn')
}

#the true cate function in our example

tau <- function(X1, Tt) {
  r1 <- sqrt((X1 - 0.25)^2 + (Tt - 0.25)^2)
  r2 <- sqrt((X1 - 0.75)^2 + (Tt - 0.75)^2)
  if (r1 < 0.25) {
    res <- sin(2 * pi * r1 / (1 / 2) + pi / 2)  + 1
  } else{
    if (r2 < 0.25) {
      res <- sin(2 * pi * r2 / (1 / 2) + pi / 2)  + 1
    } else{
      res = 0
    }
  }
  return(res)
}

#plot the true CATE function

Ttt <- seq(0,1,0.01)
V1 <- seq(0,1,0.01)
grids <- Ttt %o% V1
truth <- matrix(0,nrow = length(Ttt), ncol = length(V1))
for (i in 1:length(Ttt)){
  for (j in 1:length(V1))
  {
    truth[i,j] <- tau(V1[j],Ttt[i])
  }
}
plot_fun(Ttt,V1,truth)

#data generation

set.seed(123456)
tau <- Vectorize(tau)
n <- 5000
X <- cbind(runif(n, min = 0, max = 1))
V <- X
Tt <- runif(n, min = 0, max = 1) %>% matrix
Y <-  sin(2 * X[, 1]) + X[, 1]^2 + tau(X[, 1], Tt) + rnorm(n, 0, 0.3)

#in this code, we use V to represent the covariates X
#the degree of the b-spline function is chosen to be 2
#the numbers of basis for both treatment and covariate are both 8

degree <- 2
df_T <- 8
df_V <- 8
#the range of B-spline function in each direction is (0,1)
Trange = c(0, 1) %>% list
Vrange = c(0, 1) %>% list
V_knot_list <- list()
T_knot_list <- list()
V_knot_list[[1]] <- build_knot_X(X = X[, 1], df = df_V, degree = degree)
T_knot_list[[1]] <- build_knot_X(X = Tt[, 1], df = df_T, degree = degree)
rho <- 0.05

#train nuisance

#Obtain the \hat{\Gamma}(X_i), i = 1,2,......, for our proposed method
#For the randomized trail setting, the generalized propensity score is always a uniform distribution
#So the \hat{\Gamma}(X_i) can be calculated by a direct integration
#For the observation study, please refer to Section S2.1 for the estimation of \hat{\Gamma}(X_i)
Basis_Matrix <- Basis_M_Build(
    V %>% as.matrix,
    Tt  %>% as.matrix,
    df_V,
    df_T,
    Vrange,
    Trange,
    V_knot_list,
    T_knot_list
  )
basis_V <- Basis_Matrix$basis_V
basis_M <- Basis_Matrix$basis_M
basis_T <- Basis_Matrix$basis_T
Gamma_T <- matrix(0, nrow = n, ncol = df_T)
for (i in 1:df_T) {
  for (l in 1:n) {
    t_x <- function(Tt) {
      b = bSpline(
        x = Tt,
        knots = T_knot_list[[1]],
        degree = degree,
        Boundary.knots = Trange[[1]],
        intercept = T
      )[, i]
      b %>% as.numeric
    }
    Gamma_T[l, i] <- MASS::area(f = t_x, a = 0.001, b = 0.999)
  }
}
Gamma <- matrix(1, nrow = n, ncol = 1)
for (i in 1:1) {
  Gamma <- rowwise_kron(Gamma, basis_V[[i]])
}
Gamma <- rowwise_kron(Gamma, Gamma_T)

#Obtain the \Psi(X_i,T_i) - \hat{\Gamma}(X_i) in our proposed method
Gamma_diff <- basis_M - Gamma

#Obtain the outcome regression estimator for our proposed method
mX <- cbind(X)
m_model <- SuperLearner(
  Y = Y,
  X = cbind(X) %>% data.frame,
  family = gaussian(),
  SL.library = c("SL.polymars", "SL.nnet", "SL.lm", "SL.earth", "SL.rpart")
)
m <- predict(m_model, mX %>% data.frame, onlySL = TRUE)$pred

#obtain the coefficients of our B-spline estimator

coeff <- solve_coefficient(
    X = X,
    Tt = Tt,
    Y = Y,
    df_V = df_V,
    df_T = df_T,
    Gamma_diff = Gamma_diff,
    m_e = m,
    Vrange_list = Vrange,
    Trange_list = Trange,
    V_knot_list = V_knot_list,
    T_knot_list = T_knot_list,
    rho)

#plot the estimated CATE function

Ttt <- seq(0,1,0.01)
V1 <- seq(0,1,0.01)
tau_sample_res <- matrix(0,nrow = length(Ttt), ncol = length(V1))
for (j in 1:length(Ttt)){
  fit_model = est_tau(V_pred = V1 %>% as.matrix,
                      Tt_pred = rep(Ttt[j],length(V1)) %>% as.matrix,
                      para = coeff$phi %>% as.vector,
                      df_V = df_V,
                      df_T = df_T,
                      Vrange_list = Vrange,
                      Trange_list = Trange,
                      V_knot_list = V_knot_list,
                      T_knot_list = T_knot_list,
                      rho)
  tau_sample_res[j,] <- fit_model$tau
}
plot_fun(Ttt,V1,tau_sample_res)

