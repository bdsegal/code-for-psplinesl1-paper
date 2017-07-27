# setwd("~/L1Pspline/Bayes")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(fda)
library(Matrix)

# function for setting up splines ---------------------------------------------
ps <- function(x, norder, k, by = NULL, width = NULL, delta = 1e-7) {

  xRange <- c(floor(min(x)), ceiling(max(x)) + delta)
  if(is.null(width)) width <- diff(xRange) / 10
  basis <- create.bspline.basis(rangeval = xRange,
                                breaks = seq(xRange[1], xRange[2], width),
                                norder = norder)
  F <- eval.basis(x, basis)
  p <- basis$nbasis
  D <- diff(diag(p), diff = k + 1)

  if (is.null(by)) {

    # centering constraints
    Ct <- as.matrix(colSums(F), nrow = ncol(F))
    Q <- qr.Q(qr(Ct), complete = TRUE)[, -1]
    Fq <- F %*% Q
    Dq <- D %*% Q

    return(list(basis = basis, F = Fq, Q = Q, p = p,
                D = Dq, by = by, k = k, width = width))
  }
  else if (!is.null(by)) {
    F <- sweep(F, 1, by, `*`)

    return(list(basis = basis, F = F, Q = NULL, p = p,
                D = D, by = by, k = k, width = width))
  }
  else {
    warning("'by' must be either NULL or equal to variable name")
  }
}

# data prep -------------------------------------------------------------------
groupA <- read.csv("groupA.csv")
data <- groupA
data <- data[with(data, order(id, x)), ]
n <- table(data$id)

# get list of fixed effect smooth matrices and associated quantities
X <- list(ps(x = data$x, norder = 4, k = 1, width = 5),
          ps(x = data$x, norder = 4, k = 1, width = 5, 
                   by = data$high))

# random curves
xRange <- c(floor(min(data$x)), ceiling(max(data$x)))
basisZ <- create.bspline.basis(rangeval = xRange,
                                breaks = seq(xRange[1], xRange[2], 5),
                                norder = 4)
Zlong <- eval.basis(data$x, basisZ)
Zlist <- list()
cumn <- c(0,cumsum(n))

for (i in 1:length(n)) {
  Zlist[[i]] <- Zlong[(cumn[i]+1):cumn[i+1], ]
}
Z <- bdiag(Zlist)

Stemp <- bsplinepen(basisZ, Lfdobj=2)
Slist <- list()
Sqlist <- list()
for (i in 1:length(n)) {
  Slist[[i]] <- Stemp
}
S <- bdiag(Slist)

# rotate random effect design matrix and split into
# penalized and unpenalized parameters
eig <-eigen(S, symmetric = TRUE)
q1 <- sum(eig$values > .Machine$double.eps)
q2 <- ncol(Z) - q1
Ztilde <- Z %*% eig$vectors
Ztilde1 <- as.matrix(Ztilde[, 1:q1])
Ztilde2 <- as.matrix(Ztilde[, (q1+1):ncol(Ztilde)])

Lambda1Sqrt <- diag(sqrt(eig$values[1:q1]))
Zcheck1 = Ztilde1 %*% Lambda1Sqrt

stanDat <- list(
  n = sum(n),
  q1 = q1,
  q2 = q2,
  y = data$y,
  p1 = X[[1]]$p,
  p2 = X[[2]]$p,
  k1 = X[[1]]$k,
  k2 = X[[2]]$k,
  X1 = X[[1]]$F,
  X2 = X[[2]]$F,
  D1 = X[[1]]$D,
  D2 = X[[2]]$D,
  Zcheck1 = Zcheck1,
  Ztilde2 = Ztilde2
)

save(X, Zcheck1, Ztilde2, file = "sqrt_stanDat_norm.RData")

fitLap <- stan(file = "EDA_lap_sqrt_norm.stan", 
               data = stanDat, iter = 5000, chains = 4) 
save(fitLap, file = "EDA_lap_sqrt_norm.RData")
