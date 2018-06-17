# run on cluster
args <- commandArgs(trailingOnly = TRUE)
M <- args[1]

setwd("/home/bdsegal/L1Pspline")

library(psplinesl1)
library(mgcv)
library(reshape2)

results <- list()

for (iter in 1:100) {
  print(iter)

  # Simulate data ---------------------------------------------------------------
  xsim <- seq(0, 1, 0.01)
  basis <- create.bspline.basis(rangeval = c(0, 1),
                                  breaks = seq(0, 1, 0.075),
                                  norder = 4)
  Xsim <- eval.basis(xsim, basis)
  Xsim <- cbind(1, Xsim)

  yTrue <- 1 + (xsim > 0.2) * (xsim - 0.2)*2 - 
               (xsim > 0.4) * (xsim - 0.4)*4 +
               (xsim > 0.6) * (xsim - 0.6)*8 - 
               (xsim > 0.8) * (xsim - 0.8)*16 

  # note: see change_point_simulation.R for possible sd values
  N <- 50
  datTrue <- data.frame(x = xsim, y = yTrue, id = 0)
  b <- rnorm(N, sd = 1)

  Y <- matrix(rep(yTrue, N), ncol = N)
  Y <- Y + matrix(rnorm(length(Y), sd = 0.1), ncol = N)
  Y <- sweep(Y, 2, b, '+')
  rownames(Y) <- xsim
  colnames(Y) <- 1:N
  dat <- melt(Y, varnames = c("x", "id"), value.name = "y")

  ids <- unique(dat$id)
  ni <- length(xsim)
  run <- 20
  keepList <- NULL
  for (i in 1:length(ids)) {
    begin <- sample(size = 1, 1:ni)
    if (begin < ni / 2) {
      end <- min(ni, begin + run)
    } else {
      end <- max(0, begin - run)
    }
    keepList[[i]] <- (ni*(i-1) + begin):(ni*(i-1) + end)
  }
  keep <- do.call(c, keepList)
  newDat <- dat[keep, ]

  # thin the data
  leave <- sample(size = nrow(newDat) / 1.75, 1:nrow(newDat))
  newDat <- newDat[! 1:nrow(newDat) %in% leave, ]

  newDat <- newDat[with(newDat, order(id, x)), ]
  n <- table(newDat$id)

  # Fit models ------------------------------------------------------------------

  # l2 P-spline
  m1 <- gamm(y ~ s(x, bs = "ps", m = c(0,1), k = 21), data = newDat, 
             random = list(id = ~1))

  # l1 P-splines
  X <- list(ps(x = "x", norder = 2, k = 1, width = 0.05, data = newDat,
               center = TRUE))
  rand <- re(x = "x", id = "id", data = newDat, 
             randomCurves = FALSE)

  # fitting one path at a time
  cvOut <- cv(y = "y", X = X, rand = rand,
              id = "id",
              K = 5,
              se1 = FALSE,
              pathLength = 20,
              epsilonAbs = 1e-4,
              epsilonRel = 1e-4,
              iterMax = 1e3,
              data = newDat,
              verbose = 2
              )

  a1 <- admm(y = "y", id = "id", X = X, rand = rand,
               lambda = cvOut$smoothOpt[2:(length(X)+1)],
               tau = cvOut$smoothOpt[1],
               rho = min(max(cvOut$smoothOpt), 5),
               epsilonAbs = 1e-4,
               epsilonRel = 1e-4,
               iterMax = 1e3,
               warm = NULL,
               data = newDat,
               verbose = FALSE
               )

  # x for evaluation
  xNotDup <- !duplicated(m1$lme$data$x)
  xUniq <- m1$lme$data$x[xNotDup]
  ord <- order(xUniq)
  x <- xUniq[ord]
  dx <- diff(x, differences = 1)

  # yHat from L2 fit
  yHatL2 <- predict(m1$gam)[xNotDup][ord]

  # yHat from L1 fit
  F1 <- X[[1]]$F[xNotDup, ][ord, ]
  yHatL1 <- a1$coef$beta0 + F1 %*% a1$coef$beta[[1]]

  # results
  fit <- data.frame(x = x, yHatL1 = yHatL1, yHatL2 = yHatL2)

  results[[iter]] <- fit

}

save(results, file = paste("batch/results", M, ".RData", sep = ""))

# plot(fit$x, fit$yHatL1, type = "l", col = "red")
# lines(fit$x, fit$yHatL2)