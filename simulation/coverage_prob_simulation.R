library(psplinesl1)
library(mgcv)
library(reshape2)

iterMax <- 1000

paperPath <- "../../paper/plots"

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

N <- 50
datTrue <- data.frame(x = xsim, y = yTrue, id = 0)

cp1 <- matrix(NA, nrow = length(xsim), ncol = iterMax)
cp2 <- matrix(NA, nrow = length(xsim), ncol = iterMax)
width1 <- matrix(NA, nrow = length(xsim), ncol = iterMax)
width2 <- matrix(NA, nrow = length(xsim), ncol = iterMax)

rownames(cp1) <- xsim
rownames(cp2) <- xsim
rownames(width1) <- xsim
rownames(width2) <- xsim

for (iter in 1:iterMax) {
  print(iter)

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

  # Fit models ------------------------------------------------------------------

  try({
    # l2 P-spline
    m2 <- gamm(y ~ s(x, bs = "ps", m = c(0,1), k = 21), data = newDat, 
               random = list(id = ~1))

    # confidence bands
    B <- 10000
    beta <- coef(m2$gam)
    br <- rmvn(B, beta, m2$gam$Vp)
    X <- predict(m2$gam, newdata = datTrue, type = "lpmatrix")
    yPred <- X %*% t(br)
    alpha <- 0.05
    CI <- t(apply(yPred, 1, quantile, probs = c(alpha/2, 1-alpha/2)))

    datTrue$lower <- CI[, 1]
    datTrue$upper <- CI[, 2]

    cp2[, iter] <- with(datTrue, y <= upper & y >= lower)
    width2[, iter] <- with(datTrue, upper - lower)
  })

  try({
    # l1 P-splines
    X <- list(ps(x = "x", norder = 2, k = 1, width = 0.05, data = newDat,
                 center = TRUE))
    rand <- re(x = "x", id = "id", data = newDat, 
               randomCurves = FALSE)

    cvOut <- cv(y = "y", X = X, rand = rand,
                id = "id",
                K = 5,
                se1 = FALSE,
                pathLength = 20,
                epsilonAbs = 1e-4,
                epsilonRel = 1e-4,
                iterMax = 1e3,
                data = newDat,
                verbose = 0)

    m1 <- admm(y = "y", id = "id", X = X, rand = rand,
                 lambda = cvOut$smoothOpt[2:(length(X)+1)],
                 rho = min(max(cvOut$smoothOpt), 5),
                 epsilonAbs = 1e-4,
                 epsilonRel = 1e-4,
                 iterMax = 1e3,
                 warm = NULL,
                 data = newDat,
                 lmeUpdate = TRUE,
                 verbose = FALSE
                 )

    confInt <- ci(m1, newData = datTrue)
    confIntTrue <- merge(confInt[[1]], datTrue[, 1:2], by = "x")
    cp1[, iter] <- with(confIntTrue, y <= upper & y >= lower)
    width1[, iter] <- with(confIntTrue, upper - lower)

  })

# save(cp1, cp2, width1, with2, file = "cp_sim.RData")
}

load("cp_sim.RData")

cp1mean <- apply(cp1, 1, mean, na.rm = TRUE)
cp2mean <- apply(cp2, 1, mean, na.rm = TRUE)

cp <- data.frame(x = rep(xsim, times = 2), 
                 mean = c(cp1mean, cp2mean),
                 Penalty = rep(c("L1", "L2"), each = length(xsim)))
cp$lower <- cp$mean - 1.96 * sqrt(cp$mean * (1 - cp$mean) / iterMax)
cp$upper <- cp$mean + 1.96 * sqrt(cp$mean * (1 - cp$mean) / iterMax)

dev.new(width = 7, height = 5)
ggplot(aes(x = x, y = mean, color = Penalty, shape = Penalty), data = cp) +
  geom_point() +
  geom_line() +
  theme_bw(18) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_y_continuous(lim = c(0.6, 1)) +
  labs(y = "Coverage probability")+
  scale_color_discrete(labels = c(quote('\u2113'[1]), quote('\u2113'[2]))) +
  scale_shape_discrete(labels = c(quote('\u2113'[1]), quote('\u2113'[2])))
ggsave(file.path(paperPath, "cp.png"))
