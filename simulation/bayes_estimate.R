
library(psplinesl1)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
marDefault <- c(5, 4, 4, 2) + 0.1 

if(length(grep("bdsegal",getwd())) > 0) {
    computer <- "C:/Users/bdsegal"
} else {
  computer <- "/home/bsegal"
}
paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")

source(file.path(computer,"Dropbox/Research/psplines_L1_penalty/Rpackage/psplinesL1/R/l1mixedPackage.R"))
# source(file.path(computer,"Dropbox/Research/psplines_L1_penalty/Rpackage/psplinesL1/R/L1_CV.R"))

simData <- read.csv("simData.csv")
trueMean <- read.csv("trueMean.csv")
n <- table(simData$id)
N <- length(n)

# get list of fixed effect smooths
k <- 1
X <- list(ps(x = "x", norder = 2, k = k, width = 0.05, data = simData))
Xpoint <- list(ps(x = "x", norder = 2, k = k, width = 0.05, data = trueMean))
p <- X[[1]]$basis$nbasis

Zlist <- list()
for (i in 1:length(n)) {
  Zlist[[i]] <- rep(1, n[i])
}
Z <- as.matrix(bdiag(Zlist))

stanDat <- list(
  n = sum(n),
  N = length(n),
  p = p,
  k = k,
  y = simData$y,
  F = X[[1]]$F,
  D = X[[1]]$D,
  id = simData$id
)

# laplace prior
fitLap <- stan(file="bayes_lap.stan", data=stanDat, iter=2000, chains=4)	
save(fitLap, file = "fit_lap.RData")
rm(fitLap)

# normal prior
fitNorm <- stan(file="bayes_norm.stan", data=stanDat, iter=2000, chains=4)	
save(fitNorm, file = "fit_norm.RData")
rm(fitNorm)

# diffuse prior
fitNoPen <- stan(file="bayes_noPen.stan", data=stanDat, iter=2000, chains=4)
save(fitNoPen, file = "fit_noPen.RData")
rm(fitNoPen)


# plot results ----------------------------------------------------------------

alpha <- 0.05

# Laplace
load("fit_noPen.RData")
fit <- fitLap

traceplot(fit, pars = c("Dbeta"), inc_warmup = FALSE)
traceplot(fit, pars = c("b"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB", "sigmaLambda"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB"), inc_warmup = FALSE)
traceplot(fit, pars = c("beta0", "beta"), inc_warmup = FALSE)

print(fit)
plot(fit)

print(fit, c("sigmaB"), probs=c(0.1,0.9))

fitExt <- extract(fit, permuted = TRUE) # return a list of arrays 
beta0Hat<- mean(fitExt$beta0)
betaHat <- apply(fitExt$beta, 2, mean)
bHat <- apply(fitExt$b, 2, mean)

# Credible interval for marginal fit
curves0 <- Xpoint[[1]]$F %*% t(fitExt$beta)
curves <- curves0 + matrix(fitExt$beta0, nrow = 1)[rep(1, nrow(curves0)), ]

CI <- t(apply(curves, 1, quantile, probs = c(alpha/2, 1-alpha/2)) )

simData$yHat <- beta0Hat + X[[1]]$F %*% betaHat + Z %*% bHat
trueMean$yMarg <- beta0Hat + Xpoint[[1]]$F %*% betaHat
trueMean$lower <- CI[, 1]
trueMean$upper <- CI[, 2]

CIpoly <- data.frame(x = c(trueMean$x, rev(trueMean$x)), 
                     y = c(trueMean$lower, rev(trueMean$upper)),
                     id = 0)

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  scale_y_continuous(lim = c(-2, 4))+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_lap_CI.png"))

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_lap_point.png"))

# normal
load("fit_norm.RData")
fit <- fitNorm

traceplot(fit, pars = c("Dbeta"), inc_warmup = FALSE)
traceplot(fit, pars = c("b"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB", "sigmaLambda"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB"), inc_warmup = FALSE)
traceplot(fit, pars = c("beta0", "beta"), inc_warmup = FALSE)

print(fit)
plot(fit)

print(fit, c("sigmaB"), probs=c(0.1,0.9))

fitExt <- extract(fit, permuted = TRUE) # return a list of arrays 
beta0Hat<- mean(fitExt$beta0)
betaHat <- apply(fitExt$beta, 2, mean)
bHat <- apply(fitExt$b, 2, mean)

# Credible interval for marginal fit
curves0 <- Xpoint[[1]]$F %*% t(fitExt$beta)
curves <- curves0 + matrix(fitExt$beta0, nrow = 1)[rep(1, nrow(curves0)), ]

CI <- t(apply(curves, 1, quantile, probs = c(alpha/2, 1-alpha/2)) )

simData$yHat <- beta0Hat + X[[1]]$F %*% betaHat + Z %*% bHat
trueMean$yMarg <- beta0Hat + Xpoint[[1]]$F %*% betaHat
trueMean$lower <- CI[, 1]
trueMean$upper <- CI[, 2]

CIpoly <- data.frame(x = c(trueMean$x, rev(trueMean$x)), 
                     y = c(trueMean$lower, rev(trueMean$upper)),
                     id = 0)

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  scale_y_continuous(lim = c(-2, 4))+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_normal_CI.png"))

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_normal_point.png"))

# no prior
load("fit_noPen.RData")
fit <- fitNoPen

traceplot(fit, pars = c("Dbeta"), inc_warmup = FALSE)
traceplot(fit, pars = c("b"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB", "sigmaLambda"), inc_warmup = FALSE)
traceplot(fit, pars = c("sigmaEpsilon", "sigmaB"), inc_warmup = FALSE)
traceplot(fit, pars = c("beta0", "beta"), inc_warmup = FALSE)

print(fit)
plot(fit)

print(fit, c("sigmaB"), probs=c(0.1,0.9))

fitExt <- extract(fit, permuted = TRUE) # return a list of arrays 
beta0Hat<- mean(fitExt$beta0)
betaHat <- apply(fitExt$beta, 2, mean)
bHat <- apply(fitExt$b, 2, mean)

# Credible interval for marginal fit
curves0 <- Xpoint[[1]]$F %*% t(fitExt$beta)
curves <- curves0 + matrix(fitExt$beta0, nrow = 1)[rep(1, nrow(curves0)), ]

CI <- t(apply(curves, 1, quantile, probs = c(alpha/2, 1-alpha/2)) )

simData$yHat <- beta0Hat + X[[1]]$F %*% betaHat + Z %*% bHat
trueMean$yMarg <- beta0Hat + Xpoint[[1]]$F %*% betaHat
trueMean$lower <- CI[, 1]
trueMean$upper <- CI[, 2]

CIpoly <- data.frame(x = c(trueMean$x, rev(trueMean$x)), 
                     y = c(trueMean$lower, rev(trueMean$upper)),
                     id = 0)

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  scale_y_continuous(lim = c(-2, 4))+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_noPen_CI.png"))

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  theme_bw(38)
ggsave(file.path(paperPath,"Bayes_changePoint_noPen_point.png"))
