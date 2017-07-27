
library(fda)
library(mgcv)
library(reshape2)
library(ggplot2)

if(length(grep("bdsegal",getwd()))>0 ){
  computer <- "C:/Users/bdsegal"
} else{
  computer <- "/home/bsegal"
}
paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")
posterPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/poster/plots")

# load existing data
simData = read.csv("simData.csv")
trueMean = read.csv("trueMean.csv")

# l2 models
m1 <- gamm(y ~ s(x, bs = "ps", m = c(0,1), k = 21), data = simData, 
           random = list(id = ~1))

summary(m1$gam)
VC <- VarCorr(m1$lme)
lambdaGAMM <-as.numeric(VC[nrow(VC), 1]) / as.numeric(VC[2,1])
tauGAMM <-as.numeric(VC[nrow(VC), 1]) / as.numeric(VC[nrow(VC) - 1,1])
cbind(tauGAMM, lambdaGAMM)
#         tauGAMM lambdaGAMM
# [1,] 0.01016958  0.6052437

# sigma2epsilon, sigma2beta, sigma2b
c(as.numeric(VC[nrow(VC), 1]), as.numeric(VC[2,1]), as.numeric(VC[nrow(VC) - 1,1]))
# [1] 0.01062245 0.01755070 1.04453170

# p <- plot(m1$gam, se = 2)
# dev.new()
# rang <- range(c(p[[1]]$fit[, 1] + 2 * p[[1]]$se, p[[1]]$fit[, 1] - 2 * p[[1]]$se))
# plot(p[[1]]$x, p[[1]]$fit[, 1], ylim = rang, type = "l")
# lines(p[[1]]$x, p[[1]]$fit[, 1] + p[[1]]$se, lty = 2)
# lines(p[[1]]$x, p[[1]]$fit[, 1] - p[[1]]$se, lty = 2)

# confidence bands
B <- 10000
beta <- coef(m1$gam)
br <- rmvn(B, beta, m1$gam$Vp)
X <- predict(m1$gam, newdata = trueMean, type= "lpmatrix")
yPred <- X %*% t(br)
alpha <- 0.05
CI <- t(apply(yPred, 1, quantile, probs = c(alpha/2, 1-alpha/2)))

simData$yHat <- predict(m1$lme)
trueMean$yMarg <- predict(m1$gam, newdata = trueMean)
trueMean$lower <- CI[, 1]
trueMean$upper <- CI[, 2]

CIpoly <- data.frame(x = c(trueMean$x, rev(trueMean$x)), 
                     y = c(trueMean$lower, rev(trueMean$upper)),
                     id = 0)
                     
dev.new()
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", linetype = "solid", size = 1)+
  theme_bw(30)+
  scale_y_continuous(lim = c(-2, 4))
ggsave(file.path(paperPath, "l2_GAMM_changePoint_CI_ord2_poster.png"))
ggsave(file.path(presentPath, "l2_GAMM_changePoint_CI_ord2_poster.png"))
ggsave(file.path(posterPath, "l2_GAMM_changePoint_CI_ord2_poster.png"))

ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yMarg), data = trueMean, color = "red", size = 1)+
  theme_bw(30)
ggsave(file.path(paperPath, "l2_GAMM_changePoint_point_ord2_poster.png"))
ggsave(file.path(presentPath, "l2_GAMM_changePoint_point_ord2_poster.png"))
ggsave(file.path(posterPath, "l2_GAMM_changePoint_point_ord2_poster.png"))
