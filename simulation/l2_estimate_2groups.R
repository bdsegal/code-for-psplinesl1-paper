
library(psplinesl1)
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
data(simData2groups)
data(trueMean)
data(trueInteraction)

# l2 models
m1 <- gamm(y ~ s(x, bs = "ps", m = c(0,1), k = 21) +
               s(x, bs = "ps", m = c(0,1), k = 21, by = group),
           data = simData2groups, 
           random = list(id = ~1))

m1$gam$smooth[[2]]$knots
summary(m1$gam)
VC <- VarCorr(m1$lme)

# # sigma2epsilon, sigma2b
c(as.numeric(VC[nrow(VC), 1]), as.numeric(VC[nrow(VC) - 1,1]))
# [1] 0.009892428 1.024829289

# confidence bands
beta <- coef(m1$gam)
highGroup <- grep("group", names(beta))
newData <- data.frame(x = seq(min(simData2groups$x), max(simData2groups$x), 0.1),
                      group = 1)
X <- predict(m1$gam, newdata = newData, type= "lpmatrix")
newData$low <- as.vector(X[ ,-highGroup] %*% beta[-highGroup])
newData$diff <- as.vector(X[ ,highGroup] %*% beta[highGroup])

B <- 10000
br <- t(rmvn(B, beta, m1$gam$Vp))
predLow <- X[ ,-highGroup] %*% br[-highGroup, ]
predDiff <- X[ ,highGroup] %*% br[highGroup, ]

alpha <- 0.05
CILow <- t(apply(predLow, 1, quantile, probs = c(alpha/2, 1-alpha/2)))
CIDiff <- t(apply(predDiff, 1, quantile, probs = c(alpha/2, 1-alpha/2)))

newData$lowLower <- CILow[, 1]
newData$lowUpper <- CILow[, 2]
newData$diffLower <- CIDiff[, 1]
newData$diffUpper <- CIDiff[, 2]

# plot low vigilance
CIpoly <- data.frame(x = c(newData$x, rev(newData$x)), 
                     low = c(newData$lowLower, rev(newData$lowUpper)),
                     est = "Low vigilance")
dev.new()
ggplot(aes(x = x, y = low), data = newData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(size = 1, color = "red")+
  geom_line(aes(y = y), data = trueMean, color = "black", size = 1)+
  labs(x = "x", y = "y") +
  theme_bw(30) +
  scale_y_continuous(lim = c(-0.1, 2.5))
ggsave(file.path(paperPath, "l2_2groups_group0.png"))

# plot high - low vigilance
CIpoly <- data.frame(x = c(newData$x, rev(newData$x)), 
                     diff = c(newData$diffLower, rev(newData$diffUpper)),
                     est = "Low vigilance")
dev.new()
ggplot(aes(x = x, y = diff), data = newData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(size = 1, color = "purple")+
  labs(x = "x", y = "y") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line(aes(y = y), data = trueInteraction, color = "black", size = 1)+
  theme_bw(30) +
  scale_y_continuous(lim = c(-0.5, 1.25))
ggsave(file.path(paperPath, "l2_2groups_group1.png"))
