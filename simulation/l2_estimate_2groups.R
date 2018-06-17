
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
  geom_line(size = 1, color = "red")+
  labs(x = "x", y = "y") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_line(aes(y = y), data = trueInteraction, color = "black", size = 1)+
  theme_bw(30) +
  scale_y_continuous(lim = c(-0.5, 1.25))
ggsave(file.path(paperPath, "l2_2groups_group1.png"))








B <- 10000
br <- rmvn(B, beta, m1$gam$Vp)
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






beta <- coef(m1$gam)
highGroup <- grep("group", names(beta))

# simulate from posterior for simultaneous credible 
newData <- data.frame(x = seq(min(simData2groups$x), max(simData2groups$x), 0.1), group = 1)
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
str(CIpoly)
ggplot(aes(x = x, y = low), data = newData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(size = 1, color = "red") +
  # scale_y_continuous(lim = c(-3, 1))+
  labs(y = expression(paste("lo",g[10],"(EDA)", sep = "")), x = "Minute")+
  theme_bw(28)
# ggsave(file.path(paperPath, "eda_l2_low_fit_alt.png"))

# plot high - low vigilance
CIpoly <- data.frame(x = c(newData$x, rev(newData$x)), 
                     diff = c(newData$diffLower, rev(newData$diffUpper)),
                     est = "Low vigilance")
ggplot(aes(x = x, y = diff), data = newData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(size = 1, color = "purple")+
  scale_y_continuous(lim = c(-1, 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(y = expression(paste("lo",g[10],"(EDA)", sep = "")), x = "Minute")+
  theme_bw(28)
# ggsave(file.path(paperPath, "eda_l2_diff_fit_alt.png"))





# plot predicted curves
simData2groups$yHat <- predict(m1$lme)
dataM <- melt(simData2groups, id.vars = c("x", "id", "group"),
                    measure.vars = c("y", "yHat"))
dataM$group <- factor(dataM$group, levels = c("low", "high"))
levels(dataM$group) <- c("Low vigilance", "High vigilance")

ID <- unique(dataM$id)

gg <- ggplot() + 
      theme_bw(26) +
      scale_color_manual(guide = FALSE, "", values = c("blue", "red"), labels = c("High", "Low"))+
      scale_alpha_manual("", values = c(.3, 1), labels = c("Observed", "Predicted"))+
      scale_linetype_manual("", values = c("solid", "dashed"), labels = c("Observed", "Predicted"))+
      labs(y = expression(paste("lo", g[10], "(EDA)", sep = "")), x = "Minute")
for (i in 1:length(ID)) {
 gg <- gg + geom_line(data = dataM[which(dataM$id == ID[i]), ],
            aes(x = x, y = value,
                alpha = variable,
                color = group,
                linetype = variable))
}

dev.new(width = 9.5, height = 6.5)
gg  +
guides(alpha = guide_legend(
                 keywidth=0.35,
                 keyheight=0.35,
                 default.unit="inch"),
            linetype = guide_legend(
                 keywidth=0.35,
                 keyheight=0.35,
                 default.unit="inch")
      )+
theme(legend.position = "bottom")+
facet_wrap(~group)
# ggsave(file.path(paperPath, "l2Pred_alt.png"))


ggplot(aes(x = x, y = yHat, group = id, color = as.factor(group)), data = simData2groups) +
  geom_line(aes(y = y), color = "grey") +
  geom_line() +
  facet_wrap(~ group)