# analyze EDA data

library(fda)
library(mgcv)
library(ggplot2)
library(reshape2)

if(length(grep("bdsegal",getwd()))>0 ) {
    path <- "C:/Users/bdsegal/Documents/vigilance"
    computer <- "C:/Users/bdsegal"
} else {
	path <- "/home/bsegal/Documents/Research/data_empatica"
  computer <- "/home/bsegal"
}

paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")
posterPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/poster/plots")

marDefault <- c(5, 4, 4, 2) + 0.1 

# data prep -------------------------------------------------------------------
groupA <- read.csv(file.path(path, "groupA.csv"))
data <- groupA
data <- data[with(data, order(id, x)), ]
n <- table(data$id)

# random curves
xRange <- c(floor(min(data$x)), ceiling(max(data$x)))
basisZ <- create.bspline.basis(rangeval = xRange,
                               breaks = seq(xRange[1], xRange[2], 5),
                               norder = 4)
k <- basisZ$nbasis
Zlong <- eval.basis(data$x, basisZ)
Zlist <- list()
# Zqlist <- list()
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
Zbreve <- Z %*% eig$vectors
Zbreve1 <- as.matrix(Zbreve[, 1:q1])
Zbreve2 <- as.matrix(Zbreve[, (q1+1):ncol(Zbreve)])

Lambda1Sqrt <- diag(sqrt(eig$values[1:q1]))
Zcheck1 <- Zbreve1 %*% Lambda1Sqrt

system.time({
g1 <- gamm(y ~ s(x, bs = "ps", m = c(2, 2), k = k) + 
               s(x, bs = "ps", m = c(2, 2), k = k, by = high),
               # random = list(id = pdIdent(~Zcheck1 - 1)),
               random = list(id = pdIdent(~Zcheck1 - 1),
                             id = pdIdent(~Zbreve2 - 1)),
           method = "REML",    
           data = data)
})
#    user  system elapsed 
# 165.072   0.708 165.947 

plot(g1$gam)
summary(g1$gam)
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.07941    0.06298  -17.14   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Approximate significance of smooth terms:
#              edf Ref.df      F p-value    
# s(x)      14.160 14.160 35.537 < 2e-16 ***
# s(x):high  6.954  6.954  3.021 0.00614 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# R-sq.(adj) =  0.281   
#   Scale est. = 0.04196   n = 1700

# get indices of high vigilance subjects for plotting
beta <- coef(g1$gam)
highGroup <- grep("high", names(beta))

# simulate from posterior for simultaneous credible 
newData <- data.frame(x = seq(min(data$x), max(data$x), 0.1), high = 1)
X <- predict(g1$gam, newdata = newData, type= "lpmatrix")
newData$low <- as.vector(X[ ,-highGroup] %*% beta[-highGroup])
newData$diff <- as.vector(X[ ,highGroup] %*% beta[highGroup])

B <- 10000
br <- t(rmvn(B, beta, g1$gam$Vp))
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
ggplot(aes(x = x, y = low), data = newData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(size = 1, color = "red")+
  scale_y_continuous(lim = c(-3, 1))+
  labs(y = expression(paste("lo",g[10],"(EDA)", sep = "")), x = "Minute")+
  theme_bw(28)
ggsave(file.path(paperPath, "eda_l2_low_fit.png"))

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
ggsave(file.path(paperPath, "eda_l2_diff_fit.png"))

# statistically significant differences at 5% level
range(newData$x[newData$diffLower > 0])
# [1] 44.5 65.6

# plot predicted curves
data$yHat <- predict(g1$lme)
dataM <- melt(data, id.vars = c("x", "id", "type"),
                    measure.vars = c("y", "yHat"))
dataM$type <- factor(dataM$type, levels = c("low", "high"))
levels(dataM$type) <- c("Low vigilance", "High vigilance")

ID <- unique(dataM$id)

gg <- ggplot() + theme_bw(20) +
      scale_color_manual(guid = FALSE, "", values = c("blue", "red"), labels = c("High", "Low"))+
      scale_alpha_manual("", values = c(.3, 1), labels = c("Observed", "Predicted"))+
      scale_linetype_manual("", values = c("solid", "dashed"), labels = c("Observed", "Predicted"))+
      labs(y = expression(paste("lo", g[10], "(EDA)", sep = "")), x = "Minute")
for (i in 1:length(ID)) {
 gg <- gg + geom_line(data = dataM[which(dataM$id == ID[i]), ],
            aes(x = x, y = value,
                alpha = variable,
                color = type,
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
facet_wrap(~type)
ggsave(file.path(paperPath, "l2Pred.png"))

# MSE
with(data, mean((y - yHat)^2))
# [1] 0.03478479