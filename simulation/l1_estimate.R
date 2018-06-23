
library(psplinesl1)
library(reshape2)

if(length(grep("bdsegal",getwd()))>0 ){
  computer <- "C:/Users/bdsegal"
} else{
	computer <- "/home/bsegal"
}
paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")
posterPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/poster/plots")

# load existing data
data(simData)
data(trueMean)

# with L1 P-splines -----------------------------------------------------------

# get list of fixed effect smooths
X <- list(ps(x = "x", data = simData, norder = 2, k = 1, width = 0.05,
             center = TRUE))
rand <- re(x = "x", id = "id", data = simData, 
           randomCurves = FALSE)

# fitting one path at a time --------------------------------------------------
system.time({cvOut <- cv(y = "y", X = X,
             rand = rand,
             id = "id",
             K = 5,
             smoothInit = c(0, 0),
             lmeUpdate = FALSE,
             data = simData,
             se1 = FALSE)})
 #   user  system elapsed 
 # 14.988   0.040  15.053

dev.new(width = 9, height = 5)
cvOut
ggsave(file.path(presentPath, "cv.png"))

cvOut$smoothOpt
#        tau    lambda1 
# 0.02739983 0.16874138 

system.time({
  a1 <- admm(y = "y", id = "id", X = X, rand = rand,
             lambda = cvOut$smoothOpt[2:(length(X)+1)],
             tau = cvOut$smoothOpt[1],
             lmeUpdate = TRUE,
             rho = min(5, max(cvOut$smoothOpt)),
             centerZ = FALSE,
             data = simData,
             forCV = FALSE
             )
  })
  #  user  system elapsed 
  # 4.376   0.004   4.385 
a1$conv$iter
# 244

plot(log(a1$conv$rNorm))
plot(log(a1$conv$sNorm))

a1$fit$sigma2
# [1] 0.009289998

sum(a1$fit$residuals^2)
# [1] 4.047725

signif(a1$fit$df, 3)
#                  Overall   F1    Z
# Stein               14.3 10.0 3.29
# Restricted          14.6 10.0 3.63
# ADMM                13.6  9.0 3.63
# Ridge               22.1 17.7 3.31
# Ridge restricted    22.4 17.8 3.63

#sigma2b
a1$fit$sigma2b
# [1] 1.064899

# confidence bands
CI <- ci(model = a1, alpha = 0.05)

Fmarg <- eval.basis(trueMean$x, a1$params$X[[1]]$basis) %*% X[[1]]$Q
outOfRange <- which(trueMean$x < min(simData$x) | trueMean$x > max(simData$x))
trueMean$yHat <- a1$coef$beta0 + Fmarg %*% a1$coef$beta[[1]]
trueMean$yHat[outOfRange] <- NA

CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$lower, rev(CI[[1]]$upper)))

dev.new()              
ggplot(aes(x = x, y = y), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yHat), data = trueMean, color = "red", size = 1)+
  theme_bw(30)+
  scale_y_continuous(lim = c(-2, 4))
ggsave(file.path(paperPath,"l1_changePoint_CI_ord2_Bayes_poster.png"))
ggsave(file.path(presentPath,"l1_changePoint_CI_ord2_Bayes_poster.png"))


# example plot for package using only CI and no baseline truth
CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$lower, rev(CI[[1]]$upper)))

ggplot(aes(x = x, y = y), data = simData)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(aes(x = CI[[1]]$x, y = CI[[1]]$smooth), size = 1)+
  theme_bw(18)

# predicted curves
simData$yHat <- a1$fit$yHat

dev.new()
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yHat), data = trueMean, color = "red", size = 1)+
  theme_bw(24)
ggsave(file.path(paperPath,"l1_changePoint_point_ord2.png"))

dev.new()
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yHat), data = trueMean, color = "red", size = 1)+
  theme_bw(30)
ggsave(file.path(paperPath,"l1_changePoint_point_ord2_poster.png"))
ggsave(file.path(presentPath,"l1_changePoint_point_ord2_poster.png"))
