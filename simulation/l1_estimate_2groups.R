library(psplinesl1)
library(reshape2)

paperPath <- "../../paper/plots"
presentPath <- "../../presentation/plots"
posterPath <- "../../poster/plots"

# load existing data
data(simData2groups)
data(trueMean)
data(trueInteraction)

# with L1 P-splines -----------------------------------------------------------

# get list of fixed effect smooths
X <- list(ps(x = "x", data = simData2groups, norder = 2, k = 1, width = 0.05),
          ps(x = "x", data = simData2groups, norder = 2, k = 1, width = 0.05,
             by = "group", center = FALSE))
rand <- re(x = "x", id = "id", data = simData2groups, 
           randomCurves = FALSE)

# fitting one path at a time --------------------------------------------------
system.time({cvOut <- cv(y = "y",
             X = X,
             rand = rand,
             id = "id",
             K = 5,
             lmeUpdate = FALSE,
             data = simData2groups,
             se1 = FALSE)
})
 #   user  system elapsed 
 # 77.848   0.004  77.926 

dev.new(width = 10, height = 5)
cvOut$gg + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(presentPath, "cv_2groups.png"))

cvOut$smoothOpt
#        tau    lambda1    lambda2 
# 0.06233692 0.20944073 0.17236338 

system.time({
  a1 <- admm(y = "y", id = "id", X = X, rand = rand,
             lambda = cvOut$smoothOpt[-1],
             lmeUpdate = TRUE,
             rho = min(5, max(cvOut$smoothOpt)),
             centerZ = FALSE,
             data = simData2groups,
             forCV = FALSE
             )
  })
  #  user  system elapsed 
  # 5.564   0.016   5.588 

a1$conv$iter
# [1] 198

plot(log(a1$conv$rNorm))
plot(log(a1$conv$sNorm))

a1$fit$sigma2
# [1] 0.008952774

sum(a1$fit$residuals^2)
# [1] 7.802635

signif(a1$fit$df, 3)
#                  Overall   F1   F2    Z
# Stein               28.5 12.0  9.0 6.47
# Restricted          29.1 12.0  9.0 7.14
# ADMM                28.1 11.0  9.0 7.14
# Ridge               43.0 17.7 17.8 6.46
# Ridge restricted    45.7 18.6 19.0 7.14

#sigma2b
a1$fit$sigma2b
# [1] 1.042587

# confidence bands
CI <- ci(model = a1, alpha = 0.05)

# baseline group
Fmarg <- eval.basis(trueMean$x, a1$params$X[[1]]$basis) %*% X[[1]]$Q
outOfRange <- which(trueMean$x < min(simData2groups$x) | 
                    trueMean$x > max(simData2groups$x))
trueMean$yHat <- a1$coef$beta0 + Fmarg %*% a1$coef$beta[[1]]
trueMean$yHat[outOfRange] <- NA

CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$lower, rev(CI[[1]]$upper)))

dev.new()              
ggplot(aes(x = x, y = y), data = simData2groups)+
  geom_polygon(data = CIpoly, fill = "grey")+
  geom_line(data = trueMean, color = "black", size = 1)+
  geom_line(aes(y = yHat), data = trueMean, color = "red", size = 1)+
  theme_bw(30) +
  scale_y_continuous(lim = c(-0.1, 2.5))
ggsave(file.path(paperPath,"l1_2groups_group0.png"))
ggsave(file.path(presentPath,"l1_2groups_group0.png"))

Fmarg <- eval.basis(trueInteraction$x, a1$params$X[[2]]$basis)
outOfRange <- which(trueInteraction$x < min(simData2groups$x) | 
                    trueInteraction$x > max(simData2groups$x))
trueInteraction$yHat <- Fmarg %*% a1$coef$beta[[2]]
trueInteraction$yHat[outOfRange] <- NA

CIpoly <- data.frame(x = c(CI[[2]]$x, rev(CI[[2]]$x)), 
                     y = c(CI[[2]]$lower, rev(CI[[2]]$upper)))

dev.new()              
ggplot(aes(x = x, y = y), data = simData2groups)+
  geom_polygon(data = CIpoly, fill = "grey") +
  geom_line(data = trueInteraction, color = "black", size = 1) +
  geom_line(aes(y = yHat), data = trueInteraction, color = "purple", size = 1)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(30) +
  scale_y_continuous(lim = c(-0.5, 1.25))
ggsave(file.path(paperPath,"l1_2groups_group1.png"))
ggsave(file.path(presentPath,"l1_2groups_group1.png"))
