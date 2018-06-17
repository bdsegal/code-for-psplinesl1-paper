# analyze EDA data

library(psplinesl1)
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

# data prep -------------------------------------------------------------------
groupA <- read.csv(file.path(path, "groupA.csv"))
data <- groupA
data <- data[with(data, order(id, x)), ]
# n <- table(data$id)

# setup design matrices
X <- list(ps(x = "x", norder = 4, k = 1, data = data, width = 5),
          ps(x = "x", norder = 4, k = 1, data = data, width = 5,
             by = "high", center = FALSE))
rand <- re(x = "x", id = "id", data = data, 
            randomCurves = TRUE,  width = 5, 
            norder = 4, derivOrder = 2)

S <- diag(1, nrow(rand$S))
colnames(S) <- colnames(rand$S)
rownames(S) <- rownames(rand$S)
rand$S <- S

# fitting one path at a time --------------------------------------------------
system.time({cvOut <- cv(y = "y", X = X,
             rand = rand,
             id = "id",
             K = 3,
             data = data,
             se1 = FALSE)
             })

save(cvOut, file = "cvOut_EDA_smoothInit0.Rdata")

# To do: implement lme for random curves
system.time({
  m1 <- admm(y = "y", X = X, rand = rand,
             id = "id",
             tau = 0.05,
             lambda = c(0.09, .5),
             # tau = 450,
             # lambda = c(1, 10),
             # tau = 10,
             # lambda = c(0.5, 5),
             rho = 5,
             epsilonAbs = 1e-4,
             epsilonRel = 1e-4,
             iterMax = 1e3,
             warm = NULL,
             data = data,
             forCV = FALSE,
             centerZ = FALSE)
})
 #   user  system elapsed 
 # 11.108   0.028  11.154 

m1$conv$iter
# [1] 323

signif(m1$fit$df, 3)
#                  Overall   F1   F2   Z
# Stein                470 12.0  2.0 455
# Restricted           483 12.0  2.0 468
# ADMM                 482 11.0  2.0 468
# Ridge                471 19.9 10.6 439
# Ridge restricted     519 27.6 22.5 468

plot(log(m1$conv$rNorm))
plot(log(m1$conv$sNorm))

m1$coefs$beta0
# [1] -1.067598

#sigma2b
m1$fit$sigma2 / m1$params$tau
# [1] 0.1792059

m1$fit$sigma2
# [1] 0.008960295

m1$fit$sigma2Ridge
# [1] 0.008967313

data$yHat <- m1$fit$yHat
data$randEff <- with(m1, as.vector(params$Z %*% coefs$b))
data$bRem <- with(data, y - randEff)

# residuals
dev.new()
qplot(data$yHat, data$y - data$yHat)+
  theme_bw(22)+
  labs(x = expression(hat(y)), y = expression(paste("y - ",hat(y))))
ggsave(file.path(paperPath, "EDA_L1_resid_2.png"))

# confidence bands
CI <- ci(model = m1, alpha = 0.05)

CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$yLowerBayesQuick, rev(CI[[1]]$yUpperBayesQuick)),
                     id = 0)

CIpoint <- data.frame(x= CI[[1]]$x, y = CI[[1]]$smooth, id = 0)
            
ggplot(aes(x = x, y = y, group = id), data = CIpoly)+
  geom_polygon(fill = "grey")+
  geom_line(data = CIpoint, size = 1, color = "red")+
  theme_bw(28)+
  labs(x = "Minute", y = expression(paste("lo",g[10], "(EDA)", sep = "")))+
  scale_y_continuous(lim = c(-3,1))
ggsave(file.path(paperPath, "EDA_L1_smooth1_poster_alt.png"))
ggsave(file.path(presentPath, "EDA_L1_smooth1_poster_alt.png"))

CIpoly <- data.frame(x = c(CI[[2]]$x, rev(CI[[2]]$x)), 
                     y = c(CI[[2]]$yLowerBayesQuick, rev(CI[[2]]$yUpperBayesQuick)),
                     id = 0)
CIpoint <- data.frame(x= CI[[2]]$x, y = CI[[2]]$smooth, id = 0)
            
ggplot(aes(x = x, y = y, group = id), data = CIpoly)+
  geom_polygon(fill = "grey")+
  geom_line(data = CIpoint, size = 1, color = "purple")+
  theme_bw(28)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x = "Minute", y = expression(paste("lo",g[10], "(EDA)", sep = "")))+
  scale_y_continuous(lim = c(-1, 1))
ggsave(file.path(paperPath, "EDA_L1_smooth2_poster_alt.png"))
ggsave(file.path(presentPath, "EDA_L1_smooth2_poster_alt.png"))

# predicted curves
dataM <- melt(data, id.vars = c("x", "id", "type"),
                    measure.vars = c("y", "yHat"))
dataM$type <- factor(dataM$type, levels = c("low", "high"))
levels(dataM$type) <- c("Low vigilance", "High vigilance")

ID <- unique(dataM$id)
gg <- ggplot() + theme_bw(26) +
      # scale_color_manual("", values = c("grey", "black"), labels = c("Observed", "Predicted"))+
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
ggsave(file.path(paperPath, "EDA_L1_obs_pred_alt.png"))
ggsave(file.path(presentPath, "EDA_L1_obs_pred_alt.png"))
