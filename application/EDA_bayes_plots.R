
library(rstan)
library(fda)
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
# groupA <- read.csv("groupA.csv")
groupA <- read.csv(file.path(path, "groupA.csv"))
data <- groupA
data <- data[with(data, order(id, x)), ]
n <- table(data$id)
ord <- order(data$x)

# View results ----------------------------------------------------------------
# load("bayes_cauchy/sqrt_stanDat_cauchy.RData")
# load("bayes_cauhcy/EDA_lap_sqrt_cauchy.RData")

load("bayes_norm/sqrt_stanDat_norm.RData")
load("bayes_norm/EDA_lap_sqrt_norm.RData")

# order design matrices by increasing x to facilitate plotting
F <- list()
for(j in 1:length(X)) {
  F[[j]] <- X[[j]]$F[ord, ]
}

# WARNING: these might take a long time to load
traceplot(fitLap, pars = c("Dbeta1"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("sigmaEpsilon", "sigmaB", "sigmaB2",
          "sigmaLambda1", "logsigmaLambda2"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("beta0"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("beta1"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("beta2"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("btilde1"), inc_warmup = FALSE)
traceplot(fitLap, pars = c("btilde2"), inc_warmup = FALSE)

print(fitLap)
# plot(fitLap)

# extract chains and get posterior means
fitExt <- extract(fitLap, permuted = TRUE)
beta0Hat<- mean(fitExt$beta0)
beta1Hat <- apply(fitExt$beta1, 2, mean)
beta2Hat <- apply(fitExt$beta2, 2, mean)
btilde1Hat <- apply(fitExt$btilde1, 2, mean)
btilde2Hat <- apply(fitExt$btilde2, 2, mean)

# plot random effects
matplot(matrix(Zcheck1 %*% btilde1Hat + Ztilde2 %*% btilde2Hat, nrow = 100, ncol = 17), type = "l")

# Credible interval for marginal means ----------------------------------------
alpha <- 0.05

curves10 <- F[[1]] %*% t(fitExt$beta1)
curves1 <- curves10 + matrix(fitExt$beta0, nrow = 1)[rep(1, nrow(curves10)), ]
curves2 <- F[[2]] %*% t(fitExt$beta2)

CI1 <- t(apply(curves1, 1, quantile, probs = c(alpha/2, 1-alpha/2)) )
CI2 <- t(apply(curves2, 1, quantile, probs = c(alpha/2, 1-alpha/2)) )

# dataset ordered in the same way as design matrices
dataOrd <- data[ord, ]
dataOrd$lowPoint <- beta0Hat + F[[1]] %*% beta1Hat
dataOrd$lowLower <- CI1[, 1]
dataOrd$lowUpper <- CI1[, 2]
dataOrd$diffPoint <- F[[2]] %*% beta2Hat
dataOrd$diffLower <- CI2[, 1]
dataOrd$diffUpper <- CI2[, 2]

CIpoly <- data.frame(x = c(dataOrd$x, rev(dataOrd$x)), 
                     ci = c(dataOrd$lowLower, rev(dataOrd$lowUpper)),
                     id = 0)
ggplot(aes(x = x, y = lowPoint), data = dataOrd)+
  geom_polygon(aes(y = ci), data = CIpoly, fill = "grey")+
  geom_line(linetype = "solid", size = 1, color = "red")+
  # scale_y_continuous(lim = c(-3, 1))+
  theme_bw(28)+
  labs(y = expression(paste("lo",g[10], "(EDA)")),
       x = "Minute")+
       # title = "Low vigilance")+
  scale_y_continuous(lim = c(-3, 1))
ggsave(file.path(paperPath,"EDA_Bayes_CI_low.png"))

keep <- which(rowSums(F[[2]]) != 0)
CIpoly <- data.frame(x = c(dataOrd$x[keep], rev(dataOrd$x[keep])), 
                     ci = c(dataOrd$diffLower[keep], rev(dataOrd$diffUpper[keep])),
                     id = 0)
ggplot(aes(x = x, y = diffPoint), data = dataOrd[keep, ])+
  geom_polygon(aes(y = ci), data = CIpoly, fill = "grey")+
  geom_line(linetype = "solid", size = 1, color = "purple")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # scale_y_continuous(lim = c(-1, 1))+
  theme_bw(28)+
  labs(y = expression(paste("lo",g[10], "(EDA)")),
       x = "Minute")+
       # title = "High - low vigilance")+
  scale_y_continuous(lim = c(-1, 1))
ggsave(file.path(paperPath,"EDA_Bayes_CI_diff.png"))

# statistically significant differences at 5% level
range(dataOrd$x[dataOrd$diffLower > 0])
# [1] 42.39583 64.58750

# random curves ---------------------------------------------------------------
data$yHat <- Zcheck1 %*% btilde1Hat + Ztilde2 %*% btilde2Hat
data$yHat <- as.vector(Z %*% b)

ggplot(aes(x = x, y = y, group = id, color = factor(type)), data = data)+
  # geom_line(color = "grey")+
  geom_line(aes(y = yHat), linetype = "solid")+
  # geom_line(data = data, color = "black", size = 1)+
  # geom_line(aes(y = yMarg), data = data, color = "red", linetype = "solid", size = 1)+
  theme_bw(26)+
  labs(y = expression(paste("lo",g[10], "(EDA)")),
       x = "Minute",
       title = "Random curves")+
  scale_color_manual("Stress", values = c("blue", "red"))
ggsave(file.path(paperPath,"EDA_Bayes_rand.png"))

# predicted curves
data$yHat <- beta0Hat + as.vector(X[[1]]$F %*% beta1Hat + X[[2]]$F %*% beta2Hat + 
             Zcheck1 %*% btilde1Hat + Ztilde2 %*% btilde2Hat)

ggplot(aes(x = x, y = y, group = id), data = data)+
  geom_line(color = "grey")+
  geom_line(aes(y = yHat), color = "blue", linetype = "dashed")+
  theme_bw(30)+
  labs(y = expression(paste("lo",g[10], "(EDA)")),
       x = "Minute",
       title = "Subject-specific curves")
ggsave(file.path(paperPath,"EDA_Bayes_subj_curves.png"))

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
ggsave(file.path(paperPath, "EDA_Bayes_obs_pred.png"))
ggsave(file.path(presentPath, "EDA_Bayes_obs_pred.png"))
