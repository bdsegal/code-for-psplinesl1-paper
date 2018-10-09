library(psplinesl1)
library(reshape2)

paperPath <- "../../paper/plots"
presentPath <- "../../presentation/plots"
posterPath <- "../../poster/plots"

path <- "/home/bsegal/Documents/Research/data_empatica"

# data prep -------------------------------------------------------------------
data <- read.csv(file.path(path, "groupA.csv"))
data <- data[with(data, order(id, x)), ]

# setup design matrices
X <- list(ps(x = "x", norder = 4, k = 1, data = data, width = 5),
          ps(x = "x", norder = 4, k = 1, data = data, width = 5,
             by = "high", center = FALSE))
rand <- re(x = "x", id = "id", data = data, 
            randomCurves = TRUE,  width = 5, 
            norder = 4, derivOrder = 2)

system.time({
  m1 <- admm(y = "y", X = X, rand = rand,
             id = "id",
             tau = 450,
             lambda = c(1, 10),
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
 # 10.792   0.052  10.858

m1$conv$iter
# [1] 300

signif(m1$fit$df, 3)
#                  Overall    F1    F2   Z
# Stein                185  9.98  1.97 172
# Restricted           194 10.00  2.00 181
# ADMM                 193  9.00  2.00 181
# Ridge                196 19.50  8.00 167
# Ridge restricted     216 21.10 13.50 181

plot(log(m1$conv$rNorm))
plot(log(m1$conv$sNorm))

m1$coefs$beta0
# [1] -0.9897339

#sigma2b
m1$fit$sigma2 / m1$params$tau
# [1] 4.852319e-05

m1$fit$sigma2
# [1] 0.02183543

m1$fit$sigma2Ridge
# [1] 0.02198702

data$yHat <- m1$fit$yHat
data$randEff <- with(m1, as.vector(params$rand$Z %*% coefs$b))
data$bRem <- with(data, y - randEff)

# residuals
dev.new()
qplot(data$yHat, data$y - data$yHat)+
  theme_bw(22)+
  labs(x = expression(hat(y)), y = expression(paste("y - ",hat(y))))
ggsave(file.path(paperPath, "EDA_L1_resid_2.png"))

# confidence bands
CI <- ci(model = m1, alpha = 0.05)
plot(CI)

CIpoly <- data.frame(x = c(CI[[1]]$x, rev(CI[[1]]$x)), 
                     y = c(CI[[1]]$lower, rev(CI[[1]]$upper)),
                     id = 0)

CIpoint <- data.frame(x= CI[[1]]$x, y = CI[[1]]$smooth, id = 0)
            
ggplot(aes(x = x, y = y, group = id), data = CIpoly)+
  geom_polygon(fill = "grey")+
  geom_line(data = CIpoint, size = 1, color = "red")+
  theme_bw(28)+
  labs(x = "Minute", y = expression(paste("lo",g[10], "(EDA)", sep = "")))+
  scale_y_continuous(lim = c(-3,1))
ggsave(file.path(paperPath, "EDA_L1_smooth1_poster.png"))
ggsave(file.path(presentPath, "EDA_L1_smooth1_poster.png"))

CIpoly <- data.frame(x = c(CI[[2]]$x, rev(CI[[2]]$x)), 
                     y = c(CI[[2]]$lower, rev(CI[[2]]$upper)),
                     id = 0)
CIpoint <- data.frame(x= CI[[2]]$x, y = CI[[2]]$smooth, id = 0)
            
ggplot(aes(x = x, y = y, group = id), data = CIpoly)+
  geom_polygon(fill = "grey")+
  geom_line(data = CIpoint, size = 1, color = "purple")+
  theme_bw(28)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x = "Minute", y = expression(paste("lo",g[10], "(EDA)", sep = "")))+
  scale_y_continuous(lim = c(-1, 1))
ggsave(file.path(paperPath, "EDA_L1_smooth2_poster.png"))
ggsave(file.path(presentPath, "EDA_L1_smooth2_poster.png"))

# predicted curves
dataM <- melt(data, id.vars = c("x", "id", "type"),
                    measure.vars = c("y", "yHat"))
dataM$type <- factor(dataM$type, levels = c("low", "high"))
levels(dataM$type) <- c("Low vigilance", "High vigilance")

ID <- unique(dataM$id)
gg <- ggplot() + theme_bw(20) +
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
gg +
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
ggsave(file.path(paperPath, "EDA_L1_obs_pred.png"))
ggsave(file.path(presentPath, "EDA_L1_obs_pred.png"))

# MSE
with(data, mean((y - yHat)^2))
# [1] 0.01945408