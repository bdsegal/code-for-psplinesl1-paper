
library(reshape2)

if(length(grep("bdsegal",getwd())) > 0){
  computer <- "C:/Users/bdsegal"
} else{
  computer <- "/home/bsegal"
}
paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")
posterPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/poster/plots")

# plot for paper --------------------------------------------------------------
# must fit models a1 and m1 first in files 

# set model equal to a1 or m1
model <- a1
model <- m1

J <- length(model$params$X)
j <- 1
F <- model$params$X[[j]]$F
rho <- model$params$rho
lambda <- model$params$lambda[j]
u <- as.vector(model$params$u[[j]])
w <- as.vector(model$params$wNew[[j]])
D <- model$params$X[[j]]$D
y <- model$data$y
x <- model$data$x
b <- model$coefs$b
beta0 <- model$coefs$beta0

Hrho <- F %*% solve(crossprod(F) + rho * crossprod(D), t(F))
Hlambda <- F %*% solve(crossprod(F) + lambda * crossprod(D), t(F))

yb <- y - beta0 - rand$Z %*% b
for (l in (1:J)[-j]) {
  Fs <- model$params$X[[s]]$F
  betas <- model$coefs$beta[[s]]
  yb <- yb - Fs %*% betas
}

R <- as.vector(rho * solve(crossprod(F) + rho * crossprod(D), crossprod(D, w - u)))
Hy <- as.vector(Hlambda %*% yb)
HyR <- as.vector(Hrho %*% yb + F %*% R)
# Hr <- rho * F %*% solve(t(F) %*% F + rho * t(D) %*% D, t(D))

if (j == 1) {
  Hy <- Hy + beta0
  HyR <- HyR + beta0
}

dat <- data.frame(x = x, Hy = Hy, HyR = HyR, FR = F %*% R)
datM <- melt(dat, id.vars = "x")
levels(datM$variable) <- c("Hy", "Fbeta", "Fdelta")
datM$variable <- factor(datM$variable, levels = c("Fbeta", "Hy", "Fdelta"))

labs <- c(expression(paste("F", hat(beta), sep = "")),
          "Hy",
          expression(paste("F", delta, sep = "")))

dev.new(width = 7, height = 5)
ggplot(aes(x = x, y = value, color = variable, linetype = variable), data = datM)+
  geom_line(size = 1)+
  theme_bw(30)+
  geom_hline(yintercept = 0, color = "grey")+
  scale_color_discrete("", labels = labs)+
  scale_linetype_discrete("", labels = labs)+
  labs(y = "y")+
  guides(color = guide_legend(
                 keywidth=0.35,
                 keyheight=0.45,
                 default.unit="inch")
      )

# save plots from model = a1
ggsave(file.path(paperPath, "Hy_breakdown_sim.png"))
ggsave(file.path(presentPath, "Hy_breakdown_sim_present.png"))

# save plots from model = m1
ggsave(file.path(paperPath, "Hy_breakdown_EDA.png"))
ggsave(file.path(presentPath, "Hy_breakdown_EDA_present.png"))
