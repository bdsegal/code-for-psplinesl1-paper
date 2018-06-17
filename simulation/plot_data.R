
# library(mgcv)
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
data(simData)
data(simData2groups)
data(trueMean)
data(trueInteraction)

range(table(simData$id))
# [1]  4 14

dev.new()
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(data = trueMean, size = 1)+
  theme_bw(24)
ggsave(file.path(paperPath, "nonDiff_smalln.png"))

dev.new()
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(data = trueMean, size = 1)+
  theme_bw(28)
ggsave(file.path(posterPath, "nonDiff_smalln_poster.png"))

dev.new(width = 7, height = 6)
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  # geom_line(data = trueMean, size = 1)+
  theme_bw(28)
ggsave(file.path(presentPath, "nonDiff_smalln_noMarg.png"))

dev.new(width = 7, height = 6)
ggplot(aes(x = x, y = y, group = id), data = simData)+
  geom_point(color = "grey")+
  geom_line(color = "grey")+
  geom_line(data = trueMean, size = 1)+
  theme_bw(28)
ggsave(file.path(presentPath, "nonDiff_smalln_Marg.png"))

truth <- data.frame(x = rep(trueMean$x, times = 2), 
                    y = c(trueMean$y, trueInteraction$y),
                    group = c(rep("beta[0] + f[1](x)", length(trueMean$x)),
                              rep("f[2](x)", length(trueMean$x))))

dev.new(width = 7, height = 3.5)
ggplot(aes(x = x, y = y), data = truth) +
  geom_line() +
  theme_bw(26) +
  facet_wrap(~group, labeller = label_parsed) +
  scale_x_continuous(breaks = c(0, 0.25, .5, .75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1))
ggsave(file.path(presentPath, "nonDiff_2groups_truth.png"))

simData2groups$group <- simData2groups$group + 1

dev.new(width = 6.75, height = 4)
ggplot(aes(x = x, y = y, group = id, color = as.factor(group),
           linetype = as.factor(group)), 
       data = simData2groups) +
  geom_line() +
  geom_point() +
  theme_bw(26) +
  scale_color_discrete("Group") +
  scale_linetype_discrete("Group")
ggsave(file.path(presentPath, "nonDiff_2groups.png"))
