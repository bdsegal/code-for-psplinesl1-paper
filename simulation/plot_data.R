
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
simData = read.csv("simData.csv")
trueMean = read.csv("trueMean.csv")

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
