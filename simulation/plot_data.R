library(psplinesl1)

paperPath <- "../../paper/plots"
presentPath <- "../../presentation/plots"
posterPath <- "../../poster/plots"

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

# for paper
dev.new(width = 7, height = 3.5)
ggplot(aes(x = x, y = y), data = truth) +
  geom_line() +
  theme_bw(20) +
  facet_wrap(~group, labeller = label_parsed) +
  scale_x_continuous(breaks = c(0, 0.25, .5, .75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1))
ggsave(file.path(paperPath, "nonDiff_2groups_truth.png"))

simData2groups$group <- paste0("Group ", simData2groups$group)

dev.new(width = 9, height = 4)
ggplot(aes(x = x, y = y, group = id, color = as.factor(group)), 
       data = simData2groups) +
  geom_line() +
  geom_point() +
  theme_bw(20) +
  facet_wrap(~ group) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(guide = FALSE, values = c("red", "blue"))
ggsave(file.path(paperPath, "nonDiff_2groups_facet.png"))