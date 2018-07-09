# Plotting EDA measurements from E3 Wristband

library(ggplot2)

#TODO: Get times of tags
path <- "/home/bsegal/Documents/Research/data_empatica"
dataPath <- "/home/bsegal/Documents/Research/data_empatica/data_complete/unzipped/"
paperPath <- "/home/bsegal/Dropbox/Research/psplines_L1_penalty/paper/plots"
presentPath <- "/home/bsegal/Dropbox/Research/psplines_L1_penalty/presentation/plots"

marDefault <- c(5, 4, 4, 2) + 0.1 

# functions for processing raw data
getDataFun <- function(path, outcome, subject){

  # read data from csv file
  raw <- read.csv(paste(path, subject, "/", outcome, ".csv", sep = ""), header = FALSE)

  # get initial time stamp and sampling rate
  initialTime <- raw[1, 1]
  Hz <- raw[2, 1]

  # remove time stamp and sampling rate from the EDA data
  raw <- raw[-(1:2), ]
  
  if (is.vector(raw)){
    len <- length(raw)
  } else {
    len <- nrow(raw)
  }
  
  # create vector of seconds (in 1/Hz increments)
  sec <- (0:(len - 1)) / Hz

  # convert initial time stamp to time of day
  # t1 <- as.POSIXct(initialTime, origin = "1970-01-01")

  # add initial time to vector of seconds get a vector of time of day
  # time <- sec + t1

  # put data in a data frame for ggplot
  data <- data.frame(raw = raw, sec = sec)

  tagsNotRead <- TRUE
  try({
  tags <- read.csv(paste(path, subject, "/tags.csv", sep = ""), header = FALSE)
  tagsNotRead <- FALSE
  tags$sec <- tags[, 1] - initialTime
  })
  
  if(tagsNotRead) {
    tags <- NULL
  }

  return(list(data = data, tags = tags))
}

groupVig <- read.csv(file.path(path, "groupAndVigilance.csv"))
groupVig$idpart <- as.factor(groupVig$idpart)
groups <- c("A", "B", "C")
vigLevel <- c("Low", "High")

subjects <- list.files(dataPath)

# EDA -------------------------------------------------------------------------
EDAlist <- list()
for (i in 1:length(subjects)) {
  print(i)
  EDAlist[[i]] <- getDataFun(path = dataPath, subject = subjects[i], outcome = "EDA")
}
names(EDAlist) <- as.factor(gsub("ID", "", subjects))

# original scale
png(file.path(paperPath, "EDA_by_group.png"))
par(mfrow = c(3, 2), mar = marDefault + c(0, 1, 0, 0))
for (g in groups) {
  for (v in vigLevel) {
    sub <- groupVig[with(groupVig, which(group == g & vig == v)), ]
    subj <- which(names(EDAlist) %in% sub$idpart)
    plot(x = EDAlist[[subj[1]]]$data$sec / 60, y = EDAlist[[subj[1]]]$data$raw, type = "l",
      ylim = c(0, 7), xlim = c(0, 10000 / 60),
      main = paste(g, " ", v, ", ", length(subj), " subjects", sep = ""), 
      xlab = "Minute", ylab = "EDA",
      cex.lab = 1.9, cex.axis = 1.9, cex.main = 1.9)
    # abline(v = EDAlist[[subj[1]]]$tags$sec, lty = 3)
    for (i in 2:length(subj)) {
      points(x = EDAlist[[subj[i]]]$data$sec / 60, y = EDAlist[[subj[i]]]$data$raw,
        type = "l", cex.lab = 1.9, cex.axis = 1.9, col = i)
    }
  }
}
dev.off()

# original scale, just group A

png(filename = file.path(paperPath, "EDA_by_group_A.png"), height = 400, width = 800)
par(mfrow = c(1, 2), mar = marDefault + c(0, 1, 0, 0))
for (g in c("A")) {
  for (v in vigLevel) {
    sub <- groupVig[with(groupVig, which(group == g & vig == v)), ]
    subj <- which(names(EDAlist) %in% sub$idpart)
    plot(x = EDAlist[[subj[1]]]$data$sec / 60, y = EDAlist[[subj[1]]]$data$raw, type = "l",
      ylim = c(0, 7), xlim = c(0, 10000 / 60),
      main = paste(v, " vigilance", sep = ""), 
      xlab = "Minute", ylab = "EDA",
      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    # abline(v = EDAlist[[subj[1]]]$tags$sec, lty = 3)
    for (i in 2:length(subj)) {
      points(x = EDAlist[[subj[i]]]$data$sec / 60, y = EDAlist[[subj[i]]]$data$raw,
        type = "l", cex.lab = 1.5, cex.axis = 1.5, col = i)
    }
  }
}
dev.off()

# original scale with time marks
png(file.path(paperPath, "EDA_by_group_time_marks.png"))
par(mfrow = c(3, 2), mar = marDefault + c(0, 1, 0, 0))
for (g in groups) {
  for (v in vigLevel) {
    sub <- groupVig[with(groupVig, which(group == g & vig == v)), ]
    subj <- which(names(EDAlist) %in% sub$idpart)
    plot(x = EDAlist[[subj[1]]]$data$sec / 60, y = EDAlist[[subj[1]]]$data$raw, type = "l",
      ylim = c(0, 7), xlim = c(0, 10000 / 60),
      main = paste(g, " ", v, ", ", length(subj), " subjects", sep = ""), 
      xlab = "Minute", ylab = "EDA",
      cex.lab = 1.9, cex.axis = 1.9, cex.main = 1.9)
    # abline(v = EDAlist[[subj[1]]]$tags$sec, lty = 3)
    for (i in 2:length(subj)) {
      points(x = EDAlist[[subj[i]]]$data$sec / 60, y = EDAlist[[subj[i]]]$data$raw,
        type = "l", cex.lab = 1.9, cex.axis = 1.9, col = i)
    }
    for (i in 1:length(subj)) {
      abline(v = EDAlist[[subj[i]]]$tags$sec / 60, col = i, lty = 3)
    }
  }
}
dev.off()

# original scale - pickout reliable measurements
g <- groups[1]
v <- vigLevel[1]
par(mar = marDefault + c(0, 1, 0, 0))
sub <- groupVig[with(groupVig, which(group == g & vig == v)), ]
subj <- which(names(EDAlist) %in% sub$idpart)
for (i in 1:length(subj)) {
  plot(x = EDAlist[[subj[i]]]$data$sec / 60, y = EDAlist[[subj[i]]]$data$raw, type = "l",
       ylim = c(0, 7), xlim = c(0, 10000 / 60),
       main = names(EDAlist)[subj[i]],
       xlab = "Minute", ylab = "EDA",
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  print(names(EDAlist)[subj[i]])
  readline()
}

# A low
aLowID <- c(1587013955, 5391983104, 6016633856, 6781066752, 
          7803223040, 8471556096, 9126784000)

# maybe...
# 6781066752
# 8305564672

# original scale - pickout reliable measurements
g <- groups[1]
v <- vigLevel[2]
par(mar = marDefault + c(0, 1, 0, 0))
sub <- groupVig[with(groupVig, which(group == g & vig == v)), ]
subj <- which(names(EDAlist) %in% sub$idpart)
for (i in 1:length(subj)) {
  plot(x = EDAlist[[subj[i]]]$data$sec / 60, y = EDAlist[[subj[i]]]$data$raw, type = "l",
    ylim = c(0, 7), xlim = c(0, 10000 / 60),
    main = names(EDAlist)[subj[i]],
    xlab = "Minute", ylab = "EDA",
    cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  print(names(EDAlist)[subj[i]])
  readline()
}

# A high
aHighID <- c(1879278988, 5448769536, 5553606144, 5723965440, 
             6650021376, 7672177664, 8379824128, 8414769664,
             8641914880, 9279671296)

# maybe...
# 6484029952
# 7702754816
# 9279671296

# plot reliable measurements
aLowSubj <- which(names(EDAlist) %in% aLowID)
aHighSubj <- which(names(EDAlist) %in% aHighID)

# plot usuable measurements for group A
png(file.path(paperPath, "groupA_useable.png"))
par(mar = marDefault + c(0, 1, 0, 0))
plot(x = EDAlist[[aLowSubj[1]]]$data$sec / 60,
     y = EDAlist[[aLowSubj[1]]]$data$raw,
     type = "l",
     ylim = c(0, 8), xlim = c(0, 9000 / 60),
     main = paste("Usable measurements in group A\n", length(aHighSubj), "high and", 
                  length(aLowSubj), "low vigilance subjects"), 
     xlab = "Minute", ylab = "EDA",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, col = "red")
    for (i in 2:length(aLowSubj)) {
      points(x = EDAlist[[aLowSubj[i]]]$data$sec / 60,
             y = EDAlist[[aLowSubj[i]]]$data$raw,
             type = "l", cex.lab = 1.5, cex.axis = 1.5, col = "red")
    }
    for (i in 1:length(aHighSubj)) {
      points(x = EDAlist[[aHighSubj[i]]]$data$sec / 60,
             y = EDAlist[[aHighSubj[i]]]$data$raw,
             type = "l", cex.lab = 1.5, cex.axis = 1.5, col = "blue")
    }
legend("topright", inset = 0.02, legend = c("High", "Low"), title = "Vigilance",
      lwd = 2, col = c("blue", "red"), cex = 1.5)
dev.off()

aLowData <- NULL
for (i in 1:length(aLowSubj)) {
  aLowData <- rbind(aLowData, cbind(EDAlist[[aLowSubj[i]]]$data,
                                    id = names(EDAlist)[aLowSubj[i]]))
}

lowThinList <- list()
for (i in unique(aLowData$id)) {
  keep <- which(aLowData$id == i)
  smooth <- ksmooth(x = aLowData$sec[keep], y = aLowData$raw[keep], bandwidth = 50)

  # to do: thin data proportional to curvature
  # p <- abs(diff(ySmooth$y)) / sum(abs(diff(ySmooth$y)))
  thin <- seq(1, length(smooth$x), length = 100)
  lowThinList[[i]] <- data.frame(x = smooth$x[thin], y = smooth$y[thin], id = i)
  # plot(aLowData$sec[keep], aLowData$raw[keep], type = "l", col = "grey")
  # lines(smooth, type = "l")
  # lines(smooth$x[thin], smooth$y[thin], type = "l", col = "red")
  # legend("topright", legend = c("Raw", "Smoothed", "Thinned"),
  #        col = c("grey", "black", "red"), lwd = 2, inset = 0.02)
  # readline()
}

# plot one for demonstration in paper
i <- unique(aLowData$id)[1]
keep <- which(aLowData$id == i)
smooth <- ksmooth(x = aLowData$sec[keep] / 60, y = aLowData$raw[keep], bandwidth = 50/60)
thin <- seq(1, length(smooth$x), length = 100)

png(filename = file.path(paperPath, "thinned.png"))
par(mar = marDefault + c(0, 1, 0, 0))
plot(aLowData$sec[keep] / 60, aLowData$raw[keep], type = "l", col = "grey60",
     ylab = "EDA", xlab = "Minute",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(smooth, type = "l")
lines(smooth$x[thin], smooth$y[thin], type = "l", col = "red")
legend("topright", legend = c("Raw", "Smoothed", "Thinned"),
       col = c("grey60", "black", "red"), lwd = 2, inset = 0.02, cex = 1.5)
dev.off()

aLowThin <- do.call(rbind, lowThinList)
aLowThin$x <- aLowThin$x / 60
aLowThin$y <- log(aLowThin$y + 0.001, 10)

# A group, high vigilance 
aHighData <- NULL
for (i in 1:length(aHighSubj)) {
  aHighData <- rbind(aHighData, cbind(EDAlist[[aHighSubj[i]]]$data,
                                    id = names(EDAlist)[aHighSubj[i]]))
}

highThinList <- list()
for (i in unique(aHighData$id)) {
  keep <- which(aHighData$id == i)
  smooth <- ksmooth(x = aHighData$sec[keep], y = aHighData$raw[keep], bandwidth = 50)
  
  # to do: thin data proportional to curvature
  # p <- abs(diff(ySmooth$y)) / sum(abs(diff(ySmooth$y)))
  thin <- seq(1, length(smooth$x), length = 100)
  highThinList[[i]] <- data.frame(x = smooth$x[thin], y = smooth$y[thin], id = i)
  # plot(aHighData$sec[keep], aHighData$raw[keep], type = "l", col = "grey")
  # lines(smooth, type = "l")
  # lines(smooth$x[thin], smooth$y[thin], type = "l", col = "red")
  # legend("topright", legend = c("Raw", "Smoothed", "Thinned"),
  #        col = c("grey", "black", "red"), lwd = 2, inset = 0.02)
  # readline()
}

aHighThin <- do.call(rbind, highThinList)
aHighThin$x <- aHighThin$x / 60
aHighThin$y <- log(aHighThin$y + 0.001, 10)

# low and high groups together
groupA <- rbind(cbind(aLowThin, high = 0, type = "low"),
      cbind(aHighThin, high = 1, type = "high"))


# write to file
# write.csv(groupA, file.path(path, "groupA.csv"), row.names = FALSE)

groupA <- read.csv("groupA.csv")

groupA$type <- factor(groupA$type, levels = c("low", "high"))

levels(groupA$type) <- c("Low vigilance", "High vigilance")

# plot processed data
ggplot(aes(x = x, y = y, group = id, color = as.factor(type)), data = groupA)+
  geom_line()+
  theme_bw(20)+
  labs(x = "Minute", y = expression(paste("lo",g[10],"(EDA + 0.001)")))+
  scale_color_discrete("Vigilance", labels = c("Low", "High"))
ggsave(file.path(paperPath, "groupA_processed.png"))

# plot processed data

dev.new(width = 8, height = 5)
ggplot(aes(x = x, y = y, group = id, color = as.factor(type)), data = groupA)+
  geom_line()+
  theme_bw(18)+
  facet_wrap(~ type) +
  labs(x = "Minute", y = expression(paste("lo",g[10],"(EDA)")))+
  scale_color_manual(guide = FALSE, values = c("red", "blue"))
ggsave(file.path(paperPath, "eda_data_split.png"))

groupA <- read.csv(file.path(path, "groupA.csv"))
library(dplyr)
groupA %>% group_by(type, id) %>%
  summarize(n = n())