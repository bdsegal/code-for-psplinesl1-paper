
library(reshape2)
library(ggplot2)

if(length(grep("bdsegal",getwd())) > 0){
  computer <- "C:/Users/bdsegal"
} else{
  computer <- "/home/bsegal"
}

paperPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/paper/plots")
presentPath <- file.path(computer, "Dropbox/Research/psplines_L1_penalty/presentation/plots")

getInflect <- function(cut, dat, resFlat) {
  # function for getting inflection point statistics for different cutoff values
  numPeaks <- matrix(nrow = length(resFlat), ncol = 2)
  colnames(numPeaks) <- c(expression(paste(L[1], " penalty", sep = "")), 
                          expression(paste(L[2], " penalty", sep = "")))

  absDevTau <- matrix(nrow = length(resFlat), ncol = 2)
  colnames(absDevTau) <- c(expression(paste(L[1], " penalty", sep = "")), 
                           expression(paste(L[2], " penalty", sep = "")))

  meanAbsDevTauHat <- matrix(nrow = length(resFlat), ncol = 2)
  colnames(meanAbsDevTauHat) <- c(expression(paste(L[1], " penalty", sep = "")), 
                                  expression(paste(L[2], " penalty", sep = "")))


  absDevTauHat <- matrix(nrow = length(resFlat), ncol = 2)
  colnames(absDevTauHat) <- c(expression(paste(L[1], " penalty", sep = "")), 
                              expression(paste(L[2], " penalty", sep = "")))

  for(m in 1:1000) {

    dat <- resFlat[[m]]
    dx <- diff(dat$x)

    # L1 derivatives
    dyL1 <- diff(dat$yHatL1, differences = 1)
    d1L1 <- dyL1 / dx
    d2L1 <- diff(d1L1) / dx[-1]

    # L2 derivatives
    dyL2 <- diff(dat$yHatL2, differences = 1)
    d1L2 <- dyL2 / dx
    d2L2 <- diff(d1L2) / dx[-1]

    phi1 <- max(abs(d2L1))
    phi2 <- max(abs(d2L2))
    cutoff1 <- cut * phi1
    cutoff2 <- cut * phi2

    numPeaks[m, ] <- c(sum(abs(d2L1) >= cutoff1),
                       sum(abs(d2L2) >= cutoff2)
                      )

    absDevTau[m, ] <- c(sum(sapply(inflections, function(tau){min(abs(tau - dat$x[which(abs(d2L1) >= cutoff1)]))})),
                        sum(sapply(inflections, function(tau){min(abs(tau - dat$x[which(abs(d2L2) >= cutoff2)]))}))
                       )

    meanAbsDevTauHat[m, ] <- c(mean(sapply(dat$x[which(abs(d2L1) >= cutoff1)], function(tau){min(abs(tau - inflections))})),
                               mean(sapply(dat$x[which(abs(d2L2) >= cutoff2)], function(tau){min(abs(tau - inflections))}))
                              )

    absDevTauHat[m, ] <- c(sum(sapply(dat$x[which(abs(d2L1) >= cutoff1)], function(tau){min(abs(tau - inflections))})),
                           sum(sapply(dat$x[which(abs(d2L2) >= cutoff2)], function(tau){min(abs(tau - inflections))}))
                          )
  }

  return(list(numPeaks = numPeaks,
              absDevTau = absDevTau,
              meanAbsDevTauHat = meanAbsDevTauHat,
              absDevTauHat = absDevTauHat))
}

# true inflection points
inflections <- c(0.2, 0.4, 0.6, 0.8)

# get data
res <- list()
for (m in 1:10) {
 load(paste("batch/results", m, ".RData", sep = ""))
 res[[m]] <- results
 rm(results)
}

# flatten top level of res list
resFlat <- unlist(res, recursive = FALSE)

# get inflection point statistics for difference cutoff values
out <- getInflect(cut = 0.01, dat = data, resFlat = resFlat)
numPeaksM01 <- melt(out$numPeaks)
absDevTauHatM01 <- melt(out$absDevTauHat)
meanAbsDevTauHatM01 <- melt(out$meanAbsDevTauHat)

out <- getInflect(cut = 0.05, dat = data, resFlat = resFlat)
numPeaksM05 <- melt(out$numPeaks)
absDevTauHatM05 <- melt(out$absDevTauHat)
meanAbsDevTauHatM05 <- melt(out$meanAbsDevTauHat)

out <- getInflect(cut = 0.1, dat = data, resFlat = resFlat)
numPeaksM1 <- melt(out$numPeaks)
absDevTauHatM1 <- melt(out$absDevTauHat)
meanAbsDevTauHatM1 <- melt(out$meanAbsDevTauHat)

numPeaksM01$cut <- absDevTauHatM01$cut <- meanAbsDevTauHatM01$cut <- 0.01
numPeaksM05$cut <- absDevTauHatM05$cut <- meanAbsDevTauHatM05$cut<- 0.05
numPeaksM1$cut <- absDevTauHatM1$cut <- meanAbsDevTauHatM1$cut <- 0.1

# plot
numPeaksM <- rbind(numPeaksM01, numPeaksM05, numPeaksM1)#, numPeaksM25)#, numPeaksM5)
numPeaksM$cut <- factor(numPeaksM$cut)
levels(numPeaksM$cut) <- c(expression(paste("c = 0.01", sep = "")), 
                           expression(paste("c = 0.05", sep = "")),
                           expression(paste("c = 0.1", sep = "")))

dev.new(width = 8, height = 6)
ggplot(aes(x = value), data = numPeaksM)+
  geom_histogram(binwidth = 1)+
  facet_grid(Var2 ~ cut, labeller = label_parsed)+
  theme_bw(28)+
  labs(x = expression(n[I]), y = "Count")+
  # labs(x = "", y = "Count")+
  geom_vline(xintercept = 4, size = 1, linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
ggsave(file.path(paperPath, "numPeaks_multi_c.png"))

dev.new(width = 7.5, height = 6.5)
ggplot(aes(x = value), data = numPeaksM)+
  geom_histogram(binwidth = 1)+
  facet_grid(Var2 ~ cut, labeller = label_parsed)+
  theme_bw(30)+
  labs(x = "Number of estimated inflection points", y = "Count")+
  # labs(x = "", y = "Count")+
  geom_vline(xintercept = 4, size = 1, linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
ggsave(file.path(presentPath, "numPeaks_multi_c_poster.png"))

absDevTauM <- melt(absDevTau)
dev.new(width = 7, height = 5)
ggplot(aes(x = value), data = absDevTauM)+
  geom_histogram()+
  facet_wrap(~Var2, labeller = label_parsed)+
  theme_bw(18)+
  labs(x = "Absolute deviance", y = "Count", title = "Summed over true inflection points")
ggsave("abs_dev_tau.png")

meanAbsDevTauHatM <- rbind(meanAbsDevTauHatM01, 
                           meanAbsDevTauHatM05, 
                           meanAbsDevTauHatM1)

meanAbsDevTauHatM$cut <- factor(meanAbsDevTauHatM$cut)
levels(meanAbsDevTauHatM$cut) <- c(expression(paste("c = 0.01", sep = "")), 
                                   expression(paste("c = 0.05", sep = "")),
                                   expression(paste("c = 0.1", sep = "")))

dev.new(width = 8, height = 6)
ggplot(aes(x = value), data = meanAbsDevTauHatM)+
  geom_histogram()+
  facet_grid(Var2 ~ cut, labeller = label_parsed)+
  theme_bw(28)+
  labs(x = expression(bar(d)), y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
ggsave(file.path(paperPath, "mean_abs_dev_tauHat.png"))

dev.new(width = 7.5, height = 6.5)
ggplot(aes(x = value), data = meanAbsDevTauHatM)+
  geom_histogram()+
  facet_grid(Var2 ~ cut, labeller = label_parsed)+
  theme_bw(30)+
  labs(x = "Mean absolute deviance", y = "Count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1))
ggsave(file.path(presentPath, "mean_abs_dev_tauHat_poster.png"))
