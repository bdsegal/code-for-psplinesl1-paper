library(reshape2)
library(fda)
library(ggplot2)

# change point   

xsim <- seq(0, 1, 0.01)
yTrue <- 1 + (xsim > 0.2) * (xsim - 0.2)*2 - 
             (xsim > 0.4) * (xsim - 0.4)*4 +
             (xsim > 0.6) * (xsim - 0.6)*8 - 
             (xsim > 0.8) * (xsim - 0.8)*16 

interaction <- (xsim > 0.2) * (xsim - 0.2)*2 - 
               (xsim > 0.4) * (xsim - 0.4)*4 +
               (xsim > 0.6) * (xsim - 0.6)*2

n1 <- 50
trueMean <- data.frame(x = xsim, y = yTrue, id = 0, group = 0)
trueInteraction <- data.frame(x = xsim, y = interaction, id = 0, group = 0)

b1 <- rnorm(n1, sd = 1)

Y <- matrix(rep(yTrue, n1), ncol = n1)
Y <- Y + matrix(rnorm(length(Y), sd = 0.1), ncol = n1)
Y <- sweep(Y, 2, b1, '+')
rownames(Y) <- xsim
colnames(Y) <- 1:n1
dat <- melt(Y, varnames = c("x", "id"), value.name = "y")

ids <- unique(dat$id)
ni <- length(xsim)
run <- 20
keepList <- NULL
for (i in 1:length(ids)) {
  begin <- sample(size = 1, 1:ni)
  if (begin < ni / 2) {
    end <- min(ni, begin + run)
  } else {
    end <- max(0, begin - run)
  }
  keepList[[i]] <- (ni*(i-1) + begin):(ni*(i-1) + end)
}
keep <- do.call(c, keepList)
simData <- dat[keep, ]

# thin the data
leave <- sample(size = nrow(simData) / 1.75, 1:nrow(simData))
simData1 <- simData[! 1:nrow(simData) %in% leave, ]

# note: see change_point_simulation.R for possible sd values
n2 <- 50
b2 <- rnorm(n2, sd = 1)

Y <- matrix(rep(yTrue + interaction, n2), ncol = n2)
Y <- Y + matrix(rnorm(length(Y), sd = 0.1), ncol = n2)
Y <- sweep(Y, 2, b2, '+')
rownames(Y) <- xsim
colnames(Y) <- (n1 + 1):(n1 + n2)
dat <- melt(Y, varnames = c("x", "id"), value.name = "y")

ids <- unique(dat$id)
ni <- length(xsim)
run <- 20
keepList <- NULL
for (i in 1:length(ids)) {
  begin <- sample(size = 1, 1:ni)
  if (begin < ni / 2) {
    end <- min(ni, begin + run)
  } else {
    end <- max(0, begin - run)
  }
  keepList[[i]] <- (ni*(i-1) + begin):(ni*(i-1) + end)
}
keep <- do.call(c, keepList)
simData <- dat[keep, ]

# thin the data
leave <- sample(size = nrow(simData) / 1.75, 1:nrow(simData))
simData2 <- simData[! 1:nrow(simData) %in% leave, ]

simData1$group = 0
simData2$group = 1
simData2groups <- rbind(simData1, simData2)

write.csv(simData2groups, file = "simData2groups.csv", row.names = FALSE)
write.csv(trueInteraction, file = "trueInteraction.csv", row.names = FALSE)
