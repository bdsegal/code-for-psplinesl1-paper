library(reshape2)
library(fda)

# change point   

xsim <- seq(0, 1, 0.01)
basis <- create.bspline.basis(rangeval = c(0, 1),
                                breaks = seq(0, 1, 0.075),
                                norder = 4)
Xsim <- eval.basis(xsim, basis)
Xsim <- cbind(1, Xsim)

yTrue <- 1 + (xsim > 0.2) * (xsim - 0.2)*2 - 
             (xsim > 0.4) * (xsim - 0.4)*4 +
             (xsim > 0.6) * (xsim - 0.6)*8 - 
             (xsim > 0.8) * (xsim - 0.8)*16 

plot(xsim, yTrue, lwd = 2, type = "l")

# note: see change_point_simulation.R for possible sd values
N <- 50
trueMean <- data.frame(x = xsim, y = yTrue, id = 0)
b <- rnorm(N, sd = 1)

Y <- matrix(rep(yTrue, N), ncol = N)
Y <- Y + matrix(rnorm(length(Y), sd = 0.1), ncol = N)
Y <- sweep(Y, 2, b, '+')
rownames(Y) <- xsim
colnames(Y) <- 1:N
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
simData <- simData[! 1:nrow(simData) %in% leave, ]

write.csv(simData, file = "simData.csv", row.names = FALSE)
write.csv(trueMean, file = "trueMean.csv", row.names = FALSE)
