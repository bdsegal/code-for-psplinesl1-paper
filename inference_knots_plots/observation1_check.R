# code for confirming Observation 1
bspline <- function(x, i, m, tau) {
  n <- length(x)
  if (m>1) {
    aDen <- (tau[i+m-1] - tau[i])
    bDen <- (tau[i+m] - tau[i+1])

    if(aDen == 0) {
      a <- rep(0, n)
    } else{
      a <- (x - tau[i]) / aDen
    }

    if(bDen == 0) {
      b <- rep(0, n)
    } else{
      b <- (tau[i+m] - x) / bDen
    }

    a * bspline(x, i, m-1, tau) +  b * bspline(x, i+1, m-1, tau)

  } else {
    ((tau[i] <= x) & (x < tau[i+1]))*1
  }
}

M <- 2
x <- seq(0, 1, .05)

tau <- c(rev(seq(0, by = -0.1, length.out = M)),
         x[-c(1,length(x))],
         seq(1, by = 0.1, length.out = M)
         )

c <- length(tau) - 2*M
p <- M + c

n <- length(x)
F <- matrix(nrow = n, ncol = p)
for (j in 1:p) {
  for (i in 1:n) {
    F[i,j] <- bspline(x=x[i], i=j, m=M, tau=tau)
  }
}

matplot(x = x, F,type = "l")

sum(abs(F-diag(n)))
