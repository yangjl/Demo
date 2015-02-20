


getc <- function(tot = 10^7){
  a <- sample(1:10^9, tot)
  c <- sum(a %in% 1:1000)
  return(c)
}

n7 <- replicate(100, getc(tot=10^7))
n8 <- replicate(100, getc(tot=10^8))

res <- read.csv("~/Desktop/res.csv")

par(mfrow=c(2,1))
hist(res$n7, xlim=c(0, 20), main = expression(paste("Sampling ", 10^7, " Molecules")), xlab="" )
hist(res$n8, xlim=c(0, 200), main = expression(paste("Sampling ", 10^8, " Molecules")), xlab="" )

