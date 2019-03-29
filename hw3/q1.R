# cs208 HW 3 - Question 1
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw3/")
require(plyr); require(dplyr); require(ggplot2)
pums <- read.csv("MaPUMS5full.csv")

# a
sgn <- function(x) {     # function borrowed from class
  return(ifelse(x < 0, -1, 1))
}

rlap = function(mu=0, b=1, size=1) {     # function borrowed from class
  p <- runif(size) - 0.5
  draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
  return(draws)
}

trimmedMean <- function(x, d, n, epsilon) {
  scale <- d/(epsilon*0.9*n)
  quants <- quantile(x, c(0.05,0.95))
  x_trimmed <- x[x>quants[1] && x<quants[2]]
  mean_trimmed <- (1/(0.9*n))*sum(x_trimmed)
  mean_release <- mean_trimmed + rlap(mu=0,b=scale)
  return(mean_release)
}

# c
exponentialPercentile <- function(x, t, lower, upper, epsilon, nbins=0){
  t <- t/100
  
  if (nbins==0) {
    bins <- floor(lower):ceiling(upper)
    nbins <- length(bins)
  }
  else {
    bins <- seq(lower, upper, length.out=nbins)
  }
  
  quality <- rep(NA, nbins)
  for(i in 1:length(quality)){
    quality[i] <- 1 - abs(t/(1-t) - sum(x<bins[i])/sum(x>=bins[i]))
  }
  
  likelihoods <- exp(epsilon * quality) / 2
  probabilities <- likelihoods/sum(likelihoods)
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities) 
  
  bin_low = ifelse(sum(bins[flag==FALSE])==0, lower, max(bins[flag==FALSE]))
  bin_high = ifelse(sum(bins[flag==TRUE])==0, upper, min(bins[flag==TRUE]))
  DPrelease <- runif(n=1, min=bin_low, max=bin_high)
  
  return(list(dp=DPrelease, bins=bins, flag=flag, q=quality))
}

h = exponentialPercentile(pums$income,t=25,lower=0,upper=1000000,epsilon=0.5,nbins=1000)

save <- rep(NA, 50)
for (i in 1:length(save)){
  save[i] <- exponentialPercentile(pums$income,t=5,lower=0,upper=1000000,epsilon=0.5,nbins=100)
}
hist(save)

# d


