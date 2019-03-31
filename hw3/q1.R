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
unique_x <- 
nunique_x <- length(unique_x)
bins <- unique_x 
nbins <- nunique_x


exponentialPercentile <- function(x, t, epsilon){
  t <- t/100
  n <- length(x)
  bins <- sort(x)
  #nbins = n-1 ?

  quality <- rep(NA, n)
  for(i in 1:length(quality)){
    quality[i] <- n - abs(n*t - i)
  }
  
  likelihoods <- exp(epsilon * quality) / 2
  weighted_likelihoods <- rep(NA, length(likelihoods))
  for (i in 1:length(likelihoods)){
    weighted_likelihoods[i] <- weighted_likelihoods[i]*(bins[i]-bins[i+1]+1)
  }
  probabilities <- weighted_likelihoods/sum(weighted_likelihoods)
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities) 
  
  bin_low = max(bins[flag==FALSE])
  bin_high = min(bins[flag==TRUE])
  DPrelease <- sample(x=bin_low:bin_high, size=1)
  
  return(list(dp=DPrelease, bins=bins, flag=flag, q=quality))
}

h = exponentialPercentile(pums$income, t=25, epsilon=0.5)


# d


