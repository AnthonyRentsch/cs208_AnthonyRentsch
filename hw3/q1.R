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

trimmedMean <- function(x, d, epsilon) {
  n <- length(x)
  scale <- d/(epsilon*0.9*n)
  quants <- quantile(x, c(0.05,0.95))
  x_trimmed <- x[x>quants[1] && x<quants[2]]
  mean_trimmed <- (1/(0.9*n))*sum(x_trimmed)
  mean_release <- mean_trimmed + rlap(mu=0,b=scale,size=1)
  return(mean_release)
}

# c
exponentialPercentile <- function(x, t, epsilon){
  t <- t/100
  bins <- sort(x)
  nbins = length(bins)

  likelihoods <- rep(NA, nbins)
  for(i in 1:nbins){
    quality <- nbins - abs(nbins*t - i)
    bin_width <- ifelse((i+1) %in% 1:nbins, bins[i]-bins[i+1]+1, 1)
    likelihoods[i] <- (exp(epsilon * quality) / 2)*(bin_width)
  }
  
  probabilities <- likelihoods/sum(likelihoods)
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities) 
  
  bin_low_ind = which(flag)[1]
  bin_high_ind = bin_low_ind + 1
  DPrelease <- runif(n=1, min=bins[bin_low_ind], max=bins[bin_high_ind])
  
  return(DPrelease)
}

# d
e3trimmedMean <- function(x, epsilon, tile_low=5, tile_high=95) {
  
  n <- length(x)
  
  tile_low_value <- exponentialPercentile(x, t=tile_low, epsilon=epsilon/3)
  tile_high_value <- exponentialPercentile(x, t=tile_high, epsilon=epsilon/3)
  x_trimmed <- x[x>tile_low_value && x<tile_high_value]
  mean_trimmed <- (1/(0.9*n))*sum(x_trimmed)
  
  scale <- (3*(tile_high_value-tile_low_value))/(0.9*epsilon*n)
  mean_release <- mean_trimmed + rlap(mu=0,b=scale,size=1)
  return(mean_release)
}

# f
clip <- function(x, lower, upper){     # borrowed from class
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

rmse <- function(pred, true) {
  rmse <- sqrt(mean((pred-true)^2))
  return(rmse)
}

laplaceMeanRelease <- function(x, lower, upper, epsilon){     # borrowed from class
  n <- length(x)
  
  sensitivity <- (upper - lower)/n
  scale <- sensitivity / epsilon
  
  x.clipped <- clip(x, lower, upper)
  DPrelease <- mean(x.clipped) + rlap(mu=0, b=scale, size=1)
  
  return(DPrelease)
}

x = pums
pumas <- unique(x$puma)
num_pumas <- length(pumas)

means <- rep(NA, num_pumas)
lapmean_rels = rep(NA, num_pumas)
e3trimmedmean_rels = rep(NA, num_pumas)
rmses <- rep(NA, num_pumas)

for(i in 1:num_pumas){
  print(i)
  means[i] <- mean(x$income[x$puma==pumas[i]])
  lapmean_rels[i] <- laplaceMeanRelease(x=x$income[x$puma==pumas[i]],lower=0, upper=100000, epsilon=0.5)
  e3trimmedmean_rels[i] <- e3trimmedMean(x=x$income[x$puma==pumas[i]], epsilon=0.5)
  rmses[i] <- rmse(means[i], lapmean_rels[i])
}

plot(density(lapmean_rels), col="blue")
lines(density(e3trimmedmean_rels), col="red")
lines(density(means), col="green")

