# cs208 HW 4a - Question 1
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw4a/")
require(plyr); require(dplyr); require(ggplot2)
test_data <- read.csv("hw4testdata.csv")
test_data_matrix <- as.matrix(test_data)
pums <- read.csv("CAPUMS5full.csv")

# a
## general ##
get_conjunction_matrix <- function(x, y, negative_val=0, positive_val=1){
  n <- nrow(x)
  d <- ncol(x)
  conjunction <- matrix(NA, nrow=n, ncol=d)
  for(j in 1:d) {
    conjunction[,j] <- as.integer(x[,j]==negative_val & y==positive_val)
  }
  return(conjunction)
}

calc_pj <- function(conjunction, epsilon, type) { 
  
  if(type=="c"){
    # centralized - add Laplace noise
    scale <- 1/(epsilon*length(conjunction))
    pj <- mean(conjunction) + rlap(mu=0, b=scale, size=1) 
  }
  else if(type=="l"){
    # localized - bias correction for RR
    pj <- correction(release=conjunction, epsilon=epsilon)
  }
  else{
    print("Type not supported.")
  }
  return(pj)
}

## centralized ##
sgn <- function(x) {     # function borrowed from class
  return(ifelse(x < 0, -1, 1))
}

rlap = function(mu=0, b=1, size=1) {     # function borrowed from class
  p <- runif(size) - 0.5
  draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
  return(draws)
}

SQcentralized <- function(conjunction, t, epsilon=1) { 
  d <- ncol(conjunction)
  attributes <- c()
  for(j in 1:d) {
    cur_pj <- calc_pj(conjunction=conjunction[,j], epsilon=epsilon/d, type="c")
    if(cur_pj < t) {
      attributes <- c(attributes, j)
    }
  }
  return(attributes) 
}

## localized ##
localRelease <- function(x, values=c(0,1), epsilon){     # function borrowed from class
  draw <- runif(n=1, min=0, max=1)
  cutoff <- 1/(1+exp(epsilon))
  if(draw < cutoff){
    return_val <- values[!values %in% x]
  }
  else{
    return_val <- x
  }
  return(return_val) 
}

correction <- function(release, epsilon){     # function updated from class
  n <- length(release)
  mulitiplicative <- (exp(epsilon) + 1)/(exp(epsilon) - 1)
  additive <- -n/(exp(epsilon) - 1)
  expectation <- (sum(release)*mulitiplicative + additive)/n
  return(expectation)
}

SQlocalized <- function(conjunction, t, epsilon=1, negative_val=0, positive_val=1) { 
  n <- nrow(conjunction)
  d <- ncol(conjunction)
  new_conjunction <- matrix(NA, nrow=n, ncol=d)
  attributes <- c()
  
  # RR 
  for(j in 1:d){
    for(i in 1:n){
      new_conjunction[i,j] <- localRelease(conjunction[i,j], values=c(0,1), epsilon=epsilon/d)
    }
  }
  
  # pj calculation (bias correction happens inside here)
  for(j in 1:d) {
    cur_pj <- calc_pj(conjunction=new_conjunction[,j], epsilon=epsilon/d, type="l")
    if(cur_pj < t) {
      attributes <- c(attributes, j)
    }
  }
  return(attributes)
}

## test ##
test_conj_mat <- get_conjunction_matrix(x=test_data[,1:10], y=test_data[,11])
SQcentralized(test_conj_mat, t=0.01, epsilon=1)
SQlocalized(conjunction=test_conj_mat, t=0.01, epsilon=1)



# c
pums_x <- pums[,c("sex","married","black","asian","collegedegree","employed","militaryservice",
                  "uscitizen","disability","englishability")]
pums_y <- pums$targetted
pums_conj_mat <- get_conjunction_matrix(pums_x, pums_y)

names(pums[SQcentralized(conjunction=pums_conj_mat, t=0.0001, epsilon=1)])
names(pums[SQlocalized(conjunction=pums_conj_mat, t=0.0001, epsilon=1)])



