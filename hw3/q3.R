# cs208 HW 3 - Question 3
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw3/")
require(plyr); require(dplyr); require(ggplot2)
pums <- read.csv("MaPUMS5full.csv")

# functions

sgn <- function(x) {     # borrowed from class
  return(ifelse(x < 0, -1, 1))
}

rlap = function(mu=0, b=1, size=1) {     # borrowed from class
  p <- runif(size) - 0.5
  draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
  return(draws)
}

clip <- function(x, lower, upper){     # borrowed from class
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

bootstrap <- function(x, y=NULL, n){     # borrowed from class
  index <- sample(x=1:length(x), size=n, replace=TRUE) 
  
  if(is.null(y)){
    return(x[index])
  }else{
    return(list(x=x[index], y=y[index]))
  }
}

normalize <- function(x){     # borrowed from class
  x[x<0] <- 0
  x <- x/sum(x)
  return(x)
}

mse <- function(pred, true) {
  mse = mean((pred-true)^2)
  return(mse)
}

xyzHistogramRelease <- function(x, y, z, xlower, xupper, ylower, yupper, zlower, zupper, xnbins=0, ynbins=0, znbins=0, epsilon){

  if(xnbins==0){
    xlower <- floor(xlower)
    xupper <- ceiling(xupper)
    xbins <- xlower:(xupper+1)    
    xnbins <- length(xbins)-1
    xgranularity <- 1
    xcodebook <- xbins[1:xnbins]
  } else {
    xbins <- seq(from=xlower, to=xupper, length=xnbins+1)
    xgranularity <- (xupper-xlower)/xnbins
    xbins[xnbins+1] <-  xbins[xnbins+1] + xgranularity
    xcodebook <- xbins[1:xnbins] + 0.5*xgranularity
  }
  
  if(ynbins==0){
    ylower <- floor(ylower)
    yupper <- ceiling(yupper)
    ybins <- ylower:(yupper+1)    
    ynbins <- length(ybins)-1
    ygranularity <- 1
    ycodebook <- ybins[1:ynbins]
  } else {
    ybins <- seq(from=ylower, to=yupper, length=ynbins+1)
    ygranularity <- (yupper-ylower)/ynbins
    ybins[ynbins+1] <-  ybins[ynbins+1] + ygranularity
    ycodebook <- ybins[1:ynbins] + 0.5*ygranularity
  }
  
  if(znbins==0){
    zlower <- floor(zlower)
    zupper <- ceiling(zupper)
    zbins <- zlower:(zupper+1)    
    znbins <- length(zbins)-1
    zgranularity <- 1
    zcodebook <- zbins[1:znbins]
  } else {
    zbins <- seq(from=zlower, to=zupper, length=znbins+1)
    zgranularity <- (zupper-ylower)/znbins
    zbins[znbins+1] <-  zbins[ynbins+1] + zgranularity
    zcodebook <- zbins[1:znbins] + 0.5*zgranularity
  }
  
  x.clipped <- clip(x=x, lower=xlower, upper=xupper)
  y.clipped <- clip(x=y, lower=ylower, upper=yupper)
  z.clipped <- clip(x=z, lower=zlower, upper=zupper)
  
  sensitivity <- 2
  scale <- sensitivity / (epsilon)
  
  sensitiveValue <- DPrelease <- matrix(NA, nrow=(xnbins*ynbins*znbins), ncol=4)
  
  row_ind <- 1
  for(i in 1:xnbins){
    for(j in 1:ynbins){
      for(k in 1:znbins){
        sensitiveValue[row_ind,1] <- DPrelease[row_ind,1] <- i
        sensitiveValue[row_ind,2] <- DPrelease[row_ind,2] <- j
        sensitiveValue[row_ind,3] <- DPrelease[row_ind,3] <- k
        sensitiveValue[row_ind,4] <- sum(x.clipped >= xbins[i] & x.clipped < xbins[i+1] & y.clipped >= ybins[j] & y.clipped < ybins[j+1] & z.clipped >= zbins[k] & z.clipped < zbins[k+1])
        DPrelease[row_ind,4] <- sensitiveValue[row_ind,4] + rlap(mu=0, b=scale, size=1)
        
        row_ind = row_ind + 1
      }
    }
  }
  
  sensitiveValue[is.na(sensitiveValue)] <- 0
  DPrelease[is.na(DPrelease)] <- 0
  
  return(list(release=DPrelease, true=sensitiveValue, xcodebook=xcodebook, ycodebook=ycodebook, zcodebook=zcodebook))
}

#
# feel like I should make this matrix not a 3d array 
# but matrix with educ, age, income columns you know?
#


# run histogram release
# is 10 bins for income too few? probably
res = xyzHistogramRelease(x=pums$educ, y=pums$age, z=pums$income, 
                          xlower=0, xupper=16, ylower=18, yupper=100, zlower=0, zupper=500000, 
                          xnbins=0, ynbins=0, znbins=10, epsilon=1)

# private results

# turn these bin sizes into probs
# sample from them
# create matrix that has educ, age, income values corresponding to bins that get sampled
# turn this into a df and run lm()
# that's essentially what James did in his code, but I need to figure out:
# 1. do i create some master synthetic data, so size nrow(pums)? or smaller synthetic data?
# I feel like create one master and bootstrap from it is way to go, because no extra privacy price
# 2. do i use rmultinom or can i just use sample w/ replacement like i did?

synthetic_n <- 10000
synthetic_bin_probs <- normalize(res$release[,4])
synthetic_inds <- sample(x=1:nrow(res$release), size=synthetic_n, replace=TRUE, prob=synthetic_bin_probs)
synthetic_data <- res$release[synthetic_inds,1:3]
synthetic_data_df <- data.frame(synthetic_data)
names(synthetic_data_df) <- c("educ", "age", "income")

synthetic_data_df$educ <- plyr::mapvalues(synthetic_data_df$educ, 
                                          from=sort(unique(synthetic_data_df$educ)), 
                                          to=res$xcodebook)
synthetic_data_df$age <- plyr::mapvalues(synthetic_data_df$age, 
                                          from=sort(unique(synthetic_data_df$age)), 
                                          to=res$ycodebook)
synthetic_data_df$income <- plyr::mapvalues(synthetic_data_df$income, 
                                         from=sort(unique(synthetic_data_df$income)), 
                                         to=res$zcodebook)

synthetic_reg <- lm(income ~ age + educ, data = synthetic_data_df)
synthetic_slopes <- coef(synthetic_reg)[2:3]

# sensitive results
true_reg <- lm(income ~ age + educ, data=pums)
true_slopes <- coef(true_reg)[2:3]

# error
paste0("MSE for age coefficient: ", mse(synthetic_slopes[1], true_slopes[1]))
paste0("MSE for education coefficient: ", mse(synthetic_slopes[2], true_slopes[2]))

# bootstrap
n_sims <- 1000




