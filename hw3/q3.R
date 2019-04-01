# cs208 HW 3 - Question 3
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw3/")
require(plyr); require(dplyr); require(ggplot2)
pums <- read.csv("MaPUMS5full.csv")
run_sims_flag <- FALSE
options(scipen=999)
set.seed(208)

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
    zgranularity <- (zupper-zlower)/znbins
    zbins[znbins+1] <-  zbins[znbins+1] + zgranularity
    zcodebook <- zbins[1:znbins] + 0.5*zgranularity
  }
  
  x.clipped <- clip(x=x, lower=xlower, upper=xupper)
  y.clipped <- clip(x=y, lower=ylower, upper=yupper)
  z.clipped <- clip(x=z, lower=zlower, upper=zupper)
  
  sensitivity <- 2
  scale <- sensitivity / (epsilon)
  
  sensitiveValue <- DPrelease <- matrix(NA, nrow=xnbins*ynbins*znbins, ncol=4)
  
  row_ind <- 1
  for(i in 1:xnbins){
    for(j in 1:ynbins){
      for(k in 1:znbins){
        bin_count <- sum(x.clipped >= xbins[i] & x.clipped < xbins[i+1] & y.clipped >= ybins[j] & y.clipped < ybins[j+1] & z.clipped >= zbins[k] & z.clipped < zbins[k+1])
        release_count <- bin_count + rlap(mu=0, b=scale, size=1)
        sensitiveValue[row_ind,] <- c(i,j,k,bin_count)
        DPrelease[row_ind,] <- c(i,j,k,release_count)
        row_ind = row_ind + 1
      }
    }
  }
  
  return(list(release=DPrelease, true=sensitiveValue, xcodebook=xcodebook, ycodebook=ycodebook, zcodebook=zcodebook))
}

# clip and log income
pums$log_income_clipped <- log(clip(x=pums$income, lower=1, upper=1000000))
log_income_low <- floor(log(1))
log_income_high <- ceiling(log(1000000))

# run histogram release
start = Sys.time()
res = xyzHistogramRelease(x=pums$educ, y=pums$age, z=pums$log_income_clipped, 
                          xlower=1, xupper=16, ylower=18, yupper=100, zlower=log_income_low, zupper=log_income_high, 
                          xnbins=0, ynbins=0, znbins=10, epsilon=1)
end = Sys.time()
end-start

# private results

synthetic_n <- 10000
synthetic_bin_probs <- normalize(res$release[,4])
synthetic_inds <- sample(x=1:nrow(res$release), size=synthetic_n, replace=TRUE, prob=synthetic_bin_probs)
synthetic_data <- res$release[synthetic_inds,1:3]
synthetic_data_df <- data.frame(synthetic_data)
names(synthetic_data_df) <- c("educ", "age", "log_income")

synthetic_data_df$educ <- plyr::mapvalues(synthetic_data_df$educ, 
                                          from=sort(unique(synthetic_data_df$educ)), 
                                          to=res$xcodebook)
synthetic_data_df$age <- plyr::mapvalues(synthetic_data_df$age, 
                                          from=sort(unique(synthetic_data_df$age)), 
                                          to=res$ycodebook)
synthetic_data_df$log_income <- plyr::mapvalues(synthetic_data_df$log_income, 
                                         from=sort(unique(synthetic_data_df$log_income)), 
                                         to=res$zcodebook)

synthetic_reg <- lm(log_income ~ age + educ, data = synthetic_data_df)
synthetic_slopes <- coef(synthetic_reg)[2:3]

# sensitive results
true_reg <- lm(log_income_clipped ~ age + educ, data = pums)
true_slopes <- coef(true_reg)[2:3]

# error
paste0("MSE for age coefficient: ", mse(synthetic_slopes[1], true_slopes[1]))
paste0("MSE for education coefficient: ", mse(synthetic_slopes[2], true_slopes[2]))

# simulations to examine contribution to MSE of bias and variance
n_sims <- 10
sim_history <- matrix(NA, nrow=n_sims, ncol=2)

if(run_sims_flag){
  for(i in 1:n_sims){
    print(i)
    res_sim = xyzHistogramRelease(x=pums$educ, y=pums$age, z=pums$log_income_clipped,
                                  xlower=1, xupper=16, ylower=18, yupper=100, zlower=log_income_low, zupper=log_income_high,
                                  xnbins=0, ynbins=0, znbins=10, epsilon=1)
    synthetic_bin_probs_sim <- normalize(res_sim$release[,4])
    synthetic_inds_sim <- sample(x=1:nrow(res_sim$release), size=synthetic_n, replace=TRUE, prob=synthetic_bin_probs_sim)
    synthetic_data_sim <- res_sim$release[synthetic_inds_sim,1:3]
    synthetic_data_sim_df <- data.frame(synthetic_data_sim)
    names(synthetic_data_sim_df) <- c("educ", "age", "log_income")
    
    synthetic_data_sim_df$educ <- plyr::mapvalues(synthetic_data_sim_df$educ,
                                                  from=sort(unique(synthetic_data_sim_df$educ)),
                                                  to=res_sim$xcodebook)
    synthetic_data_sim_df$age <- plyr::mapvalues(synthetic_data_sim_df$age,
                                                 from=sort(unique(synthetic_data_sim_df$age)),
                                                 to=res_sim$ycodebook)
    synthetic_data_sim_df$log_income <- plyr::mapvalues(synthetic_data_sim_df$log_income,
                                                    from=sort(unique(synthetic_data_sim_df$log_income)),
                                                    to=res_sim$zcodebook)
    
    synthetic_reg_sim <- lm(log_income ~ age + educ, data = synthetic_data_sim_df)
    synthetic_slopes_sim <- coef(synthetic_reg_sim)[2:3]
    
    sim_history[i,1] <- synthetic_slopes_sim[1]
    sim_history[i,2] <- synthetic_slopes_sim[2]
  }
  write.csv(sim_history, "sim_history.csv")
}

sim_history <- read.csv("sim_history.csv")


# bootstrap to examine sampling error
n_boots <- 1000
boot_size <- 1000
boot_history <- matrix(NA, nrow=n_boots, ncol=2)
for(i in 1:n_boots){
  boot_inds <- bootstrap(x=1:nrow(pums), n=boot_size)
  boot_data <- pums[boot_inds,]
  boot_reg <- lm(log_income_clipped ~ age + educ, data=boot_data)
  boot_slopes <- coef(boot_reg)[2:3]
  boot_history[i,1] <- boot_slopes[1]
  boot_history[i,2] <- boot_slopes[2]
}

# print results
mse_age <- mse(pred=sim_history[,2], true=true_slopes[1])
var_age <- var(sim_history[,2])
bias_sq_age <- mse_age - var_age
sampling_mse_age <- mse(boot_history[,1], true_slopes[1])

mse_educ <- mse(pred=sim_history[,3], true=true_slopes[2])
var_educ <- var(sim_history[,3])
bias_sq_educ <- mse_educ-var_educ
sampling_mse_educ <- mse(boot_history[,2], true_slopes[2])

cat("Age","\nDP MSE: ", mse_age, "\nDP Variance: ", var_age, "\nDP Bias^2: ", bias_sq_age, "\nSampling MSE: ", sampling_mse_age)
cat("Educ","\nDP MSE: ", mse_educ, "\nDP Variance: ", var_educ, "\nDP Bias^2: ", bias_sq_educ, "\nSampling MSE: ", sampling_mse_educ)

# save results
q3_results <- data.frame("Metric"=c("DP MSE","DP Variance","DP Bias^2","Sampling MSE"),
                         "Education"=c(mse_age,var_age,bias_sq_age,sampling_mse_age),
                         "Age"=c(mse_educ,var_educ,bias_sq_educ,sampling_mse_educ))
write.csv(q3_results, "q3_results.csv")

