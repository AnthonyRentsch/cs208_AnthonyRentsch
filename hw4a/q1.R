# cs208 HW 4a - Question 1
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw4a/")
require(plyr); require(dplyr); require(ggplot2)
test_data <- read.csv("hw4testdata.csv")
test_data_matrix <- as.matrix(test_data)
pums <- read.csv("CAPUMS5full.csv")
pums$sex <- 1-pums$sex # recode sex so that we have 0=male, 1=female

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

names(pums[SQcentralized(conjunction=pums_conj_mat, t=0.01, epsilon=1)])
names(pums[SQlocalized(conjunction=pums_conj_mat, t=0.01, epsilon=1)])

bootstrap <- function(x, n){     # function updated from class
  index <- sample(x=1:nrow(x), size=n, replace=TRUE) 
  return(x[index,])
}

return_rates <- function(predictors, y_true){
  d <- ncol(predictors)
  if(d == 0){
    y_pred <- rep(0, length(y_true))
  }
  else if(d > 1){
    y_pred <- ifelse(rowSums(predictors) == d, 1, 0)
  }
  else{
    y_pred <- ifelse(predictors == d, 1, 0)
  }
  fpr <- mean(y_pred==1 & y_true==0)
  fnr <- mean(y_pred==0 & y_true==1)
  return(list(fpr=fpr, fnr=fnr))
}

ns <- round(seq(from=100, to=10000, length.out=30))
num_boots <- 10
predictors <- c("sex","married","black","asian","collegedegree","employed","militaryservice",
                "uscitizen","disability","englishability")
eps <- 1

rates <- matrix(NA, nrow=length(ns), ncol=5)
i <- 1
for(n in ns){
  central_fprs <- c()
  central_fnrs <- c()
  local_fprs <- c()
  local_fnrs <- c()
  
  # optimal t calc
  thresh_central <- -1*(10*eps/n)*log(0.2/10)*0.1
  thresh_local <- (sqrt(exp(eps))*qnorm(0.1/10)+sqrt(n))/(sqrt(n)*(1+exp(eps)))*0.1
  
  for(boot in 1:num_boots){
    
    new_data <- bootstrap(x=pums, n=n)
    new_conj_mat <- get_conjunction_matrix(new_data[,predictors], new_data$targetted)

    central_res <- names(pums[SQcentralized(conjunction=new_conj_mat, t=thresh_central, epsilon=eps)])
    local_res <- names(pums[SQlocalized(conjunction=new_conj_mat, t=thresh_local, epsilon=eps)])

    central_predictor_mat <- as.matrix(new_data[,central_res])
    local_predictor_mat <- as.matrix(new_data[,local_res])
    central_rates <- return_rates(predictors=central_predictor_mat, y_true=new_data$targetted)
    local_rates <- return_rates(predictors=local_predictor_mat, y_true=new_data$targetted)
    
    central_fprs <- c(central_fprs, central_rates$fpr)
    central_fnrs <- c(central_fnrs, central_rates$fnr)
    local_fprs <- c(local_fprs, local_rates$fpr)
    local_fnrs <- c(local_fnrs, local_rates$fnr)
  }
  
  rates[i, 1] <- n
  rates[i, 2] <- mean(central_fprs)
  rates[i, 3] <- mean(central_fnrs)
  rates[i, 4] <- mean(local_fprs)
  rates[i, 5] <- mean(local_fnrs)
  i = i + 1
}
rates_df <- as.data.frame(rates)
names(rates_df) <- c("n", "central_fpr", "central_fnr", "local_fpr", "local_fnr")

q1_plot <- ggplot(rates_df) +
  geom_line(aes(x=n, y=central_fpr, colour="Central", lty="False positive")) + 
  geom_line(aes(x=n, y=central_fnr, colour="Central", lty="False negative")) + 
  geom_line(aes(x=n, y=local_fpr, colour="Local", lty="False positive")) + 
  geom_line(aes(x=n, y=local_fnr, colour="Local", lty="False negative")) +
  scale_colour_brewer(palette="Set1") +
  theme_bw() +
  labs(x="Sample size", y="Rate") +
  theme(legend.title=element_blank())
pdf("plots/q1_plot.pdf", width=8, height=8)
q1_plot
dev.off()

