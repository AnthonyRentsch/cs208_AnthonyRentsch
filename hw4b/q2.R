# cs208 HW 4b - Question 2
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw4b/")
require(plyr); require(dplyr); require(ggplot2)
pums <- read.csv("data/MAPUMS5full.csv")
pums <- pums[c("married","educ")]

# Implementation
calcllik<-function(b,data){      # function borrowed from class           
  y<-data[,1]
  x<-data[,2]
  
  pi<- 1/(1+exp(-b[1] - b[2]*x))
  
  if(pi == 1){ pi = 0.999 }
  
  llik<-y * log(pi) + (1-y) * log(1-pi) 
  return(-llik)
}

gaussianReleaseNoise <- function(size=1, epsilon, delta, clip, num_steps, batch_size, n){
  scale <- sqrt((clip/epsilon)^2 * (num_steps*batch_size/n) * log(1/delta))
  noise <- rnorm(n=size, mean=0, sd=scale)
  return(noise)
}

clip <- function(x, lower, upper){     # function borrowed from class
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

calcgradient <- function(B, C, theta, fun){     # function borrowed from class
  dx <- 0.0001
  
  out1 <-	eval(fun(b=theta, data=B))
  out2 <- eval(fun(b=theta + c(0,dx), data=B))
  out3 <- eval(fun(b=theta + c(dx,0), data=B))
  
  Del.1 <- (out3 - out1)/dx
  Del.1 <- clip(Del.1,-C,C)
  mean.Del.1 <- mean(Del.1)
  
  
  Del.2 <- (out2 - out1)/dx
  Del.2 <- clip(Del.2,-C,C)
  mean.Del.2 <- mean(Del.2)
  
  return(c(mean.Del.1,mean.Del.2))
}

calc_class_error <- function(x, y, theta1, theta2){
  y_preds_raw <- 1 / (1 + exp(-1*(theta1 + theta2*x)))
  y_preds <- round(y_preds_raw)
  error <- 1 - mean(y_preds==y)
  return(error)
}

rmse <- function(y_pred, y_true){
  return(sqrt(mean((y_pred-y_true)^2)))
}

# Run
# non-private model
true.out <- glm(married ~ educ, family="binomial", data=pums)

# private model
N <- nrow(pums)
L <- round(sqrt(nrow(pums)))
steps <- L
new.inds <- sample(1:nrow(pums))
shuffled.pums <- pums[new.inds,]
C <- 10			 
nu <- c(0.05, 0.01)

epsilons <- seq(from=0.1, to=1, by=0.1)
num_sims <- 5

res_rmse <- matrix(NA, nrow=length(epsilons), ncol=3)
row_ind_rmse <- 1
res_error <- matrix(NA, nrow=length(epsilons)*num_sims, ncol=3)
row_ind_error <- 1

for(epsilon in epsilons){
  print(paste0("Epsilon: ", epsilon))
  
  thetas <- matrix(NA, nrow=num_sims, ncol=2)
  for(sim in 1:num_sims){
    print(paste0("Simulation: ", sim))
    
    theta <- c(0,0)  
    
    for(i in 1:steps){
      startB <- ((i-1)*L+1)
      if(i<L){
        stopB <- i*L
      }else{
        stopB <- nrow(shuffled.pums)
      }
      
      index <- sample(1:nrow(shuffled.pums), L)
      B <- shuffled.pums[startB:stopB, ]
      grads <- c(0.0, 0.0)
      for(j in 1:nrow(B)){
        grads <- grads + calcgradient(B[j,], C, theta, fun=calcllik) + gaussianReleaseNoise(size=2, epsilon=epsilon/2, delta=1e-5, clip=C, num_steps=steps, batch_size=nrow(B), n=N)
      }
      theta <- theta - (grads/nrow(B) * nu) 
    }
    
    thetas[sim, 1] <- theta[1]
    thetas[sim, 2] <- theta[2]
    
    res_error[row_ind_error, 1] <- calc_class_error(pums$educ, pums$married, theta[1], theta[2])
    res_error[row_ind_error, 2] <- sim
    res_error[row_ind_error, 3] <- epsilon
    row_ind_error <- row_ind_error + 1
  }
  res_rmse[row_ind_rmse, 1] <-  rmse(thetas[1], true.out$coef[1])
  res_rmse[row_ind_rmse, 2] <-  rmse(thetas[2], true.out$coef[2])
  res_rmse[row_ind_rmse, 3] <- epsilon
  row_ind_rmse <- row_ind_rmse + 1
}

res_rmse_df <- as.data.frame(res_rmse)
names(res_rmse_df) <- c("rmse_theta1", "rmse_theta2", "epsilon")
res_error_df <- as.data.frame(res_error)
names(res_error_df) <- c("class_error", "simulation", "epsilon")

q2_plot_rmse <- ggplot(data=res_rmse_df, aes(x=epsilon)) +
  geom_line(aes(y=rmse_theta1, colour="Intercept")) +
  geom_line(aes(y=rmse_theta2, colour="Education coefficient")) +
  labs(x="Epsilon", y="RMSE") +
  scale_color_manual(values=c("blue","red")) + 
  theme_bw() +
  theme(legend.title=element_blank())
pdf("plots/q2_plot_rmse.pdf", width=8, height=8)
q2_plot_rmse
dev.off()

q2_plot_error <- ggplot(data=res_error_df, aes(x=epsilon, y=class_error)) +
  geom_point() +
  labs(x="Epsilon", y="Classification error") +
  theme_bw() +
  theme(legend.title=element_blank())
pdf("plots/q2_plot_error.pdf", width=8, height=8)
q2_plot_error
dev.off()
