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

trimmedMean <- function(x, n, d, epsilon) {
  scale <- d/(epsilon*0.9*n)
  quants <- quantile(x, c(0.05,0.95))
  x_trimmed <- x[x>quants[1] & x<quants[2]]
  mean_trimmed <- (1/(0.9*n))*sum(x_trimmed)
  mean_release <- mean_trimmed + rlap(mu=0,b=scale,size=1)
  return(mean_release)
}

# c
exponentialPercentile <- function(x, t, d, epsilon){
  t <- t/100
  bins <- sort(c(x,d))
  nbins = length(bins)

  likelihoods <- rep(NA, nbins)
  for(i in 1:nbins){
    quality <- -1 * abs(nbins*t - i)
    if(i < nbins) { bin_width <- bins[i]-bins[i+1]+1 }
    else { bin_width <- 1 }
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
e3trimmedMean <- function(x, epsilon, n, d, tile_low=5, tile_high=95) {
  
  tile_low_value <- exponentialPercentile(x, t=tile_low, d=d, epsilon=epsilon/3)
  tile_high_value <- exponentialPercentile(x, t=tile_high, d=d, epsilon=epsilon/3)
  x_trimmed <- x[x>tile_low_value & x<tile_high_value]
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
  val <- sqrt(mean((pred-true)^2))
  return(val)
}

laplaceMeanRelease <- function(x, lower, upper, epsilon){     # borrowed from class
  n <- length(x)
  
  sensitivity <- (upper - lower)/n
  scale <- sensitivity / epsilon
  
  x.clipped <- clip(x, lower, upper)
  DPrelease <- mean(x.clipped) + rlap(mu=0, b=scale, size=1)
  
  return(DPrelease)
}

### simulations
pumas <- unique(pums$puma)
n_pumas <- length(pumas)
n_reps <- 100
means_mat <- matrix(NA, nrow=n_pumas*n_reps, ncol=6)
row_ind <- 1

for(i in 1:n_pumas){
  print(i)
  puma_region <- pumas[i]
  dat <- pums$income[pums$puma==puma_region]
  true_mean <- mean(dat)
  
  for(j in 1:n_reps){
    laplace_noise <- laplaceMeanRelease(x=dat, lower=0, upper=1000000 , epsilon=1)
    e3_trimmed <- e3trimmedMean(x=dat, epsilon=1, n=length(dat), d=1000000)
    means_mat[row_ind,1] <- puma_region
    means_mat[row_ind,2] <- true_mean
    means_mat[row_ind,3] <- laplace_noise
    means_mat[row_ind,4] <- e3_trimmed
    means_mat[row_ind,5] <- rmse(laplace_noise, true_mean)
    means_mat[row_ind,6] <- rmse(e3_trimmed, true_mean)
    row_ind <- row_ind + 1
  }
}

means_df <- data.frame(means_mat)
names(means_df) <- c("puma","true_mean","laplace_mean","exponential_mean","rmse_laplace","rmse_exponential")

### results
q1_results <- means_df %>% group_by(puma) %>% summarise(avg_lap=mean(rmse_laplace),
                                          avg_exp=mean(rmse_exponential))
write.csv(q1_results, "q1_results.csv")


q1_plot1 <- ggplot(data=means_df, aes(x=reorder(factor(puma), -laplace_mean))) + 
  geom_boxplot(aes(y=laplace_mean, shape="DP means"), outlier.shape=NA, alpha=0.7) + 
  geom_point(aes(y=true_mean, colour="True mean")) +
  scale_colour_manual(values=c("red")) +
  labs(x="PUMA region", y="Mean income", title="Laplace mechanism") + 
  theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
pdf("plots/q1_plot1.pdf", width=8, height=8)
q1_plot1
dev.off()


q1_plot2 <- ggplot(data=means_df, aes(x=reorder(factor(puma), -exponential_mean))) + 
  geom_boxplot(aes(y=exponential_mean, shape="DP means"), outlier.shape=NA, alpha=0.7) + 
  geom_point(aes(y=true_mean, colour="True mean"), alpha=0.7) +
  scale_colour_manual(values=c("red")) +
  labs(x="PUMA region", y="Mean income", title="Trimmed exponential mechanism") + 
  theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
pdf("plots/q1_plot2.pdf", width=8, height=8)
q1_plot2
dev.off()


