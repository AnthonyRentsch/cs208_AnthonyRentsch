# cs208 HW 2 - Question 3
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw2/")
require(plyr); require(dplyr); require(ggplot2)

# a
poissonDGP <- function(n){ return(rpois(n, lambda=10)) }
noisyLinearDGP <- function(x, n, alpha, beta, mu=0, sd=1) { return(beta*x + alpha + rnorm(n, mu, sd)) }

sgn <- function(x) {     # function borrowed from class
  return(ifelse(x < 0, -1, 1))
}

rlap = function(mu=0, b=1, size=1) {     # function borrowed from class
  p <- runif(size) - 0.5
  draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
  return(draws)
}

clip <- function(x, lower, upper){     # function borrowed from class
  x.clipped <- x
  x.clipped[x.clipped<lower] <- lower
  x.clipped[x.clipped>upper] <- upper
  return(x.clipped)	
}

rmse <- function(pred, true){ return(sqrt(mean((pred-true)^2))) }

laplaceClampMeanRelease <- function(x, epsilon, a=0, b=1){     # from q2
  n <- length(x)
  sensitivity <- (a - b)/n
  scale <- sensitivity / epsilon
  
  x.clipped <- clip(x, a, b)
  clipped.mean <- mean(x.clipped)
  noisy.mean <- clipped.mean + rlap(mu=0, b=scale, size=1)
  release.mean <- clip(noisy.mean, a, b)
  true.mean <- mean(x)
  
  return(list(release=release.mean, true=true.mean))
}

regressionRelease <- function(y, x, ylower=0, yupper=17, xlower=0, xupper=19, eplsilon, partition){     # augmented from class
  x <- clip(x, xlower, xupper)
  y <- clip(y, ylower, yupper)
  
  n <- length(x)
  sens.Sxy <- ((xupper-xlower)*(yupper-ylower))   
  sens.Sxx  <- ((xupper-xlower)^2)   
  
  scale.Sxy <- sens.Sxy / (epsilon*partition$Sxy)
  scale.Sxx <- sens.Sxx / (epsilon*partition$Sxx)
  
  true.beta <- sum((x - mean(x))*(y - mean(y))) / sum((x - mean(x))^2) 
  true.alpha <- mean(y) - true.beta*mean(x)
  
  release.Sxy <- sum((x - mean(x))*(y - mean(y)))  + rlap(mu=0, b=scale.Sxy, size=1)
  release.Sxx <- sum((x - mean(x))^2) + rlap(mu=0, b=scale.Sxx, size=1)  
  release.beta <- release.Sxy/release.Sxx
  
  release.x.bar <- laplaceClampMeanRelease(x, epsilon*partition$x.bar, a=xlower, b=xupper)$release
  release.y.bar <- laplaceClampMeanRelease(y, epsilon*partition$y.bar, a=ylower, b=yupper)$release
  release.alpha <- release.y.bar - release.beta*release.x.bar
  
  release.mean.sq.residuals <- mean((y - release.beta*x - release.alpha)^2)
  true.mean.sq.residuals <- mean((y - true.beta*x - true.alpha)^2)
  
  return(list(release.beta=release.beta, 
              release.alpha=release.alpha,
              true.beta=true.beta, 
              true.alpha=true.alpha,
              release.mean.sq.residuals=release.mean.sq.residuals, 
              true.mean.sq.residuals=true.mean.sq.residuals))
}
  
# get optimal upper bound for y
n = 200
epsilon = 1
b_vals = seq(from=1, to=100, by=1)
n_sims <- 100

results_y <- matrix(NA, nrow=(length(b_vals)*n_sims), ncol=4)
i = 1
for (b in b_vals){
  dat <- poissonDGP(n)
  y <- noisyLinearDGP(dat, n, alpha=1, beta=1, mu=0, sd=1)
  for (j in 1:n_sims){
    DPrelease <- laplaceClampMeanRelease(y, epsilon, a=0, b=b)
    results_y[i,1] <- b
    results_y[i,2] <- j
    results_y[i,3] <- DPrelease$release
    results_y[i,4] <- DPrelease$true
    i = i + 1
  }
}
results_y_df <- data.frame(results_y)
names(results_y_df) <- c("b", "sim", "release", "true")
avg_results_y_df <- results_y_df %>% group_by(b) %>% summarise(rmse = rmse(release, true))

q3_plot1 <- ggplot(data=avg_results_y_df, aes(x=b, y=rmse)) +
  geom_point() + geom_hline(yintercept = min(avg_results_y_df$rmse), col="red", lty=2) +
  geom_vline(xintercept = avg_results_y_df[which.min(avg_results_y_df$rmse), ]$b, col="red", lty=2) +
  labs(x="Upper bound for y", y="RMSE") + theme_bw()
pdf("plots/q3_plot1.pdf", width=8, height=8)
q3_plot1
dev.off()

# b
equal_partition <- list(Sxy=0.25, Sxx=0.25, x.bar=0.25, y.bar=0.25)
n = 1000
alpha = beta = eplison = sd = 1
n_sims = 1000

results_reg <- matrix(NA, nrow=n_sims, ncol=6)
for (i in 1:n_sims){
  x <- poissonDGP(n)
  y <- noisyLinearDGP(dat, n, alpha=1, beta=1, mu=0, sd=1)
  DPrelease <- regressionRelease(y, x, ylower=0, yupper=17, xlower=0, xupper=19, eplsilon, equal_partition)
  results_reg[i,1] <- DPrelease$release.beta
  results_reg[i,2] <- DPrelease$release.alpha
  results_reg[i,3] <- DPrelease$true.beta
  results_reg[i,4] <- DPrelease$true.alpha
  results_reg[i,5] <- DPrelease$release.mean.sq.residuals
  results_reg[i,6] <- DPrelease$true.mean.sq.residuals
}
results_reg_df <- data.frame(results_reg)
names(results_reg_df) <- c("release.beta", "release.alpha", "true.beta", 
                           "true.alpha", "release.mean.sq.residuals", "true.mean.sq.residuals")

semi.blue <- rgb(0,90,239,50,maxColorValue=255)
semi.red  <- rgb(239,90,0,200,maxColorValue=255)
q3_plot2 <- ggplot(data=results_reg_df) + 
  geom_density(aes(x=release.mean.sq.residuals, fill="DP regression"), alpha=0.4, colour=NA) +
  geom_density(aes(x=true.mean.sq.residuals, fill="Actual regression"), alpha=0.4, colour=NA) +
  labs(x="Mean Squared residuals", y="") + theme_bw() + 
  scale_fill_manual(values=c(semi.blue, semi.red)) +
  theme(legend.position="bottom", legend.title=element_blank())
pdf("plots/q3_plot2.pdf", width=8, height=8)
q3_plot2
dev.off()


  