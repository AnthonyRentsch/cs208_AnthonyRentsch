# cs208 HW 2 - Question 2
# Anthony Rentsch

# Set up
rm(list = ls())
require(plyr); require(dplyr); require(ggplot2)

# a
poissonDGP <- function(n){ return(rpois(n, lambda=10)) }

# b
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

laplaceClampMeanRelease <- function(x, epsilon, a=0, b=1){
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

rmse <- function(pred, true){ return(sqrt(mean((pred-true)^2))) }

# c
n = 200
epsilon = 0.5
b_vals = seq(from=1, to=100, by=1)
n_sims <- 100

results <- matrix(NA, nrow=(length(b_vals)*n_sims), ncol=4)
i = 1
for (b in b_vals){
  dat <- poissonDGP(n)
  for (j in 1:n_sims){
    DPrelease <- laplaceClampMeanRelease(dat, epsilon, a=0, b=b)
    results[i,1] <- b
    results[i,2] <- j
    results[i,3] <- DPrelease$release
    results[i,4] <- DPrelease$true
    i = i + 1
  }
}
results_df <- data.frame(results)
names(results_df) <- c("b", "sim", "release", "true")
avg_results_df <- results_df %>% group_by(b) %>% summarise(rmse = rmse(release, true))

q2_plot <- ggplot(data=avg_results_df, aes(x=b, y=rmse)) +
  geom_point() + geom_hline(yintercept = min(avg_results_df$rmse), col="red", lty=2) +
  geom_vline(xintercept = avg_results_df[which.min(avg_results_df$rmse), ]$b, col="red", lty=2) +
  labs(x="Upper bound", y="RMSE") + theme_bw()
pdf("plots/q2_plot.pdf", width=5, height=5)
q2_plot
dev.off()





