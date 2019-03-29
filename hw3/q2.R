# cs208 HW 3 - Question 2
# Anthony Rentsch

# Set up
rm(list = ls())
setwd("~/Desktop/Harvard/S19/cs208/cs208_AnthonyRentsch/hw3/")
require(plyr); require(dplyr); require(ggplot2)

james_update_packages <- function(packageList){
  availableRepos <- getCRANmirrors()
  flag <- availableRepos$Country=="USA" & grepl("https",availableRepos$URL)
  useRepos <- sample(availableRepos$URL[flag],1)
  
  ## install missing packages, and update if newer version available
  for(i in 1:length(packageList)){
    if (!require(packageList[i],character.only = TRUE)){
      install.packages(packageList[i], repos=useRepos)
    }
  }
  
  update.packages(ask = FALSE, dependencies = c('Suggests'), oldPkgs=packageList, repos=useRepos)
}
packagelist <- c("devtools", "jsonlite", "openssl")
devtools::install_github("privacytoolsproject/PSI-Library", ref="develop") 
library("PSIlence")

# parameters
global_epsilon = 1
global_delta = 10e-9
global_sens = 1
max_k = 100

# Laplace sd
laplaceSD <- function(epsilon) { return((sqrt(2)/epsilon)) }

# basic
basicComposition <- function(epsilon, k) { return(epsilon/k) }

# advanced
advancedComposition <- function(epsilon, k, delta) { return(epsilon/sqrt(2*k*log(1/delta))) }

# optimal
# use PSIlence:::update_parameters

# compute sds
sds <- matrix(NA, nrow=100, ncol=4)
for (k in 1:max_k) {
  epsilon_comp <- basicComposition(global_epsilon, k)
  epsilon_adv <- advancedComposition(global_epsilon, k, global_delta)
  
  # ?
  init <- rep(c(1/k, 0), k)
  params <- matrix(init, nrow=k, ncol=2, byrow=TRUE)
  inverse <- PSIlence:::update_parameters(params, hold=0, eps=global_epsilon, del=global_delta)
  epsilon_opt <- max(inverse[,1])
  
  sds[k,1] <- laplaceSD(epsilon_comp)
  sds[k,2] <- laplaceSD(epsilon_adv)
  sds[k,3] <- laplaceSD(epsilon_opt)
  sds[k,4] <- k
}

sds_df <- data.frame(sds)
names(sds_df) <- c("basic", "advanced", "optimal", "k")

q2_plot <- ggplot(sds_df) + geom_line(aes(x=k, y=basic, colour="Basic")) +
  geom_line(aes(x=k, y=advanced, colour="Advanced")) +
  geom_line(aes(x=k, y=optimal, colour="Optimal")) +
  labs(x="k", y="Standard deviation") +
  scale_colour_manual(values=c("#FF6347", "#0000FF", "#228B22")) +
  theme_bw() +
  theme(legend.title = element_blank())
pdf("plots/q2_plot.pdf", width=8, height=8)
q2_plot
dev.off()
