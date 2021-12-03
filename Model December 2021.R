############################## TOP INCOMES MODEL 3 #############################

# Solve for general eqbmn and build figures for the third version of the model

################################ INITIALISATION ################################

rm(list = ls())
gc()

library(ggplot2)
library(tidyr)

############################### PARAMETER VALUES ###############################

a <- 0.5 # share of managerial labour
tau <- 0.2 # labour market "friction"
mu <- c(0, 0) # means of the normal distributions
o <- c(1, 1) # standard deviations of the normal distributions

# Use (m[1],o[1]) as the parameters of the m distribution, (m[2],o[2]) as the 
# parameters of the s distribution 

################################ LABOUR SUPPLY #################################

# Given w and parameters, this function returns labour supply
get_Ns <- function(w, mu, o, a, tau){
  
  chi <- (w / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
  
  temp <- function(ln_z){
    
    return((1 / o[2]) * 
             dnorm((ln_z - (mu[2] + o[2] ^ 2)) / o[2]) * 
             pnorm((ln_z + log(chi) + log((1 - tau) ^ (1 / a)) - mu[1]) / o[1]))
    
  }
  
  return(exp(mu[2] + (o[2] ^ 2) / 2) * integrate(temp, -Inf, Inf)$value)
  
}

# Choose the maximum w so that chi is 1000
w_max <- 1000 ^ a * (1 - a) ^ (1 - a) * a ^ a
w <- seq(0, w_max, length.out = 100)
Ns <- rep(0, 100)

# Labour Supply values
Ns[2:100] <- sapply(w[2:100], get_Ns, mu = mu, o = o, a = a, tau = tau)

################################ LABOUR DEMAND #################################

# Given w and parameters, this function returns labour demand
get_Nd <- function(w, mu, o, a, tau){
  
  chi <- (w / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
  
  temp <- function(ln_z){
    
    return((1 / o[2]) * dnorm((ln_z - mu[2]) / o[2]) * 
             (1 - pnorm((ln_z + log(chi) - (mu[1] + o[1] ^ 2)) / o[1])))
    
  }

  return(exp(mu[1] + (o[1] ^ 2) / 2) * ((1 - a) / w) ^ (1 / a) *
           integrate(temp, -Inf, Inf)$value)
  
}

# labour demand values
Nd <- sapply(w, get_Nd, mu = mu, o = o, a = a, tau = tau)

########################### DEMAND AND SUPPLY PLOT #############################

# Combined plot of labour demand and supply
data.frame(w = w, Ns = Ns, Nd = Nd) %>% 
  ggplot() +
  geom_line(aes(x = Ns, y = w, colour = "Labour Supply")) +
  geom_line(aes(x = Nd, y = w, colour = "Labour Demand")) +
  theme_bw() +
  theme(legend.title=element_blank())

######################### SOLVE FOR EQUILIBIRUM WAGE ###########################

# Use a simple split the difference algorithm to solve for the equilibrium wage

# Target Ns - Nd. This is positive for w too high and negative for w too low

# Initial lower and upper bounds obtained from Nd and Ns vectors
low <- w[Ns - Nd < 0][which.max((Ns - Nd)[Ns - Nd < 0])]
high <- w[Ns - Nd > 0][which.min((Ns - Nd)[Ns - Nd > 0])]

exit <- 0

while(exit == 0){
  
  mid <- (high + low) / 2
  
  gap <- get_Ns(mid, mu, o, a, tau) - get_Nd(mid, mu, o, a, tau)
  
  if (gap > 0){
    
    high <- mid
    
  }
  
  if (gap < 0){
    
    low <- mid
    
  }
  
  if (abs(gap) < 0.00001){
    
    exit <- 1
    
  }
  
}

w_star <- mid

######################## FRACTION IN SELF EMPLOYMENT ###########################

chi <- (w_star / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)

temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      pnorm((ln_z + log(chi) + log((1 - tau) ^ (1 / a)) - mu[1]) / o[1]))
 
}

se_share <- 1 - integrate(temp, -Inf, Inf)$value
se_share

############################# INCOME DISTRIBUTION ##############################

# sequence of incomes
K <- 1000
y <- seq(0, 10, length.out = K)

# useful quantity
chi <- (w_star / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)

#### WAGE WORKERS ####

get_ww_dist <- function(w_star, mu, o, a, tau, chi, y){
  
  temp <- function(ln_z){
  
    return(
      (1 / o[2]) * 
        dnorm((ln_z - mu[2]) / o[2]) * 
        pnorm((ln_z + log(chi) + log((1 - tau) ^ (1 / a)) - mu[1]) / o[1]))
  
  }
  
  top <- log(y) - log(w_star * (1 - tau))
  
  return(integrate(temp, -Inf, top)$value)
  
      
}

ww_dist <- rep(0, K)
ww_dist[2:K] <- sapply(y[2:K], get_ww_dist, w_star = w_star, mu = mu, o = o,
                         a = a, tau = tau, chi = chi)

#### SELF EMPLOYED (OWN ACCOUNT) ####

get_oa_dist <- function(w_star, mu, o, a, tau, chi, y){
  
  temp_one <- function(ln_z){
    
    return(
      (1 / o[2]) * 
        dnorm((ln_z - mu[2]) / o[2]) * 
        (pnorm((log(chi) + ln_z - mu[1]) / o[1]) -
           pnorm((log(chi) + log((1 - tau) ^ (1 / a)) + ln_z - mu[1]) / o[1])))
  }
  
  temp_two <- function(ln_z){
    
    return(
      (1 / o[2]) * 
        dnorm((ln_z - mu[2]) / o[2]) * 
        (pnorm((log(chi) + (1 / a) * log(y / w_star) - ((1 - a) / a) * ln_z - mu[1]) / o[1]) -
           pnorm((log(chi) + log((1 - tau) ^ (1 / a)) + ln_z - mu[1]) / o[1])))
  }
  
  mid <- log(y / w_star)
  top <- log(y / (w_star * (1  - tau)))
  
  return(integrate(temp_one, -Inf, mid)$value + 
           integrate(temp_two, mid, top)$value)
  
  
}

oa_dist <- rep(0, K)
oa_dist[2:K] <- sapply(y[2:K], get_oa_dist, w_star = w_star, mu = mu, o = o,
                          a = a, tau = tau, chi = chi)


#### SELF EMPLOYED (ENTREPRENEURS) ####

get_ent_dist <- function(w_star, mu, o, a, tau, chi, y){
  
  temp <- function(ln_z){
    
    return(
      (1 / o[2]) * 
        dnorm((ln_z - mu[2]) / o[2]) * 
        (pnorm((log(chi) + log(y / w_star) - mu[1]) / o[1]) -
        pnorm((log(chi) + ln_z - mu[1]) / o[1])))
    }
  
  top <- log(y / w_star)
  
  return(integrate(temp, -Inf, top)$value)
  
  
}

ent_dist <- rep(0, K)
ent_dist[2:K] <- sapply(y[2:K], get_ent_dist, w_star = w_star, mu = mu, o = o,
                         a = a, tau = tau, chi = chi)

#### MAKE A PLOT ####

# Distributions
data.frame(y = y, ww = ww_dist, oa = oa_dist, ent = ent_dist) %>%
  ggplot() +
  geom_line(aes(x = y, y = ww, colour = "Wage Workers")) +
  geom_line(aes(x = y, y = oa, colour = "Own Account Workers")) +
  geom_line(aes(x = y, y = ent, colour = "Entrepreneurs")) +
  theme_bw() +
  theme(legend.title=element_blank())

# Approximate Densities
data.frame(y = y[1:(K - 1)], 
           ww = ww_dist[2:K] - ww_dist[1:(K - 1)], 
           oa = oa_dist[2:K] - oa_dist[1:(K - 1)], 
           ent = ent_dist[2:K] - ent_dist[1:(K - 1)]) %>%
  ggplot() +
  geom_line(aes(x = y, y = ww, colour = "Wage Workers")) +
  geom_line(aes(x = y, y = oa, colour = "Own Account Workers")) +
  geom_line(aes(x = y, y = ent, colour = "Entrepreneurs")) +
  theme_bw() +
  theme(legend.title=element_blank())

#### TOP 1% ####

total <- ww_dist + oa_dist + ent_dist

i <- sum(total - 0.99 < 0)

ww_top1 <- ww_dist[K] - ww_dist[i]
oa_top1 <- oa_dist[K] - oa_dist[i]
ent_top1 <- ent_dist[K] - ent_dist[i]

ww_top1 + oa_top1 + ent_top1
