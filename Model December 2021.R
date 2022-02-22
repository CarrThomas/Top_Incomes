############################## TOP INCOMES MODEL 3 #############################

# Solve for general eqbmn and build figures for the third version of the model

################################ INITIALISATION ################################

rm(list = ls())
gc()

library(ggplot2)
library(tidyr)

############################### PARAMETER VALUES ###############################

a <- 0.25 # share of managerial labour
tau <- c(0.5, 0.4) # labour market "friction" between 0 and 1, 
mu <- c(0, 0) # means of the normal distributions
o <- c(1, 1) # standard deviations of the normal distributions

# Use (m[1],o[1]) as the parameters of the m distribution, (m[2],o[2]) as the 
# parameters of the s distribution 

####################### A SPLIT THE DIFFERENCE FUNCTION ########################

# gap_func should be defined such that 
# i) function is positive below the solution
# ii) function is negative above the solution

split_the_diff <- function(gap_func, top, bot){
  
  exit <- 0
  
  while(exit == 0){
    
    mid <- (top + bot) / 2
    
    gap <- gap_func(mid)
    
    if (gap > 0){
      
      bot <- mid
      
    }
    
    if (gap <= 0){
      
      top <- mid
    }
    
    if (abs(gap) < 0.00001){
      
      exit <- 1
      
    }
  }
  
  return(mid)
  
}

########################## LABOUR SUPPLY AND DEMAND ############################

# Given w and parameters, this function returns labour supply
get_Ns <- function(w, mu, o, a, tau){
  
  chi <- (w / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
  
  temp <- function(ln_z){
    
    return((1 / o[2]) * 
             dnorm((ln_z - (mu[2] + o[2] ^ 2)) / o[2]) * 
             pnorm((ln_z + log(chi) + log((1 - tau) ^ (1 / a)) - mu[1]) / o[1]))
    
  }
  
  return((1 - tau) * exp(mu[2] + (o[2] ^ 2) / 2) * 
           integrate(temp, -Inf, Inf)$value)
  
}

# Given w and parameters, this function returns labour demand
get_Nd <- function(w, mu, o, a){
  
  chi <- (w / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
  
  temp <- function(ln_z){
    
    return((1 / o[2]) * dnorm((ln_z - mu[2]) / o[2]) * 
             (1 - pnorm((ln_z + log(chi) - (mu[1] + o[1] ^ 2)) / o[1])))
    
  }

  return(exp(mu[1] + (o[1] ^ 2) / 2) * ((1 - a) / w) ^ (1 / a) *
           integrate(temp, -Inf, Inf)$value)
  
}

# construct a sequence of w to find Ns and Nd

# lower bound chosen so that Nd = max(Ns) (for the minimum tau)
max_Ns <- (1 - min(tau)) * exp(mu[2] + o[2] / 2)

top <-  (1 - a) * (exp(mu[1] - mu[2] + (o[1] - o[2]) / 2) ^ a) /  
  ((1 - tau[1]) ^ a)
bot <- 0

gap_func <- function(w){
  
  return(get_Nd(w, mu, o, a) - (1 - min(tau)) * exp(mu[2] + o[2] / 2))
  
}

w_bot <- split_the_diff(gap_func, top, bot)

# upper bound chosen so that Ns > Nd(for the maximum tau)
top <- (1 - a) / (0.00001 ^ a) # set so that ((1 - a) / w) ^ (1 / a) = 0.00001
bot <- w_bot

gap_func <- function(w){
  
  return(exp(mu[1] + (o[1] ^ 2) / 2) * ((1 - a) / w) ^ (1 / a) - 
           get_Ns(w, mu, o, a, tau[1]))
  
}

w_top <- split_the_diff(gap_func, top, bot)

w <- seq(w_bot, w_top, length.out = 100)
Ns <-  sapply(w, get_Ns, mu = mu, o = o, a = a, tau = tau[1]) 
Nd <- sapply(w, get_Nd, mu = mu, o = o, a = a)

########################### DEMAND AND SUPPLY PLOT #############################

# Combined plot of labour demand and supply
data.frame(w = w, Ns = Ns, Nd = Nd) %>% 
  ggplot() +
  geom_line(aes(x = Ns, y = w, colour = "Labour Supply")) +
  geom_line(aes(x = Nd, y = w, colour = "Labour Demand")) +
  labs(x = "N", y = "w") +
  theme_bw() +
  theme(legend.title=element_blank())

################################ CHANGE IN TAU #################################

Ns_2 <-  sapply(w, get_Ns, mu = mu, o = o, a = a, tau = tau[2]) 

data.frame(w = w, Ns = Ns, Nd = Nd) %>% 
  ggplot() +
  geom_line(aes(x = Ns, y = w, colour = "Labour Supply")) +
  geom_line(aes(x = Nd, y = w, colour = "Labour Demand")) +
  geom_line(aes(x = Ns_2, y = w, colour = "Labour Supply 2")) +
  theme_bw() +
  theme(legend.title=element_blank())


######################### SOLVE FOR EQUILIBIRUM WAGE ###########################

# Solve for the equilibrium wage for both values of tau
w_star <- rep(0, 2)

# Target Nd - Ns. This is positive for w too low and negative for w too high

# First value of tau
low <- w[Ns - Nd < 0][which.max((Nd - Ns)[Nd - Ns > 0])]
high <- w[Ns - Nd > 0][which.min((Nd - Ns)[Nd - Ns < 0])]
gap_func <- function(w){
  
  return(get_Nd(w, mu, o, tau[1]) - get_Ns(w, mu, o, a, tau[1]))
  
}
w_star[1] <- split_the_diff(gap_func, high, low)

# Second value of tau
low <- w[Ns_2 - Nd < 0][which.max((Nd - Ns_2)[Nd - Ns_2 > 0])]
high <- w[Ns_2 - Nd > 0][which.min((Nd - Ns_2)[Nd - Ns_2 < 0])]
gap_func <- function(w){
  
  return(get_Nd(w, mu, o, tau[2]) - get_Ns(w, mu, o, a, tau[2]))
  
}
w_star[2] <- split_the_diff(gap_func, high, low)


############# SHARES OF WAGE WORKERS, OWN-ACCOUNT AND EMPLOYERS ################

ww_share <- rep(0, 2)
oa_share <- rep(0, 2)
emp_share <- rep(0, 2)

# first value of tau
chi <- (w_star[1] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      pnorm((ln_z + log(chi) + log((1 - tau[1]) ^ (1 / a)) - mu[1]) / o[1]))
 
}
ww_share[1] <- integrate(temp, -Inf, Inf)$value

temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      (pnorm((ln_z + log(chi) - mu[1]) / o[1]) 
        - pnorm((ln_z + log(chi) + log((1 - tau[1]) ^ (1 / a)) - mu[1]) / o[1])))
}
oa_share[1] <- integrate(temp, -Inf, Inf)$value

temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      (1 - pnorm((ln_z + log(chi) - mu[1]) / o[1])))
}
emp_share[1] <- integrate(temp, -Inf, Inf)$value

# second value of tau
chi <- (w_star[2] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      pnorm((ln_z + log(chi) + log((1 - tau[2]) ^ (1 / a)) - mu[1]) / o[1]))
  
}
ww_share[2] <- integrate(temp, -Inf, Inf)$value

temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      (pnorm((ln_z + log(chi) - mu[1]) / o[1]) 
       - pnorm((ln_z + log(chi) + log((1 - tau[2]) ^ (1 / a)) - mu[1]) / o[1])))
}
oa_share[2] <- integrate(temp, -Inf, Inf)$value

temp <- function(ln_z){
  
  return(
    (1 / o[2]) * 
      dnorm((ln_z - mu[2]) / o[2]) * 
      (1 - pnorm((ln_z + log(chi) - mu[1]) / o[1])))
}
emp_share[2] <- integrate(temp, -Inf, Inf)$value

ww_share
oa_share
emp_share

############################### INCOME DENSITIES ###############################

# sequence of incomes
K <- 1000
y <- seq(0, 10, length.out = K)

ww_dens <- matrix(0, K, 2)
oa_dens <- matrix(0, K, 2)
emp_dens <- matrix(0, K, 2)


# density functions
get_ww_dens <- function(y, w, tau, chi, mu, o, a){
  
  temp <- y / (w * (1 - tau))
  
  return((1 / y) * dnorm(log(temp), mu[2], o[2]) * 
    pnorm(log(temp) + log(chi * (1 - tau) ^ (1 / a)), mu[1], o[1]))
  
}

get_oa_dens <- function(y, w, tau, chi, mu, o, a){
  
  temp <- y / w
  
  temp_func <- function(z){
    
     return((1 / (a * y * z)) * 
      dnorm(a ^ (-1) * log(temp * chi ^ a / (z ^ (1 - a))), mu[1], o[1]) *
      dnorm(log(z), mu[2], o[2]))
      
  }
  
  return(integrate(temp_func, temp, temp / (1 - tau))$value)

}

get_emp_dens <- function(y, w, tau, chi, mu, o, a){
  
  temp <- y / w
  
  return((1 / y) * dnorm(log(chi * temp), mu[1], o[1]) * 
    pnorm(log(temp), mu[2], o[2]))
  
}


#### First equilibrium

chi <- (w_star[1] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)

ww_dens[2:K, 1] <- sapply(y[2:K], get_ww_dens, w = w_star[1], tau = tau[1], 
                          chi = chi, mu = mu, o = o, a = a)
oa_dens[2:K, 1] <- sapply(y[2:K], get_oa_dens, w = w_star[1], tau = tau[1], 
                          chi = chi, mu = mu, o = o, a = a)
emp_dens[2:K, 1] <- sapply(y[2:K], get_emp_dens, w = w_star[1], tau = tau[1], 
                           chi = chi, mu = mu, o = o, a = a)

data.frame(y = y, ww = ww_dens[, 1], oa = oa_dens[, 1], emp = emp_dens[, 1]) %>% 
  ggplot() +
  geom_line(aes(x = y, y = ww, colour = "Wage Workers")) +
  geom_line(aes(x = y, y = oa, colour = "Own Account Workers")) +
  geom_line(aes(x = y, y = emp, colour = "Employers")) +
  labs(x = "income", y = "") +
  theme_bw() +
  theme(legend.title=element_blank())

#### Second equilibrium

chi <- (w_star[2] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)

ww_dens[2:K, 2] <- sapply(y[2:K], get_ww_dens, w = w_star[2], tau = tau[2], 
                          chi = chi, mu = mu, o = o, a = a)
oa_dens[2:K, 2] <- sapply(y[2:K], get_oa_dens, w = w_star[2], tau = tau[2], 
                          chi = chi, mu = mu, o = o, a = a)
emp_dens[2:K, 2] <- sapply(y[2:K], get_emp_dens, w = w_star[2], tau = tau[2], 
                           chi = chi, mu = mu, o = o, a = a)

data.frame(y = y, ww = ww_dens[, 2], oa = oa_dens[, 2], emp = emp_dens[, 2]) %>% 
  ggplot() +
  geom_line(aes(x = y, y = ww, colour = "Wage Workers")) +
  geom_line(aes(x = y, y = oa, colour = "Own Account Workers")) +
  geom_line(aes(x = y, y = emp, colour = "Employers")) +
  labs(x = "income", y = "") +
  theme_bw() +
  theme(legend.title=element_blank())

#### Comparison plots for individual groups

# wage workers
data.frame(y = y, ww_1 = ww_dens[, 1], ww_2 = ww_dens[, 2]) %>%
  ggplot() + 
  geom_line(aes(x = y, y = ww_1, colour = "Wage Workers Before")) +
  geom_line(aes(x = y, y = ww_2, colour = "Wage Workers After")) +
  labs(x = "income", y = "") +
  theme_bw() +
  theme(legend.title=element_blank())

# own-account workers
data.frame(y = y, oa_1 = oa_dens[, 1], oa_2 = oa_dens[, 2]) %>%
  ggplot() + 
  geom_line(aes(x = y, y = oa_1, colour = "Non-Employers Before")) +
  geom_line(aes(x = y, y = oa_2, colour = "Non-Employers After")) +
  labs(x = "income", y = "") +
  theme_bw() +
  theme(legend.title=element_blank())

# employers
data.frame(y = y, emp_1 = emp_dens[, 1], emp_2 = emp_dens[, 2]) %>%
  ggplot() + 
  geom_line(aes(x = y, y = emp_1, colour = "Employers Before")) +
  geom_line(aes(x = y, y = emp_2, colour = "Employers After")) +
  labs(x = "income", y = "") +
  theme_bw() +
  theme(legend.title=element_blank())

#################################### TOP 1% ####################################

# Functions

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

get_emp_dist <- function(w_star, mu, o, a, tau, chi, y){
  
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

# Proceed in three steps
# i) find the overall shares of wage workers, own-account workers, and employers
# ii) find y such that the sum of the three distributions is 0.99
# iii) share in 1% can then be calculated as total share less mass below 0.99

#### first equilbrium

chi <- (w_star[1] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
ww_share <- get_ww_dist(w_star[1], mu, o, a, tau[1], chi, 10000)
oa_share <- get_oa_dist(w_star[1], mu, o, a, tau[1], chi, 10000)
emp_share <- get_emp_dist(w_star[1], mu, o, a, tau[1], chi, 10000)

ww_share + oa_share + emp_share

gap_func <- function(y){
  
  return(0.99 -  get_ww_dist(w_star[1], mu, o, a, tau[1], chi, y) -
           get_oa_dist(w_star[1], mu, o, a, tau[1], chi, y) -
           get_emp_dist(w_star[1], mu, o, a, tau[1], chi, y))
  
}

y_top1 <- split_the_diff(gap_func, 10000, 0)   

# shares
ww_top1 <- ww_share - get_ww_dist(w_star[1], mu, o, a, tau[1], chi, y_top1)
oa_top1 <- oa_share - get_oa_dist(w_star[1], mu, o, a, tau[1], chi, y_top1)
emp_top1 <- emp_share - get_emp_dist(w_star[1], mu, o, a, tau[1], chi, y_top1)

ww_top1 / (ww_top1 + oa_top1 + emp_top1)
oa_top1 / (ww_top1 + oa_top1 + emp_top1)
emp_top1 / (ww_top1 + oa_top1 + emp_top1)
(oa_top1 + emp_top1) / (ww_top1 + oa_top1 + emp_top1)

#### second equilibrium

chi <- (w_star[2] / (a ^ a * (1 - a) ^ (1 - a))) ^ (1 / a)
ww_share <- get_ww_dist(w_star[2], mu, o, a, tau[2], chi, 10000)
oa_share <- get_oa_dist(w_star[2], mu, o, a, tau[2], chi, 10000)
emp_share <- get_emp_dist(w_star[2], mu, o, a, tau[2], chi, 10000)

ww_share + oa_share + emp_share

gap_func <- function(y){
  
  return(0.99 -  get_ww_dist(w_star[2], mu, o, a, tau[2], chi, y) -
           get_oa_dist(w_star[2], mu, o, a, tau[2], chi, y) -
           get_emp_dist(w_star[2], mu, o, a, tau[2], chi, y))
  
}

y_top1 <- split_the_diff(gap_func, 10000, 0)   

# shares
ww_top1 <- ww_share - get_ww_dist(w_star[2], mu, o, a, tau[2], chi, y_top1)
oa_top1 <- oa_share - get_oa_dist(w_star[2], mu, o, a, tau[2], chi, y_top1)
emp_top1 <- emp_share - get_emp_dist(w_star[2], mu, o, a, tau[2], chi, y_top1)

ww_top1 / (ww_top1 + oa_top1 + emp_top1)
oa_top1 / (ww_top1 + oa_top1 + emp_top1)
emp_top1 / (ww_top1 + oa_top1 + emp_top1)
(oa_top1 + emp_top1) / (ww_top1 + oa_top1 + emp_top1)
