# Script for solving for general equilibrium

rm(list = ls())

# Rough outline of the algorithm

# STEP 1: obtain a `supply curve' for each industry

# i) fix a grid of values for prices
# ii) at each price, solve for the quantity supplied
# iii) interpolate to obtain an approximate supply curve

# Subroutine: solving for output conditional on a price
# i) Guess a value for ebar, the mismatch threshold
# ii) Back out the implied theta using get_theta
# iii) Calculate the implied occupational choice threshold
# iv) Using the choice threshold, find the mass of workers in the industry and the integral over Z
# v) put all the above together and calculate the value of creating a vacancy
# vi) find the gap between the value and cost of vacancy creation
# vii) perform a bisection search to find the ebar that makes this gap zero.
# viii) with ebar in hand, solve for u and output

# STEP 2: Find the equilibrium

# i) Guess prices
# ii) Using the interpolated supply curves, calculate the supplies of both goods at these prices and aggregate output
# iii) Calculate demands from aggregate output and prices
# iv) Find the gap between supply and demand in each industry
# v) If supply exceeds demand, revise prices down. If demand exceeds supply, revise prices up

# Notes on parameters

# Gamma
# The pdf of the mismatch shock is gamma * e ^ (gamma - 1)
# - For the pdf to be valid, it must be that gamma > 0
# - When 0 < gamma < 1, the pdf is decreasing and convex: mass is concentrated around low values of e
# - When gamma = 1, e is U[0,1]
# - When 1 < gamma < 2, the pdf is increasing and concave. In this region, higher gammas place more weight on high e
# - When 2 < gamma, the pdf is increasing and convex: mass is concentrated around high values of e

################################################## FUNCTIONS ##################################################

#### Functions for both industries ####

get_theta <- function(ebar, a, B, d, gamma, A){
  
  # returns theta given e_bar and parameters
  theta <- (1 - B * (1 - d)) * ebar / (B * A)
  theta <- theta / (ebar ^ (1 + gamma) / (1 + gamma) + gamma / (1 + gamma) - ebar)
  return(theta ^ (1 / (1 - a)))
  
  # theta <- (B / (1 - B * (1 - d))) * ((1 + gamma) / A)
  # theta <- theta * (1 - ebar ^ (1 / gamma)) / ebar ^ ((1 + gamma) / gamma)
  # return(theta ^ (1 / (1 - a)))
  
}

get_u <- function(theta, ebar, a, d, A, gamma){
  
  # returns u given theta, ebar, and parameters
  u <- d  + A * theta ^ (1 - a) * (1 - ebar ^ gamma)
  return (d / u)
  
  #u <- d  + A * theta ^ (1 - a) * (1 - ebar)
  #return (d / u)
  
}


#### Industry 1 functions ####

get_z_top <- function(P, ebar, rho, sigma){
  
  # returns the top threshold given P, ebar, and parameters
  z_top <- (rho * P * ebar) ^ (sigma - 1)
  z_top <- (z_top - 1) ^ (1 / (1 - sigma))
  
  # z_top <- (rho * P * (1 - ebar ^ (1 / gamma))) ^ (sigma - 1)
  # z_top <- (z_top - 1) ^ (1 / (1 - sigma))
  return(z_top)
  
}

get_EZ_1 <- function(z, Zmin, zeta){
  
  # returns expected Z for industry 1 given the threshold and parameters
  EZ <- Zmin * (zeta / (zeta - 1)  - (zeta * z ^ zeta) / (2 * zeta - 1)) * (z < 1) 
  EZ <- EZ + (Zmin / (z ^ (zeta - 1))) * (zeta / (zeta - 1)  - zeta / (2 * zeta - 1)) * (z >= 1)
  return(EZ)
  
}

get_N_1 <- function(z, Zmin, zeta){
  
  # returns N for industry 1 given the threshold and parameters
  N <- (1 - (z ^ zeta) / 2) * (z < 1) + (1 / (2 * z ^ zeta)) * (z >= 1)
  return(N)
  
}

free_entry_1 <- function(P, ebar, a, B, d, gamma, A, rho, sigma, Zmin, zeta, kappa){
  
  # calculates the gap between cost and benefits of entry for industry 1 given P, ebar, and parameters
  theta <- get_theta(ebar, a, B, d, gamma, A)
  z <- get_z_top(P, ebar, rho, sigma)
  N <- get_N_1(z, Zmin, zeta)
  EZ <- get_EZ_1(z, Zmin, zeta)
  
  gap <- A * B * theta ^ (- a) * (1 - rho) * P / (1 - B * (1 - d))
  gap <- gap * (gamma / (gamma + 1)) * (1 - ebar ^ (1 + gamma))
  gap <- gap * EZ / N
  gap <- gap - kappa
  
  return(gap)
  # postive if expected benefits exceed expected costs (ebar too high given p), negative otherwise (ebar too low)
  
}

get_Y_1 <- function(ebar, P, a, B, d, gamma, A, rho, sigma, Zmin, Zeta){
  
  # returns industry 1 output given P, ebar, and parameters
  theta <- get_theta(ebar, a, B, d, gamma, A)
  u <- get_u(theta, ebar, a, d, A, gamma)
  z <- get_z_top(P, ebar, rho, sigma)
  EZ <- get_EZ_1(z, Zmin, zeta)
  
  Y <- (1 - u) * EZ
  Y <- Y * (gamma / (1 + gamma)) * (1 - ebar ^ (1 + gamma))
  Y <- Y/ (1 - ebar ^ gamma)
  
  return(Y)
  
}

#### Industry 2 functions ####

get_z_bot <- function(P, ebar, rho, sigma){
  
  # returns the lower threshold given P, ebar, and parameters
  z_bot <- (rho * P * ebar) ^ (sigma - 1)
  z_bot <- (z_bot - 1) ^ (1 / (sigma - 1))
  return(z_bot)
  
}

get_EZ_2 <- function(z, Zmin, zeta){
  
  # returns expected Z for industry 1 given the threshold and parameters
  EZ <- Zmin * (z ^ (zeta - 1)) * (zeta / (zeta - 1)  - zeta / (2 * zeta - 1)) * (z < 1) 
  EZ <- EZ + Zmin * (zeta / (zeta - 1)  - zeta / ((2 * zeta - 1) * z ^ zeta)) * (z >= 1)
  return(EZ)
  
}

get_N_2 <- function(z, Zmin, zeta){
  
  # returns N for industry 1 given the threshold and parameters
  N <- ((z ^ zeta) / 2) * (z < 1) + (1 - 1 / (2 * z ^ zeta)) * (z >= 1)
  return(N)
  
}

free_entry_2 <- function(P, ebar, a, B, d, gamma, A, rho, sigma, Zmin, zeta, kappa){
  
  # calculates the gap between cost and benefits of entry for industry 1 given P, ebar, and parameters
  theta <- get_theta(ebar, a, B, d, gamma, A)
  z <- get_z_bot(P, ebar, rho, sigma)
  N <- get_N_2(z, Zmin, zeta)
  EZ <- get_EZ_2(z, Zmin, zeta)
  
  gap <- A * B * theta ^ (- a) * (1 - rho) * P / (1 - B * (1 - d))
  gap <- gap * (gamma / (1 + gamma)) * (1 - ebar ^ (1 + gamma))
  gap <- gap * EZ / N
  gap <- gap - kappa
  
  return(gap)
  # postive if expected benefits exceed expected costs (ebar too high given p), negative otherwise (ebar too low)
  
}

get_Y_2 <- function(ebar, P, a, B, d, gamma, A, rho, sigma, Zmin, Zeta){
           
  # returns industry 1 output given P, ebar, and parameters
  theta <- get_theta(ebar, a, B, d, gamma, A)
  u <- get_u(theta, ebar, a, d, A, gamma)
  z <- get_z_bot(P, ebar, rho, sigma)
  EZ <- get_EZ_2(z, Zmin, zeta)
  
  Y <- (1 - u) * EZ
  Y <- Y * (gamma / (1 + gamma)) * (1 - ebar ^ (1 + gamma))
  Y <- Y / (1 - ebar ^ gamma)
  
  return(Y)
  
}

#### Solvers ####

ebar_solve <- function(entry_func, P, mid, a, B, d, gamma, A, rho, sigma, Zmin, zeta, kappa){
  
  top <- 1
  bot <- 0
  solution <- 0
  count <- 0
  
  while (solution == 0){
    mid <- (top + bot) / 2
    gap <- entry_func(P, mid, a, B, d, gamma, A, rho, sigma, Zmin, zeta, kappa)
    if (abs(gap) < tol){
      solution <- 1
    }
    if (gap >= 0){
      bot <- mid
    }
    if (gap < 0){
      top <- mid
    }
    
    count <- count + 1
    if (count > 1000){
      print(c('failed', P))
      solution <- 1
    }
  }
  
  return(mid)
  
}

#### Income Distribution Functions ####

income_dist_1 <- function(rho, P, Zmin, z, ebar, y, gamma, zeta, u){
  
  min_inc <- rho * P * Zmin * max(1, z) * ebar
  inc_thres <- min_inc / ebar
  Zratio <- min_inc / (y * max(1, z))
  A1 <- gamma / (gamma + zeta)
  A2 <- gamma / (2 * gamma + (4 * zeta))
  A3 <- zeta / (gamma + zeta)
  A4 <- zeta / (gamma + 2 * zeta)
  
  # Region 1: below the minimum income
  if (y <= min_inc){
    return(0)
  }
  # Region 2: no agent has income below y for all epsilon
  if ((y > min_inc) & (y <= inc_thres)){
    part1 <- A1 * Zratio ^ zeta - A2 * z ^ zeta * Zratio ^ (2 * zeta)
    part2 <- Zratio ^ (-gamma) * (A3 - A4 * z ^ zeta) + z ^ zeta / 2 - 1
    part3 <- Zratio ^ (-gamma) * (A3 - A4) * z ^ (- zeta - gamma) - 1 / (2 * z ^ zeta)
    return((1 - u) * ebar ^ gamma *(part1 + part2 * (z < 1) + part3 * (z >= 1)) / (1 - ebar ^ gamma))
  }
  # Region 3: some agents have income below y for all epsilon
  if (y > inc_thres){
    part1 <- A1 * Zratio ^ zeta * (ebar ^ gamma - ebar ^ (-zeta)) - 
      Zratio ^ (2 * zeta) * A2 * z ^ zeta * (ebar ^ gamma - ebar ^ (- 2 * zeta))
    part2 <- (1 - ebar ^ gamma) * (1 - z ^ zeta / 2)
    part3 <- (1 - ebar ^ gamma) * (1 / (2 * z ^ zeta))
    return((1 - u) * (part1 + part2 * (z < 1) + part3 * (z >= 1)) / (1 - ebar ^ gamma))
  }
  
}

income_dist_2 <- function(rho, P, Zmin, z, ebar, y, gamma, zeta, u){
  
  min_inc <- rho * P * Zmin * max(1, 1 / z) * ebar
  inc_thres <- min_inc / ebar
  Zratio <- min_inc / (y * max(1, 1 / z))
  A1 <- gamma / (gamma + zeta)
  A2 <- gamma / (2 * gamma + (4 * zeta))
  A3 <- zeta / (gamma + zeta)
  A4 <- zeta / (gamma + 2 * zeta)
  
  # Region 1: below the minimum income
  if (y <= min_inc){
    return(0)
  }
  # Region 2: no agent has income below y for all epsilon
  if ((y > min_inc) & (y <= inc_thres)){
    part1 <- A1 * Zratio ^ zeta - A2 * z ^ (- zeta) * Zratio ^ (2 * zeta)
    part2 <- Zratio ^ (-gamma) * (A3 - A4) * z ^ (zeta + gamma) - z ^ zeta / 2
    part3 <- Zratio ^ (-gamma) * (A3 - A4 * z ^ (- zeta)) + 1 / (2 * z ^ zeta) - 1
    return((1 - u) * ebar ^ gamma *(part1 + part2 * (z < 1) + part3 * (z >= 1)) / (1 - ebar ^ gamma))
  }
  # Region 3: some agents have income below y for all epsilon
  if (y > inc_thres){
    part1 <- A1 * Zratio ^ zeta * (ebar ^ gamma - ebar ^ (-zeta)) - 
      Zratio ^ (2 * zeta) * A2 * z ^ (- zeta) * (ebar ^ gamma - ebar ^ (- 2 * zeta))
    part2 <- (1 - ebar ^ gamma) * (z ^ zeta / 2)
    part3 <- (1 - ebar ^ gamma) * (1 -  1 / (2 * z ^ zeta))
    return((1 - u) * (part1 + part2 * (z < 1) + part3 * (z >= 1)) / (1 - ebar ^ gamma))
  }
  
}

selfemp_mininc <- function(z_lower, z_upper, sigma, Zmin){
  if (z_upper > 1 & z_lower < 1){
    return(Zmin * 2 ^ (1 / (sigma - 1)))
  }
  if (z_upper <= 1){
    return(Zmin * (1 + z_upper ^ (1 - sigma)) ^ (1 / (sigma - 1)))
  }
    
  if (z_lower >= 1){
    return(Zmin * (1 + z_lower ^ (sigma - 1)) ^ (1 / (sigma - 1)))
  }
}

selfemp_lowerbound <- function(Zmin, sigma, y, z_lower){
  thres <- Zmin * (z_lower ^ (1 - sigma) + 1) ^ (1 / (sigma - 1))
  
  if (y <= thres){
    return(Zmin)
  }
  
  if (y > thres){
    return(y * (1 + z_lower ^ (1 - sigma)) ^ (1 / (1 - sigma)))
  }
  
}

selfemp_upperbound <- function(Zmin, sigma, y, z_upper){
  thres <- Zmin * z_upper / (1 + z_upper ^ (1 - sigma)) ^ (1 / (1 - sigma))
  
  if (y <= thres){
    return(y * Zmin / (Zmin ^ (1 - sigma) - y ^ (1 - sigma)) ^ (1 / (1 - sigma)))
  }
  if (y > thres){
    return(y * (1 + z_upper ^ (1 - sigma)) ^ (1 / (1 - sigma)))
  }
}

selfemp_integral <- function(Z, sigma, y, zeta, Zmin){
  int <- (Z ^ (1 - sigma) - y ^ (1 - sigma)) ^ ((1 + zeta) / (1 - sigma))
  int <- zeta ^ 2 * Zmin ^ (2 * zeta) * int
  int <- int / (Z ^ (2 * (1 + zeta)) * y ^ (1 + zeta))
  return(int)
}

selfemp_density <- function(z_upper, z_lower, Zmin, sigma, y, zeta){
  
  lower_bound <- selfemp_lowerbound(Zmin, sigma, y, z_lower)
  upper_bound <- selfemp_upperbound(Zmin, sigma, y, z_upper)
  
  dens <- integrate(selfemp_integral, lower = lower_bound, upper = upper_bound, sigma = sigma, y = y, zeta = zeta,
                    Zmin = Zmin)
  
  return(dens$value)
  
}
################################################## MODEL ##################################################

# Global Parameters
B <- 0.975
a <- 0.65
d <- 0.1
sigma <- 1 /10
rho <- 0.8
kappa <- 0.132
Zmin <- 0.36
zeta <- 2.12

# Numerical function parameters
tol <- 0.00001

#### Industry 1 Supply Function ####

A <- c(100, 100)
gamma <- c(3 / 2, 9 / 5)

# specify a sequence of prices. For each price, obtain the mismatch threshold and supply.
# With sigma < 1, the price permitted by the price index for either good is 1. 
P <- seq(0, 1, 0.001)
n <- length(P)
ind1 <- data.frame(P, "ebar" = rep(0, n), "Y" = rep(0, n), "theta" = rep(0, n), "z" = rep(0, n), 
                    "u" = rep(0, n), "N" = rep(0, n), "EZ" = rep(0, n))

ind1$Y[1] <- 0 # no need to solve at P = 0 

for (i in 2:n){
  P <- ind1$P[i]
  ind1$ebar[i] <- ebar_solve(free_entry_1, P, mid, a, B, d, gamma[1], A[1], rho, sigma, Zmin, zeta, kappa)
  ind1$Y[i] <- get_Y_1(ind1$ebar[i], P, a, B, d, gamma[1], A[1], rho, sigma, Zmin, Zeta)
  ind1$theta[i] <- get_theta(ind1$ebar[i], a, B, d, gamma[1], A[1])
  ind1$z[i] <- get_z_top(P, ind1$ebar[i], rho, sigma)
  ind1$u[i] <- get_u(ind1$theta[i], ind1$ebar[i], a, d, A[1], gamma[1])
  ind1$N[i] <- get_N_1(ind1$z[i], Zmin, zeta)
  ind1$EZ[i] <- get_EZ_1(ind1$z[i], Zmin, zeta)
}

plot(ind1$Y, ind1$P)

#### Industry 2 Supply Function ####

P <- seq(0, 1, 0.001)
n <- length(P)
ind2 <- data.frame(P, "ebar" = rep(0, n), "Y" = rep(0, n), "theta" = rep(0, n), "z" = rep(0, n), 
                   "u" = rep(0, n), "N" = rep(0, n), "EZ" = rep(0, n))

ind2$Y[1] <- 0

for (i in 2:n){
  P <- ind2$P[i]
  ind2$ebar[i] <- ebar_solve(free_entry_2, P, mid, a, B, d, gamma[2], A[2], rho, sigma, Zmin, zeta, kappa)
  ind2$Y[i] <- get_Y_2(ind2$ebar[i], P, a, B, d, gamma[2], A[2], rho, sigma, Zmin, Zeta)
  ind2$theta[i] <- get_theta(ind2$ebar[i], a, B, d, gamma[2], A[2])
  ind2$z[i] <- get_z_top(P, ind2$ebar[i], rho, sigma)
  ind2$u[i] <- get_u(ind2$theta[i], ind2$ebar[i], a, d, A[2], gamma[2])
  ind2$N[i] <- get_N_2(ind2$z[i], Zmin, zeta)
  ind2$EZ[i] <- get_EZ_2(ind2$z[i], Zmin, zeta)
  
}

plot(ind2$Y, ind2$P)

#### Solve for eqbm prices ####

eqbm <- data.frame("P" = rep(0, 2), "ebar" = rep(0, 2), "Y" = rep(0, 2), "theta" = rep(0, 2), "z" = rep(0, 2), 
                   "u" = rep(0, 2), "N" = rep(0, 2), "EZ" = rep(0, 2))

# interpolate the functions above to get supply functions, and then find prices that set excess demand to 0
P_top <- 1
P_bot <- 0
count <- 0
solution <- 0

while (solution == 0){
  eqbm$P[1] <- (P_top + P_bot) / 2
  eqbm$P[2] <- (1 - eqbm$P[1] ^ (1 - sigma)) ^ (1 / (1 - sigma))

  YS_1 <- approx(ind1$P, ind1$Y, eqbm$P[1])$y
  YS_2 <- approx(ind2$P, ind2$Y, eqbm$P[2])$y

  Y <- (YS_1 ^ ((sigma - 1) / sigma) + YS_2 ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1))

  YD_1 <- Y * eqbm$P[1] ^ (-sigma)
  YD_2 <- Y * eqbm$P[2] ^ (-sigma)

  Y1_gap <- YS_1 - YD_1
 
  if (abs(Y1_gap) < tol){
    print(c('converged', eqbm$P[1]))
    solution <- 1
  }
  if (Y1_gap >= 0){
    P_top <- eqbm$P[1]
  }
  if (Y1_gap < 0){
    P_bot <- eqbm$P[1]
  }

  count <- count + 1
  if (count > 1000){
    print(c('failed', eqbm$P[1]))
    solution <- 1
  }
}

#### Given eqbm prices, obtain other quantities ####

free_entry_funcs <- list(free_entry_1, free_entry_2)
get_z_funcs <- list(get_z_top, get_z_bot)
get_N_funcs <- list(get_N_1, get_N_2)
get_EZ_funcs <- list(get_EZ_1, get_EZ_2)

for (i in 1:2){
  eqbm$ebar[i] <- ebar_solve(free_entry_funcs[[i]], eqbm$P[i], mid, a, B, d, gamma[i], A[i], rho, sigma, Zmin, zeta, kappa)
  eqbm$Y[i] <- get_Y_1(eqbm$ebar[i], eqbm$P[i], a, B, d, gamma[i], A[i], rho, sigma, Zmin, Zeta)
  eqbm$theta[i] <- get_theta(eqbm$ebar[i], a, B, d, gamma[i], A[i])
  eqbm$z[i] <- get_z_funcs[[i]](eqbm$P[i], eqbm$ebar[i], rho, sigma)
  eqbm$u[i] <- get_u(eqbm$theta[i], eqbm$ebar[i], a, d, A[i], gamma[i])
  eqbm$N[i] <- get_N_funcs[[i]](eqbm$z[i], Zmin, zeta)
  eqbm$EZ[i] <- get_EZ_funcs[[i]](eqbm$z[i], Zmin, zeta)
}

# check that thresholds make sense
1 > (rho * eqbm$P[1] * eqbm$ebar[1]) ^ (1 - sigma) + (rho * eqbm$P[2] * eqbm$ebar[2]) ^ (1 - sigma)

############################################### INCOME DISTRIBUTION ###############################################

# Industry 1

min_inc_1 <- rho * eqbm$P[1] * Zmin * max(1, eqbm$z[1]) * eqbm$ebar[1]

temp <- seq(min_inc_1, 2, length.out = 100)
mass <- sapply(temp, income_dist_1, rho = rho, P = eqbm$P[1], Zmin = Zmin, z = eqbm$z[1], 
               ebar = eqbm$ebar[1], gamma = gamma[1], zeta = zeta, u = eqbm$u[1])
mass <- mass / eqbm$N[1] # get as a fraction of total people participating in industry 1

plot(temp, mass)

# Industry 2

min_inc_2 <- rho * eqbm$P[2] * Zmin * max(1, 1 / eqbm$z[2]) * eqbm$ebar[2]
temp_2 <- seq(min_inc_2, 2, length.out = 100)
mass_2 <- sapply(temp_2, income_dist_2, rho = rho, P = eqbm$P[2], Zmin = Zmin, z = eqbm$z[2], 
               ebar = eqbm$ebar[2], gamma = gamma[2], zeta = zeta, u = eqbm$u[2])
mass_2 <- mass_2 / eqbm$N[2] # get as a fraction of total people participating in industry 2

lines(temp_2, mass_2)

# Self Employed

min_inc_selfemp <- selfemp_mininc(eqbm$z[2], eqbm$z[1], sigma, Zmin)

temp_3 <- seq(min_inc_selfemp, 2, length.out = 100)
mass_3 <- sapply(temp_3, selfemp_density, z_upper = eqbm$z[1], z_lower = eqbm$z[2], Zmin = Zmin, sigma = sigma, 
                 zeta = zeta)
mass_3 <- mass_3 / (eqbm$N[1] + eqbm$N[2]) 

plot(temp_3, mass_3)
