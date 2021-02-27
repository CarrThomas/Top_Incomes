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


#### Functions #### 

# Functions that can be used for both industries

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


# Industry 1 functions

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

# Industry 2 functions

get_z_bot <- function(P, ebar, rho, gamma, sigma){
  
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
  
  Y <- (1 - u) * ebar * EZ
  Y <- Y * (gamma / (1 + gamma)) * (1 - ebar ^ (1 + gamma))
  Y <- Y / (1 - ebar ^ gamma)
  
  return(Y)
  
}


#### Model #####

# Global Parameters
B <- 0.975
a <- 0.65
d <- 0.1
sigma <- 1 / 2
rho <- 0.73
kappa <- 0.132
Zmin <- 0.36
zeta <- 2.12

# Numerical function parameters
tol <- 0.00001

#### Industry 1 Supply Function ####

A <- c(2.3, 2.3)
gamma <- c(3 / 2, 3 / 2)

# specify a sequence of prices. For each price, obtain the mismatch threshold and supply.
# With sigma < 1, the price permitted by the price index for either good is 1. 
P <- seq(0, 1, 0.01)
n <- length(P)
ind1 <- data.frame(P, "ebar" = rep(0, n), "Y" = rep(0, n), "theta" = rep(0, n), "z" = rep(0, n), 
                    "u" = rep(0, n), "N" = rep(0, n), "EZ" = rep(0, n))

ind1$Y[1] <- 0 # no need to solve at P = 0 

for (i in 2:n){
  P <- ind1$P[i]
  # bisection search
  top <- 1
  bot <- 0
  solution <- 0
  count <- 0
  
  while (solution == 0){
    mid <- (top + bot) / 2
    gap <- free_entry_1(P, mid, a, B, d, gamma[1], A[1], rho, sigma, Zmin, zeta, kappa)
    
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
  
  ind1$ebar[i] <- mid
  ind1$Y[i] <- get_Y_1(ind1$ebar[i], P, a, B, d, gamma[1], A[1], rho, sigma, Zmin, Zeta)
  ind1$theta[i] <- get_theta(ind1$ebar[i], a, B, d, gamma[1], A[1])
  ind1$z[i] <- get_z_top(P, ind1$ebar[i], rho, sigma)
  ind1$u[i] <- get_u(ind1$theta[i], ind1$ebar[i], a, d, A[1], gamma[1])
  ind1$N[i] <- get_N_1(ind1$z[i], Zmin, zeta)
  ind1$EZ[i] <- get_EZ_1(ind1$z[i], Zmin, zeta)
}

plot(ind1$Y, ind1$P)

# #### Industry 2 Supply Function ####
# 
# A_2 <- 20.3
# gamma_2 <- 2/3
# 
# # specify a sequence of prices. For each price, obtain the mismatch threshold and supply 
# #P_2 <- (1 - P_1 ^ (1 - sigma)) ^ (1 / (1 - sigma))
# P_2 <- seq(0, 1, 0.01)
# n <- length(P_1)
# ebar_2 <- rep(0, n)
# Y_2 <- rep(0, n)
# 
# Y_2[1] <- 0
# for (i in 2:n){
#   P <- P_2[i]
#   # bisection search
#   top <- 1
#   bot <- 0
#   solution <- 0
#   count <- 0
#   
#   while (solution == 0){
#     mid <- (top + bot) / 2
#     gap <- free_entry_2(P, mid, a, B, d, gamma_2, A_2, rho, sigma, Zmin, zeta, kappa)
#     
#     if (abs(gap) < tol){
#       solution <- 1
#     }
#     if (gap >= 0){
#       top <- mid
#     }  
#     if (gap < 0){
#       bot <- mid
#     }
#     
#     count <- count + 1
#     if (count > 1000){
#       print(c('failed', P))
#       solution <- 1
#     }
#   }
#   
#   ebar_2[i] <- mid
#   Y_2[i] <- get_Y_2(ebar_2[i], P, a, B, d, gamma_2, A_2, rho, sigma, Zmin, Zeta)
# }
# 
# #plot(Y_2, P_2)
# 
# #### Solve for eqbm prices ####
# 
# # interpolate the functions above to get supply functions, and then find prices that set excess demand to 0
# P_top <- 1
# P_bot <- 0
# count <- 0
# solution <- 0
# 
# while (solution == 0){
#   p1 <- (P_top + P_bot) / 2
#   p2 <- (1 - p1 ^ (1 - sigma)) ^ (1 / (1 - sigma))
#   
#   y1_S <- approx(P_1, Y_1, p1)$y
#   y2_S <- approx(P_2, Y_2, p2)$y
#   
#   Y <- (y1_S ^ ((sigma - 1) / sigma) + y2_S ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1))
#   
#   y1_D <- Y * p1 ^ (-sigma)
#   y2_D <- Y * p2 ^ (-sigma)
#   
#   y1_gap <- y1_S - y1_D
#   y2_gap <- y2_S - y2_D
#   
#   if (abs(y1_gap) < tol){
#     print(c('converged', p1))
#     solution <- 1
#   }
#   if (y1_gap >= 0){
#     P_top <- p1
#   }  
#   if (y1_gap < 0){
#     P_bot <- p1
#   }
#   
#   count <- count + 1
#   if (count > 1000){
#     print(c('failed', p1))
#     solution <- 1
#   }
# }
#   
# #### Given eqbm prices, obtain other quantities #### 
# 
# # bisection search for industry 1 equilibrium quantities
# top <- 1
# bot <- 0
# solution <- 0
# count <- 0
#   
# while (solution == 0){
#   mid <- (top + bot) / 2
#   gap <- free_entry_1(p1, mid, a, B, d, gamma_1, A_1, rho, sigma, Zmin, zeta, kappa)
#     
#   if (abs(gap) < tol){
#     solution <- 1
#   }
#   if (gap >= 0){
#     top <- mid
#   }  
#   if (gap < 0){
#     bot <- mid
#   }
#   
#   count <- count + 1
#   if (count > 1000){
#     print(c('failed', P))
#     solution <- 1
#   }
# }
#   
# ebar_1 <- mid
# Y_1 <- get_Y_1(ebar_1, p1, a, B, d, gamma_1, A_1, rho, sigma, Zmin, Zeta)
# theta_1 <- get_theta(ebar_1, a, B, d, gamma_1, A_1)
# z_top <- get_z_top(p1, ebar_1, rho, gamma_1, sigma)
# u_1<- get_u(theta_1, ebar_1, a, d, A_1)
# N_1 <- get_N_1(z_top, Zmin, zeta)
# EZ_1 <- get_EZ_1(z_top, Zmin, zeta)
# 
# # bisection search for industry 2 eqbm quantities
# top <- 1
# bot <- 0
# solution <- 0
# count <- 0
#   
# while (solution == 0){
#   mid <- (top + bot) / 2
#   gap <- free_entry_2(p2, mid, a, B, d, gamma_2, A_2, rho, sigma, Zmin, zeta, kappa)
#     
#   if (abs(gap) < tol){
#     solution <- 1
#   }
#   if (gap >= 0){
#     top <- mid
#   }  
#   if (gap < 0){
#     bot <- mid
#   }
#     
#   count <- count + 1
#   if (count > 1000){
#     print(c('failed', P))
#     solution <- 1
#   }
# }
#   
# ebar_2 <- mid
# Y_2 <- get_Y_2(ebar_2, p2, a, B, d, gamma_2, A_2, rho, sigma, Zmin, Zeta)
# theta_2 <- get_theta(ebar_2, a, B, d, gamma_2, A_2)
# z_bot <- get_z_bot(p2, ebar_2, rho, gamma_2, sigma)
# u_2 <- get_u(theta_2, ebar_2, a, d, A_2)
# N_2 <- get_N_2(z_bot, Zmin, zeta)
# EZ_2 <- get_EZ_2(z_bot, Zmin, zeta)
