k <- 0.25
Z_min <- 1
zeta <- 3.42
Z_bar <- zeta * Z_min / (zeta - 1)

pi <- (Z_min / (Z_bar + k)) ^ zeta
y_star <- ((pi * Z_bar ^ zeta + (1 - pi) * Z_min ^ zeta) / 0.01) ^ (1 / zeta)

top1_se <- pi * (Z_bar / y_star) ^ zeta
top1_emp <- (1 - pi) * (Z_min / y_star) ^ zeta

dy_emp <- function(y, Z_bar, Z_min, pi, zeta){
  if (y < Z_min){
    return(0)
    
  }
  else{
    return((1 - pi) * zeta * Z_min ^ zeta * y ^ (-(zeta + 1)))
  }
}

dy_se <- function(y, Z_bar, Z_min, pi, zeta){
  if (y < Z_bar){
    return(0)
  }
  else{
    return(pi * zeta * Z_bar ^ zeta * y ^ (-(zeta + 1)))
  } 
}


y <- seq(Z_min, y_star + 1, length.out = 100)
dist_emp <- sapply(y, dy_emp, Z_bar = Z_bar, Z_min = Z_min, pi = pi, zeta = zeta)
dist_se <- sapply(y, dy_se, Z_bar = Z_bar, Z_min = Z_min, pi = pi, zeta = zeta)

plot(y, dist_emp, type = 'l', col = 4)
lines(y, dist_se, col = 2)
legend(4, 2.5, c("employees", "self employed"), col = c(4, 2), lty = c(1, 1))


