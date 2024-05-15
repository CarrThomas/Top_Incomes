rm(list = ls())
gc()

library(ggplot2)
library(tidyr)

data <- read.csv("Business Dynamics NAICS 5411.csv")

# format the data
data <- data[, c(4,5,6,7,8)]
names(data) <- c("emp_size", "year", "firms", "est", "emp")
data[, 3:5] <- apply(data[, 3:5], 2, sub, pattern = ",", replacement = "")
data[, 3:5] <- apply(data[, 3:5], 2, sub, pattern = ",", replacement = "")
data[, 3:5] <- apply(data[, 3:5], 2, sub, pattern = "D", replacement = NA)
data[, 3:5] <- apply(data[, 3:5], 2, as.numeric)
data$emp_size <- sub("Firms with ", "", data$emp_size)

sizes <- unique(data$emp_size)

i <- 9

data %>% 
  dplyr::filter(emp_size == sizes[i]) %>%
  ggplot() +
  geom_line(aes(x = year, y = emp, col = sizes[i])) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.9))


temp <- dplyr::filter(data, year == 1990)
