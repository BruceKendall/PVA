library(PVA)
library(ggplot2)
timehorizon <- 10
SEGsim <- replicate(1000,simulateSEG(-0.05, 0.02, 1000, timehorizon))

