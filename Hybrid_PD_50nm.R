## set working directory
setwd("D:/PhD WUR/2nd output/PBK model")

## download library

library(poweRlaw)
library(ggplot2)
library(extrafont)
font_import()
loadfonts(device = "win")
install.packages("triangle")
library(triangle)
library(mixtools)
library(fitdistrplus)
library(rmutil)
library(msm)
library(readr)
library(sn)
library(dplyr)
library(tidyverse)
library(sigmoid)
library(aod)
library(ggplot2)
library(LaplacesDemon)
library(writexl)


##################
#### Realignment ########

### number of monte-carlo simulations
n = 10000
df <-data.frame(n = c(1:n),
                Alpha = numeric(n),
                Size = numeric(n));

############################
#### SIZE (derived the parameters from Heasi model) ########
m.alpha.w = 1.60
sd.alpha.w = 0.51
min.alpha.w = 1.01
max.alpha.w = 2.56

## Alpha (-)
df$Alpha <- rtnorm(n = n, mean = m.alpha.w, sd = sd.alpha.w, lower = min.alpha.w, upper = max.alpha.w)

## generate new size distribution (1 nm-5000 um)
set.seed = 123
n = 10000
xmin = 0.001 #0.001 um (1 nm) is the lower limit

i = 1
for(i in 1:n){
  # generate data, should be smaller than 5000 um
  X.func <- function (X){
    success <- FALSE
    while (!success){
      U = runif(1, 0, 1)
      X = xmin*(1-U)^(1/(1-df$Alpha[i]))
      success <- X < 5000} ##should be smaller than 5000 um
    return(X)
  }
  
  df$Size[i] <- X.func()
}

write_csv(df,"D:/PhD WUR/2nd output/PBK model/PSD_wardani.csv")

####plot the distribution PSD_wardani vs PSD_heasi#####
###Check the distribution from the realignment
df$log.size <- log10(df$Size)
plot(df$Alpha~df$log.size)

###
PSD_wardani <- read.csv("D:/PhD WUR/2nd output/PBK model/PSD_wardani.csv")
PSD_heasi <- read.csv("D:/PhD WUR/2nd output/PBK model/PSD_heasi.csv")

ggplot(NULL, aes(x=log10(Size))) +
  geom_density(data=PSD_wardani, color="midnightblue", fill="skyblue", alpha=0.5) +
  geom_density(data=PSD_heasi, color="red", fill="pink", alpha=0.5) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "Log10 (Particle Size Distribution)",
       y = "Density",
       color = "Legend")

### Build Hybrid P-D size selection ###

# defining parameters (slope, x0)
# calculate 90%CI (to get smallest range of Hybrid P-D size)

set.seed(123)
n <- 100000       # number of sample size
xbar <- 0.05       # midpoint (x0) 50 nm
s <- sd(PSD_wardani$Size)
s

margin <- qt(0.95, df=n-1) * s / sqrt(n)
lowerinterval <- xbar - margin
upperinterval <- xbar + margin
delta_x <- upperinterval - lowerinterval
CI <- c(lowerinterval, upperinterval)

s    # 202.6425
margin # 1.214601 (log10 scale)

lowerinterval  # Lower bound of the 90% CI  1.054052 (log10 scale)
upperinterval  # Upper bound of the 90% CI -1.104052 (log10 scale)
delta_x        # Width of the 90% CI 2.108103 (log10 scale)
CI             # Confidence interval -1.004052  1.104052 (log10 scale)

slope <- 1/log10(delta_x)

slope ## 3.08745

# Define logistic function
logistic <- function(x, k, x0) {
  1 / (1 + exp(-k * (log10(x) - log10(x0))))
}

# Assuming PSD_wardani is already loaded and has a column named 'Size'

# Create the diffusion data
diffusion <- logistic(PSD_wardani$Size, k=slope, x0=0.05)

# Define perfusion function
perfusion <- function(x) {
  1 - logistic(x, k=slope, x0=0.05)
}

# Create the perfusion data
perfusion_data_1 <- perfusion(PSD_wardani$Size)

# Create the new dataframe called PSD_human
PSD_human_1 <- data.frame(Size = PSD_wardani$Size, P_Diffusion = diffusion, P_Perfusion = perfusion_data_1)

# Print the dataframe
print(PSD_human_1)
write.csv(PSD_human_1, "Hybrid_PD_selection.csv", row.names=FALSE)

# Prepare data for plotting
plot_data <- PSD_human_1 %>%
  mutate(logSize = log10(Size))

# Plot the data
ggplot(plot_data) + 
  theme_bw() +
  geom_point(aes(x = logSize, y = P_Diffusion), color = 'blue') + 
  geom_point(aes(x = logSize, y = P_Perfusion), color = 'red') +
  geom_vline(xintercept = log10(0.05), linetype = "dashed", color = "black", linewidth = 1) + 
  scale_x_continuous(
    breaks = seq(-3, 4, by = 1),  
    labels = seq(-3, 4, by = 1)   
  ) +
  theme(
    axis.line = element_line(colour = "black", size = 1.2),
    axis.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    plot.title = element_text(size = 18, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  labs(
    x = "Log10 (particle size , in micrometer)",
    y = "Probability",
    color = "Legend"
  )


