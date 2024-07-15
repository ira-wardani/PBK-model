## set working directory
setwd("D:/PhD WUR/2nd output/PBK model")

## Install packages
install.packages("poweRlaw")
library(poweRlaw)
install.packages("ggplot2")
library(ggplot2)
install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device = "win")
install.packages("triangle")
library(triangle)
install.packages("mixtools")
library(mixtools)
install.packages("fitdistrplus")
library(fitdistrplus)
install.packages("rmutil")
library(rmutil)
install.packages("msm")
library(msm)
install.packages("readr")
library(readr)
install.packages("sn")
library(sn)
library(dplyr)
install.packages("tidyverse")
library(tidyverse)
install.packages("sigmoid")
library(sigmoid)
library(aod)
library(ggplot2)
library(LaplacesDemon)
library(writexl)


##################
#### GENERAL ########

### number of monte carlo simulations
n = 10000
df <-data.frame(n = c(1:n),
                Alpha = numeric(n),
                Size = numeric(n));

############################
#### SIZE ########
m.alpha.w = 1.60
sd.alpha.w = 0.51
min.alpha.w = 1.01
max.alpha.w = 2.56

## Alpha (-)
df$Alpha <- rtnorm(n = n, mean = m.alpha.w, sd = sd.alpha.w, lower = min.alpha.w, upper = max.alpha.w)

## generate new size distribution (0.01-5000 um)
set.seed = 123
n = 10000
xmin = 0.01 #0.01 um is the lower limit

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

write_csv(df,"D:/PhD WUR/2nd output/PBK model/PSD_human.csv")

####plot the distribution
df$log.size <- log10(df$Size)

plot(df$Alpha~df$log.size)

PSD_human <- read.csv("D:/PhD WUR/2nd output/PBK model/PSD_human.csv")
PSD_heasi <- read.csv("D:/PhD WUR/2nd output/PBK model/PSD_heasi.csv")
data_1 = PSD_human
data_2 = PSD_heasi
ggplot(NULL, aes(x=log10(Size))) +
  geom_density(data=data_1, color="midnightblue", fill="skyblue", alpha=0.5) +
  geom_density(data=data_2, color="red", fill="pink", alpha=0.5) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "Log10 (Particle Size Distribution)",
       y = "Density",
       color = "Legend")
 
PSD_heasi <- read.csv("D:/PhD WUR/2nd output/PBK model/PSD_heasi.csv")



# defining parameters (slope, x0)
# calculate 95%CI
set.seed <- 123
n <- 10000       # number of sample size
xbar <- 10        # midpoint (x0) diameter of capillary
s <- sd(PSD_human$Size)     # standard deviation of sample size
s  ### 271.2727

margin <- qt(0.975,df=n-1)*s/sqrt(n)

lowerinterval <- xbar - margin
lowerinterval  ### 4.68

upperinterval <- xbar + margin
upperinterval ## 15.32

delta_x <- upperinterval-lowerinterval
delta_x ## 10.63

CI <- c(lowerinterval,upperinterval)
CI ## 4.68-15.32

slope <- 1/log10(delta_x)
slope ## 0.97


# make the first sigmoid curve (diffusion)
y1 <- logistic(PSD_human$Size, k=0.97, x0=10)
y1 <- as.numeric(y1)

y2 = function(x) {
  1-y1
}

PSD_new <- data.frame(PSD_human$Size, y1)
PSD_new

# make the second sigmoid curve (perfusion)
perfusion = function(y1) {
  y2 <- 1 - y1
  return(y2)
}

perfusion(y1)
y2 <- perfusion(y1)
y2

## make a new dataframe
PSD_human_2 <- data.frame(PSD_human$Size, y1, y2)
PSD_human_2

colnames(PSD_human_2) <- c("Size", "P_Diffusion","P_Perfusion") 
print(PSD_human_2)

plot(y1~PSD_human_2$Size,xlab="Size", ylab="probability", 
     main="Hybrid Perfusion-Difusion conceptual model", pch=20, col = "blue")
par(new=TRUE)
plot(y2~PSD_human_2$Size,xlab="Size", ylab="probability", pch=20, col= "red")

plot_data <- data.frame(x=log10(PSD_human_2$Size), y1, y2)

ggplot(plot_data, aes(x=x, y=y1)) + 
  geom_point()+ 
  labs(title='Hybrid P-D model (0.01-5000 um)', x='Log(size)', y='Probability')

ggplot() + 
  theme_bw() +
  geom_point(data=plot_data, aes(x=x, y=y1), color='blue') + 
  geom_point(data=plot_data, aes(x=x, y=y2), color='red') +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title = "Hybrid P-D concept (0.01-5000 um)",
      x = "Log10 (Particle Size Distribution)",
       y = "Probability",
       color = "Legend")

ggplot(data=PSD_human_2, aes(x = log10(Size), y = P_Perfusion)) +
  geom_point(color = "blue", alpha = 0.5) +  # Scatter plot
  labs(title = "Size vs. Perfusion",
       x = "Size",
       y = "Perfusion") +
  theme_minimal()
