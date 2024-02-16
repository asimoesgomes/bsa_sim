# Generate all results

# A few inputs needed to generate other parameters -----
# We use pbc_spread and default_cm from the default data inputs file
library(tidyverse)
load("data/default_inputs.Rdata")

#We use reproduction rate and vaccination from OWID
vax_owid <- read_csv("data/owid-covid-data.csv") 
vax_owid <- vax_owid %>% filter(location=="World"&date<'2023-02-01') %>% arrange(date) %>% 
  mutate(people_partially_vax=people_vaccinated-people_fully_vaccinated)
period_size <- as.numeric(max(vax_owid$date)-min(vax_owid$date))+1#+1 to include day 1
new_v1 <- c(0,diff(vax_owid$people_partially_vax,lag=1))/(unique(vax_owid$population)*100)
new_v1[new_v1<0] <- 0
new_v1[is.na(new_v1)] <- 0
new_v2 <- c(0,diff(vax_owid$people_fully_vaccinated,lag=1))/(unique(vax_owid$population)*100)
new_v2[is.na(new_v2)] <- 0
R_series <- vax_owid$reproduction_rate
R_series[is.na(R_series)] <- 0

# Sensitivity analyses (global parameters to modify) ------

#Baseline efficacy
default_e <- 0.8

# Main FDF assumptions - not needed for antiviral analysis
delay_default <- 28 - 10
delay_fdf <- 84 - 10
delay_hybrid_k <- sapply(1:8, function(k){
  c(rep(delay_fdf, k), rep(delay_default, 9-k))
})
colnames(delay_hybrid_k) <- 1:8
all_k <- TRUE
default_group_seq <- TRUE

# Demographics (for comparing HIC vs LIC)
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))
pop <- hic_pop/sum(hic_pop)
default_pdeath <- ifr_hic# Changing to a single age group

#Parameters for analysis of broad spectrum antiviral
default_pbsa <- 0.35
world_pop <- 9.5e9 #Projected world population for 2050
# use_delta <- TRUE

# Case with losing immunity
kappa_default <- 0
# Case with lower supply (25% vs 100%)
default_supply_ceiling <- 1

def_labels <- list(
  "speed" = "Percentage of pop. vaccinated daily"
)

# default_speeds <- c(seq(60, 360, 10), 450, 540, 630, 730, 1460, Inf)
# main3speeds <- c(360, 180, 90)
default_speeds <- round(100/c(2, rev(seq(.05, 1, .05)), .025, .01, 0), 5) # % per day
default_speeds <- unique(c(default_speeds,
                           round(1/c(c(0.001,0.0025,0.005,0.0075,0.01,0.02),
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.2,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.4,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.6,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*2,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*4),5)))
fdf_speeds <- rev(round(100/c(.1, .25, .5, .75, 1, 2), 5))
# d1_general <- c(90, 120, 180, 360, 730, 1460)
d1_general <- 100/c(2, 1, .75, .5, .25, .1) # % per day
default_delta_value <- .0025 #for LE scenario
# le_speeds <- round(100/c(.25, .3, .35, .4, .5, 1, 2), 5)
le_speeds <- round(1/0.0025*c(1,3/4,1/2,1/4,1/3,1/8), 5)
default_speeds <- unique(c(default_speeds,le_speeds))

source("R/setup.R")
