# Generate all results

# A few inputs needed to generate other parameters -----
# We use pbc_spread and default_cm from the default data inputs file
library(tidyverse)
load("data/default_inputs.Rdata")

#Testing new contact matrix data - estimated separately for each country - not used currently
load("data/contact_all.rdata")

#Projected population at country level
countries_population_2050 <- read_csv("data/countries_population_2050.csv")

#Period of analysis
start_date <- as.Date('2020-01-30')#official announcement by WHO
end_date <- as.Date('2023-02-01')#end of last significant wave

#We use reproduction number and vaccination from OWID
vax_owid <- read_csv("data/owid-covid-data.csv") 
vax_owid_world <- vax_owid %>% filter(location=="World"&date>=start_date&date<end_date) %>% arrange(date) %>% 
  mutate(people_partially_vax=people_vaccinated-people_fully_vaccinated)
curr_pop <- vax_owid_world$population[1]

#List of continents/regions
iso_code_continents <- c("OWID_AFR","OWID_ASI","OWID_EUR","OWID_NAM","OWID_OCE","OWID_SAM")
Ncountries <- length(iso_code_continents)

vax_owid_continents <- vax_owid %>% filter(iso_code %in% iso_code_continents&date>=start_date&date<end_date) %>% arrange(date) %>%
  mutate(people_partially_vax=people_vaccinated-people_fully_vaccinated) %>% arrange(iso_code)
pop_continents <- (vax_owid_continents %>% select(iso_code,population) %>% unique())$population

if(countrylevel=="all"){
  vax_owid <- vax_owid_continents
  period_size <- as.numeric(max(vax_owid$date)-min(vax_owid$date))+1#+1 to include day 1
  
  #new_v1 and new_v2 are the daily number of new vaccinations for first-dose and second-dose, respectively
  new_v1 <- vax_owid %>% group_by(iso_code) %>% 
    mutate(people_partially_vax=c(0,diff(people_partially_vax,lag=1))/curr_pop) %>% 
    ungroup() %>%
    select(iso_code,people_partially_vax,date)
  
  #If new_v1 is negative we set it to 0
  #This happens when there's more people receiving the second-dose than receiving the first-dose, in "general-example" we check that this is still conservative
  new_v1$people_partially_vax[new_v1$people_partially_vax<0] <- 0
  new_v1$people_partially_vax[is.na(new_v1$people_partially_vax)] <- 0
  new_v1 <- as.matrix(sapply(spread(new_v1, 
                                    key="iso_code", 
                                    value="people_partially_vax") %>% select(-date), as.numeric))  
  
  
  new_v2 <- vax_owid %>% group_by(iso_code) %>% 
    mutate(people_fully_vaccinated=as.numeric(c(0,diff(people_fully_vaccinated,lag=1)))/curr_pop) %>% 
    ungroup() %>%
    select(iso_code,people_fully_vaccinated,date)
  new_v2$people_fully_vaccinated[new_v2$people_fully_vaccinated<0] <- 0
  new_v2$people_fully_vaccinated[is.na(new_v2$people_fully_vaccinated)] <- 0
  new_v2 <- as.matrix(sapply(spread(new_v2, 
                                    key="iso_code", 
                                    value="people_fully_vaccinated") %>% select(-date), as.numeric))  
} 
if(countrylevel=="OWID_WRL"){
  Ncountries <- 1
  vax_owid <- vax_owid %>% filter(iso_code==countrylevel&date>=start_date&date<end_date) %>% arrange(date) %>%
    mutate(people_partially_vax=people_vaccinated-people_fully_vaccinated)
  
  period_size <- as.numeric(max(vax_owid$date)-min(vax_owid$date))+1#+1 to include day 1
  new_v1 <- c(0,diff(vax_owid$people_partially_vax,lag=1))/(unique(vax_owid$population))
  new_v1[new_v1<0] <- 0
  new_v1[is.na(new_v1)] <- 0
  new_v2 <- c(0,diff(vax_owid$people_fully_vaccinated,lag=1))/(unique(vax_owid$population))
  new_v2[is.na(new_v2)] <- 0
}

#Time series of the reproduction number
R_series <- vax_owid_world$reproduction_rate
R_series[is.na(R_series)] <- 0

#Matrix with mobility between regions
#Just a placeholder for now 
mob <- matrix(rep(0.0005,36),ncol=6)

#Baseline vaccine efficacy - after 1 dose
default_e <- 0.8
#Baseline vaccine efficacy - after 1 dose
default_e2 <- 0.9

# Demographics (for comparing HIC vs LIC)
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
world_pop_9c <- pbc_spread[countries["World"],] %>% as.numeric()
if(countrylevel=="OWID_WRL"|countrylevel=="all"){
  ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
  world_pop <- 9.5e9 #Projected world population for 2050
} else {#We're not using this case currently
  ifr_hic <- c(0.002, 0.002, 0.006, 0.006, 0.03, 0.03, 0.08, 0.08, 0.15, 0.15, 0.60, 0.60, 2.2, 2.2, 5.1, 9.3)/100 #Version with 16 age bins
  default_cm <-  contact_all[[countrylevel]]
  world_pop <- as.numeric(countries_population_2050[countries_population_2050$code==countrylevel,]$pop_2050)
}
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))
pop <- world_pop_9c/sum(world_pop_9c)
default_pdeath <- ifr_hic

#Parameters for analysis of broad spectrum antiviral
default_pbsa <- 0.7 #share of infected taking antivirals

# Case with losing immunity
kappa_default <- 0
# Case with lower supply (25% vs 100%)
default_supply_ceiling <- 1

source("R/setup.R")
