
# Parameters for each scenario -----
Ndays <- period_size
Ngroups <- 9
pre_immunity <- 0
pre_immunity_prop <- 0

# infected0 <- c("fast" = 1.085e-05,#Using OWID cumulative cases
#                "slow" = 1.085e-05,
#                "decreasing" = 1.085e-05)

# R0 adjustment:
ev <- eigen(default_cm)$values[1]
# Pre-immunity adjustment:
ev_pi <- eigen(default_cm*(1-pre_immunity))$values[1]
r0 <- function(r) r/(5*ev_pi)
ev/ev_pi


# Two doses model -----
pars_fdf_slow <- lst(
  Nc = 13, #Number of compartments
  Ngroups, #Number of age-groups
  Ndays,#Simulation window
  y0 = y0_gen(13, Ngroups, pre_immunity, 1.085e-05),#Initial values (last argument is initial number of infected individuals)
  q = rep(r0(1.1), Ngroups),#reproduction number
  contacts = default_cm,#contact matrix
  gamma1 = rep(.2, Ngroups), #(inverse) duration of exposure period 0.2
  gammaT = rep(0.5, Ngroups), #(inverse) time from infection until taking antivirals
  gamma2 = rep(.18, Ngroups), #(inverse) duration of infectious period 0.18
  delta1 = rep(1, Ngroups),#antiviral protection against transmission
  delta2 = rep(1, Ngroups),#antiviral protection against death
  kappa1 = rep(kappa_default, Ngroups),#reinfection rate after first dose of vaccine
  kappa2 = rep(kappa_default, Ngroups),#reinfection rate after second dose of vaccine
  phi = rep(1/100, Ngroups), #reinfection rate after recovery from infection
  ta = rep(0, Ngroups),#start date of vaccinations (not used currently - using OWID time series)
  e1 = rep(default_e, Ngroups),#default 
  e2 = rep(default_e2, Ngroups),
  new_v1=new_v1,#daily vaccinations - first dose
  new_v2=new_v2,#daily vaccinations - second dose
  full_period=1:Ndays,
  pdeath = default_pdeath,#IFR (=CFR in our model)
  pbsa = default_pbsa,#antiviral take-up rate
  vrf = 1,#relative vaccine hesitancy among recovered (compared to susceptible)
  vstop = rep(.8, Ngroups), #around 80% vaccinated we slow down - not used with OWID time series
  constantrisk = 0#offset of reproduction number
)

#Parameters for model with only vaccines (not used for antivirals)
pars_le_fast <- list_modify(pars_fdf_slow,
                            e1 = 0.95, e2 = 0,
                            ta1 = rep(0, Ngroups),
                            ta2 = rep(0, Ngroups),
                            tmore1 = rep(Inf, Ngroups),
                            tmore2 = rep(Inf, Ngroups),
                            ts1 = rep(Ndays, Ngroups))

#Parameters for model with antiviral 
pars_le_covid <- list_modify(pars_fdf_slow,
                             y0 = y0_gen(20, Ngroups, pre_immunity, vax_owid_world$total_cases[1]/curr_pop),
                             Nc = 20,
                            e1 = default_e, e2 = default_e2, 
                            ta1 = rep(0, Ngroups), #Same as ta but for first-dose
                            ta2 = rep(0, Ngroups), #Same as ta but for second-dose
                            tmore1 = rep(Inf, Ngroups),#Time of increase in vax rate for first-dose
                            tmore2 = rep(Inf, Ngroups),#Time of increase in vax rate for second-dose
                            ts1 = rep(Ndays, Ngroups),
                            q=matrix(rep(array(as.numeric(lapply(R_series,r0))),Ngroups),Ndays,Ngroups))

#Parameters for model with antiviral and cross-continent dynamics
i0_cont <- ((vax_owid %>% arrange(date,iso_code))$total_cases[1:Ncountries])/curr_pop
i0_cont[is.na(i0_cont)]<-1e-20
pars_le_covid_mob <- list_modify(pars_fdf_slow,
                             y0 = y0_gen_mob(20, Ngroups, Ncountries, pre_immunity, 
                                             ii=i0_cont),
                             Nc = 20,
                             Ncountries = Ncountries,
                             mob=mob,
                             e1 = default_e, e2 = default_e2, 
                             ta1 = rep(0, Ngroups), 
                             ta2 = rep(0, Ngroups), 
                             tmore1 = rep(Inf, Ngroups),
                             tmore2 = rep(Inf, Ngroups),
                             ts1 = rep(Ndays, Ngroups),#Time when first-dose starts to get priority (not used for antiviral analysis)
                             q=matrix(rep(array(as.numeric(lapply(R_series,r0))),Ngroups),Ndays,Ngroups))

scenario_list_2v <- lst(
  "Covid epidemic" = pars_le_covid,
) %>%
  setNames(c("Covid epidemic"))

