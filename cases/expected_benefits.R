library(stats)
library(MASS)
library(EnvStats)
library(cascsim)

select <- dplyr::select

#Assumed distribution for reproduction number
R0_draw <- rtri(1000,0.5,2.8,(0.5+2.8)/2)
#Assumed distribution for CFR
#cfr_draw <- #rtbeta(n, shape1, shape2, ncp = 0, min = 0, max = 1) - not sure what to assume as shape1 and 2

exp_ben <- function(model, tbsa, R0, cfr=1, rm = FALSE) {
  d1 <- 1.5
  d2 <- 0.5
  e <- 0.8
  pars <- grab_2v_parms(model)
  
  d1 <- rep(d1,Ndays)
  d1[1:tbsa]<-1
  d2 <- rep(d2,Ndays)
  d2[1:tbsa]<-1
  
  R_weighted <- mean(R_series[R_series>0][1:7])/R0*R_series
  q <- matrix(rep(array(as.numeric(lapply(R_weighted,r0))),Ngroups),Ndays,Ngroups)
  
  #pdeath <- rep(cfr,Ngroups)
  pdeath <- default_pdeath
  
  y <- sr(list_modify(pars, pdeath=pdeath, q=q, e1 = e, delta1=d1, delta2=d2), "bsa")
  if(rm) return(y)
  main_metrics(y, pop)
}

##Iterating over sampled R0s and CFRs for reference case (50% protection available after 100 days)----


df_r0_cfr_raw.base_eff_grid <- expand_grid(tbsa=c(100,Ndays),
                                                   R0=R0_draw,
                                                   #cfr=cfr_draw,
                                                   model = "pars_le_covid") %>%
  mutate(data = pmap(list(model, tbsa, R0), function(x, o, q) data.frame(value = exp_ben(x, o, q), 
                                                                                    var = metric_nms)))  %>%
  unnest(data) %>%
  spread(var, value) %>%
  #mutate(nab_factor = factor(round(e,3),levels=round(effs,3),labels=nab_factor_mod)) %>% 
  group_by(R0)  %>%#group_by(R0, cfr)  %>%
  mutate(d_rel=1-d/max(d)) %>%
  mutate(tbsa = factor(tbsa, levels = c(100,Ndays),
                       labels = c("100 day mission","Default"))) 

df_r0_cfr_raw.base_eff_grid <- df_r0_cfr_raw.base_eff_grid %>%
  mutate(econ_loss = exp(0.74+0.46*log(d*1000/(Ndays/365)))) %>% 
  # select(R0, cfr, tbsa, i,d,d_rel,harm, econ_loss) %>%
  # gather(var, value, -R0, -cfr, -tbsa) %>%
  # mutate(raw_value = value) %>%
  # group_by(R0,cfr,var) %>%
  select(R0, tbsa, i,d,d_rel,harm, econ_loss) %>%
  gather(var, value, -R0, -tbsa) %>%
  mutate(raw_value = value) %>%
  group_by(R0,var) %>%
  mutate(ref = value[tbsa == "Default"]) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, seq(0,1.0001,0.2)
                         # labels = c("Fractional dose better",
                         #            "Full dose better")
  )) %>%
  mutate(var = factor(var, levels = c("i", "d", "econ_loss"),#,"d_rel"
                      labels = c("Infections", "Deaths","Economic Loss"))) %>%#,"Relative Deaths"
  mutate(value = round(r, 2)) %>%
  filter(var != "Economic harm")

#Expected outcomes
#With antivirals
mean((df_r0_cfr_raw.base_eff_grid %>% filter(tbsa=="100 day mission"&var=="Deaths") %>% unique())$raw_value)
#Without antivirals
mean((df_r0_cfr_raw.base_eff_grid %>% filter(tbsa!="100 day mission"&var=="Deaths") %>% unique())$raw_value)
