library(kableExtra)

#Daily infections
gg1_comp.df <- lapply(scenario_list_2v[1], function(pars) {
  ll <- list(
    "Status quo" = list_modify(pars, e1 = default_e, delta1=rep(1,Ndays), delta2=rep(1,Ndays)),
    "Immediate" = list_modify(pars, e1 = default_e, delta1=rep(1.5,Ndays), delta2=rep(0.5,Ndays)),
    "100 days" = list_modify(pars, e1 = default_e, delta1=c(rep(1,100),rep(1.5,Ndays-100)), delta2=c(rep(1,100),rep(0.5,Ndays-100))),
    "1 year" = list_modify(pars, e1 = default_e, delta1=c(rep(1,365),rep(1.5,Ndays-365)), delta2=c(rep(1,365),rep(0.5,Ndays-365)))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(var = factor(var))

#Daily number of antiviral doses required in each scenario
gg2_comp.df <- gg1_comp.df %>% mutate(value=round(default_pbsa*value*world_pop))
gg2_comp.df$month <- NA
for(scen in unique(gg2_comp.df$var)){
  if(scen=="100 days"){
    gg2_comp.df$value[gg2_comp.df$var==scen][1:100]<-0
  }
  if(scen=="1 year"){
    gg2_comp.df$value[gg2_comp.df$var==scen][1:365]<-0
  }
  gg2_comp.df$month[gg2_comp.df$var==scen] <- rep(1:ceiling(Ndays/30), each = 30)[1:Ndays]
}
gg2_comp.df <- gg2_comp.df %>% select(month,var,value) %>% group_by(var,month) %>%
  mutate(value=sum(value)) %>% ungroup() %>% unique() %>% filter(var!="Status quo")

#Cumulative deaths
gg3_comp.df <- lapply(scenario_list_2v[1], function(pars) {
  ll <- list(
    "Status quo" = list_modify(pars, e1 = default_e, delta1=rep(1,Ndays), delta2=rep(1,Ndays)),
    "Immediate" = list_modify(pars, e1 = default_e, delta1=rep(1.5,Ndays), delta2=rep(0.5,Ndays)),
    "100 days" = list_modify(pars, e1 = default_e, delta1=c(rep(1,100),rep(1.5,Ndays-100)), delta2=c(rep(1,100),rep(0.5,Ndays-100))),
    "1 year" = list_modify(pars, e1 = default_e, delta1=c(rep(1,365),rep(1.5,Ndays-365)), delta2=c(rep(1,365),rep(0.5,Ndays-365)))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"D",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(var = factor(var))

#Cumulative infections
gg4_comp.df <- lapply(scenario_list_2v[1], function(pars) {
  ll <- list(
    "Status quo" = list_modify(pars, e1 = default_e, delta1=rep(1,Ndays), delta2=rep(1,Ndays)),
    "Immediate" = list_modify(pars, e1 = default_e, delta1=rep(1.5,Ndays), delta2=rep(0.5,Ndays)),
    "100 days" = list_modify(pars, e1 = default_e, delta1=c(rep(1,100),rep(1.5,Ndays-100)), delta2=c(rep(1,100),rep(0.5,Ndays-100))),
    "1 year" = list_modify(pars, e1 = default_e, delta1=c(rep(1,365),rep(1.5,Ndays-365)), delta2=c(rep(1,365),rep(0.5,Ndays-365)))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"cumI",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(var = factor(var))

gg1_comp <- gg1_comp.df %>% select(var,value,time) %>% unique() %>%  
  ggplot(aes(x = time, y = 100*value)) +
  geom_line(aes(color=var))+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) + ylab("Current infections (% pop.)")+
lightness(scale_color_brewer(name = "Scenario", palette = "RdGy",direction = 1),scalefac(0.85))
# theme(legend.position = "none")

gg2_comp <- gg2_comp.df %>%  
  ggplot(aes(x = month, y = value)) + geom_line()+ facet_wrap(.~var, scales = "free",ncol=1) +#, linetype=var
  xlab("Time (months)") + scale_x_continuous(breaks = seq(0, ceiling(Ndays/30), 2)) + 
  ylab("Number of antiviral doses required (per month)")+theme(legend.spacing.x = unit(0.1, 'in'),
                                                               text = element_text(size=13),
                                                               axis.title=element_text(size=14),
                                                               legend.text = element_text(size=9),
                                                               legend.key.size = unit(0.4, "cm"))

write_csv(gg2_comp.df %>% rename(scenario=var,doses=value),"results/monthly_antiviral_doses.csv")

gg3_comp <- gg3_comp.df %>% select(value,time,var) %>% unique() %>%  
  ggplot(aes(x = time, y = 100*value)) +
  geom_line(aes(color=var))+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) +
  ylab("Cumulative deaths (% pop.)")+
  lightness(scale_color_brewer(name = "Scenario", palette = "RdGy",direction = 1),scalefac(0.85))

gg4_comp <- gg4_comp.df %>% select(value,time,var) %>% unique() %>%  
  ggplot(aes(x = time, y = 100*value)) + 
  geom_line(aes(color=var))+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) +
  ylab("Cumulative infections (% pop.)")+
  lightness(scale_color_brewer(name = "Scenario", palette = "RdGy",direction = 1),scalefac(0.85))

g1_comp_joint <- ggarrange(
  gg1_comp + ggtitle("Infections") +
    theme(legend.spacing.x = unit(0.1, 'in'),
          text = element_text(size=9),
          legend.text = element_text(size=9),
          legend.key.size = unit(0.4, "cm")),
  common.legend = TRUE,
  gg4_comp + ggtitle("Cumulative Infections") +
    theme(text = element_text(size=9)),
  gg3_comp + ggtitle("Cumulative Deaths") +
    theme(text = element_text(size=9)),
  ncol=1)

