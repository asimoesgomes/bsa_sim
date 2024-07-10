library(stats)
library(MASS)

select <- dplyr::select

#Defining function for sensitivity analysis of the antiviral take-up rate
sens_takeup <- function(model, d1, d2, e, tbsa, pbsa, rm = FALSE) {
  pars <- grab_2v_parms(model)
  d1 <- rep(d1,Ndays)
  d1[1:tbsa]<-1
  d2 <- rep(d2,Ndays)
  d2[1:tbsa]<-1
  y <- sr(list_modify(pars, e1 = e, delta1=d1, delta2=d2, pbsa = pbsa), "bsa")
  if(rm) return(y)
  main_metrics(y, pop)
}

#Reference antiviral efficacy
d1s <- c(1.5)
d2s <- c(0.5)

#List of antiviral take-up rates to iterate over
pbsas <- seq(0,1,0.05)

#Vaccine efficacy
base_eff <- 0.8

df_sens_takeup_raw <- expand_grid(d1 = d1s,
                                                   d2 = d2s,
                                                   e = base_eff,
                                                   tbsa=c(1,100,365),
                                                   model = "pars_le_covid",
                                                   pbsa=pbsas) %>%
  mutate(data = pmap(list(model, d1, d2, e, tbsa, pbsa), function(x,y,w, z, o, p) data.frame(value = sens_takeup(x,y,w, z, o, p), 
                                                                                    var = metric_nms)))  %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(tbsa, e, d1, d2)  %>%
  mutate(d_rel=1-d/max(d)) %>%
  mutate(tbsa = factor(tbsa, levels = c(1,100,365),
                       labels = c("Immediate","100 day mission","1 year"))) %>% 
  ungroup()

df_sens_takeup <- df_sens_takeup_raw %>%
  mutate(econ_loss = exp(0.74+0.46*log(d*1000/(Ndays/365))),
         red_bsa_infec = (d1-1)*100,
         red_bsa_death = (1-d2)*100) %>% 
  select(tbsa, pbsa, i,d,d_rel, econ_loss) %>% 
  gather(var, value, -pbsa, -tbsa) %>%
  mutate(raw_value = value) %>%
  group_by(tbsa,var) %>%
  mutate(ref = value[pbsa==0]) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(var = factor(var, levels = c("i", "d", "econ_loss"),#,"d_rel"
                      labels = c("Infections", "Deaths","Economic Loss"))) %>%#,"Relative Deaths"
  mutate(value = round(r, 2)) %>% 
  filter(var != "Economic harm")

df_sens_takeup.plot <- df_sens_takeup %>%
  ggplot(aes(x = pbsa*100, y = (1-value))) + geom_line(linewidth=1) +
  theme(legend.position = "bottom", axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_grid(var~tbsa) +
  ylab("Percentage reduction in deaths (relative to no antiviral)") +
  xlab("Share of infected individuals taking antivirals [%]")
