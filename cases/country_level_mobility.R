library(stats)
library(MASS)

select <- dplyr::select

# nab_factor <- c(0.2,0.4,0.8,1)
# speedups <- c(1,2/3,1/2,1/3,1/4)

#Vaccine efficacy after first-dose grid
baseline_effs <- c(0.8)

#Antiviral efficacy against transmission grid (efficacy is 1-d1)
d1s <- seq(1,1.8,0.1)
#Antiviral efficacy against death grid
d2s <- seq(0, 1, .25)

# curve <- read.csv('data/curve.csv',header = FALSE)
# colnames(curve) <- c("x","y")

##Interpolating curve----
# model.eff <- approxfun(curve$x,curve$y)
# model.neutr <- approxfun(curve$y,curve$x)

##Iterating over baseline efficacies----
list_plots <- vector('list', length(baseline_effs))
le.base_eff_grid <- data.frame()
for (i in 1:length(baseline_effs)){
  base_eff <- baseline_effs[i]
  
  df_efficacy_delta_raw.base_eff_grid <- expand_grid(d1 = d1s,
                                                     d2 = d2s,
                                                     e = base_eff,
                                                     tbsa=c(1,100,365),
                                                     model = "pars_le_covid_mob") %>%
    mutate(data = pmap(list(model, d1, d2, e, tbsa), function(x,y,w, z, o) data.frame(value = model_i(x,y,w, z, o, ode_model="bsa_mob"), 
                                                                                      var = metric_nms)))  %>%
    unnest(data) %>%
    spread(var, value) %>%
    group_by(tbsa, e)  %>%
    mutate(d_rel=1-d/max(d)) %>%
    mutate(tbsa = factor(tbsa, levels = c(1,100,365),
                         labels = c("Immediate","100 day mission","1 year"))) 
  
  df_efficacy_delta.base_eff_grid <- df_efficacy_delta_raw.base_eff_grid %>%
    mutate(econ_loss = exp(0.74+0.46*log(d*1000/(Ndays/365))),
           red_bsa_infec = (d1-1)*100,
           red_bsa_death = (1-d2)*100) %>% 
    select(red_bsa_infec, red_bsa_death, e, tbsa, i,d,d_rel, econ_loss) %>%
    gather(var, value, -red_bsa_infec, -red_bsa_death, -e, -tbsa) %>%
    mutate(raw_value = value) %>%
    group_by(tbsa,var) %>%
    mutate(e_ref=base_eff) %>% 
    mutate(ref = value[e == base_eff & red_bsa_infec==0 & red_bsa_death==0]) %>%
    ungroup() %>%
    mutate(r = (value/ref)) %>% 
    mutate(le_better = cut(r, seq(0,1.0001,0.2))) %>%
    mutate(var = factor(var, levels = c("i", "d", "econ_loss"),
                        labels = c("Infections", "Deaths","Economic Loss"))) %>%
    mutate(value = round(r, 2)) %>%
    mutate(e = factor(round(e,2))) %>%
    mutate(red_bsa_infec = factor(red_bsa_infec,
                                  levels = (d1s-1)*100,
                                  labels = as.character((d1s-1)*100))) %>% 
    mutate(red_bsa_death = factor(red_bsa_death,
                                  levels = (1-d2s)*100,
                                  labels = as.character((1-d2s)*100))) %>%
    filter(var != "Economic harm")
  
  df_efficacy_delta.base_eff_grid$le_better[is.na(df_efficacy_delta.base_eff_grid$le_better)]<-"(0.8,1]"
  
  le.base_eff_grid <- rbind(le.base_eff_grid,df_efficacy_delta.base_eff_grid)
  
  df_efficacy_delta.base_eff_grid.plot <- df_efficacy_delta.base_eff_grid %>%
    ggplot(aes(x = red_bsa_infec, y = red_bsa_death, fill = le_better)) + geom_tile() +
    scale_fill_manual(values = c("grey80","grey60", "grey40", "grey20", "black"),
                      name = "") +
    theme(legend.position = "bottom", axis.text.x = element_text(hjust = 1,angle = 45)) +
    facet_grid(var~tbsa) +
    ylab("Percentage reduction in probability of death") +
    xlab("Percentage reduction in probability of transmission") +
    geom_text(aes(label = format(value,3)), color = "white", size = 3)
  list_plots[[i]] <- local(print(df_efficacy_delta.base_eff_grid.plot))
}

#Matrix of results
le_continent_mob<-list_plots[[1]]

view(df_efficacy_delta.base_eff_grid %>% filter((red_bsa_infec==50&red_bsa_death==50)|(red_bsa_infec==0&red_bsa_death==0)))
