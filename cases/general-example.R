library(kableExtra)

gg1.df <- lapply(scenario_list_2v[4], function(pars) {
  ll <- list(
    "0.25%" = list_modify(pars, e1 = 0.8, delta1=rep(1,Ndays), delta2=rep(1,Ndays))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Covid epidemic"),
                           labels = c("Covid epidemic")),
         owid=vax_owid$new_cases_smoothed_per_million/10000)

gg2.df <- lapply(scenario_list_2v[4], function(pars) {
  ll <- list(
    "0.25%" = list_modify(pars, e1 = 0.8, delta1=rep(1,Ndays), delta2=rep(1,Ndays))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"cumV",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Covid epidemic"),
                           labels = c("Covid epidemic")))

gg3.df <- lapply(scenario_list_2v[4], function(pars) {
  ll <- list(
    "0.25%" = list_modify(pars, e1 = 0.8, delta1=rep(1,Ndays), delta2=rep(1,Ndays))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"D",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Covid epidemic"),
                           labels = c("Covid epidemic")),
         owid=vax_owid$total_deaths_per_million/10000)

gg4.df <- lapply(scenario_list_2v[4], function(pars) {
  ll <- list(
    "0.25%" = list_modify(pars, e1 = 0.8, delta1=rep(1,Ndays), delta2=rep(1,Ndays))) %>%
    lapply(sr, f = "bsa") %>%
    lapply(rescale_rcs, pop, merge=T) %>%
    abind::abind()
  as.data.frame(ll[,"cumI",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Covid epidemic"),
                           labels = c("Covid epidemic")),
         owid=vax_owid$total_cases_per_million/10000)

gg1 <- gg1.df %>% select(scenario,value,time,owid) %>% unique() %>%  
  # mutate(var=factor(var,levels=c("None","0.25%","0.50%","1.00%"))) %>% 
  ggplot(aes(x = time, y = 100*value)) + geom_line() + #facet_wrap(.~scenario, scales = "free") +#, linetype=var
  geom_line(aes(y=owid),linetype = "dashed",color="red")+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) + ylab("Current infections (% pop.)")
  # lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,breaks = c("None","0.25%","0.50%","1.00%")),scalefac(0.85))+
  # theme(legend.position = "none")

gg2 <- gg2.df %>% select(value,time) %>% unique() %>%  
    ggplot(aes(x = time, y = 100*value)) + geom_line(linetype = "dashed",color="red") +
    xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) +
    ylab("Cumulative vaccinations (% pop.)")
    # lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,breaks = c("None","0.25%","0.50%","1.00%")),scalefac(0.85))+
    # labs(color="Vaccinated per day")

gg3 <- gg3.df %>% select(value,time,owid) %>% unique() %>%  
  ggplot(aes(x = time, y = 100*value)) + geom_line() +
  geom_line(aes(y=owid),linetype = "dashed",color="red")+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) +
  ylab("Cumulative deaths (% pop.)")
# lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,breaks = c("None","0.25%","0.50%","1.00%")),scalefac(0.85))+
# labs(color="Vaccinated per day")

gg4 <- gg4.df %>% select(value,time,owid) %>% unique() %>%  
  ggplot(aes(x = time, y = 100*value)) + geom_line() +
  geom_line(aes(y=owid),linetype = "dashed",color="red")+
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, Ndays, 120)) +
  ylab("Cumulative infections (% pop.)")
# lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,breaks = c("None","0.25%","0.50%","1.00%")),scalefac(0.85))+
# labs(color="Vaccinated per day")

g1_joint <- ggarrange(
  gg2 + ggtitle("Vaccinations") +
    theme(legend.spacing.x = unit(0.1, 'in'),
          text = element_text(size=9),
          legend.text = element_text(size=9),
          legend.key.size = unit(0.4, "cm")),
  common.legend = TRUE,
  gg1 + ggtitle("Infections") +
  theme(text = element_text(size=9)),
  gg4 +
    theme(text = element_text(size=9)),
  gg3 + ggtitle("Deaths") +
    theme(text = element_text(size=9)),
  ncol=1)

