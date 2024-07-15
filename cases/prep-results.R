model_i <- function(model, d1, d2, e, tbsa, ode_model="bsa",rm = FALSE) {
  pars <- grab_2v_parms(model)
  d1 <- rep(d1,Ndays)
  d1[1:tbsa]<-1
  d2 <- rep(d2,Ndays)
  d2[1:tbsa]<-1
  y <- sr(list_modify(pars, e1 = e, delta1=d1, delta2=d2), ode_model)
  if(rm) return(y)
  if(ode_model=="bsa_mob"){
    y.full <- main_metrics(y[,,1,], pop)
    for(i in 2:Ncountries){
      y.full <- rbind(y.full,main_metrics(y[,,i,], pop))
    }
    rownames(y.full)<-1:dim(y.full)[1]
    y.full <- data.frame(y.full)
    y.full['cont']<-sort(iso_code_continents)
    return(y.full)
  }
  main_metrics(y, pop)
}
# 
# df_efficacy_delta_raw <- expand_grid(d1 = seq(1, 1.5, .1),#reduction in infection period/transmissibility
#                                      d2 = seq(0, 1, .2),#reduction in mortality
#                                      e = 0.8,
#                                      tbsa=1,
#                                      model = "pars_le_covid") %>%
#   mutate(data = pmap(list(model, d1, d2, e, tbsa), function(x,y,w, z, o) data.frame(value = model_i(x,y,w, z, o), 
#                                                                      var = metric_nms))) %>%
#   unnest(data) %>%
#   spread(var, value) %>%
#   group_by(model, e)  %>%
#   mutate(model = factor(model,
#                         levels = c("Covid epidemic"),
#                         labels = c("Covid epidemic")))
# 
# saveRDS(df_efficacy_delta_raw, file = "results/df_efficacy_delta_raw.rds")
