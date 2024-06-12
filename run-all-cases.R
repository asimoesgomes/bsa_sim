# list_countries <- scan("data/list_countries.txt", what="", sep="\n")

list_countries <- c('OWID_WRL')#c("ARG","BEL","BRA","CAN","COL","COM","ESP","DEU","ISR","ITA","USA","ZAF","CHN","IND","KEN")

for(countrylevel in list_countries){
  start.time <- Sys.time()
  source("project-setup.R")

  # source("cases/prep-results.R")
  # source("cases/general-example.R")
  # source("cases/general-example-and-main.R")
  source("cases/expected_benefits.R")
  # source("cases/lower_efficacy_baseline_grid.R")

  if(countrylevel=="OWID_WRL"){
    fig_folder <- "figures/"
  } else if(countrylevel=="all"){
    source("cases/country_level_mobility.R")
  } else {
    dir.create(paste0("figures/",countrylevel), showWarnings = FALSE)
    fig_folder <- paste0("figures/",countrylevel)
  }
  width <- 6.5
  source("cases/generate-figures.R")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0(countrylevel,": ",time.taken))
}
