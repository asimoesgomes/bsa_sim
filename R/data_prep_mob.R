library(tidyverse)
library(ggplot2)
library(stats)
library(stringi)

#Skewd, decimals, R1 and R2, number of countries?

##DF1 201901-201906##----
df1 <- read_delim("data/University_Of_Chicago/SM_ILLUM_UNIV_CHI_07012024A.txt")
colnames(df1)
head(df1)
table(df1$TVLDT)
# 201901  201902  201903  201904  201905  201906 
# 1762151 1632752 1753318 1756702 1823975 1902108

df1 <- df1 %>% select(ORG_CTRY,DST_CTRY,TTL_BKGS) %>% 
  rename(org_ctry=ORG_CTRY,
         dst_ctry=DST_CTRY,
         bkgs=TTL_BKGS) %>% 
  mutate(bkgs=as.numeric(bkgs))

ggplot(df1 %>% filter(bkgs<1000), aes(x=bkgs)) + geom_histogram()
min(df1$bkgs)#0.04
max(df1$bkgs)#229727.3
median(df1$bkgs)#2.52
quantile(df1$bkgs,0.95)#104.99 - EXTREMELY SKEWED

df1 <- df1 %>% group_by(org_ctry,dst_ctry) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique()

write_csv(df1,"data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024A.csv")
rm(df1)

##DF2 201907-201912##----
df2 <- read_delim("data/University_Of_Chicago/SM_ILLUM_UNIV_CHI_07012024B.txt")
colnames(df2)
head(df2)
table(df2$TVLDT)
# 201907  201908  201909  201910  201911  201912 
# 1904185 1939146 1894715 1898627 1709333 1784573

df2 <- df2 %>% select(ORG_CTRY,DST_CTRY,TTL_BKGS) %>% 
  rename(org_ctry=ORG_CTRY,
         dst_ctry=DST_CTRY,
         bkgs=TTL_BKGS) %>% 
  mutate(bkgs=as.numeric(bkgs))

ggplot(df2%>% filter(bkgs<1000), aes(x=bkgs)) + geom_histogram()
min(df2$bkgs)#0.04
max(df2$bkgs)#235483.7
median(df2$bkgs)#2.52
quantile(df2$bkgs,0.95)#106.241 - EXTREMELY SKEWED

df2 <- df2 %>% group_by(org_ctry,dst_ctry) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique()

write_csv(df2,"data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024B.csv")
rm(df2) 
##DF3 202301-202306##----
df3 <- read_delim("data/University_Of_Chicago/SM_ILLUM_UNIV_CHI_07012024C.txt")
colnames(df3)
head(df3)
table(df3$TVLDT)
# 202301  202302  202303  202304  202305  202306 
# 1752221 1665812 1914808 1834220 1919259 1977881

df3 <- df3 %>% select(ORG_CTRY,DST_CTRY,TTL_BKGS) %>% 
  rename(org_ctry=ORG_CTRY,
         dst_ctry=DST_CTRY,
         bkgs=TTL_BKGS) %>% 
  mutate(bkgs=as.numeric(bkgs))

min(df3$bkgs)#0.06
max(df3$bkgs)#142863.9
median(df3$bkgs)#2
quantile(df3$bkgs,0.95)#94.8 - EXTREMELY SKEWED
ggplot(df3 %>% filter(bkgs<1000), aes(x=bkgs)) + geom_histogram()

df3 <- df3 %>% group_by(org_ctry,dst_ctry) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique()

write.csv(df3,"data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024C.csv")
rm(df3) 
##DF4 202307-202312##----
df4 <- read_delim("data/University_Of_Chicago/SM_ILLUM_UNIV_CHI_07012024D.txt")
colnames(df4)
head(df4)
table(df4$TVLDT)
# 202307  202308  202309  202310  202311  202312 
# 1947619 1929436 1942396 1986102 1827495 1921348

df4 <- df4 %>% select(ORG_CTRY,DST_CTRY,TTL_BKGS) %>% 
  rename(org_ctry=ORG_CTRY,
         dst_ctry=DST_CTRY,
         bkgs=TTL_BKGS) %>% 
  mutate(bkgs=as.numeric(bkgs))

min(df4$bkgs)#0.03
max(df4$bkgs)#171924.4
median(df4$bkgs)#2.2
quantile(df4$bkgs,0.95)#106.08 - EXTREMELY SKEWED
ggplot(df4 %>% filter(bkgs<1000), aes(x=bkgs)) + geom_histogram()

df4 <- df4 %>% group_by(org_ctry,dst_ctry) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique()

write.csv(df4,"data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024D.csv")
rm(df4) 


##Exploring at country level 2019##----
df1 <- read_csv("data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024A.csv")
df2 <- read_csv("data/University_Of_Chicago/agg_SM_ILLUM_UNIV_CHI_07012024B.csv")

df <- rbind(df1,df2) %>% group_by(org_ctry,dst_ctry) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique()

dim(df)
dim(df %>% select(org_ctry,dst_ctry) %>% unique())

write_csv(df,"data/University_Of_Chicago/agg_ctry_2019.csv")

df.org <- df %>% select(org_ctry,bkgs) %>% 
  group_by(org_ctry) %>% mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% unique()
write_csv(df.org,"data/University_Of_Chicago/agg_orgctry_2019.csv")
df.dst <- df %>% select(dst_ctry,bkgs) %>% 
  group_by(dst_ctry) %>% mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% unique()
write_csv(df.dst,"data/University_Of_Chicago/agg_dstctry_2019.csv")

##(continuing)##----
df <- read_csv("data/University_Of_Chicago/agg_ctry_2019.csv")

ctry_dict <- read_csv("data/University_Of_Chicago/ctry_dict.csv",locale = readr::locale(encoding = "latin1")) %>% 
  mutate(country=stri_trans_general(str = country, id = "Latin-ASCII"))
head(ctry_dict)

ctry_dict <- ctry_dict %>% select(country,code2,code3) %>% 
  rename(org_ctry=code2,
         org_ctry_code3=code3,
         org_ctry_name=country)
df <- df %>% left_join(ctry_dict)

ctry_dict <- ctry_dict %>%
  rename(dst_ctry=org_ctry,
         dst_ctry_code3=org_ctry_code3,
         dst_ctry_name=org_ctry_name)
df <- df %>% left_join(ctry_dict)

df %>% filter(is.na(org_ctry_name)) %>% select(org_ctry) %>% unique()
df %>% filter(is.na(dst_ctry_name)) %>% select(dst_ctry) %>% unique()
#Only R1 and R2 are unmatched! Let's ignore those for now

owid <- read_csv("data/owid-covid-data.csv")
owid_world <- owid %>% filter(location=="World")
curr_pop <- owid_world$population[1]
owid <- owid %>% select(iso_code,location,continent) %>% unique()

iso_code_continent <- c("OWID_AFR","OWID_ASI","OWID_EUR","OWID_NAM","OWID_OCE","OWID_SAM")
continent <- c("Africa","Asia","Europe","North America","Oceania","South America")
owid_cont <- data.frame(iso_code_continent,continent)
owid <- owid %>% left_join(owid_cont) %>% filter(!is.na(continent))

owid <- owid %>% select(iso_code,iso_code_continent) %>% 
  rename(org_ctry_code3=iso_code,
         org_ctnt_code=iso_code_continent)
df <- df %>% left_join(owid)

owid <- owid %>%
  rename(dst_ctry_code3=org_ctry_code3,
         dst_ctnt_code=org_ctnt_code)
df <- df %>% left_join(owid)

##Aggregating and generating the mobility matrix##----
df.agg <- df %>% filter(!is.na(org_ctnt_code)&!is.na(dst_ctnt_code)) %>% 
  select(org_ctnt_code,dst_ctnt_code,bkgs) %>% 
  arrange(org_ctnt_code,dst_ctnt_code) %>% 
  group_by(org_ctnt_code,dst_ctnt_code) %>% 
  mutate(bkgs=sum(bkgs)) %>% 
  ungroup() %>% 
  unique() %>% 
  filter(org_ctnt_code!=dst_ctnt_code)

df.agg <- df.agg %>% pivot_wider(names_from = dst_ctnt_code, values_from = bkgs)
df.agg[is.na(df.agg)] <- 0
df.agg <- df.agg[,iso_code_continent]/(curr_pop*365)

write_csv(df.agg,"data/cont_mob_matrix.csv")
