#load libraries
library(brms)         # stan wrapper for building bayesian regression models easily in r
library(ggplot2)      # plotting and viz
library(plyr)         # legacy df manip
library(dplyr)        # variable grouping and manipulation
library(reshape)      # legacy df manipulation
library(data.table)   # legacy functions on df 
library(itsmr)        # time series
library(tidyr)        # gathering and spreading
library(lubridate)    # date handling
library(zoo)          # roll mean
library(pdp)

##############################
###  Custom ggPlot theme   ###
##############################
themeo <- ggplot2::theme_classic()+
  ggplot2::theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(margin = margin( 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(c(1, 0.2), unit = "cm")),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0) )

#read in an merge palila data
# palila UC_davis and kyle ID cross reference table
Palila_meta <- read.delim('data/palila_samples/palila_sample_metadata_MBA_UCDavis_IDs.txt')
# UC davis stable isotope data
CSSIAPalila <- read.csv('data/palila_samples/Palila_CSSIA_Aug20.csv')
# location, individual, date, metadata
addt_meta <- read.csv('data/palila_samples/Palila_specimens_feather_database.csv')

# turn the palila ID sheet in to a mergeable file
Palila_meta <- Palila_meta %>% 
  dplyr::select( SIF.Internal.ID,ID_NO,YEAR, MONTH) %>% 
  mutate(ucdavis_id = SIF.Internal.ID, acession_id = ID_NO, year = YEAR) %>% 
  dplyr::select(-SIF.Internal.ID,-ID_NO, - YEAR)

# Join meta with CSSIA data by ucdavis_id
CSSIAPalila <- CSSIAPalila %>% 
  left_join(Palila_meta, by = "ucdavis_id") %>% 
  mutate(accession_id = acession_id.y, year = year.y) %>% 
  dplyr::select(-acession_id.x,-acession_id.y,- year.y,-year.x)

# bring additional metadata from other sheet... and join again below
addt_meta <- addt_meta %>% 
  dplyr::select(ID_NO, SEX, AGE, SPEC_NOTES,SPEC_PREP, LOCATION,REGION) %>% 
  mutate(accession_id = ID_NO) %>% dplyr::select(-ID_NO)

# join
CSSIAPalila <- CSSIAPalila %>% left_join(addt_meta, by = "accession_id")

# count number of samples per year
CSSIAPalila %>% filter(value == "ave") %>% group_by(year) %>% dplyr::summarise(n())

# clean up post joining
rm(addt_meta,Palila_meta)

# Simulating additional observations given the SD defined in the lab results
data <- CSSIAPalila
dataList = split(data, data$ucdavis_id) #split raw data up based on specimen ID

#build function to generate 1000 TP estimates based on random draws from AA gaussian distributions
GetTP <- function(anID){
  
  ## Chikaraishi has the below constants for marine SIA studies
  ## b = 8.4
  ## TEF = 7.6
  
  # Nielsen TEF
  ## B is ours, empirically derived difference between trophic and source AAs in primary producers
   b =   2.1
   TEF = 5.9
  
   ## num_of_draws = 1000
   num_of_draws = 100
  
  ala.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  glu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamic acid
  leu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  pro.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
 
  thr.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine
  ser.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  phe.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
 
  TPforID <- (((((glu.est)/1)) - ((thr.est+ser.est+phe.est)/3) - b)/TEF ) + 1  # function to estimate TP via source and trophic AA with constants 
  
  return(TPforID)
}

TPforID_all = lapply(dataList, GetTP )          # apply that function to each specimen
TPforID_all = reshape2::melt(TPforID_all)       # melt that dataframe
TPforID_all$L1 <- as.factor(TPforID_all$L1)                              
setnames(TPforID_all, old=c("L1","value"), new=c("ucdavis_id", "tp"))
data<-subset(data, value == "ave")
total<-merge(data,TPforID_all, by="ucdavis_id")

# cleanup post sampling
rm(dataList,data,TPforID_all)

# random sample by year to reduce multiple year bias
total <- total %>% group_by(year) %>% dplyr::sample_n(1000) %>% ungroup()

str(total)

# is there annual phenology over the entire dataset?
ggplot(total,aes(x=MONTH, y= tp, color = year %>% as.factor()))+
  geom_point(position = "jitter", size = .05)+
  themeo

# is there annual phenology just within the year 1991 (n=20) only?
total %>%  
  filter(year == 1991) %>% 
    ggplot(aes(x=MONTH, y= tp, color = year %>% as.factor()))+
      geom_point(position = "jitter", size = .05)+
      themeo


# correct TP based on interannual phenology
mod91 <- lm(tp~MONTH, filter(total, year == 1991) ) # we have full coverage accross the year for a single year in 1991
beta <- mod91$coefficients[2] %>% as.numeric()      # extract beta coeficient from the model fit to 1991 data

# month that we wish to adjust all TP values to
new_month <- 6

# adjust all TPs to June (6 mnth) value using the beta estimated from the 1991 linear model 
total$adj_tp <- beta * (new_month) + (total$tp - beta * total$MONTH)  

# predict accross the entire year to visualize it below
prediction <- predict(mod91,interval = "confidence") %>% as.data.frame() 

# plot 1991 TP phenology with fitted model
ggplot(filter(total, year == 1991) )+
  geom_point(aes(x=MONTH, y= tp), position = "jitter", size = .05)+
  geom_ribbon(aes(x=MONTH, ymin = prediction$lwr, ymax = prediction$upr), alpha = .5)+
  geom_line(aes(x=MONTH,y=prediction$fit))+
  themeo+
  labs(title = "Linear model to correct annual phenology, based on 1991 full year coverage")

## what would the region of practical equivalence (ROPE) be for trophic position changes
# total %>% 
#   group_by(year) %>% 
#   summarise(mean_tp = mean(tp),
#             sd_tp = sd(tp),
#             range = abs(range(tp)[1] - range(tp)[2]),
#             num_of_months = length(unique(MONTH)) )
# 
# total %>% 
#   summarise(mean_tp = mean(tp),
#             sd_tp = sd(tp),
#             range = abs(range(tp)[1] - range(tp)[2]),
#             num_of_months = length(unique(MONTH)) )
# 
# max(prediction$fit) - min(prediction$fit)

# plot the corrected 
ggplot(total)+
  geom_violin(aes(x=year, y= tp, group = year),fill = "#2b8cbe",outlier.shape = NA, width = 7)+
  geom_violin(aes(x=year, y= adj_tp, group = year),color = "#e34a33",fill = "#e34a33", alpha = .4,outlier.shape = NA, width = 5)+
  labs(title = "Raw TL by year vs annual phenology adjusted, red = adjusted, blue = raw")+
  themeo

# assign TP to be the adjusted TP
total$tp <- total$adj_tp; total$adj_tp <- NULL

# cleanup uneeded objects
rm(mod91,beta,new_month,prediction, GetTP)

#########################################################################################################
##### resample subsamples and draw multiple loess models from which an ensemble model is summarized ####
#########################################################################################################
null_df <- NULL                                                                      # empty object to chuck results in to:  


  for( i in 1:100){
  ## for( i in 1:1000){
  new_df <- ddply(total,.(year),function(x) x[sample(nrow(x),1),])                   # sample value from each available year
  tp_est <- loess(tp ~ year , new_df,span = 1, surface = "direct")                   # fit loess model through the series
  tp_predict <- as.data.frame(predict(tp_est,data.frame(year = seq(1890, 2015,1))))  # predict for all years
  # fit loess through the data
  round <- data.frame(tp_predict = as.vector(tp_predict), year = seq(1890, 2015,1) , sim_num = i )
  null_df <- rbind(null_df,round)
}

colnames(null_df)[1] <- "prediction"

# get quantiles from the resampled models
null_df_quant <- null_df %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(upper = quantile(prediction,.025,na.rm = T),   # get 2.5%, 50%, 95% quantile
            median = quantile(prediction,.5,na.rm = T) ,
            lower =  quantile(prediction,.975,na.rm = T))

# put the dates all in the same format
null_df_quant$year <- paste0("1","/","1","/",as.character(null_df_quant$year)) %>% mdy()
total$year <- paste0("1","/","1","/",as.character(total$year)) %>% mdy()

## RF/ICE looking plot with multiple runs, no mean
ggplot()+
  geom_line(data = null_df, aes(x=year,y=prediction,group = sim_num), size = .1)+
  #geom_boxplot(data = total,aes(x=year, y= tp, group = year),outlier.shape = NA)+
  themeo


# b/w version of what ends up becoming Figure 1 time series of palila TP
tp <- ggplot(null_df_quant)+
  geom_boxplot(data = total,aes(x=year, y= tp, group = year),outlier.shape = NA)+
  geom_ribbon( aes(x=year,ymin=lower,ymax=upper), size = .1, alpha = .5)+
  geom_line( aes(x=year,y = median), size = .1)+
  #geom_violin(data = total,aes(x=year, y= tp, group = year), width = 20)+
  #scale_y_continuous(limits = c(1,2))+
  #scale_x_date(expand = c(0,0), breaks =  c(as.Date(0)))+
  scale_x_date(expand = c(0,0))+
  labs(x=NULL)+
  themeo
tp

# what ends up becoming Figure 1 time series of palila TP
ggplot(null_df_quant)+
  geom_ribbon( aes(x=year,ymin=lower,ymax=upper), size = .1, alpha = .5)+
  geom_line( aes(x=year,y = median), color = "#ffde17", size = 1.5)+
  #geom_violin(data = total,aes(x=year, y= tp, group = year), width = 20)+
  #scale_y_continuous(limits = c(1.5,2.5))+
  scale_x_date(expand = c(0,0))+
  labs(x=NULL)+
  themeo

# plotting the traverse that Palila takes across constant TL of prey items
# could be interesting supplemental figure that shows 3 major forage item groups

source_csv <- read.csv("./data/mixing model/source_data.csv") # estimated prey item TL, dev in prey_TL.R saved for mixing model as .csv
source_csv_plot <- rbind(source_csv,source_csv)  
source_csv_plot$year <- NULL
source_csv_plot$year[1:(nrow(source_csv_plot)/2)] <- '1/1/1890'  %>% mdy()
source_csv_plot$year[(nrow(source_csv_plot)/2 + 1):nrow(source_csv_plot)] <- '1/1/2015'  %>% mdy()
source_csv_plot$ymin <- source_csv_plot$MeanTL - source_csv_plot$SDTL*1.96 
source_csv_plot$ymax <- source_csv_plot$MeanTL + source_csv_plot$SDTL*1.96 
source_csv_plot$year <- as.Date(source_csv_plot$year)

tp+
  geom_hline(data = source_csv_plot, aes(yintercept = MeanTL, color = spp), alpha = .7, size = 2)+
  scale_x_date(expand = c(0,0))+
  scale_color_manual(values = c("#faa61a","#008745","#3a60ac"))+
  scale_fill_manual(values = c("#faa61a","#008745","#3a60ac"))+
  geom_ribbon(data = source_csv_plot, aes(ymin = ymin, ymax = ymax, x = year , fill = spp), alpha = .7)+
  geom_ribbon( aes(x=year,ymin=lower,ymax=upper), size = .1, alpha = .5)+
  geom_line( aes(x=year,y = median), size = .1)
  

# Create time series of TP in all Palila through time for mixing model
# Consider breaking up by sex later
 mix_csv <- null_df %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(
    TL = mean(prediction),
    sd = sd(prediction),
    spp = "PALI")
# write.csv(mix_csv,"data/mixing model/mixture_data.csv")

#cleanup
rm(new_df,null_df,i,tp_est, tp_predict,round,mix_csv, mod91, new_month,prediction, GetTP, source_csv_plot,source_csv)


# Validating what amino acids are source vs trophic
CSSIAPalila_sub <- CSSIAPalila %>% 
  filter(value == "ave") %>% 
  mutate(avg_source = (ser + thr + phe) / 3)

tall_cssia <- CSSIAPalila_sub %>% gather(key = AA, value = AA_val, ala,asp,glu,gly,leu,phe,pro,val,ser,thr,avg_source )
str(tall_cssia)

time<- ggplot(tall_cssia,aes(x=year,y=AA_val,group = AA,color = AA))+
  geom_point(show.legend = F, size = .5)+
  geom_smooth(show.legend = F, size = .5)+
  scale_color_brewer(palette = "Spectral")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  themeo+
  labs(title = "Amino acid N15 through time")

space <- ggplot(tall_cssia,aes(x=AA_val,group = AA,fill = AA, color = AA))+
  geom_density(alpha = .5)+
  #geom_histogram(alpha = .5, bins = 100)+
  scale_fill_brewer(palette = "Spectral")+
  scale_color_brewer(palette = "Spectral")+
  scale_x_continuous(expand = c(0,0))+
  #scale_y_continuous(expand = c(0,0), limits = c(0,.53))+
  themeo

gridExtra::grid.arrange(time,space, layout_matrix = c(1,1,2) %>% as.matrix())
rm(time,space,CSSIAPalila_sub, tall_cssia)


# drought index attempts
SPEI <- read.csv('./data/environ_covars/spei_-155.25_19.75.csv',header = T, skip = 1)
SPEI$year <- seq.Date(from = mdy('01/01/1901'), by = "month" , length.out = nrow(SPEI)) #seq(1:nrow(SPEI))
str(SPEI)

# units is months since 1901-01
# range: January 1901, December 2011
drought <- ggplot(SPEI)+
  geom_line(aes(x = year, y = SPEI01),size = .2, color = "light gray")+
  geom_area(aes(x = year, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
  geom_hline(yintercept = 0)+
  themeo
drought

## evapotranspiration plot over the whole dataset
drought <- ggplot(SPEI)+
  geom_line(aes(x = ((year)), y = SPEI01),size = .1, color = "gray")+
  geom_bar(stat = "identity", aes(x = ((year)), y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
  #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
  #scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
  scale_fill_manual(values = c("tomato3","cadetblue3"))+
  geom_hline(yintercept = 0)+
  #scale_x_date(expand = c(0,0), breaks = c(as.Date(0)))+
  scale_x_date(expand = c(0,0))+
  labs(x = NULL)+
  themeo
drought

# temp
# two cells over relavant areas
# temp <- read.table('https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.01/ge/grid05/cells/N17.5W157.5/tmp/N19.75W155.75.tmp.txt', header = T, skip = 4)
# closer to mauna kea/loa


## Modified script to reflect new data URL link through Google Earth kml 5degree grid
## lines 339-352 curently not working

temp2 <- read.table('https://crudata.uea.ac.uk/cru/data/crutem/ge/crutem4-2019-12/N17.5W157.5_data.txt', header = T, skip = 8)
## print(temp2)

# modify temp df
colnames(temp2) <- c("year","month","degC")
temp2$date <- paste0(as.character(temp2$month),"/","1","/",as.character(temp2$year)) %>% mdy()
temp2$rollmean<- rollmean(temp2$degC,k = 12, na.pad = T, align = "right")

temp <- ggplot(temp2)+
  geom_line(aes(date,degC),size = .2, color = "gray")+
  geom_line(aes(date, rollmean))+
  scale_y_continuous(expand = c(-0.3,0))+
  #scale_x_date(expand = c(0,0), breaks = c(as.Date(0)))+
  scale_x_date(expand = c(0,0))+
  labs(x = NULL)+
  themeo

gridExtra::grid.arrange(tp,drought,temp)

 # library(gganimate)
 # ggplot(temp2,aes(x = month, y = degC, group = year, color = year))+
 #  geom_line()+themeo+scale_x_continuous(expand = c(0,0))+scale_color_distiller(palette = "RdBu")+
 #  transition_time(year) +shadow_mark(past = T) +ease_aes('linear')

# merging predictors and response variables in to dataframe. 
temp2$year <- temp2$date
SPEI <- SPEI[,c("year","SPEI36")]

empty_yrs <- data.frame(year = seq.Date(from = mdy('01/01/1901'), by = "month" , length.out = nrow(SPEI)))

str(total$year)
str(null_df_quant$year)
str(SPEI$year)
str(temp2$year)

raw_tp <- data.frame(year = total$year,raw_tp = total$tp)

mod_df <- empty_yrs %>% 
  full_join(null_df_quant, by = "year") %>% 
  full_join(temp2,by = "year") %>% 
  full_join(SPEI, by = "year") %>% 
  full_join(raw_tp, by = "year") %>% 
  arrange(year) %>% 
  
  mutate(median = na.approx(median, na.rm = F),
         upper = na.approx(upper, na.rm = F),
         lower = na.approx(lower, na.rm = F))


# generating an indice of caterpillar parasitism rate
intro_date <- mdy('01/01/1950')    # introduction date
known_date1 <- mdy('01/01/1997')   # date(s) for known data

# what years had no parasititism
zero_index <- which(empty_yrs$year %in% c(first(empty_yrs$year),intro_date), arr.ind=TRUE)
# what years were after the introduction
parasite_years_index <- which(empty_yrs$year %in% c(intro_date, last(empty_yrs$year)), arr.ind=TRUE)

# build df based on known introduction date and known rates of parasitism
parasite_data <- data.frame( year = empty_yrs %>% filter(year %in% c(intro_date, known_date1)),
                             parasitism = c(0, .40))

# fit a model to it
para_mod <- lm(parasitism ~ year, data = parasite_data)

# predict accross all years after introduction
sim_range <- empty_yrs %>% filter(year > intro_date)
pred_parasitism <- predict(para_mod, newdata =  sim_range) %>% as.data.frame()

# logistic/sigmoid alternative (needs additoinal estimates compared to lm)
y <- c(parasite_data$parasitism,.4, .47, .51)              # ratios drawn from lit, with possible final parasitation estimate
x <- c(as.numeric(parasite_data$year),9862, 10592,20089)   # relative years, but in numerical format derived via as.numeric( mdy('01/01/2025'))
log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))            # fit logistic model via nls
plot(as.Date(x),y)
# plot the model prediction across all years and into the future
lines(as.Date(min(x):max(x)), predict(log.ss, data.frame(x=min(x):max(x))), col="red")
pred_parasitism <- predict(log.ss,  newdata = data.frame( x= as.numeric(sim_range$year) ) ) %>% as.data.frame()


# bind the pre introduction rate (zero) to the model post introduction
p_rate <- c(
  rep(0,zero_index[[2]]),   # pre introduction
  pred_parasitism$.         # post introduction model
)




wasps <- data.frame( year = empty_yrs, parasitism = p_rate)
plot(wasps$year,wasps$parasitism, ylim = c(0,1), type = "l")

mod_df <- mod_df %>% 
  full_join(wasps, by = "year") %>% 
  arrange(year)

para_gg <- ggplot(mod_df)+
  geom_line(aes(x = year,
                y = parasitism), size = 3, color = "orange2" )+
  labs(x = NULL, y = "parasitism rate")+
  scale_x_date(limits = c(min(temp2$date), max(temp2$date)), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1))+
  themeo

gridExtra::grid.arrange(tp,drought,temp,para_gg, ncol = 1)


# plot to check 
par(mfrow=c(2,3))
plot(mod_df$year, mod_df$raw_tp)
plot(mod_df$year, mod_df$median)
plot(mod_df$year, mod_df$rollmean)
plot(mod_df$year, mod_df$SPEI36)
plot(mod_df$year, mod_df$parasitism)

##### setting up model building dataframe #####
mod_df <- mod_df[complete.cases(mod_df), ]


####################################################################################
     #    Fig 4 - model predicting TP given various environmental drivers    #
####################################################################################

# load librares
library(tidyverse)
library(tidybayes)    # easy means to investigate brms built models in a tidy framework
library(modelr)       # tidy brms dependency
library(strengejacke)

# predicition function for brm pdp plot
pred_fun <- function(object, newdata) {
  mean(predict(object, newdata)[1], na.rm = TRUE)
}

mod_df <- mod_df %>% # center on 0 and sd = 1 of predictors
  mutate(
    SPEI36 = scale(SPEI36),
    rollmean = scale(rollmean),
    parasitism = scale(parasitism)
    #raw_tp = scale(raw_tp)
  )

str(mod_df)

# Fit a standard linear regression with independent error
# potential models of interest
formula1 <- median ~ rollmean + SPEI36 + rollmean*SPEI36
#formula2 <- raw_tp | trunc(lb = 1) ~ rollmean + SPEI36 + rollmean*SPEI36
#formula2 <- raw_tp | trunc(lb = 1) ~ rollmean + SPEI36 + parasitism #+ rollmean*SPEI36 + parasitism*SPEI36
formula2 <- raw_tp | trunc(lb = 1) ~ rollmean + SPEI36 + parasitism + rollmean*SPEI36 + parasitism*SPEI36
formula3 <- median ~ rollmean + SPEI36 
formula4 <- raw_tp ~ rollmean + SPEI36 

# concatenate in to a list for the loop build
formula_list <- list(formula1,formula2,formula3,formula4)

# single model selection
each_formula <- 2

# for(each_formula in 1:length(formula_list)){

get_prior(formula_list[[each_formula]], data = mod_df, family = gaussian())  

data_f <- sample_frac(mod_df,.9) %>% as.data.frame()

data.brms = brm(formula_list[[each_formula]], 
                data = data_f, 
                
                iter = 10000, 
                warmup = 500, 
                chains = 4,
                thin = 5, 
                
                #cores = 5,
                #control = list(adapt_delta =.9, max_treedepth = 15),
                
                prior = c(prior(normal(0, 10), class = "Intercept"),
                          prior(normal(0, 1), class = "b")#,
                          #prior(normal(0, 1), class = "sigma")
                ),
                sample_prior = T,
                silent = T,
                #family = lognormal()
                family = gaussian()
)

data.brms <- add_waic(data.brms)
print(data.brms$waic)

# Useful diagnostics
# summary(data.brms)
# plot(marginal_effects(data.brms), points = T)
# plot(data.brms)
# pairs(data.brms)
# bayes_R2(data.brms)
# launch_shinystan(data.brms)

# according to a "small effect" Cohen 1988
#We might want to say that a regression coecient x on predictor x is practically equivalent to zero if a
#change across the \main range of x" produces only a negligible change in the predicted value ^y.
equi_test(x = data.brms, rope = c(-sd(mod_df$raw_tp)*.2,sd(mod_df$raw_tp)*.2), out = "plot")+theme_bw()+coord_flip()

library("bayesplot")
color_scheme_set("mix-blue-red")
coef_plot <- mcmc_areas(data.brms %>% as.array(), pars = c("b_rollmean", "b_SPEI36", "b_rollmean:SPEI36", "b_SPEI36:parasitism", "sigma"))
color_scheme_set("viridis")
tracer <- mcmc_trace(data.brms %>% as.array(), pars = c("b_rollmean", "b_SPEI36", "b_rollmean:SPEI36", "b_SPEI36:parasitism", "sigma"), 
                     facet_args = list(ncol = 1, strip.position = "left"))

grid.arrange(coef_plot, tracer, ncol = 2)


# posterior plots
plot(hypothesis(data.brms, "rollmean > 0"))
plot(hypothesis(data.brms, "SPEI36 > 0"))
plot(hypothesis(data.brms, "rollmean:SPEI36 > 0"))

library(RColorBrewer)

# margional effects plots with prediction intervals

SPEI36_bayes <- mod_df %>%
  data_grid(SPEI36 = seq_range(SPEI36, n = 10), rollmean = seq_range(rollmean, n = 10), parasitism = seq_range(parasitism, n = 10)) %>%
  tidybayes::add_predicted_draws(data.brms) %>%
  ggplot(aes(x = SPEI36, y = raw_tp)) +
  tidybayes::stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "black", show.legend = F) +
  geom_jitter(data = mod_df, size = .1, width = .02) +
  scale_x_continuous(expand = c(0,0))+
  themeo+
  scale_fill_brewer()

temp_bayes <- mod_df %>%
  data_grid(SPEI36 = seq_range(SPEI36, n = 10), rollmean = seq_range(rollmean, n = 10), parasitism = seq_range(parasitism, n = 10)) %>%
  tidybayes::add_predicted_draws(data.brms) %>%
  ggplot(aes(x = rollmean, y = raw_tp)) +
  tidybayes::stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "black", show.legend = F)+
  geom_jitter(data = mod_df, size = .1, width = .02)+
  scale_x_continuous(expand = c(0,0))+
  themeo+
  scale_fill_brewer()

parasitism_bayes  <-  mod_df %>%
  data_grid(SPEI36 = seq_range(SPEI36, n = 10), rollmean = seq_range(rollmean, n = 10), parasitism = seq_range(parasitism, n = 10)) %>%
  tidybayes::add_predicted_draws(data.brms) %>%
  ggplot(aes(x = parasitism, y = raw_tp)) +
  tidybayes::stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "black", show.legend = F) +
  geom_jitter(data = mod_df, size = .1, width = .02)+ 
  scale_x_continuous(expand = c(0,0))+
  themeo+
  scale_fill_brewer()

gridExtra::grid.arrange(SPEI36_bayes,parasitism_bayes,temp_bayes)



# 
# 
# 
# # quick randomforest
# mod_rf <- randomForest::randomForest(median ~ rollmean + SPEI36, data = mod_df)
# mod_rf <- randomForest::randomForest(raw_tp ~ rollmean + SPEI36, data = mod_df)
# mod_rf <- randomForest::randomForest(raw_tp ~ rollmean + SPEI36 + parasitism, data = mod_df)
# 
# 
# mod_rf
# 
# partial(mod_rf, pred.var = "rollmean", plot = T)
# partial(mod_rf, pred.var = "SPEI36", plot = T)
# partial(mod_rf, pred.var =  c("SPEI36", "rollmean"), plot = T, chull = T)
# 
# 
# df <- partial(mod_rf, pred.var =  c("SPEI36", "rollmean"), chull = T, grid.resolution = 10)
# 
# ggplot()+
#   geom_tile(data = df,aes(x = SPEI36, y = rollmean, fill = yhat), inherit.aes = F)+
#   geom_point(data = mod_df,aes(x = SPEI36, y = rollmean, fill = raw_tp), pch = 21, size = 3)+
#   scale_fill_distiller(palette = "Spectral", direction = 1)+
#   themeo
# 
# 
# data_z <- reshape2::acast(df, SPEI36~rollmean, value.var = "yhat")
# plotly::plot_ly(z = data_z,  type = "surface", colors = "Spectral")
# 
# plotPartial(partial(mod_rf, pred.var =  c("SPEI36", "rollmean")), levelplot = F, zlab = "median", colorkey = TRUE, 
#             screen = list(z = 120, x = -60)  )
# 
# 
# 
# 
# # Gradient boosting exploration
# # Load required packages
# library(xgboost)
# 
# # Load the data
# X <- mod_df[,c("rollmean","SPEI36")]
# #y <- mod_df[,"median"]
# y <- mod_df[,"raw_tp"]
# 
# pima.xgb <- lm( raw_tp ~ rollmean + SPEI36, data = mod_df)
# pima.xgb <- lm( raw_tp ~ rollmean + SPEI36 +  rollmean*SPEI36, data = mod_df)
# 
# # Parameters for XGBoost model
# param.list <- list(max_depth = 5, eta = 0.01, objective = "reg:linear")
# 
# # Fit an XGBoost model
# set.seed(101)
# pima.xgb <- xgb.train(params = param.list, 
#                       data = xgb.DMatrix(data.matrix(X), 
#                                          label = y), 
#                       nrounds = 500)
# 
# pima.xgb <- xgb.cv(params = param.list, 
#                       data = xgb.DMatrix(data.matrix(X), label = y), 
#                    nfold = 5,
#                    nthread = 4,
#                    nrounds = 500)
# pima.xgb
# 
# xgb.plot.importance(xgb.importance(model = pima.xgb))
# xgb.plot.deepness(model = pima.xgb)
# 
# # Partial of diabetes test result on glucoe
# partial(pima.xgb, pred.var = "rollmean", plot = TRUE, train = X, rug = T)
# partial(pima.xgb, pred.var = "SPEI36", plot = TRUE, train = X)
# partial(pima.xgb, pred.var =  c("SPEI36", "rollmean"), plot = T, chull = T, train = X, progress = "text")
# 
# plotPartial(partial(pima.xgb, pred.var =  c("SPEI36", "rollmean"), train = X, chull = T, progress = "text"), levelplot = FALSE, colorkey = TRUE, 
#             screen = list(z = 120, x = -60))
# 
# 
# 
# #plotly 3D surface
# library(plotly)
# 
# surface <- partial(mod_rf, pred.var =  c("SPEI36", "rollmean"), train = X, progress = "text", chull = T, grid.resolution = 40)
# surface <- partial(pima.xgb, pred.var =  c("SPEI36", "rollmean"), train = X, progress = "text", chull = T, grid.resolution = 40)
# str(surface)
# 
# library(reshape2)
# data_z <- acast(surface, SPEI36~rollmean, value.var = "yhat")
# #plot_ly(z = data_z,  type = "surface", colors = "Spectral")
# str(data_z)
# 
# x <- attr(data_z,"dimnames")[1] %>% as.data.frame(); x <- as.vector(x[,1]); x <- as.numeric(x)
# y <- attr(data_z,"dimnames")[2] %>% as.data.frame(); y <- as.vector(y[,1]); y <- as.numeric(y)
# 
# # Color palette (100 colors)
# col.pal<-colorRampPalette(c("yellow", "dark green"))
# colors<-col.pal(100)
# # height of facets
# z = data_z
# z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# # Range of the facet center on a 100-scale (number of colors)
# z.facet.range<-cut(z.facet.center, 100)
# 
# colors = colors[z.facet.range]
# 
# dev.new()
# seq360 <- rbind(seq(-180,180,by = .25),seq(-180,180,by = .25))
# 
# for(round in 1:720){
#   
# persp( x = x, y, z = data_z, 
#        col = colors,
#        #col = "green",
#        shade = 0.4,
#        expand = 0.4,
#        phi = 30 , #theta = -60
#        lwd = .1,
#        zlab = "trophic level",
#        xlab = "drought index",
#        ylab = "mean air temperature",
#        box = T
#        ,theta = seq360[round]
#        )
# 
#   }
# 
# 
# 
# ggplot(mod_df,aes(x=year %>% as.numeric()))+
#   geom_line(aes(y = SPEI36),size = .1, color = "gray")+
#   geom_bar(stat = "identity", aes( y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
#   #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
#   scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
#   geom_hline(yintercept = 0)+
#   scale_x_time(expand = c(0,0))+
#   themeo
# 
# drought<-ggplot(mod_df,aes(x=year))+
#   geom_line(aes(y = SPEI36),size =1, color = "gray")+
#   #geom_bar(stat = "identity", aes( y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
#   #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
#   scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
#   geom_hline(yintercept = 0)+
#   scale_x_date(expand = c(0,0))+
#   themeo
# 
# temp<-ggplot(mod_df,aes(x=year ))+
#   geom_line(aes( y = rollmean),size =1, color = "gray")+
#   #geom_bar(stat = "identity", aes( y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
#   #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
#   #scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
#   #geom_hline(yintercept = 0)+
#   scale_x_date(expand = c(0,0))+
#   themeo
# 
# raw_tp <-ggplot(mod_df,aes(x=year))+
#   geom_line(aes( y = raw_tp ),size =1, color = "gray")+
#   #geom_bar(stat = "identity", aes( y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
#   #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
#   #scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
#   #geom_hline(yintercept = 0)+
#   scale_x_date(expand = c(0,0))+
#   themeo
# 
# smooth_tp<-ggplot(mod_df,aes(x=year))+
#   geom_line(aes( y = median ),size =1, color = "gray")+
#   #geom_bar(stat = "identity", aes( y = SPEI36, fill = ifelse(SPEI36>0,"red","green")), show.legend = F)+
#   #geom_area(aes(x = DateSeq, y = SPEI48), size = .2, fill = "dark gray", color = "black")+
#   #scale_fill_manual(values = c("#4d9221","#d8b365") %>% rev())+
#   #geom_hline(yintercept = 0)+
#   scale_x_date(expand = c(0,0))+
#   themeo
# 
# grid.arrange(smooth_tp,raw_tp,drought,temp, nrow = 1)
# 
# 
# str(mod_df)
# 












  