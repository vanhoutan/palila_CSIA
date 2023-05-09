#### this script builds the plots for Figs 1-2 in the main text 
#### regarding deriving Beta values from producers
#### and transparently reporting AA changes over known TLs in prey 


library(ggplot2)      # plotting and viz
library(plyr)         # legacy df manip
library(dplyr)        # variable grouping and manipulation
library(reshape)      # legacy df manipulation
library(data.table)   # legacy functions on df 
library(tidyr)        # gathering and spreading
library(zoo)          # roll mean
library(ggthemes)     # helpful ggplot themes
library(forcats)


#### Custom ggPlot theme
themeKV <-theme_few()+
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(c(1, 0.2), unit = "cm")),
        axis.ticks.length=unit(-0.15, "cm"),element_line(colour = "black", linewidth=.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0))


#### read in palila prey data, output from UCDavis SIL
# setwd("/Users/kylevanhoutan/palila_CSIA/")
forage_raw <- read.csv('data/palila_prey/palila_prey.csv')
# subset for just producer data (what's needed to calculate β)
forage_TL1 <- subset(forage_raw, TL == 1)

# Generate random variates given the norm distrib params (ave, sd) given in the lab results
data <- forage_TL1
dataList = split(data, data$ucdavis_id) #split raw data up based on lab specimen ID

#### build functions to generate 100 estimates of each AA for each sample
#### use these to generate 100 estimates of Beta 
#### for each of the 6 formulations of Beta in Besser et al 2022, Figure 5 (DOI: 10.1111/1365-2745.13853)


# set # of random variates drawn from norm distrib, outside of 
num_draws = 100

# start with Beta formulation #1 
GetBeta1 <- function(anID){
#  num_of_draws = 100
  
ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 

B_glu.phe <- (glu.est) - (phe.est)  # first of 6 functions to calculate a formulation of β
return(B_glu.phe)
}

Beta1 = lapply(dataList, GetBeta1)        # apply the 1st Beta calc function to each specimen
Beta1 = reshape2::melt(Beta1)             # melt the dataframe
Beta1$L1 <- as.factor(Beta1$L1)                              
setnames(Beta1, old=c("L1","value"), new=c("ucdavis_id", "Glu-Phe"))    # replace Beta value column header with Beta formula 
Beta1 <- Beta1[, c(2, 1)]    # reorder colummns


# run Beta function 2/6
GetBeta2 <- function(anID){

  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 

B_pro.phe <- (pro.est) - (phe.est)  # 2/6 functions to calculate a formulation of β
  return(B_pro.phe)
}

Beta2 = lapply(dataList, GetBeta2)        # apply the 1st Beta calc function to each specimen
Beta2 = reshape2::melt(Beta2)             # melt the dataframe
Beta2$L1 <- as.factor(Beta2$L1)                              
setnames(Beta2, old=c("L1","value"), new=c("ucdavis_id", "Pro-Phe"))    # replace Beta value column header with Beta formula 


# run Beta function 3/6
GetBeta3 <- function(anID){
  
  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 
  
  B_thr.phe <- (thr.est) - (phe.est)  # 3/6 functions to calculate a formulation of β
  return(B_thr.phe)
}

Beta3 = lapply(dataList, GetBeta3)        # apply the 1st Beta calc function to each specimen
Beta3 = reshape2::melt(Beta3)             # melt the dataframe
Beta3$L1 <- as.factor(Beta3$L1)                              
setnames(Beta3, old=c("L1","value"), new=c("ucdavis_id", "Thr-Phe"))    # replace Beta value column header with Beta formula 



# run Beta function 4/6
GetBeta4 <- function(anID){
  
  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 
  
  B_glu.lys <- (glu.est) - (lys.est)  # 4/6 functions to calculate a formulation of β
  return(B_glu.lys)
}

Beta4 = lapply(dataList, GetBeta4)        # apply the 1st Beta calc function to each specimen
Beta4 = reshape2::melt(Beta4)             # melt the dataframe
Beta4$L1 <- as.factor(Beta4$L1)                              
setnames(Beta4, old=c("L1","value"), new=c("ucdavis_id", "Glu-Lys"))    # replace Beta value column header with Beta formula 


# run Beta function 5/6
GetBeta5 <- function(anID){
  
  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 
  
  B_pro.lys <- (pro.est) - (lys.est)  # 5/6 functions to calculate a formulation of β
  return(B_pro.lys)
}

Beta5 = lapply(dataList, GetBeta5)        # apply the 1st Beta calc function to each specimen
Beta5 = reshape2::melt(Beta5)             # melt the dataframe
Beta5$L1 <- as.factor(Beta5$L1)                              
setnames(Beta5, old=c("L1","value"), new=c("ucdavis_id", "Pro-Lys"))    # replace Beta value column header with Beta formula 


# run Beta function 6/6
GetBeta6 <- function(anID){
  
  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 
  
  B_thr.lys <- (thr.est) - (lys.est)  # 6/6 functions to calculate a formulation of β
  return(B_thr.lys)
}

Beta6 = lapply(dataList, GetBeta6)        # apply the 1st Beta calc function to each specimen
Beta6 = reshape2::melt(Beta6)             # melt the dataframe
Beta6$L1 <- as.factor(Beta6$L1)                              
setnames(Beta6, old=c("L1","value"), new=c("ucdavis_id", "Thr-Lys"))    # replace Beta value column header with Beta formula 

# bind all the Beta formulations into one single frame
Beta_tot <- cbind(Beta1, Beta2$`Pro-Phe`, Beta3$`Thr-Phe`, Beta4$`Glu-Lys`, Beta5$`Pro-Lys`, Beta6$`Thr-Lys`)
names(Beta_tot)    # check on col names, prob should rename
setnames(Beta_tot, old=c("Beta2$`Pro-Phe`", "Beta3$`Thr-Phe`", "Beta4$`Glu-Lys`", "Beta5$`Pro-Lys`", "Beta6$`Thr-Lys`"), 
         new=c("Pro-Phe", "Thr-Phe", "Glu-Lys", "Pro-Lys", "Thr-Lys"))   # first col is already correctly named from cbind
# what I have been waiting for!!!
# reshape df from wide to long using gather()
Beta_gather <- gather(Beta_tot, key="beta", value="value", 2:7)

# bring in metadata from original data sheet
# remove the duplicate rows for SD
data<-subset(data, value == "ave")
names(data)    # check on col names for cleaning final df
# remove unnecessary rows for the AAs, etc
data_sm<-subset(data, select = -c(value, ala, asp, glu, gly, lys, leu, phe, pro, val, ser, thr))
Beta_total<-merge(data_sm,Beta_gather, by="ucdavis_id")

#### make a boxplot to compare results of different Beta formulations


ggplot(Beta_total, aes(x = photo, y = value)) +
  themeKV +
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.ticks.length = unit(-.15, "cm"), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(-2,-2,-2,-2)),
        axis.text.y = element_text(hjust = 1, margin = margin(10, 10, 10, 10))) +
  geom_point(alpha=0.02, color="black", position="jitter", shape = 16, size = 3) +  
  geom_boxplot(alpha=0, colour = "black", linewidth = 0.25) +
  scale_y_continuous(limits = c(-29,10),
                     breaks = c(-30, -25, -20, -15, -10, -5, 0, 5, 10)) +
  ylab("β (‰)") +
  facet_wrap(~beta, ncol=6)


