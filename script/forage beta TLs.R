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
library(scales)
library(RColorBrewer) # pretty colors
library(ggdist)
library(colorspace)
library(ragg)
library(grid)
library(png)
library(patchwork) 

# Custom ggPlot theme
themeKV <- theme_few()+
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(colour = "black", margin = margin(0.2, unit = "cm")),
        axis.text.y = element_text(colour = "black", margin = margin(c(1, 0.2), unit = "cm")),
        axis.ticks.x = element_line(colour = "black"), axis.ticks.y = element_line(colour = "black"),
        axis.ticks.length=unit(-0.15, "cm"),element_line(colour = "black", linewidth=.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0, size=8))


#### Determine which trophic and source AA combinations are suitable 
#### from the CSIA-AA d15N data in this terrestrial setting

# read in palila prey data, output from UC Davis SIL
# setwd("/Users/kylevanhoutan/palila_CSIA/")
forage_raw <- read.csv('data/palila_prey/palila_prey.csv')
# subset for just producer data (what's needed to calculate β)
forage_TL1 <- subset(forage_raw, TL == 1)

# Generate random variates given the norm distrib params (ave, sd) given in the lab results
data <- forage_TL1
dataList = split(data, data$ucdavis_id) #split raw data up based on lab specimen ID

#### build functions to generate X estimates of each AA for each sample
#### use these to generate X estimates of Beta 
#### for each of the 6 formulations of Beta in Besser et al 2022, Figure 5 (DOI: 10.1111/1365-2745.13853)

# set # of random variates drawn from norm distrib, outside of the function to avoid repetition
num_draws = 37

# start with Beta formulation #1 
GetBeta1 <- function(anID){
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
B_glu.phe <- (glu.est) - (phe.est)  # first of 6 functions to calculate a formulation of β
return(B_glu.phe)
}
Beta1 = lapply(dataList, GetBeta1)        # apply the 1st Beta calc function to each specimen
Beta1 = reshape2::melt(Beta1)             # melt the dataframe
Beta1$L1 <- as.factor(Beta1$L1)                              
setnames(Beta1, old=c("L1","value"), new=c("ucdavis_id", "Glu-Phe"))    # replace Beta value column header with Beta formula 
Beta1 <- Beta1[, c(2, 1)]    # reorder columns

# run Beta function 2/6
GetBeta2 <- function(anID){
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
B_pro.phe <- (pro.est) - (phe.est)  # 2/6 functions to calculate a formulation of β
  return(B_pro.phe)
}
Beta2 = lapply(dataList, GetBeta2)        # apply the 1st Beta calc function to each specimen
Beta2 = reshape2::melt(Beta2)             # melt the dataframe
Beta2$L1 <- as.factor(Beta2$L1)                              
setnames(Beta2, old=c("L1","value"), new=c("ucdavis_id", "Pro-Phe"))    # replace Beta value column header with Beta formula 

# run Beta function 3/6
GetBeta3 <- function(anID){
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine 
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  B_thr.phe <- (thr.est) - (phe.est)  # 3/6 functions to calculate a formulation of β
  return(B_thr.phe)
}
Beta3 = lapply(dataList, GetBeta3)        # apply the 1st Beta calc function to each specimen
Beta3 = reshape2::melt(Beta3)             # melt the dataframe
Beta3$L1 <- as.factor(Beta3$L1)                              
setnames(Beta3, old=c("L1","value"), new=c("ucdavis_id", "Thr-Phe"))    # replace Beta value column header with Beta formula 

# run Beta function 4/6
GetBeta4 <- function(anID){
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  B_glu.lys <- (glu.est) - (lys.est)  # 4/6 functions to calculate a formulation of β
  return(B_glu.lys)
}
Beta4 = lapply(dataList, GetBeta4)        # apply the 1st Beta calc function to each specimen
Beta4 = reshape2::melt(Beta4)             # melt the dataframe
Beta4$L1 <- as.factor(Beta4$L1)                              
setnames(Beta4, old=c("L1","value"), new=c("ucdavis_id", "Glu-Lys"))    # replace Beta value column header with Beta formula 
# check the average value, evenly weighted between mamane and naio, should ~1.52
mean(Beta4$`Glu-Lys`) # 1.499065

# run Beta function 5/6
GetBeta5 <- function(anID){
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  B_pro.lys <- (pro.est) - (lys.est)  # 5/6 functions to calculate a formulation of β
  return(B_pro.lys)
}
Beta5 = lapply(dataList, GetBeta5)        # apply the 1st Beta calc function to each specimen
Beta5 = reshape2::melt(Beta5)             # melt the dataframe
Beta5$L1 <- as.factor(Beta5$L1)                              
setnames(Beta5, old=c("L1","value"), new=c("ucdavis_id", "Pro-Lys"))    # replace Beta value column header with Beta formula 

# run Beta function 6/6
GetBeta6 <- function(anID){
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # lysine
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
names(data)    # check on col names for orderly cleaning of final df
# remove unnecessary cols for the AAs, etc
data_sm<-subset(data, select = -c(value, ala, asp, glu, gly, lys, leu, phe, pro, val, ser, thr))
# finally here... so good
Beta_total<-merge(data_sm,Beta_gather, by="ucdavis_id")

#### make a box plot to compare results of different Beta formulations
# similar form to Fig 5 from Besser et al 2022
# but scale free y axes to emphasize whatever structure is there
p1 <- ggplot(Beta_total, aes(x = photo, y = value, fill = photo)) +
  themeKV +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position="bottom",
        axis.text.y = element_text(size = 7, margin = margin(c(1, 0.1), unit = "cm")),
        axis.title.y = element_text(size = 10),
        legend.key.height = unit(0.5, 'cm'), # shrink the native height of legend
        legend.text = element_text(size=7), # reduce font size on legend
        panel.spacing.x = unit(0, "lines")) + # move panels closer horiz
  geom_point(alpha=0.05, color="black", position="jitter", shape = 16, size = 3) +  
  geom_boxplot(alpha=0.9, colour = "black", linewidth = 0.25, outlier.color = NA) +
  scale_fill_manual(values=c("#5e4fa2", "#66c2a5")) +
  scale_y_continuous(breaks= pretty_breaks()) +
  ylab("Beta (d15N %)") + # this needs to read β(δ15N ‰) but it wont print to PDF
  facet_wrap(~beta, ncol=6, scales = "free_y")
  #facet_wrap(~beta, ncol=6)



#### Do a new analysis to pull random variates from AA distributions
#### to display and help determine which are intuitive trophic and source AAs
#### which will lead to calculating TEF (i.e., TDF)
#### group results by TL (all results from TL=1, TL=2...)
#### make box plot, facet by AA 

# first read in palila prey data, again, output from UCDavis SIL
# setwd("/Users/kylevanhoutan/palila_CSIA/")
forage_TL123 <- read.csv('data/palila_prey/palila_prey.csv')

# Generate random variates given the norm distrib params (ave, sd) in the lab results
data2 <- forage_TL123
data2List = split(data2, data2$ucdavis_id) #split raw data up based on lab specimen ID

#### build functions to generate X estimates of each AA for each sample in the forage data
#### use these to generate X estimates of each AA d15N value 
#### at each trophic level

# set # of random variates drawn from norm distrib, outside of the function
num_draws = 37

# Alanine, first AA function, 1 of 11
GetAla <- function(anID){
  ala.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  return(ala.est)
}
Ala = lapply(data2List, GetAla)        # apply function to the first AA using data from each specimen
Ala = reshape2::melt(Ala)             # melt the df
Ala$L1 <- as.factor(Ala$L1)
setnames(Ala, old=c("L1","value"), new=c("ucdavis_id", "Ala"))    # replace column header names 
Ala <- Ala[, c(2, 1)]    # reorder columns

# AA function 2 of 11, Asp
GetAsp <- function(anID){
  asp.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "asp")) , sd = as.numeric(subset(anID, value == "sd", select = "asp"))) # asparagine
  return(asp.est)
}
Asp = lapply(data2List, GetAsp)
Asp = reshape2::melt(Asp)
Asp$L1 <- as.factor(Asp$L1)
setnames(Asp, old=c("L1","value"), new=c("ucdavis_id", "Asp"))

# AA function 3 of 11, Glu
GetGlu <- function(anID){
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamine
  return(glu.est)
}
Glu = lapply(data2List, GetGlu)
Glu = reshape2::melt(Glu)
Glu$L1 <- as.factor(Glu$L1)
setnames(Glu, old=c("L1","value"), new=c("ucdavis_id", "Glu"))

# AA function 4 of 11, Gly
GetGly <- function(anID){
  gly.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "gly")) , sd = as.numeric(subset(anID, value == "sd", select = "gly"))) # glycine
  return(gly.est)
}
Gly = lapply(data2List, GetGly)
Gly = reshape2::melt(Gly)
Gly$L1 <- as.factor(Gly$L1)
setnames(Gly, old=c("L1","value"), new=c("ucdavis_id", "Gly"))

# AA function 5 of 11, Lys
GetLys <- function(anID){
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys"))) # Lysine
  return(lys.est)
}
Lys = lapply(data2List, GetLys)
Lys = reshape2::melt(Lys)
Lys$L1 <- as.factor(Lys$L1)
setnames(Lys, old=c("L1","value"), new=c("ucdavis_id", "Lys"))

# AA function 6 of 11, Leu
GetLeu <- function(anID){
  leu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # Leucine
  return(leu.est)
}
Leu = lapply(data2List, GetLeu)
Leu = reshape2::melt(Leu)
Leu$L1 <- as.factor(Leu$L1)
setnames(Leu, old=c("L1","value"), new=c("ucdavis_id", "Leu"))

# AA function 7 of 11, Phe
GetPhe <- function(anID){
  phe.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # phenylalanine
  return(phe.est)
}
Phe = lapply(data2List, GetPhe)
Phe = reshape2::melt(Phe)
Phe$L1 <- as.factor(Phe$L1)
setnames(Phe, old=c("L1","value"), new=c("ucdavis_id", "Phe"))

# AA function 8 of 11, Pro
GetPro <- function(anID){
  pro.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # proline
  return(pro.est)
}
Pro = lapply(data2List, GetPro)
Pro = reshape2::melt(Pro)
Pro$L1 <- as.factor(Pro$L1)
setnames(Pro, old=c("L1","value"), new=c("ucdavis_id", "Pro"))

# AA function 9 of 11, Val
GetVal <- function(anID){
  val.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # valine
  return(val.est)
}
Val = lapply(data2List, GetVal)
Val = reshape2::melt(Val)
Val$L1 <- as.factor(Val$L1)
setnames(Val, old=c("L1","value"), new=c("ucdavis_id", "Val"))

# AA function 10 of 11, Ser
GetSer <- function(anID){
  ser.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # serine
  return(ser.est)
}
Ser = lapply(data2List, GetSer)
Ser = reshape2::melt(Ser)
Ser$L1 <- as.factor(Ser$L1)
setnames(Ser, old=c("L1","value"), new=c("ucdavis_id", "Ser"))

# AA function 11 of 11, Thr
GetThr <- function(anID){
  thr.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # threonine
  return(thr.est)
}
Thr = lapply(data2List, GetThr)
Thr = reshape2::melt(Thr)
Thr$L1 <- as.factor(Thr$L1)
setnames(Thr, old=c("L1","value"), new=c("ucdavis_id", "Thr"))


#### bind all the AA random deviates into one single frame
AA_tot <- cbind(Ala, Asp$Asp,	Glu$Glu, Gly$Gly, Lys$Lys, Leu$Leu,	Phe$Phe, Pro$Pro,	Val$Val, Ser$Ser,	Thr$Thr)
names(AA_tot)    # check on col names, prob should rename
setnames(AA_tot, old=c("Asp$Asp", "Glu$Glu", "Gly$Gly", "Lys$Lys", "Leu$Leu",	"Phe$Phe", "Pro$Pro",	"Val$Val", "Ser$Ser",	"Thr$Thr"), 
         new=c("Asp", "Glu", "Gly", "Lys", "Leu",	"Phe", "Pro",	"Val", "Ser",	"Thr"))  
# reshape df from wide to long using gather()
AA_gather <- gather(AA_tot, key="AA", value="value", 2:12)
# bring in metadata from original data sheet
# remove the duplicate rows for SD
data2<-subset(data2, value == "ave")
names(data2)    # check on col names for orderly cleaning of final df
# remove unnecessary cols for the AAs, etc
data2_sm<-subset(data2, select = -c(photo, value, ala, asp, glu, gly, lys, leu, phe, pro, val, ser, thr))
# finally here... so good
AA_total<-merge(data2_sm,AA_gather, by="ucdavis_id")

#### make a box plot to compare results of AAs across trophic levels
p2 <- ggplot(AA_total, aes(x = TLc, y = value, fill = AA)) +
  themeKV + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),) + 
  geom_point(alpha=0.05, color="black", position="jitter", shape = 16, size = 3) +  
  geom_boxplot(alpha=0.75, colour = "black", linewidth = 0.25, outlier.color = NA) +
  scale_y_continuous(breaks= pretty_breaks()) +
  scale_fill_brewer(palette = "Spectral")+
  ylab("d15N (%)") + # this needs to read δ15N (‰) but it wont print to PDF
  facet_wrap(~AA, ncol=6, scales = "free_y")
  # facet_wrap(~AA, ncol=6)



#### Make a final analysis of the AA data to calculate TEF
#### Perform a random draw from the AA param data, calculate TEF for TL=2, Tl=3
#### Then make a raincloud style plot of the TEF distribs at each TL

# repeated from early line above just to make sure! 
forage_raw <- read.csv('data/palila_prey/palila_prey.csv')
# subset for just consumer TL=2 data to calculate TEF2
forage_TL2 <- subset(forage_raw, TL == 2)
# Generate random variates given the norm distrib params (ave, sd) given in the lab results
data <- forage_TL2
dataList = split(data, data$ucdavis_id) #split raw data up based on lab specimen ID

# build functions to generate X estimates of each AA for each sample
# use these to generate X estimates of TEF, by first estimating TEF at TL=2 and TL=3 
# TEF is defined by Nielsen et al 2015 (https://doi.org/10.1007/s00442-015-3305-7)
# set # of random variates drawn from norm distrib, outside of the function to avoid repetition
num_draws = 100
beta = 1.9063 # this is AAtrp - AAsrc at TL=1, where we pair Glu and Lys
# importantly, this is weighted 90% for mamane, and 10% for naio
# as per palila foraging data observations

# start with TEF at TL=2 
GetTEF2 <- function(anID){
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) 
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys")))
  TEF2 <- (glu.est) - (lys.est) - beta 
  return(TEF2)
}
TEF2 = lapply(dataList, GetTEF2)      # apply the function to each specimen
TEF2 = reshape2::melt(TEF2)           # melt the df
TEF2$L1 <- as.factor(TEF2$L1)                              
setnames(TEF2, old=c("L1","value"), new=c("ucdavis_id", "TEF_TL2"))    # replace TEF value column header with TEF formulation 
TEF2 <- TEF2[, c(2, 1)]    # reorder columns
# check data, should average to ~ 6.4
mean(TEF2$TEF_TL2) # 6.419575

# repeat above but for TL=3 to calculate TEF3
# subset SIL data for TL=3 
forage_TL3 <- subset(forage_raw, TL == 3)
# Generate random variates given the norm distrib params (ave, sd) given in the lab results
data <- forage_TL3
dataList = split(data, data$ucdavis_id) #split raw data up based on lab specimen ID

# then get TEF at TL=3 
GetTEF3 <- function(anID){
  glu.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) 
  lys.est<- rnorm(num_draws, mean = as.numeric(subset(anID, value == "ave", select = "lys")) , sd = as.numeric(subset(anID, value == "sd", select = "lys")))
  TEF3 <- (glu.est) - (lys.est) - beta 
  return(TEF3)
}
TEF3 = lapply(dataList, GetTEF3)      # apply the function to each specimen
TEF3 = reshape2::melt(TEF3)           # melt the df
TEF3$L1 <- as.factor(TEF3$L1)                              
setnames(TEF3, old=c("L1","value"), new=c("ucdavis_id", "TEF_TL3"))    # replace TEF value column header with TEF formulation 
TEF3 <- TEF3[, c(2, 1)]    # reorder columns
# check data, should average to ~ 7.6
mean(TEF3$TEF_TL3) # 7.586594

#### now bind all the TEF formulations into one single frame
TEF_tot <- cbind(TEF2, TEF3$TEF_TL3)
names(TEF_tot)    # check on col names, prob should rename
setnames(TEF_tot, old=c("TEF_TL2", "TEF3$TEF_TL3"), 
         new=c("TL2", "TL3"))   
# reshape df from wide to long using gather()
TEF_gather <- gather(TEF_tot, key="TEF", value="value", 2:3)
# no need for metadata additions from original data
# calculate global TEF from the palila forage data, should ~ 7.0
mean(TEF_gather$value) # 7.036752


#### make a raincloud style plot of the TEF values 
#### across both trophic levels: TL=2, TL=3

# read in the silhouettes for the plot
sil1 <- readPNG("/Users/kylevanhoutan/palila_CSIA/images/larvae.png", native = TRUE)
sil2 <- readPNG("/Users/kylevanhoutan/palila_CSIA/images/spider.png", native = TRUE)
img1 <- grid::rasterGrob(sil1, interpolate = TRUE)
img2 <- grid::rasterGrob(sil2, interpolate = TRUE)

#make the plot
p3 <- ggplot(TEF_gather, aes(x = value, y = TEF, fill = TEF)) + 
  themeKV + 
  theme(legend.position = "none", #ggdist is already grouping the TL categories
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10),) + 
  stat_dots(quantiles = 100, side = "bottom", color = NA, alpha = 0.8, height = 0.6) + # quantiles (e.g., "quantiles = 100") controls size of dots
  stat_halfeye(side = "top", alpha = 0.75, adjust = .5, height = 1) + # adjust regulates bandwidth/smoothness
  scale_fill_manual(values = c("#3288bd", "#990033")) +
  scale_color_manual(values = c("#3288bd", "#990033")) +
  scale_x_continuous(breaks = seq(2, 12, by = 1)) +
  xlab("TEF (d15N %)") +
  ylab("trophic level") +
  stat_summary(geom = "text", fontface = "bold", size = 4, vjust = -1.5,
               fun = "median", aes(label = round(after_stat(x), 2),
               color = TEF , color = after_scale(darken(color, 0.5)))) +
  inset_element(p = img1, left = 0.075, bottom = 0.3, right = 0.21, top = 0.4) + # insert the silhouettes
  inset_element(p = img2, left = 0.05, bottom = 0.68, right = 0.25, top = 0.88)


layout <- "
AAAAAAA
BBBBBCC"
p2 + p1 + p3 +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a') # add panel labels a, b, c... etc


# END
