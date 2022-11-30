library(tidyverse)
library(data.table)
library(ggbeeswarm)
library(ggrepel)

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


# Read in prey item CSIA
data <- read.csv('data/palila_prey/palila_prey.csv')

# Beta calculation
# beta = difference between trophic acid and source acids in primary producers
# our primary producers are NAIO and MAMANE
data %>% 
  filter(value == "ave") %>% 
  filter(spp %in% c("NAIO","MAMANE")) %>% 
  select(
    # trophic AA
    glu, 
    # source AA
    phe,ser,thr, spp, ucdavis_id) %>% 
  group_by(spp) %>% 
  mutate(beta = glu - (phe+thr+ser) ) %>% 
  mutate(beta_mu = mean(beta) )


# #### # #### # #### # #### # #### 
# #### prey TL Estimation ####
# #### # #### # #### # #### # #### 
dataList = split(data, data$ucdavis_id) #split raw data up based on specimen ID

# Simulating additional observations given the SD defined in the lab results
# build function to generate 1000 TP estimates based on random draws from AA gaussian distributions
GetTP <- function(anID){
  
  num_of_draws = 1000
  
  # Chikaraishi standard is b = 8.4 TEF = 7.6
  # b =   2.1
  # TEF = 5.9
  
  # b below is derived from methods in Chikaraishi et al 2009, using difference between trp and src AAs in primary producers
  # Nielsen et al 2015 described multiple ways to derive TEF based on multiple AA data
  # we followed the guidance of both studies
  b =   1.6   
  TEF = 5.3
    
  # trophic amino acids
  ala.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
  glu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamic acid
  leu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
  pro.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # pro
  val.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # v
  
  # source amino acids
  thr.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # ser
  ser.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # ser
  phe.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # ser
  
  TPforID <- (((((glu.est)/1)) - ((thr.est+ser.est+phe.est)/3) - b)/TEF ) + 1  # function to estimate TP via source and trophic AA with constants 
  
  return(TPforID)
}


TPforID_all = lapply(dataList, GetTP )          # apply that function to each specimen
TPforID_all = reshape2::melt(TPforID_all)       # melt that dataframe
TPforID_all$L1 <- as.factor(TPforID_all$L1)                              
setnames(TPforID_all, old=c("L1","value"), new=c("ucdavis_id", "tp"))
data<-subset(data, value == "ave")
total<-merge(data,TPforID_all, by="ucdavis_id")

# jitterplot of TL estimates
ggplot(total, aes(x = forcats::fct_reorder(spp, tp), y = tp))+
  geom_point(size = .2, position = "jitter")+
  themeo 

#####################################
# Plotting amino acid density plots #
#####################################

density<-ggplot(data = total, aes(group = spp, x = tp, fill = spp, color = spp))+
  #geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = 'dark green', alpha = .5)+
  
  geom_density(show.legend = F, alpha = .7, color = "black")+
  #geom_rug(alpha = .05,show.legend = F)+
  geom_hline(yintercept = 0, color = "black", size = .5)+
  geom_label_repel(data = total %>% dplyr::group_by(spp) %>% dplyr::summarise(median = median(tp)) , aes( x = median,y = 3, label = spp), color = "black",show.legend = F)+
  
  scale_fill_brewer(palette = "Spectral")+
  #scale_fill_manual(values = c("#908837","#8c5b2c","#cdb768","#8aa13e","#ab6129","#8f7752"))+
  #scale_color_manual(values = c("#908837","#8c5b2c","#cdb768","#8aa13e","#ab6129","#8f7752"))+
  themeo+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  xlab("trophic position")

density


# density<-ggplot(data = total, aes(group = spp, x = tp, fill = spp), color = "black")+
#   geom_density(aes(x = tp + .01),show.legend = F, kernel = "gaussian", alpha = .5, fill = "black", color = NA)+
#   geom_density(show.legend = F, kernel = "gaussian", alpha = .7)+
#   #geom_rug(aes(color = spp),alpha = .05,show.legend = F)+
#   geom_hline(yintercept = 0, color = "black", size = .5)+
#   geom_label_repel(data = total %>% dplyr::group_by(spp) %>% dplyr::summarise(median = median(tp)) , aes( x = median,y = 3, label = spp), color = "black",show.legend = F)+
#   #scale_fill_brewer(palette = "Dark2")+
#   scale_fill_manual(values = c("#908837","#8c5b2c","#cdb768","#8aa13e","#ab6129","#8f7752"))+
#   scale_color_manual(values = c("#908837","#8c5b2c","#cdb768","#8aa13e","#ab6129","#8f7752"))+
#   themeo
# density


source_csv <- total %>% dplyr::group_by(spp) %>% dplyr::summarise(MeanTL = mean(tp), SDTL = sd(tp), n = 1000)
#write.csv(source_csv, "./data/mixing model/source_data.csv", row.names = F)

# discr data
disc_data <- data.frame( Sources = unique(source_csv$spp), MeanTL = 1, SDTL = 0)
# write.csv(disc_data, "./data/mixing model/discrimination_data.csv", row.names = F)


# # merge NAIO/MAMANE, CYDIA/FOLCAT, ARTHO/SPIDER
CYDIA_ARTHO_FOLCAT <- source_csv %>% filter(spp %in% c('CYDIA', 'ARTHO', "FOLCAT")) %>% summarise(spp = 'CYDIA_ARTHO_FOLCAT', MeanTL = mean(MeanTL), SDTL = mean(SDTL), n = 1000)
MAMANE_NAIO <- source_csv %>% filter(spp %in% c('MAMANE', 'NAIO')) %>% summarise(spp = 'MAMANE_NAIO', MeanTL = mean(MeanTL), SDTL = mean(SDTL), n = 1000)
SPIDER <- source_csv %>% filter(spp %in% c('SPIDER')) %>% summarise(spp = 'SPIDER', MeanTL = mean(MeanTL), SDTL = mean(SDTL), n = 1000)

# group_source <- rbind(CYDIA_ARTHO_FOLCAT,MAMANE_NAIO,SPIDER)
# #write.csv(group_source, "./data/mixing model/source_data.csv", row.names = F)
#
# # discr data
# disc_data <- data.frame( Sources = unique(group_source$spp), MeanTL = 1, SDTL = 0)
# # write.csv(disc_data, "./data/mixing model/discrimination_data.csv", row.names = F)
#


########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  
# Script utilized to develop Figure 2, fully run Fig2_part1.R prior to running this script. Do not clear the environment.
########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  

library(gridExtra)

data <- read.csv('data/palila_prey/palila_prey.csv')
# Simulating additional observations given the SD defined in the lab results
dataList = split(data, data$ucdavis_id) #split raw data up based on specimen ID
  
  #build function to generate 1000 TP estimates based on random draws from AA gaussian distributions
  GetTP <- function(anID){
    
    num_of_draws = 1000

    # trophic amino acids
    ala.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ala")) , sd = as.numeric(subset(anID, value == "sd", select = "ala"))) # alanine
    glu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "glu")) , sd = as.numeric(subset(anID, value == "sd", select = "glu"))) # glutamic acid
    leu.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "leu")) , sd = as.numeric(subset(anID, value == "sd", select = "leu"))) # leucine
    pro.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "pro")) , sd = as.numeric(subset(anID, value == "sd", select = "pro"))) # pro
    val.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "val")) , sd = as.numeric(subset(anID, value == "sd", select = "val"))) # v
    
    # source amino acids
    thr.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "thr")) , sd = as.numeric(subset(anID, value == "sd", select = "thr"))) # ser
    ser.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "ser")) , sd = as.numeric(subset(anID, value == "sd", select = "ser"))) # ser
    phe.est<- rnorm(num_of_draws, mean = as.numeric(subset(anID, value == "ave", select = "phe")) , sd = as.numeric(subset(anID, value == "sd", select = "phe"))) # ser
    
    TPforID <- cbind(ala.est,glu.est,leu.est,pro.est,val.est,thr.est,ser.est,phe.est)

    return(TPforID)
  }


TPforID_all = lapply(dataList, GetTP )          # apply that function to each specimen
TPforID_all = reshape2::melt(TPforID_all)       # melt that dataframe
TPforID_all$L1 <- as.factor(TPforID_all$L1)                              
setnames(TPforID_all, old=c("L1","value","Var2","Var1"), new=c("ucdavis_id", "aa_val","aa_name","draw_num"))
data <- data %>% filter(value == "ave") %>% select(ucdavis_id,spp)
total<-merge(data,TPforID_all, by="ucdavis_id")


# source acids
source_AAs <- total %>% 
  filter(aa_name %in% c("phe.est","thr.est","ser.est")) %>% 
  spread(aa_name,aa_val) %>% 
  mutate(source_mean = (thr.est + ser.est + phe.est) /3 ) %>% 
  ggplot() +
  geom_density(aes(x = source_mean, fill = spp), alpha = .7, show.legend = F)+
  geom_hline(yintercept = 0, color = "black", size = .5)+
  scale_fill_brewer(palette = "Spectral")+
  themeo+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  xlab(expression(paste(delta^15,"N"["Source"]) ) )


# trophic acids
trophic_AAs <- total %>% 
  filter(aa_name %in% c( "glu.est")) %>% 
  spread(aa_name,aa_val) %>% 
  mutate(trophicAA_mean = (glu.est)) %>% 
  ggplot() +
  geom_density(aes(x = trophicAA_mean, fill = spp), alpha = .7, show.legend = F)+
  geom_hline(yintercept = 0, color = "black", size = .5)+
  scale_fill_brewer(palette = "Spectral")+
  scale_color_brewer(palette = "Spectral")+
  themeo+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  xlab(expression(paste(delta^15,"N"["Trophic"]) ) )

grid.arrange(source_AAs, trophic_AAs)

# formula build
formula_t<-ggplot()+
  geom_blank()+
  annotate('text', x = 0, y = 0, 
           label =  expression(italic(paste("trophic position = ", frac(paste( delta^15,"N"["Trophic"], " - " ,delta^15,"N"["Source"], " - ", beta), "TEF"  ), " + 1"  ) )  ),
           parse = TRUE,size=7)+
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5))+
  theme_void()
formula_t

## OSM plot for individual AA SIA values for each prey item
# density is pulled from the prey_TL.R script
grid.arrange(source_AAs, trophic_AAs, density)
grid.arrange(source_AAs, trophic_AAs, formula_t, density, ncol = 1)

total %>% 
  ggplot(aes(x = aa_val, fill = aa_name))+
    geom_histogram(color = "black", size = .25) +
    scale_fill_brewer(palette = "Spectral")+
    facet_grid(aa_name ~ spp)



          # "ARTHO"  "CYDIA"  "FOLCAT" "MAMANE" "NAIO"   "SPIDER"
TP_prior <- c(1,      10,        8,         74,      5,      2)

## plot of dietary proportion prior bar stack plot, Figure 2j
prop_plot <- data.frame(prop = c(1,10,8,74,5,2), spp = c("ARTHO","CYDIA", "FOLCAT", "MAMANE", "NAIO" , "SPIDER")) %>% 
  ggplot()+
  geom_bar(aes(x = 1, fill = fct_reorder(spp, prop, .desc = T), y = prop), stat = "identity",show.legend = F)+
  #geom_bar(aes(x = 1, fill = spp, y = prop), stat = "identity")+
  
  scale_fill_manual(values = c("#E6F598", "#FC8D59", "#FEE08B", "#99D594","#3288BD", "#D53E4F"))+
  xlab("prey item")+
  themeo+
  scale_y_continuous(expand = c(0,0), breaks = c(0,100))+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        plot.margin = unit(c(.75,.75,.75,.75), "cm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("proportional prior")
prop_plot 

grid.arrange(source_AAs, trophic_AAs, formula_t, density, prop_plot,
             
             
             layout_matrix = as.matrix( rbind(
               c(1,1,1,1,1,5),
               c(2,2,2,2,2,5),
               c(3,3,3,3,3,5),
               c(4,4,4,4,4,5)
                  )
             )
)





