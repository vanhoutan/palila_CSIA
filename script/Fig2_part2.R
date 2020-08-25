
########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  
# Script utilized to develop Figure 2, fully run Fig2_part1.R prior to running this script. Do not clear the environment.
########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  

library(gridExtra)

data <- read.csv('/Users/tgagne/Palila_2018/data/palila_prey/palila_prey.csv')
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

# prop prior plot
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





