library(tidyverse)


# read in CSIA data from both prey and palila
prey_data <- read.csv('/Users/tgagne/Palila_2018/data/palila_prey/palila_prey.csv')
CSSIAPalila <- read.csv('/Users/tgagne/Palila_2018/data/palila_samples/Palila_CSSIA_Aug20.csv')

all_aa <- merge(prey_data,CSSIAPalila, all = T)                                     # merge them
all_aa <- all_aa %>% filter(value == "ave")                                         # only pull out the ave
all_aa[,c("ucdavis_id", "acession_id", "year", "location", "ile", "lys")] <- NULL   # get rid of vars we wont need

# bulid a list of all possible combinations of reasonable TEFs and betas
combos <- expand.grid(beta = seq(from = -.5, to = 10, length.out = 10),TEF = seq(from = 4, to = 12, length.out = 10))

# repeat this list of combos as many times as there are samples in the CSIA
combos <- combos[rep(seq_len(nrow(combos)), each=nrow(all_aa)),]

# repeat the CSIA sample df enough so that we have enough samples as we have combinations
all_aa_long <- NULL
for(i in 1:length(unique(combos$beta)) - 1){
  all_aa_long <- rbind(all_aa_long,all_aa)
}

# join together, we will have a row for every sample and every combination of TEF/beta that we're testing
# bind these TEF/beta combos to the samples
all_aa_long <- cbind(all_aa_long,combos)

# calculate TP from all these combinations
all_aa_long <- all_aa_long %>% mutate(tp = ( ( (ala + glu + leu + pro + val)/5 ) - ((thr+ser+phe)/3 ) - beta) / TEF + 1 )


all_aa_long %>% 
  filter(tp < 4 & tp > 0) %>% 
  ggplot() +
  geom_tile(aes(x = TEF, y = beta, fill = tp))+
  geom_text(aes(x = TEF, y = beta, label = round(tp,2)))+
  scale_fill_distiller(palette = "Spectral")+
  facet_wrap(~spp)

# find out beta and TEF where Manmane and NAIO have a mean TP of 1
TEF_BETA <- all_aa_long %>% 
  # want to know when Manmane and Naio
  filter(spp == "MAMANE" | spp ==  "NAIO") %>% 
  
  group_by(beta, TEF) %>% 
  mutate(mean_tp = mean(tp)) %>% 
  ungroup() %>% 
  
  
  # have mean TP of around 1
  filter( mean_tp > .5 & mean_tp < 1.3) 


betas <- unique(TEF_BETA$beta) # list of betas meeting conditions above
TEFS <- unique(TEF_BETA$TEF) # same for TEFs


all_aa_long %>% 
  filter(spp == "palila") %>% 
  filter(beta == betas[1],
         TEFS == TEFS[1])







all_aa_long %>% 
  filter(tp < 2.5 & tp > 0) %>% 
  ggplot() +
  geom_point(aes(x = TEF, y = beta, fill = tp),  position = "jitter", size = 5, pch = 22, color = "black")+
  #geom_text(aes(x = TEF, y = beta, label = round(tp,0)), position = "jitter")+
  scale_fill_distiller(palette = "Spectral")+
  facet_wrap(~spp)+
  theme_bw()
