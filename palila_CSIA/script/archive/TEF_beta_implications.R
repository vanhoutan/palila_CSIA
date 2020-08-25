library(tidyverse)

data <- read.csv('/Users/tgagne/Palila_2018/data/palila_prey/palila_prey.csv') %>% filter(value == "ave") %>% select(ucdavis_id,spp,glu,phe)

str(data)

data %>% 
  mutate(beta = glu - phe)


current_beta = .4 #from chikaraishi 2010 C3 plants

data %>% filter(spp %in% c("NAIO")) %>% mutate(beta = glu - phe) %>%  summarise(mean(beta))


# Understanding TEF and beta impacts on TL
TL <- function(N15Glu,N15Phe,beta,TEF) {
  
  ( (N15Glu - N15Phe - beta)  / TEF ) + 1 
  
  }


Beta_seq <- seq(-5,10, by = 1)
TEF_seq <- seq(1,10, by = 1)

all_combos <- crossing(Beta_seq, TEF_seq)

all_combos$TL <- TL(N15Glu = 2.19, N15Phe = 14.88, beta = all_combos$Beta_seq, TEF =  all_combos$TEF_seq)


ggplot(all_combos,aes(Beta_seq,TEF_seq, fill = TL))+
  geom_tile()+
  coord_fixed()

library(reshape2)
library(plotly)
data_z <- acast(all_combos, Beta_seq~TEF_seq, value.var = "TL")
data_z <- ifelse( data_z == Inf, 0, data_z)
data_z <- ifelse( data_z == -Inf, 0, data_z)

plot_ly(z = data_z,  type = "surface", colors = "Spectral")
str(data_z)

ggplot(all_combos,aes(x = Beta_seq, y = TL, color = TEF_seq))+
  geom_point(size = 5)+
  themeo















