# exploratory script where we considered using individual AA nitrogen fractionation as individual tracers in a potential mixing model
# becuase we were able to confidently differentiate between grouped diet items, we ultimately decided not to pursue this further. 

data <- read.csv('/Users/tgagne/Palila_2018/data/palila_prey/palila_prey.csv')


data %>% filter(value == "ave") %>% 
  ggplot(aes(x = spp, y = phe - glu))+
  geom_point()

id = 1
ass = 1

plot(data$asp, data$gly, col = data$spp)
points(CSSIAPalila$asp, CSSIAPalila$gly, pch = 2)

palila_data = CSSIAPalila #%>% filter(value == "ave")

palila_data <- palila_data[,c(1,2,4,5:15)]

colnames(data)[8] <- "ile"

data <- rbind(palila_data,data)

ggplot(data %>% filter(value == "ave"),aes(x = glu - phe, y = thr, color = spp))+
  geom_point()



# Loop and resample

id_box <- NULL
davis_IDs <- levels(as.factor(data$ucdavis_id))
acids <- colnames(data)[4:14]
# for each ID
for( id in 1:length(davis_IDs)) {
  
  d_id <- filter(data, ucdavis_id == davis_IDs[id])
  
  acid_box <- NULL
  # for each amino acid
  for(AA in 1:length(acids)){
    
    ass_one <- select(d_id, acids[AA])
    ass_samp <- rnorm(n = 100, mean = ass_one[1,], sd = ass_one[2,]) %>% as.data.frame()
    ass_samp$run <- seq(1,length(ass_samp[,1]))
    ass_samp$acid <- acids[AA]
    ass_samp$ucdavis_id <- davis_IDs[id]
    acid_box <- rbind(ass_samp,acid_box)
    
  }
 
  id_box <- rbind(acid_box, id_box)
   
}


colnames(id_box)[1] <- "value"
#id_box$id <- seq(1,length(id_box$acid))

id_box <- left_join(data %>% filter(value == "ave") %>% select(ucdavis_id,spp), id_box, by = "ucdavis_id") 
acid_cors <- id_box %>% spread(acid, value) 

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
my_cols <- RColorBrewer::brewer.pal(7,"Dark2") 
upper.panel<-function(x, y){
  points(x,y, pch = 1, cex = .5,col = my_cols[acid_cors$spp])
}
# Create the plots
pairs(acid_cors[,4:14], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)







#####

#read in an merge palila data
# palila UC_davis and kyle ID cross reference table
Palila_meta <- read.delim('/Users/tgagne/Palila_2018/data/palila_samples/palila_sample_metadata_MBA_UCDavis_IDs.txt')
# UC davis stable isotope data
CSSIAPalila <- read.csv('/Users/tgagne/Palila_2018/data/palila_samples/Palila_CSSIA_Aug20.csv')
# location, individual, date, metadata
addt_meta <- read.csv('/Users/tgagne/Palila_2018/data/palila_samples/Palila_specimens_feather_database.csv')

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

# bring additional metadate from other sheet... and join again below
addt_meta <- addt_meta %>% 
  dplyr::select(ID_NO, SEX, AGE, SPEC_NOTES,SPEC_PREP, LOCATION,REGION) %>% 
  mutate(accession_id = ID_NO) %>% dplyr::select(-ID_NO)

# join
CSSIAPalila <- CSSIAPalila %>% left_join(addt_meta, by = "accession_id")

# clean up post joining
rm(addt_meta,Palila_meta)

# Simulating additional observations given the SD defined in the lab results
data <- CSSIAPalila

str(data)

id_box <- NULL

#id = 1;AA = 1

davis_IDs <- levels(as.factor(data$ucdavis_id))
acids <- c("ala","glu","leu","pro","val","thr","ser","phe")#colnames(data)[5:15]
# for each ID
for( id in 1:length(davis_IDs)) {
  
  d_id <- filter(data, ucdavis_id == davis_IDs[id])
  
  acid_box <- NULL
  # for each amino acid
  for(AA in 1:length(acids)){
    
    ass_one <- select(d_id, acids[AA])
    
    ass_samp <- rnorm(n = 100, mean = ass_one[1,], sd = ass_one[2,]) %>% as.data.frame()
    
    ass_samp$run <- seq(1,length(ass_samp[,1]))
    ass_samp$acid <- acids[AA]
    ass_samp$ucdavis_id <- davis_IDs[id]
    
    acid_box <- rbind(ass_samp,acid_box)
    
  }
  
  id_box <- rbind(acid_box, id_box)
  
}


colnames(id_box)[1] <- "value"
str(id_box)
#id_box$id <- seq(1,length(id_box$acid))

id_box <- left_join(data %>% filter(value == "ave") %>% select(ucdavis_id,spp,year), id_box, by = "ucdavis_id") 

acid_cors <- id_box %>% spread(acid, value) 

str(acid_cors)


panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
my_cols <- colorRampPalette(brewer.pal(3,"Spectral"))(length((acid_cors$year)))

upper.panel<-function(x, y){
  points(x,y, pch = 1, cex = .5 , col = my_cols[acid_cors$year])
}
# Create the plots
pairs(acid_cors[,5:12], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

#plot prey item AAs over Palila AAs





