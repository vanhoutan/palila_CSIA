library(data.table)
library(MixSIAR)
library(dplyr)
library(plyr)
library(ggplot2)
library(R2jags)

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

sppx<- "PALI"
x = 1

#for(x in 1:length(sppx)){
  mix<-read.csv('./data/mixing model/mixture_data.csv')
  str(mix)
  mix<-subset(mix, spp == sppx[x])
  write.csv(mix,file = "./data/mixing model/mixture_data_spp.csv")
  
  #mix data = i.e. consumer
  mix<- load_mix_data(filename = "./data/mixing model/mixture_data_spp.csv",iso_names = "TL",factors = "spp",
                      fac_random = TRUE,fac_nested = NULL,cont_effects = "year")
  #source data #new omma and caran TL
  source <- load_source_data( filename = "./data/mixing model/source_data.csv",
                              source_factors = NULL,conc_dep = FALSE,data_type = "means", mix)
  #widen priors of source distributions
  source$S_SIG[,1] <- source$S_SIG[,1]
  

  #discrimination data
  discr <- load_discr_data(filename = "./data/mixing model/discrimination_data.csv",mix)
  #discr$mu[,1] <- rep( .9,  6)
  #discr$sig2[,1] <- rep(0.5,6)
 
  #one dimensional isospace plot
  plot_data(filename="isospace_plot",
      plot_save_pdf=FALSE,
      plot_save_png=FALSE,
      mix,source,discr)
  
  #uninformative prior = alpha.prior = 1
  #construct informative prior from stomach_proportional_plot_apr13.xlsx

             # "ARTHO"  "CYDIA"  "FOLCAT" "MAMANE" "NAIO"   "SPIDER"
  TP_prior <- c(1,      10,        8,         74,      5,      2)
  
               # "CYDIA_ARTHO_FOLCAT" "MAMANE_NAIO" "SPIDER" 
 TP_prior <- c(19,      82,        2     )
  
  
  # Generate alpha hyperparameters scaling sum(alpha)=n.sources
  TP_prior <- TP_prior*length(TP_prior)/sum(TP_prior)
  #TP_prior <- rep(1,6)
  #TP_prior <- TP_prior*10
  
  # the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
  TP_prior[which(TP_prior==0)] <- 0.01
  
  TP_prior
  
  plot_prior(alpha.prior = TP_prior, source,
           plot_save_pdf=FALSE,
         plot_save_png=FALSE)
  #write jags model
  model_filename <- "MixSIAR_model.txt"

  
  # working error structure
  resid_err <- FALSE
  process_err <- TRUE
  
  write_JAGS_model(model_filename, resid_err,process_err, mix, source)
  #run model
  jags.1 <- run_model(run="short",mix,source,discr,model_filename,
                      alpha.prior = TP_prior,resid_err,process_err)
  
  plot_continuous_var(jags.1, mix, source, output_options)
  
  rm(jags.1,model_filename)

  ##########################
  ##########################
  #R2jags::attach.jags(jags.1)
  
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  f1 = 1
  ce = 1
  
  fac.lab <- mix$FAC[[1]]$labels[f1]
  label <- mix$cont_effects[ce]
  cont <- mix$CE[[ce]]
  ilr.cont <- get(paste("ilr.cont",ce,sep=""))
  
  get_high <- function(x){return(quantile(x,.975))}
  get_low <- function(x){return(quantile(x,.025))}
  
  n.plot = 200
  chain.len = dim(p.global)[1]
  Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
  ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
  ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
  ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
  ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
  for(src in 1:n.sources-1){
    for(i in 1:n.plot){
      ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i] + ilr.fac1[,f1,src]
      ilr.low[i,src] <- get_low(ilr.plot[i,src,])
      ilr.median[i,src] <- median(ilr.plot[i,src,])
      ilr.high[i,src] <- get_high(ilr.plot[i,src,])
    }
  }
  
  # Transform regression lines from ILR-space to p-space
  e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
  for(i in 1:(n.sources-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
  }
  cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
  p.median <- array(data=NA,dim=c(n.plot, n.sources))
  cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
  p.low <- array(data=NA,dim=c(n.plot, n.sources))
  cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
  tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
  p.high <- array(data=NA,dim=c(n.plot, n.sources))
  for(i in 1:n.plot){
    for(j in 1:(n.sources-1)){
      cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
      cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
      cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
    }
    for(src in 1:n.sources){
      tmp.p.med[i,src] <- prod(cross.med[i,src,]);
      tmp.p.low[i,src] <- prod(cross.low[i,src,]);
      tmp.p.high[i,src] <- prod(cross.high[i,src,]);
    }
    for(src in 1:n.sources){
      p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
      p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
      p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
    }
  }
  colnames(p.median) <- source_names
  
  Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
  df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
  colnames(df) <- c("source","median","x","low","high")
  
  str(df)
  
  # medians <- data.frame(cont,apply(p.ind,c(2,3),median))
  # colnames(medians) <- c("cont",source_names)
  # medians <- melt(medians,id="cont")
  
  # Plot of Diet vs. Cont effect
  # Page 370 in Francis et al (2011)
  ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
    ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
    ggplot2::labs(title = fac.lab) +
    ggplot2::ylab("Diet Proportion") +
    ggplot2::xlab(label) +

    scale_color_manual(values = c("#faa61a","#008745","#3a60ac"))+
    scale_fill_manual(values = c("#faa61a","#008745","#3a60ac"))+
    scale_x_continuous(expand = c(0,0))+
    themeo+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=.5),
          legend.title=element_blank(),
          strip.text=element_text(hjust=0),
          legend.position = c(0.8,0.5))
  
  # proportional plot through time with 95% prediction intervals shown
  ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
    geom_area(aes(x = x, y = median, fill = source, group = source), position = 'stack', show.legend = F, alpha = .85)+
    geom_area(aes(x = x, y = low, fill = source, group = source), fill = NA, color = "black", position = 'stack', lty = "dashed", size = .25)+
    geom_area(aes(x = x, y = high, fill = source, group = source), fill = NA, color = "black", position = 'stack', lty = "dashed", size = .25)+
    
    annotate("text", x = 1920, y = .9, label = "CATERPILLAR/ARTHROPOD")+
    annotate("text", x = 1980, y = .7, label = "MAMANE/NAIO")+
    annotate("text", x = 1900, y = .1, label = "SPIDER")+
   
    ggplot2::labs(title = fac.lab) +
    ggplot2::ylab("proportion of diet") +
    ggplot2::xlab(label) +
    scale_color_manual(values = c("#faa61a","#008745","#3a60ac"))+
    scale_fill_manual(values = c("#faa61a","#008745","#3a60ac"))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    
    themeo+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=.5),
          legend.title=element_blank(),
          strip.text=element_text(hjust=0),
          legend.position = c(0.7,0.5),
          legend.background =  element_rect(colour = "black", fill="light gray", size=.5),
          axis.text = element_text(color = "black"))
  
  
# plot of the ratio of plant based items to non-plant based items through time  
df %>% 
  select(x, source, median) %>%  
  spread(key = source, value = median) %>% 
  mutate(feed_ration = MAMANE_NAIO / (CYDIA_ARTHO_FOLCAT + SPIDER) ) %>% 
  ggplot() +
  geom_line(aes(x = x, y = feed_ration))+
  scale_x_continuous(expand = c(0,0))+
  themeo+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.title=element_blank(),
        strip.text=element_text(hjust=0),
        legend.position = c(0.8,0.5),
        axis.text = element_text(color = "black"))

  

  