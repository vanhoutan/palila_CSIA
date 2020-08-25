
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
levels = 50
ppoints(levels)

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

# Alternatives for other formulas
# mod_df %>%
#   data_grid(SPEI36 = seq_range(SPEI36, n = 101), rollmean = seq_range(rollmean, n = 101)) %>%
#   tidybayes::add_predicted_draws(data.brms) %>%
#   ggplot(aes(x = SPEI36, y = median)) +
#   tidybayes::stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "#08519C") +
#   geom_point(data = mod_df, size = .1, position = "jitter") +
#   scale_x_continuous(expand = c(0,0))+
#   themeo+
#   scale_fill_brewer()
# 
# mod_df %>%
#   data_grid(SPEI36 = seq_range(SPEI36, n = 101), rollmean = seq_range(rollmean, n = 101)) %>%
#   tidybayes::add_predicted_draws(data.brms) %>%
#   ggplot(aes(x = rollmean, y = median)) +
#   tidybayes::stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "#08519C") +
#   geom_point(data = mod_df, size = .1, position = "jitter") +
#   scale_x_continuous(expand = c(0,0))+
#   themeo+
#   scale_fill_brewer()


# checking residuals
qqnorm(mod_df$raw_tp, pch = 1, frame = FALSE)
qqline(mod_df$raw_tp, col = "steelblue", lwd = 2)

preds_brms <- predict(data.brms, newdata = mod_df) %>% data.frame()
mod_df$modeled <- preds_brms$Estimate
str(preds_brms)
preds_brms$raw_tp <- mod_df$raw_tp
preds_brms$year <- year(mod_df$year)

ggplot(mod_df) +
  geom_point(aes(x = rollmean, y = raw_tp, color = as.factor(year)))

ggplot()+
  geom_point(aes(x=preds_brms[,1],y=mod_df$raw_tp), size = .2)+
  geom_abline(intercept = 0, slope = 1)+
  scale_x_continuous(limits = c(min(mod_df$raw_tp),max(mod_df$raw_tp)))
  
ggplot(preds_brms)+
  geom_errorbar(aes(x = seq_len(nrow(mod_df)), ymin = Q2.5, ymax = Q97.5), color = "dark gray")+
  geom_point(aes(x = seq_len(nrow(mod_df)), y = raw_tp, fill = as.factor(year)), shape = 21)+
  geom_point(aes(x = seq_len(nrow(mod_df)), y = Estimate), size = .2, color = "red")+
  scale_fill_brewer(palette = 'Dark2')
  
ggplot(preds_brms)+
  geom_errorbarh(aes(xmax = Q97.5, xmin = Q2.5, y = mod_df$raw_tp), alpha = .5, size = .1)+
  geom_point(aes(x=Estimate,y=mod_df$raw_tp, color = as.factor(mod_df$date)),size = .02)+
  geom_abline(intercept = 0, slope = 1, color = "red")+
  scale_x_continuous(limits = c(min(mod_df$raw_tp),max(mod_df$raw_tp)))+
  themeo

ggplot()+
  geom_histogram(aes(x = mod_df$raw_tp))+
  geom_histogram(aes(x = preds_brms$Estimate), fill = "blue4", alpha = .7)

plot(preds_brms$Estimate, mod_df$raw_tp - preds_brms$Estimate)
abline(h = 0, col = "red")

ggplot(preds_brms,aes(Estimate, mod_df$raw_tp - Estimate)) +
  geom_errorbarh(aes(xmax = Q97.5, xmin = Q2.5), alpha = .5, size = .1, color = "dark gray")+
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 1.5, shape = 22, fill = "white") +
  themeo


###############################################
# comparing and visualizing marginal effects ##
###############################################

# build dataframe that covers all possible estimates
surf_data <- mod_df %>%  data_grid(SPEI36 = seq_range(SPEI36, n = 101), rollmean = seq_range(rollmean, n = 101))
surf_data <- mod_df %>%  data_grid(SPEI36 = seq_range(SPEI36, n = 101), rollmean = seq_range(rollmean, n = 101), parasitism = seq_range(parasitism, n = 101))
surf_data <- data.frame(SPEI36 = mod_df$SPEI36, rollmean = mod_df$rollmean, parasitism = mod_df$parasitism)

## this takes a while, several hours potentially depending on grid.resolution
surface <- pdp::partial(data.brms, pred.var =  c("SPEI36", "rollmean"),train = surf_data, progress = "text", chull = T, grid.resolution = 5)
surface <- pdp::partial(data.brms, pred.var =  c("SPEI36", "parasitism"), train = surf_data, progress = "text", chull = F, grid.resolution = 5)

#surface <- pdp::partial(data.brms, pred.var =  c("SPEI36", "rollmean"),train = mod_df, progress = "text", chull = T, grid.resolution = 10)
ggplot()+
  geom_tile(data = surface, aes(SPEI36,parasitism,fill = yhat)) +
  geom_point(data = mod_df, aes(SPEI36, parasitism, fill = raw_tp), shape = 21, size = 3)+
  scale_fill_distiller(palette = 'Spectral')

# estimate convex hull of the observed data
x_chull <- mod_df$SPEI36[chull(x = mod_df$SPEI36, y = mod_df$rollmean)]
y_chull <- mod_df$rollmean[chull(x = mod_df$SPEI36, y = mod_df$rollmean)]
z_chull <- mod_df$modeled[chull(x = mod_df$SPEI36, y = mod_df$rollmean)]

# plot surface with convex hull
ggplot()+
  geom_tile(data = surface, aes(SPEI36,parasitism,fill = yhat)) +
  geom_point(data = mod_df, aes(SPEI36, parasitism, fill = raw_tp), shape = 21, size = 3) +
  geom_polygon(aes(x = x_chull, y = y_chull), color = "black", fill = NA, size = 1)+
  scale_fill_distiller(palette = 'Spectral')+
  theme_bw()

ggplot(surface)+
  geom_point(aes(x = parasitism, y = yhat, fill = SPEI36), shape = 21, size = 20)

ggplot(surface)+
  geom_point(aes(x = SPEI36, y = yhat, fill = parasitism), shape = 21, size = 20)+
  scale_y_continuous(limits = c(1,3))

ggplot(surface)+
  geom_hex(aes(x = SPEI36, y = parasitism, fill = yhat), stat = "identity", size = 10)+
  coord_fixed()


# get data in to a list for persp plotting
data_z <- reshape2::acast(surface, SPEI36~rollmean, value.var = "yhat")
data_z <- reshape2::acast(surface, SPEI36~parasitism, value.var = "yhat")

#plotly::plot_ly(z = data_z,  type = "surface", colors = "Spectral")
str(data_z)

# pull numeric vectors for persp
x <- attr(data_z,"dimnames")[1] %>% as.data.frame(); x <- as.vector(x[,1]); x <- as.numeric(x)
y <- attr(data_z,"dimnames")[2] %>% as.data.frame(); y <- as.vector(y[,1]); y <- as.numeric(y)

# Color palette (100 colors)
col.pal<-colorRampPalette(c("yellow","yellow", "dark green"))
colors<-col.pal(100)
# height of facets
z = data_z
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

colors = colors[z.facet.range]

dev.new()

seq360 <- rbind(seq(-180,180,by = .25),seq(-180,180,by = .25)) 

for(i in 1:10){
  
for(round in 1:720){
  
  pmatz<-persp( x = x, y, z = data_z, 
         col = colors,
         #col = "green",
         shade = 0.4,
         expand = 0.4,
         #phi = 30 , theta = -60,
         phi = 20 , #theta = -60,
         
         lwd = .1,
         zlab = "trophic level",
         xlab = "drought index",
         ylab = "mean air temperature",
         box = T,
         zlim = c(min(data_z, na.rm = TRUE) - sd(data_z, na.rm = T), max(data_z, na.rm = TRUE) + sd(data_z, na.rm = T)),
         #axes = T,
         #nticks = 10
         theta = seq360[round]
         
  )
  #polygon(trans3d(x = x_chull, y = y_chull, z = z_chull, pmat = pmatz), 
          #density = 5, 
   #       lwd = 3, border = "dark gray")
  
  
}

}


# playing with overplotting training data
points(
  trans3d(
    x = as.vector(mod_df$SPEI36), 
    y = as.vector(mod_df$rollmean), 
    z = as.vector(mod_df$raw_tp),pmat = ptz), 
  col = RColorBrewer::brewer.pal(11,name="Spectral")[as.numeric(cut(mod_df$raw_tp,breaks = 10))], 
  cex = .25
  )









