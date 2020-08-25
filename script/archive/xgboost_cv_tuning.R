library(xgboost)
# Matrix for xgb: dtrain and dtest, "label" is the dependent variable
dtrain <- xgb.DMatrix(data.matrix(X), label = y)


#dtest <- xgb.DMatrix(X_test, label = Y_test)

best_param <- list()
best_seednumber <- 1234
best_rmse <- Inf
best_rmse_index <- 0

set.seed(123)
for (iter in 1:20) {
  param <- list(objective = "reg:linear",
                eval_metric = "rmse",
                #max_depth = sample(3:10, 1),
                #eta = runif(1, .01, .1), # Learning rate, default: 0.3
                subsample = .6, #runif(1, .6, .9),
                #colsample_bytree = runif(1, .5, .8), 
                #min_child_weight = sample(1:40, 1),
                #max_delta_step = sample(1:10, 1)
                max_depth = sample(3:10, 1), 
                eta = runif(1, .01, .1)
  )
  cv.nround <-  500
  cv.nfold <-  5 # 5-fold cross-validation
  seed.number  <-  sample.int(10000, 1) # set seed for the cv
  set.seed(seed.number)
  mdcv <- xgb.cv(data = dtrain, params = param,  
                 nfold = cv.nfold, nrounds = cv.nround,
                 verbose = T, early_stopping_rounds = 8, maximize = FALSE)
  
  min_rmse_index  <-  mdcv$best_iteration
  min_rmse <-  mdcv$evaluation_log[min_rmse_index]$test_rmse_mean
  
  if (min_rmse < best_rmse) {
    best_rmse <- min_rmse
    best_rmse_index <- min_rmse_index
    best_seednumber <- seed.number
    best_param <- param
  }
}

# The best index (min_rmse_index) is the best "nround" in the model
nround = best_rmse_index
set.seed(best_seednumber)
xg_mod <- xgboost(data = dtrain, params = best_param, nround = nround, verbose = T)

#plotly 3D surface
library(plotly)

surface <- partial(xg_mod, pred.var =  c("SPEI36", "rollmean"), train = X, progress = "text", chull = T, grid.resolution = 40)
str(surface)

library(reshape2)
data_z <- acast(surface, SPEI36~rollmean, value.var = "yhat")
plot_ly(z = data_z,  type = "surface", colors = "Spectral")
str(data_z)

x <- attr(data_z,"dimnames")[1] %>% as.data.frame(); x <- as.vector(x[,1]); x <- as.numeric(x)
y <- attr(data_z,"dimnames")[2] %>% as.data.frame(); y <- as.vector(y[,1]); y <- as.numeric(y)

# Color palette (100 colors)
col.pal<-colorRampPalette(c("yellow", "dark green"))
colors<-col.pal(100)
# height of facets
z = data_z
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

colors = colors[z.facet.range]

dev.new()
seq360 <- rbind(seq(-180,180,by = 1),seq(-180,180,by = 1))
#for(round in 1:720){

persp( x = x, y, z = data_z, 
       col = colors,
       #col = "green",
       shade = 0.4,
       expand = 0.4,
       phi = 30 , theta = -60
       ,lwd = .1,
       zlab = "trophic level",
       xlab = "drought index",
       ylab = "mean air temperature",
       box = T
       #,theta = seq360[round]
)



# Check error in testing data
yhat_xg <- predict(xg_mod, dtest)
(MSE_xgb <- mean((yhat_xg - Y_test)^2))