# correct TP based on interannual phenology
mod91 <- lm(tp~MONTH, filter(total, year == 1991) ) # we have full coverage accross the year for a single year in 1991
beta <- mod91$coefficients[2] %>% as.numeric()      # extract beta coeficient from the model fit to 1991 data

str(total)

# As you dont provide any data, here some simulated data:
x <- filter(total, year == 1991)$MONTH
y <- filter(total, year == 1991)$tp

# Plot of the noised data (y) and the noiseless line (y_clean) which you want to approximate through your model:
plot(x, y, ylim = c(-2,4))

# Model estimation and adding of the fitted values to the previous plot:
model <- nls(y~a*sin(b*x+c), 
             start = list(a=1,b=1.5, c=1),
             control = list(maxiter = 5000))

points(x, fitted(model), col = "red")

legend("topright", col=c("green", "red"), legend = c("\"true\" curve", "fitted line"), bty = "n", lty = 1)





# circular kernel regression
library(NPCirc)

dir <- filter(total, year == 1991)$MONTH
vel <- filter(total, year == 1991)$tp

estLL <- kern.reg.circ.lin(x = dir, y = vel, method = "LL")
estNW <- kern.reg.circ.lin(dir, vel, method = "NW")

res <- plot(estNW, plot.type = "circle", points.plot = TRUE,
      labels = c("JAN", "FEB", "MAR", "APR", "MAY", "JUNE", "JULY", "AUG", "SEP", "OCT", "NOV", "DEC"),
      label.pos = seq(0, 11 * pi/6, by = pi/6),
      zero = pi/2, clockwise = TRUE, lwd = 2, line.col = 2, main = "")

lines(estLL, plot.type = "circle", plot.info = res, lwd = 2,
       line.col = 3)


library(circular)
# Generate a data set of dependent circular variables.

#x <- circular(runif(50, 0, 2*pi))
#y <- atan2(0.15*cos(x) + 0.25*sin(x), 0.35*sin(x)) + rvonmises(n=50, mu=circular(0), kappa=5)
#x <- circular(total$MONTH, units = "hours")
#x <- circular(total$MONTH)

#x = c(filter(total, year == 1991)$MONTH,rep(13,length(filter(total,year==1991 & MONTH == 1)$tp)))
#y = c(filter(total, year == 1991)$tp, filter(total,year==1991 & MONTH == 1)$tp)

x = c(filter(total, year == 1991)$MONTH,12) + -1
y = c(filter(total, year == 1991)$tp, 2.1)



#y <- total$tp

# Fit a circular-circular regression model.
circ.lm <- lm.circular(y, x, order=1)

# Obtain a crude plot of the data and fitted regression line.
plot.default(x, y, lwd = .25 )

#circ.lm$fitted[circ.lm$fitted>pi] <- circ.lm$fitted[circ.lm$fitted>pi] - 2*pi 

points.default(x[order(x)], circ.lm$fitted[order(x)], type='l', lwd = 3, col = "red4")

#points.default(x = c(0,x[order(x)]), y = c(circ.lm$fitted[order(x)][length(x)],circ.lm$fitted[order(x)]), type='l', lwd = 3, col = "red4" )

lm_all_year <- predict(mod91, newdata = data.frame(MONTH = seq(0,11, by = 1) ) )

lines(x= seq(0,11, by = 1),y=lm_all_year, col= "green3", lwd = 3)

#lines(x= filter(total, year == 1991)$MONTH,y=prediction$fit, col= "green3", lwd = 3)

#plot(circ.lm$fitted)
#plot(x[order(x)], circ.lm$fitted[order(x)], type='circle', lwd = 3, col = "red4")



ggplot()+
  geom_point(aes(x = x, y = y), size = 2, shape = 21,  alpha = .5, color = "white")+
  
  # lm
  geom_line(aes(x= seq(0,11, by = 1),y=lm_all_year), col = "green3", size = 3)+
  #geom_line(aes(x= filter(total, year == 1991)$MONTH,y=prediction$fit), col = "green3", size = 3)+
  
  #geom_smooth(aes(x = total$MONTH, y = total$tp))+
  #geom_line(aes(x = c(0,x[order(x)]), y = c(circ.lm$fitted[order(x)][length(x)],circ.lm$fitted[order(x)]) ), color = "red2", size = 3)+
  #geom_line(aes(x = c(x[order(x)]), y = c(circ.lm$fitted[order(x)]) ), color = "red2", size = 3)+
  #geom_line(aes(x = c(x,0), y = c(circ.lm$fitted,circ.lm$fitted[order(x)][length(x)]) ), color = "red2", size = 3)+
  geom_line(aes(x = x, y = circ.lm$fitted), color = "red2", size = 3)+
  
  coord_polar()+
  #scale_x_continuous(breaks = seq(1,12,1))+
  theme_black()

order(unique(x))
unique( circ.lm$fitted[order(x)] )

#  how to do we interpret these coefs?

library(CircStats)

# Generate a data set of dependent circular variables.
data1 <- runif(50, 0, 2*pi)
data2 <- atan2(0.15*cos(data1) + 0.25*sin(data1), 0.35*sin(data1)) + rvm(50, 0, 5)


data1 <- c(filter(total, year == 1991)$MONTH,13)
data2 <- c(filter(total, year == 1991)$tp, 2.1)


# Fit a circular regression model.
circ.lm <- circ.reg(data1, data2, order=1)
# Obtain a crude plot a data and fitted regression line.
plot(data1, data2)
circ.lm$fitted[circ.lm$fitted>pi] <- circ.lm$fitted[circ.lm$fitted>pi] - 2*pi 

points(data1[order(data1)], circ.lm$fitted[order(data1)], type='l')



