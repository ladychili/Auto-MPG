col=c("orange", "darkgreen",  "darkorchid", "darkblue", "brown")

rss <-  function(model) {
  return(sum(sapply(residuals(model), function(x) { x^2 })))
}

myplot <- function(pred, data, degree){
  y <- mpg
  x <- switch(pred,
             d = displacement,
             h = horsepower,
             w = weight,
             a = acceleration)
  
  rss <- numeric(degree)
  aic <- numeric(degree)
  cverr <- numeric(degree)
  
  for (i in 1:degree) {
    
    mod <- glm(y ~ poly(x, degree = i), na.action = 'drop')
    rss[i] <- rss(mod)
    aic[i] <- AIC(mod)
    cverr[i] <- switch (pred,
                        d = cv.glm(data, glm(mpg ~ poly(displacement, degree = i)), K = 10)$delta[1],
                        h = cv.glm(data, glm(mpg ~ poly(horsepower, degree = i)), K = 10)$delta[1],
                        w = cv.glm(data, glm(mpg ~ poly(weight, degree = i)), K = 10)$delta[1],
                        a = cv.glm(data, glm(mpg ~ poly(acceleration, degree = i)), K = 10)$delta[1])
  }
  
  par(mfrow=c(1,3))
  # RSS
  plot(1:degree, rss, type = "o", pch = 20, col = "darkblue", 
       xlab = "Degree of Polynomial", main = "RSS")
  
  abline(v = which(rss == min(rss)), col = "orange")
  # AIC
  plot(1:degree, aic, type = "o", pch = 20, col = "darkblue", 
       xlab = "Degree of Polynomial", main = "AIC")
  abline(v = which(aic == min(aic)), col = "orange")
  # k-fold CV 
  plot(1:degree, cverr, type = "o", pch = 20, col = "darkblue", 
       xlab = "Degree of Polynomial", main = "10-fold CV",
       ylab = "MSE")
  abline(v = which(cverr == min(cverr)), col = "orange")
  par(mfrow=c(1,1))
  
  cat("Min RSS", min(rss), "at degree =",which(rss == min(rss)), "\n")
  cat("Min AIC", min(aic), "at degree =",which(aic == min(aic)), "\n")
  cat("Min CV error", min(cverr), "at degree =",which(cverr == min(cverr)), "\n")
}


mybinsmooth <- function(pred, binlength=0, knotsdef=NULL, output=1, ploto=1, opt=0){

  y <- mpg[train]
  x <- switch(pred,
              d = displacement[train],
              h = horsepower[train],
              w = weight[train],
              a = acceleration[train])
  
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  #Devide data into bins
  if(is.vector(knotsdef)) bins = knotsdef
  else bins = ceiling(length(x) / binlength)
  #Create Design Matrix without intercept
  DM <- matrix(1,length(x),bins)
  #Set all elements not corresponding to region j equal 0
  for(i in 1:bins)
  {
    if(i==1) { xstart = 1 }
    if(i>1) { xstart = (i-1)*binlength+1 }
    xend = min(xstart + binlength-1, length(x))
    binelements <- xstart:xend
    elements <- 1:length(x)
    elements[binelements] <- 0
    DM[elements,i] <- 0
  }
  
  #Perform Linear Regreesion
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
  q <- bins
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)   
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Elements per bin: ", binlength, "\n")
    cat("Number of bins: ", bins, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")

    
    #Graphic 
    yl <- "mpg"
    xl <- switch(pred,
                d = "displacement",
                h = "horsepower",
                w = "weight",
                a = "acceleration")
    if(ploto==1) plot(x,y, col='grey', pch=20, cex=2,
                      main = paste("Bin Smooth, size of bins =", binlength), xlab = xl, ylab = yl)
    j<-1
    for(i in 1:length(coef(reg)))
    {
      if(i>1) lines(c(x[xend],x[xend]), c(as.numeric(coef(reg)[i-1]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xstart = j
      if(i>1) lines(c(x[xend],x[xstart]), c(as.numeric(coef(reg)[i]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xend = min(j+binlength-1, length(x))
      lines(c(x[xstart],x[xend]), rep(as.numeric(coef(reg)[i]), 2), col="red", lwd=2)
      j<-j+binlength
      
    }
  }
  
  if(opt==1) return(c(R2adj,aic,rss))    
}

pickBin <- function(pred, iter, FUN, ...){
  
  R2adj <- numeric(iter)
  aic <- numeric(iter)
  rss <- numeric(iter)
  
  for(i in 1:iter) 
  {
    p<-i
    back <- FUN(pred, p, output=0, opt=1, ...)
    R2adj[i]<-back[1]
    aic[i]<- back[2]
    rss[i] <- back[3]
  }
  
  #!!!!!!!!!!!!!!!The next two lines are not straight forward, but easy to avoid the problem of NaN 
  #One could simply delete NaN, but then the order of the array would be chaos
  R2adj[!complete.cases(R2adj)] <- max(complete.cases(R2adj))-0.1  #!!!!!!!!!!!!!!!
  if(length(back)==2) aic[!complete.cases(aic)] <- max(aic)+10 #!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!
  
  par(mfrow=c(1,2))

  plot(1:iter, rss, type = "l", col = "darkblue", lwd=2,
       xlab = "binlength (size of bin)", main = "Minimized RSS")
  abline(v = which(rss == min(rss)), col = "orange")
  plot(1:iter, aic, type = "l", col = "darkblue", lwd=2,
       xlab = "binlength (size of bin)", main = "Minimized AIC")
  abline(v = which(aic == min(aic[aic>0])), col = "orange")
  par(mfrow=c(1,1))
  
  cat("Highest adj. R-squared: ", max(R2adj), "at binlength =", which(R2adj == max(R2adj)), "\n")
  cat("Min RSS", min(rss), "at binlength =",which(rss == min(rss)), "\n")
  cat("Min AIC", min(aic[aic>0]), "at binlength =",which(aic == min(aic[aic>0])), "\n")
}


#--------------------------------------------------

# Bin Smooth - ?size of bin

# B-spline: knots? degree?


library(splines)
library(tidyverse)
library(boot)

dat <- read_csv("../data/auto-data.csv")[,-1]
dat <- na.omit(dat)
attach(dat)


# Linear model - varying degrees ------------------------------------------

# displacement
set.seed(5059)
myplot(pred ="d", dat, 10)

lm.d <- glm(mpg ~ poly(displacement, degree = 10))
plot(mpg ~ displacement, col='grey', pch=20, cex=2,
     main = "Linear Model, degree = 10")
lines(x = min(displacement):max(displacement), 
      y = predict(lm.d, data.frame(displacement=min(displacement):max(displacement))), 
      lwd=2, lty=1, col='darkgreen')
rss(lm.d)
AIC(lm.d)
cv.glm(data = dat, lm.d,K = 10)$delta[1]

# horsepower

set.seed(5059)
myplot(pred ="h", dat, 10)

lm.h <- glm(mpg ~ poly(horsepower, degree = 8))
plot(mpg ~ horsepower, col='grey', pch=20, cex=2,
     main = "Linear Model, degree = 8")
lines(x = min(horsepower):max(horsepower), 
      y = predict(lm.h, data.frame(horsepower=min(horsepower):max(horsepower))), 
      lwd=2, lty=1, col='darkgreen')
rss(lm.h)
AIC(lm.h)
cv.glm(data = dat, lm.h, K = 10)$delta[1]


# weight
set.seed(5059)
myplot(pred = "w", dat, 10)

lm.w <- glm(mpg ~ poly(weight, degree = 2))
plot(mpg ~ weight, col='grey', pch=20, cex=2,
     main = "Linear Model, degree = 2")
lines(x = min(weight):max(weight), 
      y = predict(lm.w, data.frame(weight=min(weight):max(weight))), 
      lwd=2, lty=1, col='darkgreen')
rss(lm.w)
AIC(lm.w)
cv.glm(data = dat, lm.w, K = 10)$delta[1]

# accelerate
set.seed(5059)
myplot(pred = "a", dat, 10)

lm.a <- glm(mpg ~ poly(acceleration, degree = 4))
plot(mpg ~ acceleration, col='grey', pch=20, cex=2,
     main = "Linear Model, degree = 4")
lines(x = min(acceleration):max(acceleration), 
      y = predict(lm.a, data.frame(acceleration=min(acceleration):max(acceleration))), 
      lwd=2, lty=1, col='darkgreen')
rss(lm.a)
AIC(lm.a)
cv.glm(data = dat, lm.a, K = 10)$delta[1]


# Bin smooth - varying knots ----------------------------------------------
set.seed(5059)
train <- sample(nrow(dat), nrow(dat)/2)


# displacement
pickBin("d",200, mybinsmooth)
mybinsmooth("d",binlength = 11)

# horsepower
pickBin("h",200, mybinsmooth)
mybinsmooth("d",binlength = 4)

# weight
pickBin("w", 200, mybinsmooth)
mybinsmooth("w", binlength = 47)

# acceleration
pickBin("a", 200, mybinsmooth)
mybinsmooth("a", binlength = 3)


# B-spline - varying knots and degrees ------------------------------------
par(mfrow = c(1,1))
plot(mpg~displacement, col='grey', pch=20, cex=2)
abline(v = seq(1,455,25))
## displacement

## default - no knots, d = 3
bs.disp <- lm(mpg ~ bs(displacement))
rss(bs.disp) # 7530.594
AIC(bs.disp) # 2309.705
lines(min(displacement):max(displacement), predict(bs.disp, data.frame(displacement=min(displacement):max(displacement))), lwd=2, lty=2, col='green')

## quantile knots 0.25, 0.5, 0.75, d = 3
bs.disp.1 <- lm(mpg ~ bs(displacement, knots = quantile(displacement, c(0.25, 0.5, 0.75)), degree = 3))
rss(bs.disp.1) # 6786.116
AIC(bs.disp.1) # 2274.276
lines(min(displacement):max(displacement), predict(bs.disp.1, data.frame(displacement=min(displacement):max(displacement))), lwd=2, lty=2, col='blue')

## knots every 25 values, d = 1
bs.disp.2 <- lm(mpg ~ bs(displacement, knots = seq(1,455,30), degree = 3))
rss(bs.disp.2) # 6669.542
AIC(bs.disp.2) # 2289.379
lines(min(displacement):max(displacement), predict(bs.disp.2, data.frame(displacement=min(displacement):max(displacement))), lwd=2, lty=2, col='red')
