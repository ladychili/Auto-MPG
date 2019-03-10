## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE,
                      include = TRUE,
                      message = FALSE,
                      warning = FALSE,
                      results='hide', 
                      fig.align = 'center',
                      fig.height=3, fig.width=7)

library(tidyverse)
library(knitr)
library(kableExtra)
library(boot)
library(splines)

dat <- read_csv("../data/auto-data.csv")[,-1]
dat <- na.omit(dat)
attach(dat)

## ----load----------------------------------------------------------------
# load("Rscript/env3.RData")
source("func.R")
set.seed(5059)
train <- sample(nrow(dat), nrow(dat)/2)

## ----lm-disp, fig.cap="Candidate Linear Models for Displacement"---------
set.seed(5059) # to fix CV error
myplot(pred ="d", dat, 15)

## ----lm-d-best, results='asis', fig.cap="Selected Linear Models for Displacement", fig.height=4, fig.width=7----
lm.d <- glm(mpg ~ poly(displacement, degree = 10))
plot(mpg ~ displacement, col='grey', pch=20, cex=2, cex.main = 1,
     main = "Linear Model, degree = 10")
lines(x = min(displacement):max(displacement), 
      y = predict(lm.d, data.frame(displacement=min(displacement):max(displacement))), 
      lwd=2, lty=1, col='darkgreen')

cat("RSS: ", rss(lm.d), "\n")
cat("AIC: ", AIC(lm.d), "\n")
set.seed(5059)
lm.d.cv <- cv.glm(data = dat, lm.d,K = 10)
cat("MSE: ", lm.d.cv$delta[1], "\n")

## ----lm-h, fig.cap="Candidate Linear Models for Horsepower"--------------
set.seed(5059) # to fix CV error
myplot(pred ="d", dat, 15)

## ----lm-h-best,  results='asis', fig.cap="Selected Linear Models for Horsepower", fig.height=4, fig.width=7----
lm.h <- glm(mpg ~ poly(horsepower, degree = 8))
plot(mpg ~ horsepower, col='grey', pch=20, cex=2, cex.main = 1,
     main = "Linear Model, degree = 8")
lines(x = min(horsepower):max(horsepower), 
      y = predict(lm.h, data.frame(horsepower=min(horsepower):max(horsepower))), 
      lwd=2, lty=1, col='darkgreen')

cat("RSS: ", rss(lm.h), "\n")
cat("AIC: ", AIC(lm.h), "\n")
set.seed(5059)
lm.h.cv <- cv.glm(data = dat, lm.h,K = 10)
cat("MSE: ", lm.h.cv$delta[1], "\n")

## ----lm-w, fig.cap="Candidate Linear Models for Weight"------------------
set.seed(5059) # to fix CV error
myplot(pred ="w", dat, 15)

## ----lm-w-best,  results='asis', fig.cap="Selected Linear Models for Weight", fig.height=4, fig.width=7----
lm.w <- glm(mpg ~ poly(weight, degree = 2))
plot(mpg ~ weight, col='grey', pch=20, cex=2, cex.main = 1,
     main = "Linear Model, degree = 2")
lines(x = min(weight):max(weight), 
      y = predict(lm.w, data.frame(weight=min(weight):max(weight))), 
      lwd=2, lty=1, col='darkgreen')

cat("RSS: ", rss(lm.w), "\n")
cat("AIC: ", AIC(lm.w), "\n")
set.seed(5059)
lm.w.cv <- cv.glm(data = dat, lm.w,K = 10)
cat("MSE: ", lm.w.cv$delta[1], "\n")

## ----lm-a, fig.cap="Candidate Linear Models for Acceleration"------------
set.seed(5059) # to fix CV error
myplot(pred ="a", dat, 15)

## ----lm-a-best,  results='asis', fig.cap="Selected Linear Models for Acceleration", fig.height=4, fig.width=7----
lm.a <- glm(mpg ~ poly(acceleration, degree = 4))
plot(mpg ~ acceleration, col='grey', pch=20, cex=2, cex.main = 1,
     main = "Linear Model, degree = 4")
lines(x = min(acceleration):max(acceleration), 
      y = predict(lm.a, data.frame(acceleration=min(acceleration):max(acceleration))), 
      lwd=2, lty=1, col='darkgreen')

cat("RSS: ", rss(lm.a), "\n")
cat("AIC: ", AIC(lm.a), "\n")
set.seed(5059)
lm.a.cv <- cv.glm(data = dat, lm.a,K = 10)
cat("MSE: ", lm.a.cv$delta[1], "\n")

## ----lm-sum, results='hold'----------------------------------------------
lmlst <- list(lm.d, lm.h, lm.w, lm.a)
lmcv.lst <- list(lm.d.cv, lm.h.cv, lm.w.cv, lm.a.cv)
lm.rss <- lapply(lmlst, rss)
lm.aic <- lapply(lmlst, AIC)
lm.mse <- lapply(lmcv.lst, function(x) {'$' (x, 'delta') [1]})

sum <- data.frame(RSS = unlist(lm.rss), 
                  AIC = unlist(lm.aic),
                  MSE = unlist(lm.mse))
rownames(sum) <- c("Displacement", "Horsepower", "Weight", "Acceleration")

kable(sum, caption = "Summary of Model Fitting Measurement for the Selected Linear Model",
      digits = 2, booktabs = T, longtable = T)

## ----bin-d, fig.cap="Candidate Bin Smooth Models for Displacement"-------
pickBin("d",200, mybinsmooth)

## ----bin-d-best, results='asis', fig.cap="Selected Bin Smooth Models for Displacement", fig.height=4, fig.width=7----
bin.d <- mybinsmooth("d",binlength = 10, opt = 1)

## ----bin-h, fig.cap="Candidate Bin Smooth Models for Horsepower"---------
pickBin("h",200, mybinsmooth)


## ----bin-h-best, results='asis', fig.cap="Selected Bin Smooth Models for Horsepower", fig.height=4, fig.width=7----
bin.h <- mybinsmooth("d",binlength = 20, opt = 1)

## ----bin-w, fig.cap="Candidate Bin Smooth Models for Weight"-------------
pickBin("w",200, mybinsmooth)

## ----bin-w-best, results='asis', fig.cap="Selected Bin Smooth Models for Weight", fig.height=4, fig.width=7----
bin.w <- mybinsmooth("w", binlength = 20, opt = 1)

## ----bin-a, fig.cap="Candidate Bin Smooth Models for Acceleration"-------
pickBin("a",200, mybinsmooth)


## ----bin-a-best, results='asis', fig.cap="Selected Bin Smooth Models for Acceleration", fig.height=4, fig.width=7----
bin.a <- mybinsmooth("a", binlength = 25, opt = 1)

## ----bin-sum, results='asis'---------------------------------------------
binsum <- as.data.frame(rbind(bin.d, bin.h, bin.w, bin.a))
rownames(binsum) <- c("Displacement", "Horsepower", "Weight", "Acceleration")
colnames(binsum) <- c("Adjusted R2", "AIC", "RSS", "MSE")
binsum <- binsum[c( "RSS", "AIC", "MSE","Adjusted R2")]
kable(binsum, caption = "Summary of Model Fitting Measurement for the Selected Bin Smooth Model",
      digits = 2, booktabs = T, longtable = T)

## ----bs-d, fig.cap="Candidate B-spline Models for Displacement"----------
ppickbs("d",nknots = 10, degree = 3)

## ----bs-d-opt, results='asis', fig.cap="Selected B-spline Models for Displacement", fig.height=4, fig.width=7----
bs.d <- mybspline("d", 7, 3)

## ----bs-h, fig.cap="Candidate B-spline Models for Horsepower"------------
ppickbs("h",nknots = 10, degree = 3)

## ----bs-h-opt, results='asis', fig.cap="Selected B-spline Models for Horsepower", fig.height=4, fig.width=7----
bs.h <- mybspline("h", 5,3)

## ----bs-w, fig.cap="Candidate B-spline Models for Weight"----------------
ppickbs("w",nknots = 10, degree = 3)

## ----bs-w-opt, results='asis', fig.cap="Selected B-spline Models for Weight", fig.height=4, fig.width=7----
bs.w <- mybspline("w", 4, 3)

## ----bs-a, fig.cap="Candidate B-spline Models for Acceleration"----------
ppickbs("a",nknots = 10, degree = 3)

## ----bs-a-opt, results='asis', fig.cap="Selected B-spline Models for Acceleration", fig.height=4, fig.width=7----
bs.a <- mybspline("a",4,2)

## ----bs-sum, results='asis'----------------------------------------------
bssum <- as.data.frame(rbind(bs.d, bs.h, bs.w, bs.a))
rownames(bssum) <- c("Displacement", "Horsepower", "Weight", "Acceleration")
colnames(bssum) <- c("Adjusted R2", "AIC", "RSS", "MSE")
bssum <- bssum[c( "RSS", "AIC", "MSE","Adjusted R2")]
kable(bssum, caption = "Summary of Model Fitting Measurement for the Selected B-spline Model",
      digits = 2, booktabs= TRUE, longtable = TRUE)

