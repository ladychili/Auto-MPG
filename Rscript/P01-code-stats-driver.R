######################################
#------------------------------------#
#ID5059 - Assignment 1           #
#Driver code                         #
#------------------------------------#
######################################

#-------------------------------------
#libs
#-------------------------------------
library(splines)
library(psych)
library(car)
library(mgcv)

#-------------------------------------
#Load functions!!
#-------------------------------------
source('P01-code-stats-functions.R')

#-------------------------------------
#simulated data 
#-------------------------------------
#Input Arguments:
#       n - number of data points
#       noise - amount of noises added to the model
#       opt   - 1: returns adj R-squared, 0: returns nothing
data <- simdata(200,20, 1)
x <- data[,1]
y <- data[,2]



#-------------------------------------
#data
#-------------------------------------
data = read.csv("../data/auto-data.csv", sep=",", header=T)[,-1]
attach(data)

#MPG
  y <- mpg
#displacement
  x <- displacement
#horsepower
  x <- horsepower
# weight
  x <- weight
#acceleration
  x <-acceleration

#-------------------------------------
#Descriptive Statistics
#-------------------------------------
pairs(mpg~displacement + horsepower + weight + acceleration, col="darkblue", pch=20, main="Scatterplot Matrices")
cor(data[,-c(2,7,8,9)])
describe(data[,-c(2,7,8,9)])
#histograms
  hist(mpg, col="red", main="Histogram of mpg")
  qqPlot(mpg, col="darkblue", main="QQ-Plot of mpg", pch=20)
  hist(log(mpg))
  qqPlot(log(mpg))


#Scale
x<-x-min(x)
x<-x/max(x)
plot(x,y, col="darkblue", pch=20)

#-------------------------------------
#Linear Regression
#-------------------------------------
#The function linreg performs a linear regression
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
linreg(x, y, ploto=1)
polreg(x, y, 1, ploto =1)

#-------------------------------------
#Polynomial Regression
#-------------------------------------
#Input Arguments polreg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       p - nth order polynomial
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
p <- highestAdjR2(x,y, iter=25, polreg)
polreg(x,y,10, ploto=1)
#----------------------------


#-------------------------------------
#Binsmooth
#-------------------------------------
#Input Arguments binsmoothREG(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       binlength - size of bins
#       knotsdef - NULL
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
p <- highestAdjR2(x,y,iter=200, binsmoothREG)
binsmoothREG(x,y,60, ploto=1)


#-------------------------------------
#Truncated Spline Regression
#-------------------------------------
#Input Arguments truncReg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       knotsdef - user defined position of knots in [0,1]
#       degree - degree of basis
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
p <- highestAdjR2(x, y, iter=25, truncReg, degree=3)
truncReg(x,y,knots=24, knotsdef=c(), degree=3, ploto=1)

#-------------------------------------
#Cubic Spline Regression
#-------------------------------------
#Input Arguments csreg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       knotsdef - user defined position of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
p <- highestAdjR2(x, y, iter=25, csreg)
csreg(x, y,knots=24, knotsdef=c(), ploto=1)


#-------------------------------------
#Penalized Cubic Spline Regression
#-------------------------------------
#Input Arguments csreg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       lambda - smoothing parameter
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
GCV(x,y,10)
prscsreg(x, y, 10, knotsdef=c(), 0.01105733  , ploto=1, indvploto=1)

#-------------------------------------
#Bspline Regression
#-------------------------------------
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       knotsdef - user defined position of knots in [0,1]
#       degree - degree of basis
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
p <- highestAdjR2(x, y,iter=25, bsplinereg, degree=3)
bsplinereg(x, y, knots=9, knotsdef=c(), degree=3, ploto=1)

# knots = 7, RSS = 6818, AdjR2 = 0.7115, AIC = 2284.179
# knots = 9, RSS = 6724, AdjR2 = 0.7141, AIC = 2282.617

#-------------------------------------
#Legend
#-------------------------------------
legend(x="topleft", legend=c("linear", "polynomial", "binsmooth", "truncated", "cubic", "bspline"), pch=rep(c("-"), 6), col=c("orange", "darkgreen", "red",  "darkorchid", "green", "brown"), lwd=3)

legend(x="topright", legend=c("piecewise polynomials", "true underlying function"), pch=rep(c("-"), 2), col=c("red", "black"), lwd=3)


#-------------------------------------
#Bias-Variance Tradeoff 
#-------------------------------------
#Tradeoff for Bspline model
#Input Arguments bsplinereg(...):
#       n - number of data points
#       noise - noise added to the underlying model
#       knots - number of knots in [0,1]
#       sims - number of simulations

VarBiasTrad(n=200, noise=20, knots=3, sims=500)

################################################################
################################################################
################################################################
################################################################
#-------------------------------------
#Natural cubic splines with Spline package 
#-------------------------------------
kn <- 2   #Number of Knots
knotpos <- 1:kn/(kn+1)
NCsplinemod = lm(y~ns(x, knots=knotpos, intercept = T))
cat("RSS: ", sum(sapply(residuals(NCsplinemod), function(x) { x^2 })), "\n")
plot(x,y)
xvals <- seq(0,1,length=1000)
lines(xvals, predict(NCsplinemod, newd = data.frame(x=xvals) ),lwd=2, lty=2, col='red')

#-------------------------------------
#Bspline with Spline package 
#-------------------------------------
kn <- 2   #Number of Knots
knotpos <- 1:kn/(kn+1)
BSsplinemod = lm(y~bs(x, degree=3, knots=knotpos, Boundary.knots=c(0,1), intercept = T))
cat("RSS: ", sum(sapply(residuals(BSsplinemod), function(x) { x^2 })), "\n")
xvals <- seq(0,1,length=1000)
lines(xvals, predict(BSsplinemod, newd = data.frame(x=xvals) ),lwd=2, lty=2, col='red')

#-------------------------------------
#Smooth splines 
#-------------------------------------
Splinemode = smooth.spline(x=x,y=y, nknots=10, cv=F)
Splinemode
cat("RSS: ", sum(sapply(residuals(Splinemode), function(x) { x^2 })), "\n")
xvals <- seq(0,1,length=100)
lines(xvals, predict(Splinemode, newd = data.frame(x=xvals) ),lwd=2, lty=2, col='red')

#-------------------------------------
#Smooth splines GAM
#-------------------------------------
gamMod = gam(y ~ 1 + s(x, bs="ps", k=4, fx=F), knots=list(c(0.5)), method="GCV.Cp")
gamMod
gamMod$sp
cat("RSS: ", sum(sapply(residuals(gamMod), function(x) { x^2 })), "\n")
xvals <- seq(0,1,length=1000)
lines(xvals, predict(gamMod, data.frame(x=xvals), type="response"),lwd=2, lty=2, col='red')
plot(gamMod, residuals=T)
gam.check(gamMod)
