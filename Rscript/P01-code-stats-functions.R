######################################
#------------------------------------#
#ID5059 - Assignment 1           #
#Functions                           #
#------------------------------------#
######################################

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

linreg <- function(x, y, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Create a Design Matrix DM
  n <- length(x)
  q <- 2
  DM = matrix(1,n,q)
  DM[,2] <- x
  
  #Perform regression 
  # beta <- solve(t(DM)%*%DM)%*%t(DM)%*%y
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
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
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphic
    xp <- 0:100/100
    n <- length(xp)
    DMp = matrix(1,n,q)
    DMp[,2] <- xp
    
    if(ploto==1) par(mfrow=c(1,2))
    if(ploto==1) matplot(xp, DMp, type="l", lwd=2, main="Individual functions")
    if(ploto==1) plot(x,y, main="Linear regression", pch=20, col="darkblue")
    lines(xp, DMp%*%coef(reg), col="orange", type="l", lwd=2)
  }
  
  if(opt==1) {return(c(R2adj,aic))}   
}

#-------------------------------------
#Polynomial Regression
#-------------------------------------
#The function polreg performs a polynomial regression with p order polynomials 
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       p - nth order polynomial 
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
polreg <- function(x,y,p, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
    x<-x-min(x)
    x<-x/max(x)
  #Create a Design Matrix DM
    n <- length(x)
    q <- p + 1
    DM = matrix(1,n,q)
    DM[,2] <- x
    if(q>2) for(i in 3:q) { DM[,i] <- x**(i-1)}
  #Perform regression 
    reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
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
      cat("Number of polyn.: ", p, "\n")
      cat("RSS: ", rss, "\n")
      cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
      cat("R-squared: ", R2, "\n")
      cat("Adjusted R-squared: ", R2adj, "\n")
      cat("AIC: ", aic, "\n")
      #cat("Coefficents: \n")
      #print(coef(reg))
      #print(summary(reg))
      #print(anova(reg))
      
    #Graphic
      xp <- 0:100/100
      n <- length(xp)
      DMp = matrix(1,n,q)
      DMp[,2] <- xp
      if(q>2) for(i in 3:q) { DMp[,i] <- xp**(i-1)}
      
      if(ploto==1) par(mfrow=c(1,2))
      if(ploto==1) matplot(xp, DMp, type="l", lwd=2, main="Individual polynomial functions")
      if(ploto==1) plot(x,y, main="Polynomial regression", pch=20, col="darkblue")
      lines(xp, DMp%*%coef(reg), col="darkgreen", type="l", lwd=2)
  }
  
  if(opt==1) {return(c(R2adj,aic))}   
}



#-------------------------------------
#Binsmooth 
#-------------------------------------
#The function binsmoothREG performs a binsmooth regression with a user defined binlength
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       binlength - amount of x values per bin
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot



binsmoothREG <- function(x, y, binlength=0, knotsdef=NULL, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
    x<-x-min(x)
    x<-x/max(x)
  #Sort x values in ascending order
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
        #cat("Coefficents: \n")
        #print(coef(reg))
        #print(summary(reg))
        #print(anova(reg))
      
      #Graphic 
        if(ploto==1) plot(x,y, main="Binsmooth regression", pch=20, col="darkblue")
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
    
    if(opt==1) return(c(R2adj,aic))    
}





#-------------------------------------
#Trunceted Regression splines with fixed knots
#-------------------------------------
#The function truncREGDM performs a truncated regression with user defined knots and degree of basis
#Input Arguments truncReg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       knotsdef - user defined position of knots in [0,1]
#       degree - degree of basis
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

truncREGDM <- function(x, knotpos, degree)
{
  
  #Create desing matrix:
  n<-length(x)
  q<-length(knotpos)+2
  DM <- matrix(0,n,q)
  DM[,1] <- 1
  DM[,2] <- x
  
  #Set x values not corresponding to the bin equal 0
  for(i in 3:q)
  {
    DM[,i] <- x - knotpos[i-2]
    elements <- 1:length(x) 
    xstart = length(x[x<knotpos[i-2]])
    DM[(1:xstart-1),i] <- 0                
  }
  
  DM[,2:q] <- DM[,2:q]**degree
  return(DM)
}  


truncReg <- function(x, y, knots, knotsdef=NULL, degree=3, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
    x<-x-min(x)
    x<-x/max(x)
  #Sort x values in ascending order
    y <- y[order(x)]
    x <- sort(x)
    n <- length(x)
  
  #Calculate knot postions
    if(length(knotsdef)>0) knotpos <- knotsdef
    else  knotpos <- 1:knots / (knots+1)
   
    
  #Create Design Matrix
   DM <- truncREGDM(x, knotpos, degree)
    
  #Perform regression
    reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
   
    q <- length(knotpos)+2
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
      cat("Knot positions: ", knotpos,  "\n")
      cat("Degree: ", degree, "\n")
      cat("RSS: ", rss, "\n")
      cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
      cat("R-squared: ", R2, "\n")
      cat("Adjusted R-squared: ", R2adj, "\n")
      cat("AIC: ", aic , "\n")
      #cat("Coefficents: \n")
      #print(coef(reg))
      #print(summary(reg))
      #print(anova(reg))
      
      
      
      #Graphic 
        xp <- 0:100/100
        DMp <- truncREGDM(xp, knotpos, degree)
      
        if(ploto==1)par(mfrow=c(1,2))
        if(ploto==1) matplot(xp, DMp, type="l", lwd=2, main="Individual spline functions")
        if(ploto==1) plot(x,y, main="Truncated Spline Regression", pch=20, col="darkblue")
        
        points(xp,DMp%*%coef(reg), type="l", lwd=2, col="darkorchid")
        
      }
    if(opt==1) return(c(R2adj, aic)) 
}

#-------------------------------------
#Penalized cubic splines and regression splines, Wood 2006
#The following code is heavily  based on the Book "Generlized Additive Models" by Simon N. Wood (2006)
#-------------------------------------
#The function csreg performs a cubic spline regression with user defined knots
#Input Arguments:
#Input Arguments csreg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

#The function prscsreg performs a penalized cubic spline regression with user defined knots
#This is a straight forward extension of the function csreg
#Input Arguments:
#Input Arguments prscsreg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot



#R(x,z) for cubic spline on [0,1]
rk <- function(x,z)
{
  ((z-0.5)**2 - 1/12) * ((x-0.5)**2 - 1/12)/4 - ((abs(x-z)-0.5)**4 - (abs(x-z)-0.5)**2 / 2+7/240) / 24 
}

#Set up the penalized regression spline penalty matrix
spl.S <- function(xk)
{
  q<- length(xk)+2
  S<-matrix(0,q,q)
  S[3:q,3:q] <- outer(xk,xk,FUN=rk)
  
  return(S)
  
}

#Set up model matrix for cubic penalized regression spline
spl.X <- function(x,xk)
{
  q<-length(xk)+2
  n<-length(x)
  X<-matrix(1,n,q)
  X[,2]<-x
  X[,3:q]<-outer(x,xk,FUN=rk)
  return(X)
}

#Function for calculation the square root of a matrix
mat.sqrt <- function(S)
{
  d<-eigen(S,symmetric=T)
  rS<-d$vectors%*%diag(d$values**0.5)%*%t(d$vectors)
  return(rS)
}




csreg <- function(x, y, knots, knotsdef=NULL, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
    x<-x-min(x)
    x<-x/max(x)
  #Sort x values in ascending order
    y <- y[order(x)]
    x <- sort(x)
    n <- length(x)
  
  #Calculate knot postions
    if(length(knotsdef)>0) knotpos <- knotsdef
    else  knotpos <- 1:knots / (knots+1)
  
  #Create Design Matrix
    DM<-spl.X(x,knotpos)
   
  #Perfrom regression
    reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
    q <- knots
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
      cat("Number of knots = ", knots, "\n")
      cat("Knot positions = ", knotpos, "\n")
      cat("RSS: ", rss, "\n")
      cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
      cat("R-squared: ", R2, "\n")
      cat("Adjusted R-squared: ", R2adj, "\n")
      cat("AIC: ", aic , "\n")
      #cat("Coefficents: \n")
      #print(coef(reg))
      #print(summary(reg))
      #print(anova(reg))
      
      #Graphics
      
      #Values for prediction
        xp <- 0:100/100
        DMp <- spl.X(xp,knotpos)
      
      if(ploto==1)par(mfrow=c(1,2))
      if(ploto==1) matplot(xp, (DMp), type="l", lwd=2, main="Individual spline functions", ylim=c(-0.05,0.05))
      if(ploto==1) plot(x,y, main="Cubic Spline Regression", pch=20, col="darkblue")
      
      points(xp,DMp%*%coef(reg), type="l", lwd=2, col="green")
      
    }
    if(opt==1) return(c(R2adj, aic))     
}


#Function for performing penalized regression 
prs.fit <- function(x,y,knotpos,lambda)
{
  q<-length(knotpos)+2
  n<-length(x)
  DM <- rbind(spl.X(x,knotpos), mat.sqrt(spl.S(knotpos))*sqrt(lambda))
  y[(n+1):(n+q)] <- 0
  reg <- lm(y ~ 0 + DM)
  return(reg)
}


prscsreg <- function(x,y, knots, knotsdef=NULL, lambda, output=1, ploto=1, opt=0, indvploto=0)
{
  #Scale data to [0,1]
    x<-x-min(x)
    x<-x/max(x)
  #Sort x values in ascending order
    y <- y[order(x)]
    x <- sort(x)
    n <- length(x)
  
  ##Calculate knot postions
    if(knots == 0) knotpos <- NULL
    if(knots != 0) knotpos <- 1:knots / (knots+1)
    if(length(knotsdef)>0) knotpos <- knotsdef
  
  #Perform penalized regression
    reg <- prs.fit(x,y,knotpos,lambda)
  
  
  #Calculate goodness of fit measures
    q <- length(knotpos)
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
      cat("Number of knots = ", knots, "\n")
      cat("RSS: ", rss, "\n")
      cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
      cat("R-squared: ", R2, "\n")
      #cat("Adjusted R-squared: ", R2adj, "\n")
      cat("AIC: ", aic , "\n")
      #cat("Coefficents: \n")
      #print(coef(reg))
      #print(summary(reg))
      #print(anova(reg))
      
      #Graphics
      
      #Values for prediction
      xp <- 0:100/100
      DMp <- spl.X(xp,knotpos)
      
      if(indvploto==1)par(mfrow=c(1,2))
      if(indvploto==1) matplot(xp, (DMp), type="l", lwd=2, main="Individual spline functions", ylim=c(0,0.05))
      if(ploto==1) plot(x,y, main="Penalized Cubic Spline Regression", pch=20, col="darkblue")
      
      points(xp,DMp%*%coef(reg), type="l", lwd=2, col="brown")
      
    }
if(opt==1) return(c(R2adj, aic))       
}

#Calculate gerealized cross validation
GCV <- function(x,y,knots)
{
  x<-x-min(x)
  x<-x/max(x)
  xk <- 1:knots/(knots+1)
  
  lambda <- 1e-9
  n <- length(x)
  V<- rep(0,50)
  lambdai <- rep(0,60)
  
  for(i in 1:60)
  {
    mod <-prs.fit(x,y,xk, lambda)
    trA <- sum(influence(mod)$hat[1:n])
    rss <-  sum(sapply(residuals(mod), function(x) { x^2 }))
    rss <- sum((y-fitted(mod)[1:n])^2)
    V[i]<- n*rss / ((n-trA)**2)
    
    lambdai[i]<-lambda
    lambda <- lambda*1.5
    
  }
  
  #Summary
  minV <- min(V)
  posi <- which(V == minV)
  optlambda <- lambdai[posi]
  cat("Lowest GCV score is: ", minV, "\n")
  cat("Optimal smoothing parameter lamdda: ", optlambda, "\n")
  
  plot(1:60, log(V), type="l", main="Generalized Cross Validation Score", xlab="i")
  
}



#-------------------------------------
#Bspline regression 
#The following code is heavily based on code of
#Dr. Samiran Sinha (Department of Statistics, Texas A&M)
#http://www.stat.tamu.edu/~sinha/research/note1.pdf 
#And Jeffrey S. Racine: A PRIMER ON REGRESSION SPLINES
#-------------------------------------
#The function bsplinereg performs a bspline regression with user defined knots
#Input Arguments:
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

#Calculate basis (rekursiv)
basis <- function(x, degree, i, knots) 
{ 
    
  if(degree == 0)
  { B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0) 
  } else { 
    if((knots[degree+i] - knots[i]) == 0) 
    { alpha1 <- 0 
    } else { 
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i]) } 
    
    if((knots[i+degree+1] - knots[i+1]) == 0) 
    { alpha2 <- 0 
    } else { alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1]) } 
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots) 
  } 
  
  return(B) 
}

#Create bspline Desin Matrix
bspline <- function(x, degree, knotpos) 
{ 
  #Number of basis
    K <- length(knotpos) + degree + 1 
  #Number of observations
    n <- length(x)
  #Set Boundary knots
    Boundary.knots = c(0,1)
  #create new vector with knot positons 
    knotpos <- c(rep(Boundary.knots[1], (degree+1)), knotpos, rep(Boundary.knots[2], (degree+1))) 

    
  #Create design matrix
    DM <- matrix(0,n,K) 
    for(j in 1:K) DM[,j] <- basis(x, degree, j, knotpos) 
    if(any(x == Boundary.knots[2])) DM[x == Boundary.knots[2], K] <- 1 
  #Return DM  
    return(DM) 
}


bsplinereg <- function(x, y, knots=0, knotsdef=NULL,  degree, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Sort x values in ascending order
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  
  #Calculate knot postions
  if(knots == 0) knotpos <- NULL
  if(knots != 0) knotpos <- 1:knots / (knots+1)
  if(length(knotsdef)>0) knotpos <- knotsdef
  #Create Design Matrix 
    DM <- bspline(x, degree, knotpos) 

  #Perform penalized regression
    reg <- lm(y ~ 0 + DM)
    print(summary(reg))
  
  #Calculate goodness of fit measures
    q <- length(knotpos) + degree + 1
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
    cat("Number of knots = ", knots, "\n")
    cat("Knot positions = ", knotpos, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic , "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphics
    
    #Values for prediction
      xp <- 0:100/100
      DMp <- bspline(xp, degree, knotpos)
    
      if(ploto==1)par(mfrow=c(1,2))
      if(ploto==1) matplot(xp, (DMp), type="l", lwd=2, main="Individual spline functions")
      if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
      if(ploto==1) plot(x,y, main="BSpline Regression", pch=20, col="darkblue")
      
      points(xp,DMp%*%coef(reg), type="l", lwd=2, col="brown")
      if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
  }
  if(opt==1) return(c(R2adj, aic))       
}


#-------------------------------------
#Highest adjusted R-squared
#-------------------------------------
#The function highestAdjR2 is determining the highest adjusted R-squared and AIC
#Input Arguments:
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       iter - number of iterations
#       FUN - Function which should be analyzed 
#       ...  optional FUN specif arguments
highestAdjR2 <- function(x, y, iter, FUN, ...)
{

  R2adj <- numeric(iter)
  aic <- numeric(iter)
  for(i in 1:iter) 
  {
    p<-i
    back <- FUN(x,y,p, output=0, opt=1, ...)
    
    if(length(back)==2) {R2adj[i]<-back[1]; aic[i]<-back[2]}
    if(length(back)==1) {R2adj[i]<-back}
  }
  #!!!!!!!!!!!!!!!The next two lines are not straight forward, but easy to avoid the problem of NaN 
  #One could simply delete NaN, but then the order of the array would be chaos
    R2adj[!complete.cases(R2adj)] <- max(complete.cases(R2adj))-0.1  #!!!!!!!!!!!!!!!
    if(length(back)==2) aic[!complete.cases(aic)] <- max(aic)+10 #!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!
  if(length(back)==2) par(mfrow=c(1,2))
  plot(1:iter,R2adj, type="l", main="Maximized adj. R-squared", col="darkblue", pch=20, lwd=2, xlab="i", ylab="adj. R-squared")
  if(length(back)==2)  plot(1:iter,aic, type="l", main="Minimized AIC", col="darkblue", pch=20, lwd=2, xlab="i", ylab="AIC")
  cat("Highest adj. R-squared: ", max(R2adj), "\n")
  cat("i: ", which(R2adj == max(R2adj)), "\n")
  return ((1:iter)[R2adj==max(R2adj)])
}

#-------------------------------------
#Simulated data
#-------------------------------------
#Input Arguments:
#       n - number of data points
#       noise - amount of noises added to the model
#       opt   - 1: returns adj R-squared, 0: returns nothing
simdata <- function(n, noise, opt=0)
{
  x <- seq(-3,3,by=0.001)
  x <- sample(x=x, size=n, replace=T)
  noise = rnorm(length(x), 0, noise)
  y <- I(x**7) - 14*I(x**5) + 49*I(x**3) - 36*x + noise
  x<-x-min(x)
  x<-x/max(x)
  
  #Underlying function
    xt <- seq(-3,3,by=0.001)
    yt <- I(xt**7) - 14*I(xt**5) + 49*I(xt**3) - 36*xt
    xt<-xt-min(xt)
    xt<-xt/max(xt)

  if(opt==1) {
  #Plot  
    xt<-xt-min(xt)
    xt<-xt/max(xt)
    plot(x,y, col="darkblue", pch=20)
    lines(xt, yt, col="black", lwd=2, type="l")
  }
    
    data <- cbind(x,y)
  return(data)
}


#-------------------------------------
#Bias-Variance Tradeoff 
#-------------------------------------
#Tradeoff for Bspline model
#Input Arguments bsplinereg(...):
#       n - number of data points
#       noise - noise added to the underlying model
#       knots - number of knots in [0,1]
#       sims - number of simulations

VarBiasTrad <- function(n, noise, knots, sims) 
{
  #Calculate Knotpostions
  knotpos <- 1:knots/(knots+1)
  #Set xvals for prediction    
  xvals <- seq(0,1,length=1000)
  #Initalize vector for mean prediction   
  meanpredict <- numeric(1000)
  
  #Graphic I
  par(mfrow=c(1,1))
  data <- simdata(n, noise, 0)
  plot(data[,1],data[,2], main="Bias-Variance Tradeoff ", col="white", ylab="y", xlab="x")    
  
  print("okay")
  
  #Start simulation    
  for(i in 1:sims) {
    #Generate data
      data <- simdata(n, noise, 0)
      y <- data[,2]
      x <- data[,1]
      x<-x-min(x)
      x<-x/max(x)
    #Esimate Bspline model
      BSsplinemod = lm(y~bs(x, degree=3, knots=knotpos, Boundary.knots=c(0,1), intercept = T))
      
    #Plot model
      lines(xvals, predict(BSsplinemod, newd = data.frame(x=xvals) ),lwd=0.5, lty=1, col='lightgray')
      
    #Add predictions to meanpredict
      meanpredict <- meanpredict + predict(BSsplinemod, newd = data.frame(x=xvals))
  }
  
  
  #Graphic II
  legend(x="topright", legend=c("underlying function", "fitted models for simulated data", "mean prediction of fitted models"), pch=rep(c("-"), 3), col=c("black", "gray", "red"), lwd=3)
  #Underlying function
    xt <- seq(-3,3,by=0.001)
    yt <- I(xt**7) - 14*I(xt**5) + 49*I(xt**3) - 36*xt
    xt<-xt-min(xt)
    xt<-xt/max(xt)
    lines(xt, yt, col="black", lwd=2, type="l")
  

  #Graphic III
  meanpredict <- meanpredict/sims
  lines(xvals, meanpredict ,lwd=3, col='red')
  
}


