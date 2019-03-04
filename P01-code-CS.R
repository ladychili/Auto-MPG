##### Code for 1st ID 5059 Practical    ########
##### Tom Kelsey                        ########

data = read.table("data/auto-mpg.data", sep=" ", header=T)

attach(data)

#### convert to comma separated values for use in Excel
write.csv(data,"data/auto-data.csv")

############################################################################
###########################################################################

###bin-smooths###


binsmooth = function(x, y, binlen) {
    rss = 0
    plot(x, y)
    prev = NA
    xsorted = sort(x, na.last=NA)
    bins = ceiling(length(xsorted) / binlen)
    for (i in 1:bins) {
        xstart = xsorted[(i-1)*binlen+1]
        xend = xsorted[min(i*binlen, length(xsorted))]
        yval = mean(y[x>=xstart&x<=xend], na.rm=T)
        if(!is.na(prev)) {
            lines(c(xstart, xstart), c(prev, yval))
        }
        if(i < bins) {
            xendd = xsorted[i*binlen+1]
            lines(c(xstart, xendd), rep.int(yval, 2))
        } else {
            lines(c(xstart, xend), rep.int(yval, 2))
        }
        prev = yval
        rss = rss + sum(sapply(y[x>=xstart&x<=xend], function(x) { if(is.na(x)) { 0 } else { (x - yval)^2 } }))
    }
    return(rss)
}

## run the binsmooth function
binsmooth(displacement, mpg, 10)


## return the RSS for a model
rss = function(model) {
    return(sum(sapply(residuals(model), function(x) { x^2 })))
}


## Linear model
plot(mpg~displacement)
model = lm(mpg~displacement)
abline(coef(model))
rss(model)

## Load the library for basis functions
library(splines)

## no knots; default spline degree
splinemodel = lm(mpg~bs(displacement))
rss(splinemodel)
lines(min(displacement):max(displacement), predict(splinemodel, data.frame(displacement=min(displacement):max(displacement))))

## new plot (not meant to be pretty, just show options)
plot(mpg~displacement, bg='red', pch=23, cex=2)

## knots at 25%, 50% & 75% quartiles;  spline degree 3 (cubic spline)
quantile(displacement)
splinemodel2 = lm(mpg~bs(displacement,knots = c(104.25,148.5,262),degree=3))
rss(splinemodel2)
lines(min(displacement):max(displacement), predict(splinemodel2, data.frame(displacement=min(displacement):max(displacement))),lwd=2, lty=2, col='green')

## new plot 
plot(mpg~displacement, pch=3, cex=2)

## knots every 25 x values;  spline degree 1 (piecewise linear)
splinemodel3 = lm(mpg~bs(displacement,knots = seq(1,455,25),degree=1))
rss(splinemodel3)
lines(min(displacement):max(displacement), predict(splinemodel3, data.frame(displacement=min(displacement):max(displacement))),lwd=2, lty=2, col='blue')


## What does this code do?
x = matrix(displacement)
y = matrix(mpg)
solve(t(x) %*% x) %*% t(x) %*% y

###### Cathy's Answer:
###### to identify \hat{beta} using input X and response y of training data - the solution that minimise RSS of beta
###### L01, slides 18 - 20

