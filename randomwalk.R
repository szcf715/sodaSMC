## main effect

## intialization
n = 5
init.x <- function() sample(c(0, 1), 5, replace = TRUE, prob = c(0.7, 0.3))
## m samples
m = 100
x0 = sapply(1:m, function(i) init.x())
#y0 = sapply(1:m, function(i) init.x())

## generate x1
## calculate weights

## TODO
EBIC = calc_BIC(x0)

w = EBIC/sum(EBIC)

## resampling

