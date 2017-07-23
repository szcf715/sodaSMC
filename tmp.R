## POS = c()

## intialization
n = 5
## the number of samples
m = 100 
## init.x <- function() sample(c(0, 1), 5, replace = TRUE, prob = c(0.7, 0.3))
## m samples
## m = 100
#x0 = sapply(1:m, function(i) init.x())
#y0 = sapply(1:m, function(i) init.x())

x0 = sapply(1:m, function(x) sample(n, 1))
## POS = c(POS, x0)
origin = 1:n
x1 = sapply(1:m, function(i) origin[-x0[i]][sample(n-1, 1)])
weight1 = calc_BIC(x1)

## resampling
id = sample(1:m, m, replace = TRUE, prob = weight1)
x = cbind(x0, x1)[id, ]

x2 = sapply(1:m, function(i) origin[-x[i,]][sample(n-2, 1)])
weight2 = calc_BIC(x2)

id = sample(1:m, m, replace = TRUE, prob = weight2)
x = cbind(x0, x1, x2)[id, ]
