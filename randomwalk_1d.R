library(MASS)
library(nnet)

calc_lda_BIC = function(xx, yy, cur_set, D, K, debug=F, gam=0)
{
  N = nrow(xx);
  D = ncol(xx);
  K = max(yy);
  d = length(cur_set);  
  
  ll = 0;
  if (d == 0)
  {
    p = numeric(K);
    for(i in 1:N)
    {
      p[yy[i]] = p[yy[i]] + 1.0/N;
    }
    for(k in 1:K)
    {
      ll = ll + sum(yy==k)*log(p[k]);
    }
    BIC = -2*ll + (K-1)*(log(N) + 2*gam*log(D));
    return(BIC);
  } 
  else 
  {
    lgt = multinom(yy ~ as.matrix(xx[,cur_set]), family = "binomial", trace=F);
    BIC = lgt$deviance + (K-1)*(1+d)*(log(N) + 2*gam*log(D));
    return(BIC);
  }
}


#
#      xx: explanatory variables
#      yy: response variable
# cur_set: current set of selected variables
#   debug: if shows debug information
#     gam: gamma in EBIC
#   terms: selected linear and interaction terms
#
calc_BIC = function(xx, yy, terms, debug=F, gam=0)
{
  N = length(yy);
  D = ncol(xx);
  K = max(yy);
  d = length(terms);  
  
  ll = 0;
  if (d == 0)
  {
    p = numeric(K);
    for(i in 1:N)
    {
      p[yy[i]] = p[yy[i]] + 1.0/N;
    }
    for(k in 1:K)
    {
      ll = ll + sum(yy==k)*log(p[k]);
    }
    BIC = (K-1) * (log(N) + 2*gam*log(D));
    BIC = BIC - 2*ll;
    return(BIC);
  } 
  else
  {
    #pmatrix = create_pmatrix_from_terms(xx, terms);
    #pmatrix = xx[,terms] # NOT WORK
    pmatrix = as.matrix(xx[,terms])
    # as.factor 效果一样
    #cat(is.factor(yy))
    #yy[yy == 2] <- 10
    
    #lgtt = multinom(as.factor(yy) ~ pmatrix, family = "multinomial", trace=F);
    lgt = multinom(yy ~ pmatrix, family = "multinomial", trace=F);
    BIC = lgt$deviance;
    #BICt = lgtt$deviance;
    #cat(paste0("lgtt: ",BICt, "lgt: ", BIC,"\n"))
    BIC = BIC + (K-1)*(1+ncol(pmatrix))*(log(N) + 2*gam*log(D)); # BIC with quadratic penalty
    
    return(BIC);
  }
}


## dataset
#data = genDataset(100, 1)
data = genDataset2(100)
xdata = data$X
ydata = data$Y
## main effect

## resampling

## 
k = 0
n = 50
m = 50
origin = 1:n
## initialization
x = sapply(1:m, function(x) sample(n, 1))
x = as.matrix(x)
while(TRUE)
{
  if (k >= 2)
    break
  if (k != 0)
  {
  xx = sapply(1:m, function(i) origin[-x[i,]][sample(n-k, 1)])
  ## before sampling
  xx0 = xx
  x0 = x
  x = cbind(x, xx)
  # weight = calc_lda_BIC(xx)
  # weight = calc_BIC(xdata[xx], ydata)
  }
  weight = sapply(1:m, function(i) calc_BIC(xdata, ydata, x[i, ]))
  
  ## scale 
  weight = weight/sum(weight)
  cat(weight)
  ## resampling according to the weights
  id = sample(1:m, m, replace = TRUE, prob = weight)
  ## after sampling
  ## x = cbind(x0, xx0[id])
  x = x[id, ]
  if (!is.matrix(x))
    x = as.matrix(x)
  k = k + 1
}


