library(MASS)
library(nnet)


# create predictor matrix from terms
create_pmatrix_from_terms = function(xx, terms)
{
  nt = length(terms);
  nr = nrow(xx);
  pmatrix = matrix(0, nr, 0);
  
  if (nt > 0)
    for(it in 1:nt)
    {
      term = terms[it];
      if (grepl("*",term,fixed=T))
      {
        splits = strsplit(term,"*",fixed=T)[[1]];
        id1 = as.numeric(splits[1]);
        id2 = as.numeric(splits[2]);
        pmatrix = cbind(pmatrix, term=xx[,id1]*xx[,id2]);
      }
      else
      {
        id  = as.numeric(term);
        pmatrix = cbind(pmatrix, term=xx[,id]);
      }
    }
  return(pmatrix);
}


calc_lda_BIC = function(xx, yy, cur_set, debug=F, gam=0)
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
    pmatrix = create_pmatrix_from_terms(xx, terms);
    #pmatrix = xx[,terms] # NOT WORK
    #pmatrix = as.matrix(xx[,terms])
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
data = genDataset(100, 3)
#data = genDataset2(100)
xdata = data$X
ydata = data$Y
## main effect

## resampling

## 
k = 0
n = 50
m = 100
# origin = 1:n
origin = expand.grid(0:n, 0:n)
nc = nrow(origin)
## initialization
#x = sapply(1:m, function(x) sample(n, 1))
#x = rep(c(1:n), m/n)
#x = as.matrix(x)

## interaction
x = vector(mode = "list", m)
for (i in 1:m)
{
  x[[i]] = matrix(c(sample(0:n, 1), sample(0:n, 1)),
                    1, 2, byrow = T)
}
bic.min = 1e10
while(TRUE)
{
  if (k >= 6)
    break
  if (k != 0)
  {
    ## xx = sapply(1:m, function(i) origin[-x[i,]][sample(n-k, 1)])
    ## before sampling
    for (i in 1:m)
    {
      tmp = x[[i]]
      xx = origin[sample(1:nc, 1),]
      flag = colSums(apply(tmp, 1, function(x) x == xx))
      flagflag = sum(flag == 2)
      while(flagflag != 0)
      {
        xx = origin[sample(0:n, 1),]
        flag = colSums(apply(tmp, 1, function(x) x == xx))
        flagflag = sum(flag == 2)
      }
      x[[i]] = rbind(tmp, as.matrix(xx))
    }
    ## xx0 = xx
    ## x0 = x
    # weight = calc_lda_BIC(xx)
    # weight = calc_BIC(xdata[xx], ydata)
  }
  ## weight = sapply(1:m, function(i) calc_BIC(xdata, ydata, x[i, ]))
  weight = c()
  for (i in 1:m)
  {
    tmp = x[[i]]
    terms = c()
    for (j in 1:nrow(tmp))
    {
      if(tmp[j,1] == 0)
      {
        if (tmp[j, 2] == 0)
          break ## TODO
        else
          terms = c(terms, tmp[j, 2])
      }
      else
      {
        if (tmp[j, 2] == 0)
          terms = c(terms, tmp[j, 1])
        else
          terms = c(terms, paste0(tmp[j, 1],'*', tmp[j, 2]))
      }
    }
    weight = c(weight, calc_BIC(xdata, ydata, terms))
  }
  bic.min.tmp = min(weight)
  if (bic.min > bic.min.tmp)
  {
    bic.min = bic.min.tmp
  }
  else
  {
    cat("Terminal!")
  #  x = x0
  #  break
  }
  weight = 1/order(weight)
  ## scale 
  weight = weight/sum(weight)
  cat(weight)
  ## resampling according to the weights
  id = sample(1:m, m, replace = TRUE, prob = weight)
  ## after sampling
  ## x = cbind(x0, xx0[id])
  x = x[id]
#  if (!is.matrix(x))
#    x = as.matrix(x)
  k = k + 1
}
