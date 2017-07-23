# simulation
mu1 = c(0.5, 0, 0)
mu2 = c(-0.5, 0, 0)
omega = matrix(c(-0.6, -0.35, 0,
                 -0.35, 0, -0.35,
                 0, -0.35, -0.60), 3, 3)
omega1 = diag(1, 3, 3) - omega
omega2 = diag(1, 3, 3) + omega

# ## class 1
# df1 = mvrnorm(100, mu1, solve(omega1))
# ## class 0
# df0 = mvrnorm(100, mu2, solve(omega2))
# 
# df = as.data.frame(rbind(df1, df0))
# df$class = c(rep(1,100),rep(0,100))

QX <- function(x)
{
  x1 = x[1]
  x2 = x[2]
  x3 = x[3]
  score = 1.627+x1-0.6*x1*x1-0.6*x3*x3-0.7*x1*x2-0.7*x2*x3
  if (score > 0)
    class = 2
  else
    class = 1
  res = list(score = score, class = class)
  return(res)
}
# trueclass = df[,4]
# res = apply(df[,1:3], 1, function(x) QX(x)$class)
# 
# ## figure 1 left
# plot(df[1:100,1],df[1:100,2],xlim = c(-6,6),ylim = c(-4,4))
# points(df[101:200,1],df[101:200,2], col = "red")

## generate dataset
genDataset <- function(n)
{
  ## class 1
  df1 = mvrnorm(n, mu1, solve(omega1))
  ## class 0
  df0 = mvrnorm(n, mu2, solve(omega2))
  df = as.data.frame(rbind(df1, df0))
  
  trueclass = c(rep(2,n),rep(1,n))
  res = list(df = df, trueclass = trueclass)
  return(res)
}

example <- function(df, which.example)
{
  n = nrow(df)
  Xkl = sample(3, 2)
  xk = df[Xkl[1]]
  xl = df[Xkl[2]]
  if (which.example == 1)
  {
    bj0 = runif(n, -1, 1)
    bj1 = runif(n, -1, 1)
    bj2 = runif(n, -1, 1)
    epsilon = rnorm(n, 0, 2)
    xj = bj0 + bj1*xk + bj2*xl + epsilon   
  }
  else if (which.example == 2)
  {
    bj0 = runif(n, -1, 1)
    bj1 = runif(n, -1, 1)
    bj2 = runif(n, -1, 1)
    bj3 = runif(n, -1, 1)
    bj4 = runif(n, -1, 1)
    epsilon = rnorm(n, 0, 5)
    xj = bj0 + bj1*xk + bj2*xl + bj3*xk^2 + bj4*xl^2 + epsilon   
  }
  else if (which.example == 3)
  {
    bj1 = runif(n, -1, 1)
    bj2 = runif(n, -1, 1)
    epsilon = rnorm(n, 0, 1)
    xj = bj1*xk + bj2*xl + abs(xk)*epsilon   
  }
  return(xj)
}
