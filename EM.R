################# Missing Packages ##################
if(!require(ggplot2)){
  install.packages("ggplot2") 
  library(ggplot2)
}

if(!require(mvtnorm)){
  install.packages("mvtnorm")
  library(mvtnorm)
}
#####################################################
set.seed(42)

# Finds AIC values
aic <- function(data, class, theta)
{
  loss <- 0
  
  for(c in 1:class)
  {
    mu.temp  <- theta[(cl+1+2*(c-1)):(cl+2+2*(c-1))]
    sig.temp <- matrix(theta[(3*cl+1+4*(c-1)):(3*cl+4+4*(c-1))],byrow=TRUE,ncol=2,nrow=2) 
    loss <- loss +  theta[c]*dmvnorm(data, mu.temp, sig.temp)
  }
  
    fin <- - 2*sum(log(loss)) + 2*(6*class - 1) # No. of params  = 6*C -  1
    return(fin)
}

loglike <- function(data, class, theta)
{
  loss <- 0
  
    for(c in 1:class)
  {
    mu.temp  <- theta[(cl+1+2*(c-1)):(cl+2+2*(c-1))]
    sig.temp <- matrix(theta[(3*cl+1+4*(c-1)):(3*cl+4+4*(c-1))],byrow=TRUE,ncol=2,nrow=2)
    loss <- loss +  theta[c]*dmvnorm(data, mu.temp, sig.temp)
  }
  
    fin <-  - sum(log(loss))
    return(fin)
}


GMMforcl <- function(dat1,cl,flag=0)
{
  n = nrow(dat1)
  
  #k-means for initialisation
  dat1.kmeans <- kmeans(dat1,cl,nstart=10)
  dat1.kmeans.cluster <- dat1.kmeans$cluster
  dat1.df <- data.frame(x = dat1, cluster = dat1.kmeans.cluster)
  
  print(ggplot(data = dat1.df,aes(y = x.Y, x = x.X, color = factor(cluster))) +
          geom_point() + xlab("X") + ylab("Y") + scale_color_discrete(name = "Cluster") +
          ggtitle("K-means Clustering"))
  
  
  mu.init = as.vector(t(dat1.kmeans$centers))
  sig2.init <- rep(c(var(dat1[,1]), cov(dat1[,1],dat1[,2]), cov(dat1[,2], dat1[,1]), var(dat1[,2])),cl)
  
  dat1.kmeans$pi <- dat1.kmeans$size/sum(dat1.kmeans$size)
  pi.init <- dat1.kmeans$pi
  
  
  tol <- 1e-6
  itr <- 0
  diff <- 1000
  theta = c(pi.init,mu.init,sig2.init)
  current <- theta
  
  Exp <- matrix(0, nrow = n, ncol = cl) #Expectation Matrix
  
  while(diff > tol)
  {
    itr <- itr + 1
    
    #E step
    for(c in 1:cl)
    {
      mu.temp  <- current[(cl+1+2*(c-1)):(cl+2+2*(c-1))] #mu,sigma for class c
      sig.temp <- matrix(current[(3*cl+1+4*(c-1)):(3*cl+4+4*(c-1))],byrow=TRUE,ncol=2,nrow=2) 
      Exp[,c] <- current[c]*dmvnorm(dat1, mean=mu.temp,sigma=sig.temp)
    }
    
    Exp <- Exp/(rowSums(Exp))
    
    # M-step 
    
    for(c in 1 :cl)
    {
      theta[c] = mean(Exp[,c])
      theta[(cl+1+2*(c-1)):(cl+2+2*(c-1))] = colSums(Exp[,c]*dat1)/sum(Exp[,c])
      theta[(3*cl+(4*(c-1)+1)):(3*cl+ 4*(c-1) +4)] = cov.wt(dat1, Exp[,c])$cov
    }
    
    diff <- max(abs(theta-current)) 
    current <- theta
  }
  
  print(paste("No. of iterations = ", itr))
  print(paste("No. of classes = ", cl))
  
  if(flag==1)
  ExpGlobal <<- Exp #for final clustering
  
  return(theta)
}


# 5-fold Cross-validation
dat <- read.table("http://home.iitk.ac.in/~dootika/assets/course/MixG_data/170187.txt",
                  header = F, col.names = c("X", "Y"))
n = nrow(dat)

permutation <- sample(1:n, replace = FALSE)
K <- 5

test.index <- split(permutation, rep(1:K, length = n, each = n/K)) #where to test

classes <- 2:6
CV.Loglike <- numeric(length = length(classes)+1)
CV.AIC <- numeric(length = length(classes)+1)

for(cl in 2:6)
{
  temp3 <- 0
  temp4 <- 0
 
  for(k in 1:K)
{
  dat.train <- dat[-test.index[[k]], ] #split
  dat.test <- dat[test.index[[k]], ]
  dat.fit <- GMMforcl(dat.train,cl)
  temp3 <- temp3 + loglike(data = dat.test, class = cl, theta = dat.fit)
}
  
dat.fit <- GMMforcl(dat,cl)
temp4 <- aic(data = dat, class = cl, theta=dat.fit)

CV.Loglike[cl] <- temp3/n 
CV.AIC[cl] <- temp4/n
}

CV.Loglike <- CV.Loglike[2:6]
CV.AIC <- CV.AIC[2:6]
Class.Final <- which.min(CV.Loglike) + 1

###### Final Fit #######
ExpGlobal <- matrix(0,nrow=n,ncol=Class.Final) #final cluster expectation matrix
dat.fit <- GMMforcl(dat,Class.Final,flag=1)


mu.est <- matrix(0,nrow=Class.Final,ncol=2) 
sig.est <- list()
cl <- Class.Final
for(i in 1:Class.Final)
{
  sig.est[[i]] = matrix(dat.fit[(3*cl+1+4*(i-1)):(3*cl+4+4*(i-1))],byrow=TRUE,ncol=2,nrow=2)
  mu.est[i,] = dat.fit[(cl+1+2*(i-1)):(cl+2+2*(i-1))]
}
mix.est <- dat.fit[1:Class.Final]

####### Clustering ########

cluster <- numeric(length=n)
for(i in 1:n)
{
  cluster[i] <- which.max(ExpGlobal[i,])
}
dat$cluster <- cluster

print(ggplot(data = dat,aes(y = Y, x = X, color = factor(cluster))) +
        geom_point() + xlab("X") + ylab("Y") + scale_color_discrete(name = "Cluster") +
        ggtitle("Final Clusters"))
