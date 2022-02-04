# Gaussian example

N <- 5000
n <- 500
B <- 100
statboot <- rep(0,B)
alpha <- -n/2*log(2*pi)-n/2
alpha


## Known variance

N <- 5000
n <- 400
B <- 500
statboot <- rep(0,B)
alpha <- -n/2*log(2*pi)-n/2
alpha/n
moy <- ampli <- cover <- 0

tp <- txtProgressBar(min = 1, max = N, style = 3, char = "*")
for (i in 1:N)
{
  x <- rnorm(n)
  xbar <- mean(x)
  for (j in 1:B)
  {
    xstar <- sample(x,n,rep=TRUE)
    statboot[j] <- -1/2*sum((xstar-mean(xstar))^2)+1/2*sum((x-xbar)^2)
  }
  interva <- -n/2*log(2*pi)-1/2*sum((x-xbar)^2)-
    quantile(statboot,c(0.95,0.05),names=FALSE)
  moy <- moy+sum(interva)/2
  ampli <- ampli+diff(interva)
  cover <- cover + (alpha >= interva[1] & alpha <= interva[2])
  setTxtProgressBar(tp, i)
}
moy/(n*N)
ampli/(n*N)
cover/N

# Unknwon variance

N <- 5000
n <- 400
B <- 500
statboot <- rep(0,B)
alpha <- -1/2*log(2*pi)-1/2

moy <- ampli <- cover <- 0

tp <- txtProgressBar(min = 1, max = N, style = 3, char = "*")
for (i in 1:N)
{
  x <- rnorm(n)
  sn2 <- var(x)*(n-1)/n
  for (j in 1:B)
  {
    xstar <- sample(x,n,rep=TRUE)
    sn2star <- var(xstar)*(n-1)/n
    statboot[j] <- -1/2*log(sn2star)+1/2*log(sn2)
  }
  interva <- -1/2*(log(2*pi)+1)-1/2*log(sn2)-
    quantile(statboot,c(0.95,0.05),names=FALSE)
  moy <- moy+sum(interva)/2
  ampli <- ampli+diff(interva)
  cover <- cover + (alpha >= interva[1] & alpha <= interva[2])
  setTxtProgressBar(tp, i)
}
moy/N
ampli/N
cover/N

## Alternative 1

tp <- txtProgressBar(min = 1, max = N, style = 3, char = "*")
moy <- ampli <- cover <- 0
for (i in 1:N)
{
  x <- rnorm(n)
  xbar <- mean(x)
  for (j in 1:B)
  {
    xstar <- sample(x,n,rep=TRUE)
    statboot[j] <- -n/2*mean(xstar)^2/2+n/2*xbar^2
  }
  interva <- -n/2*(log(2*pi)+1+xbar^2)-
    quantile(statboot,c(0.95,0.05),names=FALSE)
  moy <- moy+sum(interva)/2
  ampli <- ampli+diff(interva)
  cover <- cover + (alpha >= interva[1] & alpha <= interva[2])
  setTxtProgressBar(tp, i)
}
moy/N
ampli/N
cover/N

## Alternative 2

tp <- txtProgressBar(min = 1, max = N, style = 3, char = "*")
moy <- ampli <- cover <- 0
for (i in 1:N)
{
  x <- rnorm(n)
  xbar <- mean(x)
  for (j in 1:B)
  {
    xstar <- sample(x,n,rep=TRUE)
    statboot[j] <- -n/2*log(2*pi)-1/2*sum((xstar-mean(xstar))^2)-1/2
  }
  interva <- quantile(statboot,c(0.05,0.95),names=FALSE)
  moy <- moy+sum(interva)/2
  ampli <- ampli+diff(interva)
  cover <- cover + (alpha >= interva[1] & alpha <= interva[2])
  setTxtProgressBar(tp, i)
}
moy/N
ampli/N
cover/N

## Monte Carlo

tp <- txtProgressBar(min = 1, max = N, style = 3, char = "*")
stat <- rep(0,N)
for (i in 1:N)
{
  x <- rnorm(n)
  stat[i] <- -n/2*log(2*pi)-1/2*sum((x-mean(x))^2)
  setTxtProgressBar(tp, i)
}

## Regression model

n <- 1000
sdres <- 1
X <- matrix(rnorm(n*2),n,2)
niter <- 100000
bic1 <- bic2 <- rep(0,niter)
Xm <- solve(t(X)%*%X)
pb <- txtProgressBar(min = 1, max = niter, style = 3)

for (i in 1:niter)
{
  y <- 0.1*X[,1]+0.1*X[,2]+rnorm(n,0,sdres)
  ml <- Xm%*%t(X)%*%y
  bic1[i] <- n*log(t(y-X%*%ml)%*%(y-X%*%ml)/n)+2*log(n)
  ml <- sum(X[,2]*y)/sum(X[,2]^2)  
  bic2[i] <- n*log(sum((y-X[,2]*ml)^2)/n)+log(n)
  setTxtProgressBar(pb, i)
}
res <- sum(bic1>bic2)/niter
res

niter <- 1000
best <- rep(2,niter)
B <- 1000
compteur1 <- rep(0,niter)
compteur2 <- rep(0,niter)
pb <- txtProgressBar(min = 1, max = niter, style = 3)
for (i in 1:niter)
{
  y <- 0.1*X[,1]+0.1*X[,2]+rnorm(n,0,sdres)
  ml <- Xm%*%t(X)%*%y
  bic1 <- n*log(t(y-X%*%ml)%*%(y-X%*%ml)/n)+2*log(n)
  ml <- sum(X[,2]*y)/sum(X[,2]^2)  
  bic2 <- n*log(sum((y-X[,2]*ml)^2)/n)+log(n)
  bool <- FALSE
  if (bic1<bic2) { bool <- TRUE ; best[i] <- 1 }
  for (l in 1:B)
  {
    indi <- sample(1:n,n,rep=TRUE)
    yB <- y[indi]
    XB <- X[indi,]
    XmB <- solve(t(XB)%*%XB)
    ml <- XmB%*%t(XB)%*%yB
    bic.true <- n*log(t(yB-XB%*%ml)%*%(yB-XB%*%ml)/n)+
      2*log(n)
    ml <- sum(XB[,2]*yB)/sum(XB[,2]^2)  
    bic.sub <- n*log(sum((yB-XB[,2]*ml)^2)/n)+log(n)
    if (bool==TRUE) compteur1[i] <- (bic.sub<bic.true)+compteur1[i] else
      compteur1[i] <- (bic.sub>bic.true)+compteur1[i]
    compteur2[i] <- (bic.sub<bic.true)+compteur2[i]
    setTxtProgressBar(pb, i)
  }
}
compteur1 <- compteur1/B
compteur2 <- compteur2/B
median(compteur1)
median(compteur2)

sum(compteur1>0.1)
sum(best[compteur1<0.1]==2)

sum(compteur2>0.1)
sum(best[compteur2<0.1]==2)

# boxplot(compteur1/B,res)
# boxplot(compteur1[best==1]/B,res)
# boxplot(compteur1[best==2]/B,res)
# boxplot(compteur2/B,res)


sum(diff.bic<=0)
boxplot(diff.bic)
median(diff.bic)

niter <- 100
B <- 1000
diff.bicB <- matrix(0,niter,B)

for (i in 1:niter)
{
  y <- 0.1*X[,1]+0.1*X[,2]+rnorm(n,0,sdres)
  for (l in 1:B)
  {
    indi <- sample(1:n,n,rep=TRUE)
    yB <- y[indi]
    XB <- X[indi,]
    XmB <- solve(t(XB)%*%XB)
    ml <- XmB%*%t(XB)%*%yB
    bic.true <- n*log(t(yB-XB%*%ml)%*%(yB-XB%*%ml)/n)+
      2*log(n)
    ml <- sum(XB[,2]*yB)/sum(XB[,2]^2)  
    bic.sub <- n*log(sum((yB-XB[,2]*ml)^2)/n)+log(n)
    diff.bicB[i,l] <- bic.true-bic.sub
  }
}

resB <- apply(diff.bicB,1,median)
boxplot(resB,median(diff.bic))
median(resB)
median(diff.bic)

exampleB <- rep(0,B)
y <- 0.1*X[,1]+0.1*X[,2]+rnorm(n,0,sdres)
for (l in 1:B)
{
  indi <- sample(1:n,n,rep=TRUE)
  yB <- y[indi]
  XB <- X[indi,]
  XmB <- solve(t(XB)%*%XB)
  ml <- XmB%*%t(XB)%*%yB
  bic.true <- n*log(t(yB-XB%*%ml)%*%(yB-XB%*%ml)/n)+
    2*log(n)
  ml <- sum(XB[,2]*yB)/sum(XB[,2]^2)  
  bic.sub <- n*log(sum((yB-XB[,2]*ml)^2)/n)+log(n)
  exampleB[l] <- bic.true-bic.sub
}


ml <- Xm%*%t(X)%*%y
bic.true <- n*log(t(y-X%*%ml)%*%(y-X%*%ml)/n)+2*log(n)
ml <- sum(X[,2]*y)/sum(X[,2]^2)  
bic.sub <- n*log(sum((y-X[,2]*ml)^2)/n)+log(n)
diff.example <- bic.true-bic.sub

f1 <- ecdf(exampleB)
f(diff.example)

f2 <- ecdf(diff.bic)
f2(diff.example)


