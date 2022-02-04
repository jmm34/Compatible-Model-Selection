---
title: "Confidence interval on the cross entropy"
author: "Jean-Michel Marin"
institute: "CNRS -- IMAG (Montpellier, France)"
date: February 2022
output: binb::metropolis
classoption: "aspectratio=169,12pt"
header-includes:  \usepackage{xcolor}
---

## Goal

-   True distribution $F$ (density $f$)

-   We observed an iid sample $(y_1,\ldots,y_n)$ from $F$

-   We considered the parametric model $Q_\theta$ (density $q(y|\theta)$)

-   The cross-entropy is $$
    \beta=\max_\theta\mathbb{E}\left[\log(q(y|\theta)\right]
    =\max_\theta\log(q(y|\theta))f(y)\text{d}y
    $$

\textcolor{red}{\large \bf Construct a confidence interval on $\beta$}

## Estimate of $\beta$

- Clearly, $\beta=\gamma(F)$

- $\hat F_n$ the empirical distribution associated to $(y_1,\ldots,y_n)$

- $\hat\theta_n$ the MLE of $\theta$

- An estimate of $\beta$ is 
$$
\hat\beta_n=\gamma(\hat F_n)=\frac{1}{n}\sum_{i=1}^n\log(q(y_i|\hat\theta_n))
$$

\textcolor{red}{\large \bf We propose to bootstrap the statistics $S=\hat\beta_n-\beta$}

## A gaussian toy example $F=\mathcal{N}(0,1)$ ; $Q_\theta=\mathcal{N}(\theta,1)$

In such a case 

- $\beta=-1/2\log(2\pi)-1/2$

- $\hat\theta_n=\bar y_n$ ; $\hat\beta_n=-1/2\log(2\pi)-1/(2n)\sum_{i=1}^n(y_i-\bar y_n)^2$
  
- The bootstrap version of $S=\hat\beta_n-\beta=\gamma(\hat F_n)-\gamma(F)$ is
$$
S^*=1/(2n)\sum_{i=1}^n(y_i-\bar y_n)^2-1/(2n)\sum_{i=1}^n(y_i^*-\bar y^*_n)^2
=1/2s_n^2-1/2s_n^*
$$
where $(y_1^*,\ldots,y_n^*)$ is a bootstrap sample from $(y_1,\ldots,y_n)$

## Numerical experiments for the gaussian toy example $F=\mathcal{N}(0,1)$ ; $Q_\theta=\mathcal{N}(\theta,1)$

\scriptsize 

```R
beta <- -1/2*log(2*pi)-1/2
N <- 1000 ; n <- 500 ; B <- 500 ; alpha <- 0.1 
statboot <- rep(0,B) ; cover <- 0
for (i in 1:N) {
  x <- rnorm(n) ; sn2 <- var(x)*(n-1)/n 
  for (j in 1:B) {
    xstar <- sample(x,n,rep=TRUE) ; sn2star <- var(xstar)*(n-1)/n
    statboot[j] <- sn2/2-sn2star/2 }
  interva <- -1/2*log(2*pi)-sn2/2-quantile(statboot,c(1-alpha/2,alpha/2),
  names=FALSE)
  cover <- cover + (beta >= interva[1] & beta <= interva[2]) }
cover/N
```
\textcolor{red}{\large \bf we get $0.903$}

## A gaussian example with unknown variance ; $F=\mathcal{N}(0,1)$ ; $Q_{(\mu,\sigma^2)}=\mathcal{N}(\mu,\sigma^2)$

In such a case 

- $\beta=-1/2\log(2\pi)-1/2$

- $\hat\theta_n=\left(\bar y_n,s_n^2\right)$ ; $\hat\beta_n=-1/2(\log(2\pi)+1)-1/2\log(s_n^2)$
  
- The bootstrap version of $S=\hat\beta_n-\beta=\gamma(\hat F_n)-\gamma(F)$ is
$$
S^*=1/2\log(s_n)-1/2\log(s_n^*)
$$

## Numerical experiments for the gaussian toy example $F=\mathcal{N}(0,1)$ ; $Q_{(\mu,\sigma^2)}=\mathcal{N}(\mu,\sigma^2)$

\scriptsize 

```R
beta <- -1/2*log(2*pi)-1/2
N <- 1000 ; n <- 500 ; B <- 500 ; alpha <- 0.1 
statboot <- rep(0,B) ; cover <- 0
for (i in 1:N) {
  x <- rnorm(n) ; sn2 <- var(x)*(n-1)/n 
  for (j in 1:B) {
    xstar <- sample(x,n,rep=TRUE) ; sn2star <- var(xstar)*(n-1)/n
    statboot[j] <- log(sn2)/2-log(sn2star)/2 }
  interva <- -1/2*log(2*pi)-sn2/2-quantile(statboot,c(1-alpha/2,alpha/2),
  names=FALSE)
  cover <- cover + (beta >= interva[1] & beta <= interva[2]) }
cover/N
```
\textcolor{red}{\large \bf we get $0.896$}
