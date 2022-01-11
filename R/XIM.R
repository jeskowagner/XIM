
# This code is a direct copy from Zhexiao Lin and Fang Han's code found
# here: https://sites.stat.washington.edu/people/fanghan/XIMCOR.R
# It is the implementation of their work found here: https://arxiv.org/abs/2108.06828v1
#
# Packaging by Jesko Wagner

#' @export
#' @title XIMcalculate
#'
#' @description This function computes the revised Chatterjee’s rank correlation `$\xi_{n,M}$`.
#' @param xvec Vector
#' @param yvec Vector
#' @param M the number of right nearest neighbors.
#' @return the value of `$\xi_{n,M}$`.
#' 
#' @author Zhexiao Lin and Fang Han
#' @export
XIMcalculate<-function(xvec, yvec, M){
  n <- length(xvec)
  xrank <- rank(xvec, ties.method = "random")
  yrank <- rank(yvec, ties.method = "random")
  ord <- order(xrank)
  yrank <- yrank[ord]
  coef.sum <- 0
  for (m in 1:M){
    coef.sum.temp <- sum(pmin(yrank[1:(n-m)], yrank[(m+1):n]))
    coef.sum.temp <- coef.sum.temp + sum(yrank[(n-m+1):n])
    coef.sum <- coef.sum + coef.sum.temp
  }
  coef.value <- -2+6*coef.sum/((n+1)*(n*M+M*(M+1)/4))
  return(coef.value)
}

#' @export
#' @title XIMstat
#'
#' @description This function computes the test statistics `$\xi_{n,M}^{\pm}$`.
#' @param xvec Vector
#' @param yvec Vector
#' @param M the number of right nearest neighbors.
#' @return the value of `$\xi_{n,M}^{\pm}$`.
#' 
#' @author Zhexiao Lin and Fang Han
#' @export
XIMstat<-function(xvec, yvec, M){
  n <- length(xvec)
  xrank <- rank(xvec, ties.method = "random")
  yrank <- rank(yvec, ties.method = "random")
  ord <- order(xrank)
  yrank <- yrank[ord]
  yrank.m <- n+1-yrank 
  coef.sum <- 0
  coef.sum.m <- 0
  for (m in 1:M){
    coef.sum.temp <- sum(pmin(yrank[1:(n-m)], yrank[(m+1):n]))
    coef.sum.temp <- coef.sum.temp + sum(yrank[(n-m+1):n])
    coef.sum <- coef.sum + coef.sum.temp
    coef.sum.temp.m <- sum(pmin(yrank.m[1:(n-m)], yrank.m[(m+1):n]))
    coef.sum.temp.m <- coef.sum.temp.m + sum(yrank.m[(n-m+1):n])
    coef.sum.m <- coef.sum.m + coef.sum.temp.m
  }
  coef.stat <- -2+6*max(coef.sum,coef.sum.m)/((n+1)*(n*M+M*(M+1)/4))
  return(coef.stat)
}

#' @title XIMsim
#'
#' @description This function computes the simulation statistics `$\xi_{n,M}^{\pm(b)}$`.
#' @param n number of samples
#' @param M the number of right nearest neighbors.
#' @param B the number of simulation.
#' @return the value of `$\xi_{n,M}^{\pm}$`.
#' 
#' @author Zhexiao Lin and Fang Han
#' @export
XIMsim<-function(n, M, B){
  coef.sim<-rep(0,B)
  for (b in 1:B){
    yrank <- sample(1:n,n)
    yrank.m <- n+1-yrank 
    coef.sum <- 0
    coef.sum.m <- 0
    for (m in 1:M){
      coef.sum.temp <- sum(pmin(yrank[1:(n-m)], yrank[(m+1):n]))
      coef.sum.temp <- coef.sum.temp + sum(yrank[(n-m+1):n])
      coef.sum <- coef.sum + coef.sum.temp
      coef.sum.temp.m <- sum(pmin(yrank.m[1:(n-m)], yrank.m[(m+1):n]))
      coef.sum.temp.m <- coef.sum.temp.m + sum(yrank.m[(n-m+1):n])
      coef.sum.m <- coef.sum.m + coef.sum.temp.m
    }
    coef.stat <- -2+6*max(coef.sum,coef.sum.m)/((n+1)*(n*M+M*(M+1)/4))
    coef.sim[b]<-coef.stat
  }
  return(coef.sim)
}

#' @title XIMtestT
#'
#' @description This function performs the simulation based test using test statistics and simulation statistics.
#' @param XIMstat is the value of `$\xi_{n,M}^{\pm}$`
#' @param XIMsim is the vector of `$\xi_{n,M}^{\pm(b)}$`
#' @param alpha is the signiﬁcance level.
#' @return a logical value, TRUE reject, FALSE accept.
#' 
#' @author Zhexiao Lin and Fang Han
#' @export
XIMtestT<-function(XIMstat, XIMsim, alpha = 0.05){
  B <- length(XIMsim)
  coef.test <- (1+sum(XIMstat<=XIMsim))/(1+B)
  return(coef.test <= alpha)
}

#' @title XIMtest
#'
#' @description This function performs the simulation based test directly.
#' @param xvec Vector
#' @param yvec Vector
#' @param M is the number of right nearest neighbors
#' @param B is the number of simulation
#' @param alpha is the signiﬁcance level.
#' @return a logical value, TRUE reject, FALSE accept.
#' 
#' @author Zhexiao Lin and Fang Han
#' @export
XIMtest<-function(xvec, yvec, M, B, alpha = 0.05){
  n <- length(xvec)
  coef.sim <- XIMsim(n,M,B)
  coef.stat <- XIMstat(xvec,yvec,M)
  return(XIMtestT(coef.stat, coef.sim, alpha))
}