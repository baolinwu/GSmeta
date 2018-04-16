##########################
##
##    RE2C 
##    New RE2 more powerful than RE2
##   
#####  Buhm Han & Cue Lee

### 2/18/15: v1: launch
### 2/29/15: v2: Speed up using eigen decomp. Separated calculating low/high thresholds. (Assume we always input these)
### 2/20/15: v3: I try to make p-value continuous. (not stacked at 1) 
### 9/4/15: v4: with null test function.. I will try empirical distribution one more time.
### 9/8/15: v5: fast implementation of using empirical table.
### 8/20/16: v6: set Shet.high at Stat1
### 11/18/16: v7: account correlation of shared controls
### 12/9/16: v8: update the default tlow values under GWASthreshold of 5*10^-8
### 1/11/17: v9: change the p-value estimation to sum_{x in SFE bin} ( Prob(Shet > max(S-x, Thres(x)))*Prob(x) )
### 3/20/17: v10: update to check dependency of mvtnorm package before run RE2C

### RE2C calculate the parameter tlow; To estimate tlow, we need 1. Han and Eskin pvalue table, 2.number of studies, 3.target GWAS threshold (default 5*10^-8)
## I recommend to calculate tlow for every available study number which is ranged from 2 to 50 (later 100) and run RE2h
## provided default tlow with conventional GWASthreshold of  5*10^-8

## RE2C uses R. package names <mtvnorm>
## suppressWarnings(pkgTest("mvtnorm"))

#' Improved random effects meta-analysis (RE2C) accounting for correlations and heterogeneous effects for GWAS
#'
#' RE model accounts for fixed effects (FE) and heterogeneous effects in meta-analyzing potentially dependent GWAS. 
#' This implementation is based on the original software posted at http://software.buhmhan.com/RE2C with slight modification.
#' @param beta   effect size estimates from n studies
#' @param stders  standard errors for the effct size estimates
#' @param sigma   the correlation matrix of effect size estimates
#' @return 
#' \describe{
#'  \item{RE2p}{ test p-value for overall effects (FE + RE) }
#'  \item{RE2Cp}{ test p-value for heterogenenous effect conditional on the FE }
#'  \item{LRT}{ LRT statistics for the FE and  RE in additiona to the FE  }
#' }
#' @export
#' @references
#' Lin,D.Y. and Sullivan,P.F. (2009) Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects. Am J Hum Genet 85, 862–872.
#' 
#' Han,B. and Eskin,E. (2011) Random-Effects Model Aimed at Discovering Associations in Meta-Analysis of Genome-wide Association Studies. The American Journal of Human Genetics 88, 586–598.
#'
#' Lee,C.H., Eskin,E., Han,B. (2017) Increasing the power of meta-analysis of genome-wide association studies to detect heterogeneous effects. Bioinformatics 33, i379–i388. 
## RE2h main
RE2C <- function(beta, stders, sigma) {
  ## load RE2C tabulated tables
  data(RE2Cdata, package="GSmeta")
  nstudy = length(beta)
  NR = test_n(nstudy)
  ## correction parameters
  RE2Cor = RE2Cor.list[[min(nstudy-1,6)]]
  tau2prob_cor = tau2prob_corlist[[min(nstudy-1,6)]]
  ## run correction.function(cor)
  correction.list = correction.function(sigma, RE2Cor=RE2Cor,tau2prob_cor=tau2prob_cor,raw.tau2.zero.prob=raw.tau2.zero.prob)
  RE2C_corr <- correction.list[[1]]
  tau2.zero.prob <- correction.list[[2]]
  ##
  Vb = t(sigma*stders)*stders
  out = RE2Cfc(beta,stders, sigma=Vb, stat2_cdf=as.numeric(RE2C_table[nstudy-1,])*RE2C_corr, FEt.lows=as.numeric(FEtlow_table[nstudy-1,]),
               NR=NR, FEs=FEs,FEprobs=FEprobs,tau2.zero.prob=tau2.zero.prob)
  RE2Cout<-out[[1]]
  stat1 <- RE2Cout[1]
  stat2 <- RE2Cout[2]
  RE2Cp_joint <- RE2Cout[3]
  RE2Cp <- RE2Cout[4]
  ## not calibrated test: commented out
  ## Qout<-out[[2]]
  ## Q <- Qout[1]
  ## Qp <- pchisq(q=Q,df=(nstudy-1),ncp=0,lower.tail=F)
  ## Isq <- floor(max(100*(Q-(nstudy-1))/Q,0))
  ## Q = round(Q,2)
  ## write(c(rsid, nstudy, LSs, LSse, LSp, Isq, Q, Qp, stat1, stat2, RE2Cp_cond, RE2Cp),file=outputf,ncolumns = 12,append=T)
  return( list(RE2p=RE2Cp_joint, RE2Cp = RE2Cp, LRT=c(f=stat1, h=stat2) ) )
}

## Internal function to be called to run RE2C
RE2Cfc <- function(beta, stders, sigma, stat2_cdf = FALSE, FEt.lows = FALSE,
                   NR=NR, FEs=FEs,FEprobs=FEprobs,tau2.zero.prob=tau2.zero.prob){
  ## n study numbers
  ## stopifnot (length(beta) == length(stders))
  n=length(beta)
  vars <- stders ** 2
  ws <- 1/vars 
  sigmainv <- solve(sigma)
  beta_rv <- matrix(beta,nrow=1,ncol=n)
  ones <- matrix(rep(1,n),nrow=1,ncol=n)
  sumW <- sum(ws)
  sumW2 <- sum(ws ** 2)
  meanBeta <- (ones %*% sigmainv %*% t(beta_rv)) / (ones %*% sigmainv %*% t(ones))
  Q <- (beta_rv - rep(meanBeta,n))** 2 %*% ws
  meanW <- mean(ws)
  Sw2 <- 1/(n-1) * (sumW2 - n*meanW*meanW)
  U = (n-1)*(meanW - Sw2/(n*meanW))
  tau2 <- max( (Q-(n-1))/U, 0 )
  
  ##-----------------------------------------------
  ## Eigen-decomposition optimization (EMMA- style)
  ##-----------------------------------------------
  K <- sigma
  eig.L <- eigen(K,symmetric=TRUE)
  L.values <- eig.L$values
  L.vectors <- eig.L$vectors
  S <- diag(n)-matrix(1/n,n,n)
  eig.R <- eigen (S%*%K%*%S,symmetric=TRUE)
  R.values <- eig.R$values[1:(n-1)]
  R.vectors <- eig.R$vectors[,1:(n-1)]
  etas <- crossprod(R.vectors,t(beta_rv))
  etasqs <- etas^2
  xis <- L.values
  lambdas <- R.values
  
  mle.tau2 <- NR(tau2,n,xis,etasqs,lambdas)
  Hinv <- solve(sigma+mle.tau2*diag(n))
  mle.beta <- (ones %*% Hinv %*% t(beta_rv)) / (ones %*% Hinv %*% t(ones))
  mle.ll <- -LL.fun(mle.tau2,n,xis,etasqs,lambdas)
  
  null.ll = -ll.fun(c(0,0),beta_rv,sigma)
  lrt.mle <- -2*(null.ll-mle.ll)
  
  stat1 = -2*(null.ll+ll.fun(c(meanBeta,0),beta_rv,sigma))
  stat2 = -2*(-ll.fun(c(meanBeta,0),beta_rv,sigma)+ll.fun(c(mle.beta,mle.tau2),beta_rv,sigma))
  
  p.RE2_cond <- tau2.zero.prob[n-1]*pchisq(lrt.mle,1,lower.tail=F)+(1-tau2.zero.prob[n-1])*pchisq(lrt.mle,2,lower.tail=F)
  
  ##----------------------------------------------
  ## Given two stats, let's calculate p-value
  ##----------------------------------------------
  stat=stat1+stat2
  approx_stat1=min(stat1,49.974)
  if(c(FEt.lows)[floor(approx_stat1*20)+1]>stat2) {
    ## return( list( c(stat1,stat2,p.RE2_cond,1), c(Q) ) )
    return( list( c(stat1,stat2,p.RE2_cond,1) ) )
  }
  if(stat>50){
    p.RE2C <- RE2C_ext(stat,stat2_cdf=stat2_cdf,tmax=FEt.lows[1000],n,
                       FEs_ext=FEs_ext,tau2.zero.prob=tau2.zero.prob)
    ## return( list( c(stat1,stat2,p.RE2_cond,p.RE2C), c(Q) ) )
    return( list( c(stat1,stat2,p.RE2_cond,p.RE2C) ) )
  }
  p.RE2C = ((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(FEt.lows,(rep(stat,1000)-FEs)),1,max)))+1])%*%FEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*FEt.lows[1000]+1)]*pchisq(50, 1, lower.tail=F)
  ## return( list( c(stat1,stat2,p.RE2_cond,p.RE2C), c(Q) ) )
  return( list( c(stat1,stat2,p.RE2_cond,p.RE2C) ) )
}

RE2C_ext <- function(stat,stat2_cdf,tmax,n,
                     FEs_ext=FEs_ext,tau2.zero.prob=tau2.zero.prob){
  extra<-stat-50
  modFEs<-seq(extra+0.25,stat,0.25)
  modFEprobs<-cal_FEprobs(modFEs)
  p.RE2C=(((1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*(apply(cbind(rep(tmax,200),(rep(50,200)-FEs_ext)),1,max)))+1])%*%modFEprobs + (1-tau2.zero.prob[n-1])*stat2_cdf[floor(20*tmax+1)]*pchisq(stat, 1, lower.tail=F))
  return(p.RE2C)
}

  
## RE2C uses tabulated matrices to correct stat2_cdf and FEt.lows
## correction factors
correction.function <- function(correlationMatrix, 
                                RE2Cor=RE2Cor,tau2prob_cor=tau2prob_cor,raw.tau2.zero.prob=raw.tau2.zero.prob){
  mean_cor <- mean(correlationMatrix[lower.tri(correlationMatrix)])
  int_corr <- floor(mean_cor*10)
  float_corr <- (mean_cor*10) - int_corr
  corr_coeff <- int_corr+1
  aset <- c(corr_coeff,corr_coeff+1)
  ## diff(matrix(x)) If x is a matrix then the difference operations are carried out on each column separately.
  RE2C_corr<- RE2Cor[corr_coeff,] + diff(RE2Cor[aset,])*float_corr
  tau2.zero.prob_corr <- tau2prob_cor[corr_coeff]+ diff(tau2prob_cor[aset])*float_corr
  tau2.zero.prob <- raw.tau2.zero.prob*tau2.zero.prob_corr
  return(list(RE2C_corr,tau2.zero.prob))
}
ll.fun <- function(x,beta,sigma){
  -sum(mvtnorm::dmvnorm(beta, mean=rep(x[1],length(beta)),sigma=sigma+diag(x[2],length(beta)),log=T))
}
LL.fun <- function(x,n,xis,etasqs,lambdas) {
  0.5*(n*log(2*pi)+sum(log(xis+x))+sum(etasqs/(lambdas+x)))
}
test_n <- function(nstudy){ 
  if (nstudy == 2| nstudy == 3){
    NR_R <- function(tau2,n,xis,etasqs,lambdas){
      LL.fun <- function(x) {
        0.5*(n*log(2*pi)+sum(log(xis+rep(x,n)))+sum(etasqs/(lambdas+rep(x,n-1))))
      }
      optim.rst <- optim(par=c(tau2), fn=LL.fun, method = "L-BFGS-B", lower = 0, upper=Inf)
      return(optim.rst$par[1])
    }
  } else {
    NR_R <- function(x,n,xis,etasqs,lambdas){
      init <- x 
      i = 1
      while(abs(converge <- ( -0.5 * sum( 1/(xis + rep(x,n)) ) + 0.5 * sum (  etasqs / (lambdas + rep(x,n-1))^2 ) )) > 10^-8 )   {
        newx <- x - (converge / ( 0.5 * sum(  1 / (xis + rep(x,n))^2) - sum ( etasqs / (lambdas + rep(x,n-1))^3 ))  )
        x <- newx
        i = i+1
        if(i%%100==0) {
          ## warning("NR failed to converge. Use optim to estimate tau^2")
          LL.fun <- function(x) {
            0.5*(n*log(2*pi)+sum(log(xis+rep(x,n)))+sum(etasqs/(lambdas+rep(x,n-1))))
          }
          optim.rst <- optim(par=c(x), fn=LL.fun, method = "L-BFGS-B", lower = 0, upper=Inf)
          return(optim.rst$par[1])
          break
        }
      }
      if(abs(init-x) > 10^4) {
        return(init)
      } else if (x < 0) {
        return(0)
      } else {
        return(x)
      }
    }
  }
  NR <- compiler::cmpfun(NR_R,options = list(optimize=2,suppressAll=T))
  remove(NR_R)
  return(NR)
}

cal_FEprobs <- function(FEs){
  pchisqs<-pchisq(FEs,df=1,ncp=0,lower.tail=F)
  inst=1
  FEprobs=NULL
  for (i in 1:length(FEs)){
    FEprobs <- c(FEprobs,inst-pchisqs[i])
    inst=pchisqs[i]
  }
  return(FEprobs)
}
