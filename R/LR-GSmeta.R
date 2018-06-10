## log likelihood function for MA of independent GWAS
## UX ~ N(mu*U1, diag(D))
LLma <- function(UX, U1, D){
  mu = muf = sum(U1*UX/D)/sum(U1^2/D)
  Qf = sum(U1*UX/D)^2/sum(U1^2/D);   ## Qf = sum(UX^2/D) - sum((UX-U1*muf)^2/D)
  tau2m = max(UX^2-D)
  if(tau2m<=0){
    return( list(Qf=Qf,Qh=0, pars=c(muf=muf,mu=mu,tau2=0)) )
  }
  Lf0 = function(tau2){
    mu1 = sum(U1*UX/(D+tau2))/sum(U1^2/(D+tau2))
    sum( log(D+tau2) ) + sum( (UX-U1*mu1)^2/(D+tau2) )
  }
  ans = optimize(Lf0, c(0,tau2m) )
  tau2 = ans$min; Qh = Lf0(0)-ans$obj
  if(Qh<=0){
    tau2 = 0;  Qh=0
  } else{
    mu = sum(U1*UX/(D+tau2))/sum(U1^2/(D+tau2))
  }
  ##
  return( list(Qf=Qf,Qh=Qh, pars=c(muf=muf,mu=mu,tau2=tau2)) )
}
## Monter Carlo estimation of theta: fast calc with analytical profile scores!
MCtha <- function(U1,D, B=2500, FMC=TRUE){
  ## sim
  K = length(U1); Sx = sqrt(D)
  UX = matrix(rnorm(K*B), K,B)*Sx
  muf = colSums(U1*UX/D)/sum(U1^2/D)
  f0d = sum(1/D) - colSums( (UX-outer(U1,muf))^2/D^2 )
  if(FMC){
    tha = mean(f0d>=0)
    return(tha)
  }
  id0 = which(f0d>=0); B0 = length(id0)
  lrs = rep(1, B0)
  UX0 = UX[,id0,drop=FALSE]
  Lf = sum( log(D) ) + colSums( (UX0-outer(U1,muf[id0]))^2/D )
  ##
  for(i in 1:B0){
    Ui = UX[,id0[i]]
    tau2m = max(Ui^2-D)
    if(tau2m<=0){
      lrs[i] = 0
    } else{
      fc0 = function(tau2){
        mui = sum(U1*Ui/(D+tau2))/sum(U1^2/(D+tau2))
        sum( log(D+tau2) ) + sum( (Ui-U1*mui)^2/(D+tau2) )
      }
      L1 = optimize(fc0, c(0, tau2m) )$obj
      if(Lf[i]<=L1) lrs[i] = 0
    }
  }
  tha = sum(lrs<=0)/B
  return(tha)
}


#' A fixed effects (FE) meta-analysis (MA) of multiple GWAS with potentially dependent samples
#'
#' Under FE MA, we assume the same effect sizes across all studies. They can be efficiently estimated following the Lin-Sullivan approach.
#' @param betas   effect size estimates from K studies
#' @param stders  standard errors for the effct size estimates
#' @param sigma   the correlation matrix of effect size estimates. Default to NULL for independent studies.
#' @return
#' \describe{
#'  \item{betam}{ estimated mean effect size }
#'  \item{Vm}{ variance of betam }
#'  \item{p.value}{ p-value for testing the mean effect size }
#' }
#' @export
#' @references
#' Lin,D.Y. and Sullivan,P.F. (2009) Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects. Am J Hum Genet 85, 862–872.
#'
#' Han,B. and Eskin,E. (2011) Random-Effects Model Aimed at Discovering Associations in Meta-Analysis of Genome-wide Association Studies. The American Journal of Human Genetics 88, 586–598.
#'
#' Wu,B. and Zhao,H. (2018) Powerful random effects modeling for meta-analysis of genome-wide association studies.
FEma <- function(betas, stders, sigma=NULL){
  K = length(betas)
  if(is.null(sigma)) sigma = diag(K)
  Vi = solve( t(sigma*stders)*stders )
  Vm = 1/sum(Vi)
  betam = sum(Vi*betas)*Vm
  pval = pchisq(betam^2/Vm,1, lower=FALSE)
  return(list(betam=betam, Vm=Vm, p.value=pval) )
}

#' Random effects (RE) meta-analysis (MA) of multiple (potentially dependent) GWAS
#'
#' Under RE MA, we assume heterogeneous effect sizes across all studies. This is potentially relevant
#' when meta-analyzing cross-disease GWAS with overlapped samples.
#'
#' @param betas   effect size estimates from K studies
#' @param Vb      covariance matrix of effect size estimates
#' @param B       number of Monte Carlo samples to estimate the null distribution. Default to 2500. See Wu (2018) reference.
#' @return
#' \describe{
#'  \item{p.value}{ p-values for FE, RE, RE conditional on FE }
#'  \item{theta}{ proportion parameter in the chi-square mixture dist }
#'  \item{Q}{ LRT statistics for the FE and RE }
#'  \item{pars}{ estimated parameters (muf,mu,tau2) for FE mean, RE mean/variance parameters  }
#' }
#' @export
#' @references
#' Lin,D.Y. and Sullivan,P.F. (2009) Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects. Am J Hum Genet 85, 862–872.
#'
#' Han,B. and Eskin,E. (2011) Random-Effects Model Aimed at Discovering Associations in Meta-Analysis of Genome-wide Association Studies. The American Journal of Human Genetics 88, 586–598.
#'
#' Lee,C.H., Eskin,E., Han,B. (2017) Increasing the power of meta-analysis of genome-wide association studies to detect heterogeneous effects. Bioinformatics 33, i379–i388.
#'
#' Wu,B. and Zhao,H. (2018) Powerful random effects modeling for meta-analysis of genome-wide association studies.
REma <- function(betas, Vb, B=2500){
  es = eigen(Vb, sym=TRUE); U = es$vec
  UX = colSums(U*betas); U1 = colSums(U); D = es$val
  ans = LLma(UX,U1,D)
  tha = MCtha(U1,D, B)
  ##
  Qf = ans$Qf; Qv = Qf+ans$Qh
  pvalf = pchisq(Qf,1,lower=FALSE)
  pvalv = pPc2(Qv,tha)
  pvalc = 1;   if(pvalv<=pvalf) pvalc = pRE2C(pvalv, tha)
  p.value = c(pvalf,pvalv,pvalc);  names(p.value) = c('Pf','Pr','Pc')
  Q=c(Qf,Qv); names(Q) = c('Qf','Qr')
  return(list(p.value=p.value,theta=tha, Q=Q, pars=ans$pars) )
}


#' Random effects (RE) meta-analysis (MA) of multiple independent GWAS
#'
#' Under RE MA, we assume heterogeneous effect sizes across all studies. This is potentially relevant
#' when meta-analyzing cross-disease GWAS. This function is designed and optimized for independent GWAS.
#'
#' @param betas   effect size estimates from K studies
#' @param Vb      variances of individual effect size estimates
#' @param B       number of Monte Carlo samples to estimate the null distribution. Default to 2500. See Wu (2018) reference.
#' @return
#' \describe{
#'  \item{p.value}{ p-values for FE, RE, RE conditional on FE }
#'  \item{theta}{ proportion parameter in the chi-square mixture dist }
#'  \item{Q}{ LRT statistics for the FE and RE }
#'  \item{pars}{ estimated parameters (muf,mu,tau2) for FE mean, RE mean/variance parameters  }
#' }
#' @export
#' @references
#' Lin,D.Y. and Sullivan,P.F. (2009) Meta-Analysis of Genome-wide Association Studies with Overlapping Subjects. Am J Hum Genet 85, 862–872.
#'
#' Han,B. and Eskin,E. (2011) Random-Effects Model Aimed at Discovering Associations in Meta-Analysis of Genome-wide Association Studies. The American Journal of Human Genetics 88, 586–598.
#'
#' Lee,C.H., Eskin,E., Han,B. (2017) Increasing the power of meta-analysis of genome-wide association studies to detect heterogeneous effects. Bioinformatics 33, i379–i388.
#'
#' Wu,B. and Zhao,H. (2018) Powerful random effects modeling for meta-analysis of genome-wide association studies.
REmai <- function(betas, Vb, B=2500){
  UX = betas; U1 = rep(1, length(betas)); D = Vb
  ans = LLma(UX,U1,D)
  tha = MCtha(U1,D, B)
  ##
  Qf = ans$Qf; Qv = Qf+ans$Qh
  pvalf = pchisq(Qf,1,lower=FALSE)
  pvalv = pPc2(Qv,tha)
  pvalc = 1;  if(pvalv<=pvalf) pvalc = pRE2C(pvalv, tha)
  p.value = c(pvalf,pvalv,pvalc);  names(p.value) = c('Pf','Pr','Pc')
  Q=c(Qf,Qv); names(Q) = c('Qf','Qr')
  return(list(p.value=p.value, theta=tha, Q=Q, pars=ans$pars) )
}


## internal func
pPc2 <- function(x,tha) tha*pchisq(x,1,lower=FALSE) + (1-tha)*pchisq(x,2,lower=FALSE)
qPc2 <- function(p0, tha){
  f0 = function(x) pPc2(x,tha) - p0
  q1 = qchisq(p0,1,lower=FALSE); q2 = qchisq(p0,2,lower=FALSE)
  uniroot(f0, c(q1,q2), tol=q1*1e-4)$root
}
## Pc: Pv*I(Pv<=Pf) + I(Pv>Pf)
pRE2C <- function(p0, tha){
  q0v = qPc2(p0,tha); p0f = pchisq(q0v,1,lower=FALSE)
  f1int = function(xx){
    ans = rep(0, length(xx))
    ans[xx>0] = sapply(xx[xx>0], function(x){
      x1 = qchisq(x,1,lower=FALSE)
      qx = qPc2(x,tha)
      pchisq(qx-x1,1)
    })
    ans
  }
  p1 = integrate(f1int, 0,p0f, rel.tol=p0*1e-4)$val
  f2int = function(xx){
    ans = rep(0, length(xx))
    ans[xx>0] = sapply(xx[xx>0], function(x){
      x1 = qchisq(x,1,lower=FALSE)
      qx = qPc2(x,tha)
      pchisq(qx-x1,1) - pchisq(q0v-x1,1)
    })
    ans
  }
  p2 = integrate(f2int, p0f,p0, rel.tol=p0*1e-4)$val
  p0 - (p1 + p2)*(1-tha) - tha*p0f
}

qRE2C <- function(p0, tha){
  psi0 = PrPc1(tha)
  if(p0>=psi0) return(1)
  n0 = round(1/p0)
  f0 = function(x)  ( pRE2C(x,tha) - p0 )*n0
  q1 = p0; q2 = 2*p0
  while(f0(q2)<0) q2=q2*2
  uniroot(f0, c(q1,q2), tol=q1*1e-4)$root
}

###
PrPc1 <- function(tha){
  f0 = function(xx){
    sapply(xx, function(x){
      if((x==0)|(x==1)) return(0)
      xf = qchisq(x,1, lower=FALSE)
      xv = qPc2(x,tha)
      (1-tha)*pchisq(xv-xf,1, lower=FALSE)
    })
  }
  integrate(f0, 0,1)$val
}


