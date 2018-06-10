# GSmeta
 - An R package for various useful methods for meta-analysis with main applications to GWAS.
 - Reference
    - Wu,B. and Zhao,H. (2018) Powerful random effects modeling for meta-analysis of genome-wide association studies. *tech rep*.
 - Sample R codes
```R
 ## install the package
 devtools::install_github('baolinwu/GSmeta')
 library(GSmeta)
 ## simulate data
 K = 5; rho = 0.2
 V = diag(K)*(1-rho) + rho
 U0 = rnorm(K)*sqrt(1-rho) + rnorm(1)*sqrt(rho)
 U1 = U0 + rnorm(K)*1.2
 REma(U0,V)
 REma(U1,V)
```


