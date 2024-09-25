\name{fcov.crr}
\alias{fcov.crr}

\title{Checking the functional form assumption for Fine & Gray model}
       
       \description{
Diagnostics for the linear functional form of specific covariate(s) of the Fine & Gray model. Similarly to Li (2015), we extend the Lin's approximation method to take into account competing events. We also provide an adaptation of Liu's approximation method (2008) for Fine & Gray models. P-values are derived for supremum KS test statistics.}
       
\usage{
\method{fcov}{crr}(model, fstatus, ftime, cov1, cencode=0, failcode=1,  type.test=c("Lin"), R=1000,  
plots=min(R,50),  seed=NULL, variable=NULL,...)
       }
       \arguments{
       \item{model}{Model object (\code{crr}).}
       
       \item{fstatus}{Vector corresponding to the failures of the n patients.}
       
       \item{ftime}{Vector corresponding to the failure times.}
       
       \item{cov1}{Matrix whose columns consist in the components of the p covariates.}
       
       \item{cencode}{Censoring code.}
       
       \item{failcode}{Interest event code. All the failures differing from
       cencode and faildcode are considered as competing
       events.}
       
       \item{type.test}{Type of approximation. Values are "Lin" or "Liu". Default is "Lin".}
       
       \item{R}{Generation number used for Monte-Carlo
           simulations. This is also an output argument.}
       
       \item{plots}{Realizations number of Monte-Carlo simulations to save for use in the
               plot-routine.}
       
       \item{seed}{Random seed.}
       
       \item{variable}{Vector corresponding to the labels of each covariate. This is also an output argument.}
       
       \item{...}{additional arguments.}
       }
       
       \value{
  Returns an object of class 'scproc'. The main items of this object are :
    
  \item{obs}{l x p matrix of unique covariates values. l is the maximum number of unique observation between the p covariates.}
  
  \item{W}{The process \eqn{U(\widehat{\beta},z)} adapted to Fine & Gray model.}
  
  \item{What}{The simulated limiting processes for the R-plots first Monte-Carlo realizations.}
  
  \item{sdw}{Standard error over covariates values of What.}
  
  \item{cvalues}{R x p matrix whose components are the supremum of the standardized What process for each Monte-Carlo realization. The quantiles of this output argument are used to calculate the prediction bands in the
plot-routine.}
  
  \item{KS}{Vector of the p rejection probabilities using KS type statistic.}
}  
       
       \examples{
require(cmprsk)

#Simulating survival data with competing events
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(200),nrow=200)

# Fine & Gray regression
fit.crr <- crr(ftime,fstatus,cov)

#Checking the covariates functional form assumption
fcov(model=fit.crr, ftime=ftime,fstatus=fstatus,cov1=cov)
       }
       \author{
       Patrick Sfumato and Jean-Marie Boher.
       }
       
\references{
 Li J, Scheike TH and Zhang MJ (2015). "Checking Fine & Gray subditribution hazards model with cumulative sums of residuals." Lifetime Data Analysis, 21(2), 197-217.
 
 Lin DY, Wei JL and Ying Z (1993).\emph{Checking the Cox model with cumulative sums of martingale-based residuals.} Biometrika, 80(3), 557-572.

 Liu M, Lu W and Shao (2008). \emph{A Monte Carlo approach for change-point detection in the Cox proportional hazards model.} Statistics in Medecine, 27(19), 3894-3909. } 