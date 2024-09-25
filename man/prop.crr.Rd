\name{prop.crr}
\alias{prop.crr}

\title{Checking proportional subdistribution hazard assumption for Fine & Gray model}
       
       \description{
Diagnostics for the proportional subdistribution hazards of specific covariate(s) of the Fine & Gray model. Similarly to Li (2015),  the limiting null distribution of the score process is approximated by extending Lin's method (1993) to take into account competing events. An adaptation of Liu's approximation method (2008) for Fine & Gray model is also provided. P-values are derived for KS, CvM and AD statistics.}
       
\usage{
  \method{prop}{crr}(model, fstatus, ftime, cov1, cencode=0, failcode=1, variable=NULL, 
  type.test=c("Liu"), R=1000, plots=min(R,50), seed=NULL,...)
}
\arguments{
       \item{model}{Model object (\code{crr}).}
       
       \item{fstatus}{Vector corresponding to the failures of the n patients.}
       
       \item{ftime}{Vector corresponding to the failure times.}
       
       \item{cov1}{matrix whose columns consist in the components of the p covariates.}
       
       \item{cencode}{Censoring code.}
       
       \item{failcode}{Interest event code. All the failures differing from
       cencode and failcode are considered as competing
       events.}
       
       \item{type.test}{Type of approximation. Values are "Lin" or "Liu". Default is "Liu".}
       
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
  
  \item{obs}{m x p matrix of unique times. m is the length of unique times.}

  \item{W}{The process \eqn{U(\widehat{\beta},t)} adapted to Fine & Gray model.}

  \item{What}{The simulated limiting processes for the R-plots first Monte-Carlo realizations.}

  \item{sdw}{Standard error over time of What.}

  \item{cvalues}{R x p matrix whose components are the supremum of the standardized What process for each Monte-Carlo realization . The quantiles of this output argument are used to calculate the prediction bands in the
plot-routine.}

  \item{KS}{Vector of the p rejection probabilities using KS type statistic.}
  
  \item{CvM}{Vector of the p rejection probabilities using CvM type statistic.}
  
  \item{AD}{Vector of the p rejection probabilities using AD type statistic.}
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

#Checking the proportional subdistribution hazards assumption
prop(model=fit.crr, ftime=ftime,fstatus=fstatus,cov1=cov)
}
\author{
  Patrick Sfumato and Jean-Marie Boher.
}

\references{
 Li J, Scheike TH and Zhang MJ (2015). "Checking Fine & Gray subditribution hazards model with cumulative sums of residuals." Lifetime Data Analysis, 21(2), 197-217.
 
 Lin DY, Wei JL and Ying Z (1993).\emph{Checking the Cox model with cumulative sums of martingale-based residuals.} Biometrika, 80(3), 557-572.

 Liu M, Lu W and Shao (2008). \emph{A Monte Carlo approach for change-point detection in the Cox proportional hazards model.} Statistics in Medecine, 27(19), 3894-3909.}