\name{fcov.coxph}
\alias{fcov.coxph}
\alias{fcov}
\title{Checking the functional form assumption for Cox model}
       
       \description{
Diagnostics for the linear functional form of specific covariate(s) of the Cox model. The method is based on Lin's approximation of the randomly pertubed sum process (1993). Another approximation method recommended by Liu (2008) is also proposed. P-values are derived for supremum KS test statistics.}

\usage{
\method{fcov}{coxph}(model, variable=NULL, type.test=c("Lin"), R=1000, plots=min(R,50), 
seed=NULL,...)
}
\arguments{
  \item{model}{Model object (\code{coxph}).}
  
  \item{variable}{Vector corresponding to the labels of each covariate. This is also an output argument.}
  
  \item{type.test}{Type of approximation. Values are "Lin" or "Liu". Default is "Lin".}
  
  \item{R}{Generation number used for Monte-Carlo
           simulations. This is also an output argument.}
  
  \item{plots}{Realizations number of Monte-Carlo simulations to save for use in the
               plot-routine.}
  
  \item{seed}{Random seed.}
  
  \item{...}{additional arguments.}
}

\details{
A fcov.cph function accepting models of cph type (rms package) is also included in the package, and uses the same arguments as coxph. Contrary to fcov.coxph, the use of 'x=TRUE' option for covariates design is now required in cph model. Strata and cluster options are not supported.}

\value{
  Returns an object of class 'scproc'. The main items of this object are :
    
  \item{obs}{l x p matrix of unique covariates values for functional form. l is the maximum number of unique observations between the p covariates.}
  
  \item{W}{The process \eqn{U(\widehat{\beta},z)}.}
  
  \item{What}{The simulated limiting processes for the R-plots first Monte-Carlo realizations.}
  
  \item{sdw}{Standard error over covariates values of What.}
  
  \item{cvalues}{R x p matrix whose components are the supremum of the standardized What process for each Monte-Carlo realization . The quantiles of this output argument are used to calculate the prediction bands in the
plot-routine.}
  
  \item{KS}{Vector of the p rejection probabilities using KS type statistic.}
}

\examples{
require(survival)
  
#Simulating survival data
simcox <- function(n=100, seed=1) {
  if (!is.null(seed))
  set.seed(seed)
  time<-rexp(n); cen<-2*rexp(n);
  status<-(time<cen);
  time[status==0]<-cen[status==0];
  X<-matrix(rnorm(2*n),n,2)
  return(data.frame(time=time, status=status, X))
}
n <- 100; d <- simcox(n); 
  
#Fitting Cox model
fit.coxph <- coxph(Surv(time,status)~ X1 + X2, data=d)
  
#Checking the covariates functional form assumption
fcov(fit.coxph) 
}
\author{
  Patrick Sfumato and Jean-Marie Boher.
}
\references{
  Lin DY, Wei JL and Ying Z (1993).\emph{Checking the Cox model with cumulative sums of martingale-based residuals.} Biometrika, 80(3), 557-572.
  
  Liu M, Lu W and Shao (2008). \emph{A Monte Carlo approach for change-point detection in the Cox proportional hazards model.} Statistics in Medecine, 27(19), 3894-3909.  
}
