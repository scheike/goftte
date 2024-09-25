\name{prop.coxph}
\alias{prop.coxph}
\alias{prop}
\title{Checking proportional hazards assumption for Cox model}

\description{
Diagnostics for proportional hazards of specific covariate(s) of the Cox model. The limiting null distribution of the score process is approximated using either Lin's method (1993) or Liu's (2008). P-values are derived for KS, CvM and AD statistics. }
\usage{
\method{prop}{coxph}(model, variable=NULL, type.test=c("Lin"), R=1000, plots=min(R,50), seed=NULL,...)
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

      - The same type of function was proposed in gof package (cumres.coxph) at the difference that the autors only used Lin's method to approximate the limiting null distribution of the score process.
      
      - A prop.cph function accepting models of cph type (rms package) is also included in the package, and uses the same arguments as prop.coxph. Contrary to prop.coxph, the use of 'x=TRUE' option for covariates design is now required in cph model. Strata and cluster options are not supported.
      }

\value{
  Returns an object of class 'scproc'. The main items of this object are :
  
  \item{obs}{m x p matrix of unique times. m is the length of unique times.}

  \item{W}{The process \eqn{U(\widehat{\beta},t)}.}

  \item{What}{The simulated limiting processes for the R-plots first Monte-Carlo realizations.}

  \item{sdw}{Standard error over time of What.}

  \item{cvalues}{R x p matrix whose components are the supremum of the standardized What process for each Monte-Carlo realization . The quantiles of this output argument are used to calculate the prediction bands in the
plot-routine.}

  \item{KS}{Vector of the p rejection probabilities using KS type statistic.}
  
  \item{CvM}{Vector of the p rejection probabilities using CvM type statistic.}
  
  \item{AD}{Vector of the p rejection probabilities using AD type statistic.}
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
  
#Fitting Cox's model
fit.coxph <- coxph(Surv(time,status)~ X1 + X2, data=d)
  
#Checking the proportional hazards assumption
prop(fit.coxph) 
}
\author{
  Patrick Sfumato and Jean-Marie Boher.
}
\references{

  Holst KK (2014). \emph{gof: Model-diagnostics based on cumulative residuals.} R package version 0.9.1, URL \url{https://CRAN.R-project.org/package=gof}.
  
  Lin DY, Wei JL and Ying Z (1993).\emph{Checking the Cox model with cumulative sums of martingale-based residuals.} Biometrika, 80(3), 557-572.
  
  Liu M, Lu W and Shao (2008). \emph{A Monte Carlo approach for change-point detection in the Cox proportional hazards model.} Statistics in Medecine, 27(19), 3894-3909.  
}
