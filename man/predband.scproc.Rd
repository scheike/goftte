\name{predband.scproc}
\alias{predband.scproc}
\alias{predband}

\title{Prediction bands for 'scproc' object}
       
       \description{
       Calculates the prediction bands under the null hypothesis.
       }
       
\usage{
\method{predband}{scproc}(x, idx=1:length(x$variable), level=0.95, cval=NULL, ...)
       }
       \arguments{
       \item{x}{An object of class \code{scproc}.}
       
       \item{idx}{The index number of covariates whose the prediction band is calculated.}
       
       \item{level}{The required level for prediction bands.}
       
       \item{cval}{Parameter used to calculate the symmetric prediction bands defined by the standard error multiplied by cval.}
       
       \item{...}{additional arguments.}
       }
       
\value{ 
  \item{t}{Ordered values of the variable used : unique times for proportionality and unique observations of covariates for functional form.}

  \item{yu}{Upper simultaneous confidence limit.}
}

\details{This function is derived from confint.cumres of gof package.} 
       
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
k<-fcov(model=fit.crr, ftime=ftime,fstatus=fstatus,cov1=cov)

#Prediction bands
predband(k)
}
       \author{
       Patrick Sfumato and Jean-Marie Boher.
       }
       
\references{
   Holst KK (2014). \emph{gof: Model-diagnostics based on cumulative residuals.} R package version 0.9.1, URL \url{https://CRAN.R-project.org/package=gof}.
}
