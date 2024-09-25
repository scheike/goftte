\name{print.scproc}
\alias{print.scproc}

\title{Printing method for 'scproc' objects}

\description{
  Produces a printed summary of assumption diagnostics.
}

\usage{
\method{print}{scproc}(x, idx=NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{scproc}.} 
   \item{idx}{The index number of covariates whose print is required.}  
  \item{...}{additional arguments.}
}

\details{This function is derived from print.cumres of gof package.} 

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

#Summary of covariates functional form assumption diagnostics
print(k)
}
\author{
  Patrick Sfumato and Jean-Marie Boher.
}

\references{
   Holst KK (2014). \emph{gof: Model-diagnostics based on cumulative residuals.} R package version 0.9.1, URL \url{https://CRAN.R-project.org/package=gof}.
}
