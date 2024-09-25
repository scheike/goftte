\name{plot.scproc}
\alias{plot.scproc}

\title{Graphical model diagnostic}

\description{
  Produce a graphical model diagnostic, with the possibility of adding the prediction bands calculated by \command{\link{predband.scproc}}.
}

\usage{
\method{plot}{scproc}(x, idx=1:length(x$variable), col=c("grey"), ci=FALSE,
                      col.ci="darkblue", col.alpha=0.3, lty.ci=0, level=0.95,
                      legend=c("type1","type2","none"), xlab=NULL, ylab=NULL,
                      ylim=NULL, xlim=NULL, title=NULL, cex.lab=1,cex.main=1,...)
}
\arguments{
  \item{x}{An object of class \code{scproc}.}
  
  \item{idx}{The index number of covariates whose the plot is required.}
  
  \item{col}{Boolean argument controlling for the plot of prediction bands. Values are TRUE or FALSE.}
  
  \item{ci}{Color of the processes realizations.}
  
  \item{col.ci}{Color of the prediction bands.}
  
  \item{col.alpha}{Color intensity of the prediction bands.}
  
  \item{lty.ci}{Line thickness of the prediction bands border.}
  
  \item{level}{The required level for prediction bands.}
  
  \item{legend}{Argument controling the type of legend on the plot. "type1" for p-values, "type2" for usual type of legends and "none" for no legend.}
  
  \item{xlab}{Label of the x-axis.}
  
  \item{ylab}{Label of the y-axis.}
  
  \item{xlim}{Vector of two components which indicates the extremum values of x-axis.}
  
  \item{ylim}{Vector of two components which indicates the extremum values of y-axis.}
  
  \item{title}{Vector of length idx whose components are titles of plots.}
  
  \item{cex.lab}{Magnification of x and y labels relative to cex.}
  
  \item{cex.main}{Magnification of titles relative to cex.}
  
  \item{...}{additional arguments.}
}

\details{This function is derived from plot.cumres of gof package.} 

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

#Plotting a graphical model diagnostic
plot(k)
}
\author{
  Patrick Sfumato and Jean-Marie Boher.
}

\references{
  Holst KK (2014). \emph{gof: Model-diagnostics based on cumulative residuals.} R package version 0.9.1, URL \url{https://CRAN.R-project.org/package=gof}.
}
