\docType{package}
\name{goftte-package}
\alias{goftte}
\alias{goftte-package}
\title{Omnibus tests to check the proportional regression models assumptions}
\description{

Let's denote \eqn{Z_i=(Z_{i,k})_{1\leq k \leq p}} a matrix of \eqn{p} individual regression covariates, and \eqn{\widehat{\beta}=(\widehat{\beta}_{k})_{1\leq k \leq p}} the vector of estimated regression coefficients assuming a Cox (1972) proportional hazards (PH) regression model or a Fine & Gray (1999) proportional subdistribution hazards (PSH) regression model. 
The package implements objective diagnotic tools to assess the key modeling assumptions using cumulative sums of model residuals of the form :

\deqn{U(\widehat{\beta},t,z) = \sum_{i=1}^{n} \int_{0}^{t} f(Z_{i})1\{Z_{i} \le z\}d\widehat{M}_{i}(u).}
Here \eqn{f(.)} denotes a smooth function, \eqn{z} the covariates values and \eqn{\widehat{M}_i(t)} the model residuals, so-called martingales residuals in Cox regression settings.  

To test for PH or PSH assumptions, Lin & al. (1993) and Li & al. (2015) proposed test statistics obtained by letting \eqn{f(Z_{i})=Z_{i}} and \eqn{z {\longrightarrow} \infty}. In that case, the test statistics reduce to the supremum value of the score process \eqn{U(\beta,t)}  evaluated at \eqn{\beta=\widehat{\beta}} (Kolmogorov-Smirnov type test statistic). With respect to covariate \eqn{k}, the statistic is written as follows : \deqn{KS=\sup_{t}\mid \mid U_k(\widehat{\beta},t)\mid \mid}
For checking the functional form of a covariate included in the model, same type of statistics were adopted by considering the process \eqn{U(\widehat{\beta},z)} obtained by letting \eqn{f(Z_{i})=1} and \eqn{t {\longrightarrow} \infty}.

Null limiting distribution of cumulative sum process are obtained by drawing \eqn{R} independent realizations of an asymtpotically equivalent gaussian process (Monte-Carlo simulations). The significance level of KS test statistics were derived from the approximated null distribution. 

Our package can be viewed as an extension of the gof package (Holst, 2014) allowing to check the  functional form assumptions of covariates in the  Cox model and to check Fine and Gray model assumptions. Besides this, it proposes goodness-of-fit tests for PH and PSH based on Cramer-Von Mises (CvM) and Anderson-Darling (AD) type test statistics as recommanded by Kvaloy (2004) :  

\deqn{CvM=\int_{0}^{\infty} U_k(\widehat{\beta},t)^2 dq(t)}    
\deqn{AD=\int_{0}^{\infty} \frac{U_k(\widehat{\beta},t)^2}{q(t)(1-q(t))} dq(t)}
Where \eqn{q(t)=\frac{I_{k,k}(\widehat{\beta},t)}{I_{k,k}(\widehat{\beta},\infty)}}, and \eqn{I(\widehat{\beta},t)} denotes the consistent sample estimate of the score process covariance matrix.}
 
\author{
  Patrick Sfumato and Jean-Marie Boher}
\references{
 
 Cox DR (1972). \emph{Regression models and life-tables.} Journal of the Royal Statistical Society. Series B (Methodological), 34(2), 187-220.
 
 Fine JP and Gray RJ (1999). \emph{A Proportional Hzards Model for the Subdistribution of Competing Risk.} American Statistical Association, 94(446), 496-509.

 Lin DY, Wei JL and Ying Z (1993). \emph{Checking the Cox model with cumulative sums of martingale-based residuals.} Biometrika, 80(3), 557-572.
 
 Li J, Scheike TH and Zhang MJ (2015). \emph{Checking Fine & Gray subditribution hazards model with cumulative sums of residuals.} Lifetime Data Analysis, 21(2), 197-217.
 
 Holst KK (2014). \emph{gof: Model-diagnostics based on cumulative residuals.} R package version 0.9.1, URL \url{https://CRAN.R-project.org/package=gof}.
 
 Kvaloy JT and Neef LR (2004).  \emph{Tests for the Proportional Intensity Assumption Based on the Score Process.} Lifetime Data Analysis, 10(2), 139-157.}
\keyword{package}

