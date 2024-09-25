##' @export
`scproc` <-
  function(model,...) UseMethod("scproc")

##' @S3method print cumres
print.scproc <- function(x,idx=NULL,...) {
  cat("\n");
  cat("Rejection p-values associated to ", x$type.test,"'s approximation for ", x$assumption,"\n", sep="");
  cat("---\n");
  N=1:length(x$variable)
  if (is.null(idx)==FALSE){N=idx}
  for (i in N) {
    if (is.null(x$KS[i])==FALSE) {if(x$KS[i]>=0.001){cat("Kolmogorov-Smirnov-test: p-value=", x$KS[i], "\n", sep="")}
                                  if(x$KS[i]<0.001){cat("Kolmogorov-Smirnov-test: p-value<0.001","\n", sep="")}}
    if (is.null(x$CvM[i])==FALSE) {if(x$CvM[i]>=0.001){cat("Cramer-von-Mises-test: p-value=", x$CvM[i], "\n", sep="")}
                                if(x$CvM[i]<0.001){cat("Cramer-von-Mises-test: p-value<0.001","\n", sep="")}}
    if (is.null(x$AD[i])==FALSE) {if(x$AD[i]>=0.001){cat("Anderson-Darling-test: p-value=", x$AD[i], "\n", sep="")}
                               if(x$AD[i]<0.001){cat("Anderson-Darling-test: p-value<0.001","\n", sep="")}}
    cat("Based on ", x$R, " realizations. Cumulated residuals ordered by ", x$variable[i], "-variable.\n", sep="")
    cat("---\n");
  }
  invisible(x)
}

##' @S3method summary cumres
summary.scproc <- function(object,...) print(object,...)
