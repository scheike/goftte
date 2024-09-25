##' @export
`predband` <-
  function(x,...) UseMethod("predband")

predband.scproc <- function(x, idx=1:length(x$variable), level=0.95, cval=NULL, ...) {
  t <- c(); yu <- c()
  for (i in idx) {
    if (is.null(cval))
      cval <- quantile(x$cvalues[,,i], level)
    t <- cbind(t,x$obs[,,i])
    yu <- cbind(yu,cval*x$sd[,,i]);
  }
  return(list(t=t,yu=yu));
}
