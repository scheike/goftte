##' @export
`scproc` <-
  function(model,...) UseMethod("scproc")

plot.scproc <- function(x, idx=1:length(x$variable),
                        col=c("grey"),
                        ci=FALSE,
                        col.ci="darkblue", col.alpha=0.3, lty.ci=0, level=0.95,
                        legend=c("type1","type2","none"), xlab=NULL, ylab=NULL,
                        ylim=NULL,
                        xlim=NULL,
                        title=NULL,
                        cex.lab=1,
                        cex.main=1,
                        ...) {
  
  newylab <- is.null(ylab)
  newxlab <- is.null(xlab)
  if (is.null(title)==TRUE){
    new.title<-x$variable
  }
  
  if ((is.null(title)==FALSE)&(length(idx)!=length(title))){
    stop("Title must have the same length as the number of plots desired")
  }
  
  if ((is.null(title)==FALSE)&(length(idx)==length(title))){
    new.title<-title
  }
  
  J=1
  for (i in idx) {
    legendtxt <- c(); legendpch <- c(); legendcol <- c(); legendlty <- c(); legendlwd <- c(); legendcex <- c()
    if (is.null(ylim)) {
      ylim. <- max(abs(range(x$W[, , i])))*2*c(-1,1)
    }
    if (is.null(xlim)) {
      idx.na.x<-which(is.na(x$obs[, , i])==FALSE)
      xlim. <- c(min(x$obs[, , i][idx.na.x]),max(x$obs[, , i][idx.na.x]))
    }
    
    if (is.null(xlim)==FALSE) {
      xlim. <- xlim
    }
    
    if (is.null(ylim)==FALSE) {
      ylim. <- ylim
    }
    
    ## Observed process
    main <- ""
    if (newxlab) {
      xlab <- x$variable[i]; 
      if (x$type=="prop") {
        main <- xlab; xlab <- "Time";}
      if(x$type=="fcov"){main <- xlab; xlab <- substitute(p,list(p=x$variable[i]))}
      
    }
    if (newylab) {
      ylab <- substitute(expression(W[p](x)),list(p=x$variable[i]))
    }
    legendtxt <- c(legendtxt, "Observed"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,1); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,2); legendcex <- c(legendcex,1);
    x0 <- na.omit(x$obs[, , i])
    
  
    sampleproc <- function() {
      ## Sample processes
      if (col!="none" && !is.null(col)) {
        for (k in 1:ncol(x$What[, , i])) {
          lines(x$What[, , i][,k][1:length(x0)] ~ x0, type="s", col=col, lwd=1)
        }; lines(x$W[, , i][1:length(x0)] ~ x0, type="s", lwd=2)
      }
    }
    
      if (col!="none" && !is.null(col)) {
      legendtxt <- c(legendtxt, "MC sample"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,col); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,1); legendcex <- c(legendcex,1);}

    pband <- function() {
      ## Prediction bandds
      if ( (ci[1]!="none" && !is.null(ci) && ci[1]!=0) || (ci==TRUE) ) {
        if (ci[1]=="pointwise")
          myCI <- predband(x,idx=i,cval=qnorm(1-(1-level)/2))
        else
          myCI <- predband(x,idx=i,level=level)      
        mystepf <- with(myCI, stepfun(t[1:length(x0)],c(0,yu[1:length(x0)])));
        t <- c();
        epsilon <- 1e-9
        for (k in 1:length(x0)) {
          t <- c(t, myCI$t[k]-epsilon, myCI$t[k])
        }; t <- t[-1]
        yu <- mystepf(t)
        
        lines(yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
        lines(-yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
        tt <- c(t, rev(t))
        yy <- c(yu, rev(-yu))
        polygon(tt,yy, col=col.trans, lty=0)      
      }
    }
    
    if ( (ci[1]!="none" && !is.null(ci) && ci[1]!=0) || (ci==TRUE) ) {
      if (col.alpha==0)
        col.trans <- col.ci
      else 
        col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
    legendtxt <- c(legendtxt, "95% prediction band"); legendpch <- c(legendpch,15); legendcol <- c(legendcol,col.trans); legendlty <- c(legendlty,0); legendlwd <- c(legendlwd,0); legendcex <- c(legendcex,2);}
      
        
    with(x, plot(W[, , i][1:length(x0)] ~ x0, type="n", lwd=2, ylab=ylab, ylim=ylim.,xlim=xlim.,xlab=xlab,cex.lab=cex.lab,cex.main=cex.main,
                 main=ifelse(is.null(title)==TRUE,new.title[i],new.title[J])));    
    if (col.alpha==0) {
      with(x, lines(W[, , i][1:length(x0)] ~ x0, type="s", lwd=2));
      sampleproc()
      pband()
      
    } else {
      with(x, lines(W[, , i][1:length(x0)] ~ x0, type="s", lwd=2));
      sampleproc()      
      pband()
     
      
    }
    
    if (!is.null(legend) && legend[1]!="none" && (legend!=F)) {
      if (legend[1]=="type1"){
        legend("topright", c(
          if (is.null(x$KS)==FALSE) {ifelse(x$KS[i]>=0.001,paste("KS-test: p=",x$KS[i],sep=""),paste("KS-test: p<0.001",sep=""))},
          if (is.null(x$CvM)==FALSE) {ifelse(x$CvM[i]>=0.001,paste("CvM-test: p=",x$CvM[i],sep=""),paste("CvM-test: p<0.001",sep=""))},
          if (is.null(x$AD)==FALSE) {ifelse(x$AD[i]>=0.001,paste("AD-test: p=",x$AD[i],sep=""),paste("AD-test: p<0.001",sep=""))})
      
          , bg="white")}
      else
        {legend("topright", legendtxt, lty=legendlty, pch=legendpch, col=legendcol, lwd=legendlwd, pt.cex=legendcex, bg="white")}
    }
    ylim. <- NULL
    xlim. <- NULL
    J=J+1
  }
  invisible(x)
}


