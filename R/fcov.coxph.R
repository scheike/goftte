##' @export
`fcov` <-
  function(model,...) UseMethod("fcov")

`fcov.coxph` <- function(model,
                         variable=NULL,
                         type.test=c("Lin"),
                         R=1000, plots=min(R,50), 
                         seed=NULL,
                           ...) {
  if(length(type.test)>1)
    stop("Enter a test both")
  
  type.test.num=0
  if(type.test=="Lin"){type.test.num=1}
  if(type.test=="Liu"){type.test.num=2}
  if(type.test.num==0)
    stop("Enter a valid name for the test Lin or Liu")
  
  if (is.null(seed)!=TRUE){set.seed(seed)}
  seed=round(runif(1,1,1e9))
  
  mt <- model.frame(model)
  Y <- model.extract(mt, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  if (attr(Y, "type") == "right") {
    time <- Y[, "time"]; 
    status <- Y[, "status"]
  } else stop("Expected right-censored data.");
  X <- na.omit(model.matrix(model))
  
  
  ot <- order(time); # order in time, status=1 first for ties
  time <- time[ot]; status <- status[ot]
  X <- X[ot,,drop=FALSE]
  n <- length(time)
  nd <- sum(status)
  nc <- sum(status==0)
  p <- ncol(X)
  index.dtimes <- (1:n)[status==1]
  dtimes <- time[index.dtimes]
  
  index.censtimes <- (1:n)[status==0]
  censtimes<- time[index.censtimes]
  
  idxtime=which(time==time)
  otime<-cbind(time,idxtime)
  otime<-otime[!duplicated(otime[,1]),] 
  index.otime<-otime[,2]
  otime<-otime[,1];
  m=length(index.otime)
  
  #Ties 
  
  if ((n>m)&(model$method!="breslow"))
    warning("In case of ties, use breslow method in coxph")
  
  beta <- coef(model)
  if(any(is.na(beta))) stop("Over-parametrized model")
  
  if (is.null(variable)==FALSE){
    if (length(variable)!=p) stop("Variables names must have same length than number of variables in model")}
  
  if (is.null(variable)==TRUE){
    variable <- na.omit(colnames(X))}
  
    
  variable <- unique(variable)
  UsedData <- X[,na.omit(match(variable, colnames(X))),drop=FALSE]
  
  myvars <- variable
  
  myvars.idx <- 1:p
  
  #forme de la covariable
  
  ncov<-dim(X)[2]
  if (is.null(ncov)==TRUE){ncov=1}
  
  if (is.null(ncov)==FALSE)
    l=NULL
  index.oX=NULL
  X.sort<-matrix(NA,dim(X)[1],dim(X)[2])
  for (i in 1:ncov){
    X.sort[,i]=sort(X[,i])
    idx.X=which(X.sort[,i]==X.sort[,i])
    oX<-cbind(X.sort[,i],idx.X)
    oX<-oX[!duplicated(oX[,1]),] 
    index.oX.prep<-oX[,2]
    l0=length(index.oX.prep)
    
    index.oX.prep<-c(index.oX.prep,matrix(0,1,n-l0))
    l=c(l,l0)
    index.oX=rbind(index.oX,index.oX.prep)}
  
  if (is.null(ncov)==TRUE){
    X.sort=sort(X)
    idx.X=which(X.sort==X.sort)
    oX<-cbind(X.sort,idx.X)
    oX<-oX[!duplicated(oX),] 
    index.oX.prep<-oX[,2]
    l=length(index.oX)
    index.oX=index.oX.prep}
  
    output <- .C("Wfcov",
                 R=as.integer(R),
                 n=as.integer(n),
                 m=as.integer(m),
                 nd=as.integer(nd),
                 nc=as.integer(nc),
                 p=as.integer(p),
                 l_data=as.integer(l),
                 beta_data=as.double(beta),
                 time_data=as.double(time),
                 index_otime_data=as.integer(index.otime-1),
                 index_dtimes_data=as.integer(index.dtimes-1),
                 index_censtimes_data=as.integer(index.censtimes-1),
                 X_data=as.double(X), # nxp
                 seed=seed,
                 index_ox_data=as.integer(index.oX-1),
                 X_data_sort=as.double(X.sort),
                 Mt_data=as.double(as.numeric(n)),
                 plotnum=as.integer(plots),
                 type_test_num=as.integer(type.test.num),
                 KS=as.double(numeric(p)),
                 Wsd=as.double(numeric(p*max(l))),
                 cvalues=as.double(numeric(p*R)),
                 Ws=as.double(numeric(p*max(l)*plots)),
                 W=as.double(numeric(p*max(l))),
                 pkg="goftte"
    )
  
  UsedVars <- W <- Wsd <- What <- KS <- CvM <- AS <- allcvalues <- x <- mytype <- c()
  
  mytype <- "fcov"
  KS=output$KS
  W=array(output$W, dim=c(max(l),1,p))
  What=array(output$Ws, dim=c(max(l),plots,p))
  allcvalues=array(output$cvalues,dim=c(R,1,p))
  Wsd=array(output$Wsd,dim=c(max(l),1,p))
  x=array(0,dim=c(max(l),1,p))
  for(i in 1:p)
    x[,,i]=c(unique(X.sort[,i]),rep(NA,max(l)-length(unique(X.sort[,i]))))
  
  res <- list(W=W, What=What,
              obs=x, 
              KS=KS, Wsd=Wsd,
              cvalues=allcvalues, variable=myvars,
              R=R, sd=Wsd, type=mytype, model="coxph",type.test=type.test,assumption="covariate(s) functional form assumption")
  class(res) <- "scproc"
  res
}