##' @export
`prop` <-
  function(model,...) UseMethod("prop")

`prop.coxph` <- function(model,
         variable=NULL,
         type.test=c("Lin"),
         R=1000, plots=min(R,50), seed=NULL,
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


  ot <- order(time);
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
  
    output <- .C("coxscoreW",
                 R=as.integer(R),
                 n=as.integer(n),
                 m=as.integer(m),
                 nd=as.integer(nd),
                 nc=as.integer(nc),
                 p=as.integer(p),
                 seed=seed,
                 beta_data=as.double(beta),
                 time_data=as.double(time),
                 index_otime_data=as.integer(index.otime-1),
                 index_dtimes_data=as.integer(index.dtimes-1),
                 index_censtimes_data=as.integer(index.censtimes-1),
                 X_data=as.double(X), # nxp
                 Mt_data=as.double(as.numeric(n)),
                 plotnum=as.integer(plots),
                 type_test_num=as.integer(type.test.num),
                 KS=as.double(numeric(p)),
                 CvM=as.double(numeric(p)),
                 AS=as.double(numeric(p)),
                 Wsd=as.double(numeric(p*m)),
                 cvalues=as.double(numeric(p*R)),
                 Ws=as.double(numeric(p*m*plots)),
                 W=as.double(numeric(p*m)),
                 pkg="goftte"
    )

UsedVars <- W <- Wsd <- What <- KS <- CvM <- AS <- allcvalues <- x <- mytype <- c()
  mytype <- "prop"
  KS=output$KS
  CvM=output$CvM
  AS=output$AS
  W=array(output$W, dim=c(m,1,p))
  What=array(output$Ws, dim=c(m,plots,p))
  allcvalues=array(output$cvalues,dim=c(R,1,p))
  Wsd=array(output$Wsd,dim=c(m,1,p))
  x=array(rep(otime,p),dim=c(m,1,p))

res <- list(W=W, What=What,
            obs=x, 
            KS=KS, CvM=CvM, AD=AS,
            cvalues=allcvalues, variable=myvars,
            R=R, sd=Wsd, type=mytype, model="coxph", type.test=type.test,assumption="proportional hazards assumption")
class(res) <- "scproc"
res
}