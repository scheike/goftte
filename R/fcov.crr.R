##' @export
`fcov` <-
  function(model,...) UseMethod("fcov")

                `fcov.crr` <-  function(model, 
                               fstatus, 
                               ftime,
                               cov1,
                               cencode=0,
                               failcode=1, 
                               type.test=c("Lin"),
                               R=1000, 
                               plots=min(R,50),
                               seed=NULL,
                               variable=NULL,
                               ...){
  if(length(type.test)>1)
    stop("Enter a test both")
  
  type.test.num=0
  if(type.test=="Lin"){type.test.num=1}
  if(type.test=="Liu"){type.test.num=2}
  if(type.test.num==0)
    stop("Enter a valid name for the test (Lin or Liu)")
  
  if (is.null(seed)!=TRUE){set.seed(seed)}
  seed=round(runif(1,1,1e9))
  
  if (is.null(dim(cov1)[2])==TRUE){n.cov=length(cov1)
                                   m.cov=1}
  if (is.null(dim(cov1))!=TRUE){n.cov=dim(cov1)[1]
                                m.cov=dim(cov1)[2]}
  
  
  if (m.cov>1){
    idx.na.prep=matrix(0,n.cov,m.cov)
    idx.na.count=NULL
    for (i in 1:n.cov){
      for (j in 1:m.cov){
        if(is.na(cov1[i,j])==TRUE){ 
          idx.na.prep[i,j]<-i}
      }
      idx.na.count[i]<-sum(idx.na.prep[i,])
    }
    
    idx.na<-which(idx.na.count==0)
    cov1<-cov1[idx.na,]}
  if ((m.cov==1)&(is.null(dim(cov1)[2])==TRUE))
  {idx.na<-which(is.na(cov1)==FALSE)
   cov1<-cov1[idx.na]}
  
  if ((m.cov==1)&(is.null(dim(cov1)[2])==FALSE))
  {idx.na<-which(is.na(cov1)==FALSE)
   temp=colnames(cov1)
   if (is.null(temp)==TRUE){temp=c("X 1")}
   cov1<-data.frame(cov1[idx.na,])
   colnames(cov1)=temp}
  
  
  if ((m.cov>1)|(is.null(dim(cov1)[2])!=TRUE)){
    m.X.prep=0
    for (j in 1:m.cov){
      if ((length(unique(cov1[,j]))<=2)|(is.numeric(cov1[,j])==TRUE)){
        m.X=m.X.prep+1
        m.X.prep=m.X}
      if ((length(unique(cov1[,j]))>2)&(is.numeric(cov1[,j])==FALSE)){
        m.X=m.X.prep+length(unique(cov1[,j]))-1
        m.X.prep=m.X}}
    
    X=NULL
    names=NULL
    for (j in 1:m.cov){
      if ((length(unique(cov1[,j]))<=2)|(is.numeric(cov1[,j])==TRUE)){
        m.mat=model.matrix(~cov1[,j])[,-1]
        #m.mat=cov1[,j]
        X=cbind(X,m.mat)
        if (is.null(colnames(cov1)[j])==TRUE){
          new.names<-paste(c("X"),j)
          new.names<-gsub(" ","",new.names)}
        if (is.null(colnames(cov1)[j])!=TRUE){
          new.names<-colnames(cov1)[j]
          new.names<-gsub(" ","",new.names)}
        names<-c(names,new.names)}
      
      if ((length(unique(cov1[,j]))>2)&(is.numeric(cov1[,j])==FALSE)){
        m.mat.bis<-model.matrix(~cov1[,j])[,-1]
        #m.mat.bis<-cov1[,j]
        levels<-levels(factor(cov1[,j]))[-1]
        if (is.null(colnames(cov1)[j])==TRUE){
          new.names<-paste(c("X"),j,levels)
          new.names<-gsub(" ","",new.names)}
        if (is.null(colnames(cov1)[j])!=TRUE){
          new.names<-paste(colnames(cov1)[j],levels)
          new.names<-gsub(" ","",new.names)}
        names<-c(names,new.names)
        X=cbind(X,m.mat.bis)}}
    colnames(X)=names
  }
  
  if ((m.cov==1)&(is.null(dim(cov1)[2])==TRUE)){
    m.X.prep=0
    for (j in 1:m.cov){
      if ((length(unique(cov1))<=2)|(is.numeric(cov1)==TRUE)){
        m.X=m.X.prep+1
        m.X.prep=m.X}
      if ((length(unique(cov1))>2)&(is.numeric(cov1)==FALSE)){
        m.X=m.X.prep+length(unique(cov1))-1
        m.X.prep=m.X}}
    
    X=NULL
    names=NULL
    for (j in 1:m.cov){
      if ((length(unique(cov1))<=2)|(is.numeric(cov1)==TRUE)){
        m.mat=model.matrix(~cov1)[,-1]
        X=cbind(X,m.mat)
        names<-c(names,paste(c("X"),j))}
      names<-gsub(" ","",names)
      if ((length(unique(cov1))>2)&(is.numeric(cov1)==FALSE)){
        m.mat.bis<-model.matrix(~cov1)[,-1]
        levels<-levels(factor(cov1))[-1]
        names<-c(names,paste(c("X"),j,levels))
        names<-gsub(" ","",names)
        X=cbind(X,m.mat.bis)}}
  }
  
  data.na<-data.frame(ftime=ftime[idx.na], fstatus=fstatus[idx.na], X)
  data<-na.omit(data.na)
  ftime=data$ftime
  fstatus=data$fstatus
  X=as.matrix(data[,3:((3+m.X)-1)])
  
  ot <- order(ftime);
  time <- ftime[ot]; 
  
  status <- fstatus[ot]
  
  X <- X[ot,,drop=FALSE]
  n <- length(time)
  
  ncom <- sum((status!=failcode)&(status!=cencode))
  nd <- sum(status==failcode)
  nc <- sum(status==cencode)
  
  if (m.X!=length(model$coef))
    stop("Number of variables must be the same as in model")
  p <- m.X
  
  index.censtimes <- (1:n)[status==cencode]
  censtimes<- time[index.censtimes]
  
  index.dtimes <- (1:n)[status==failcode]
  dtimes <- time[index.dtimes]
  
  index.comptimes <- (1:n)[(status!=failcode)&(status!=cencode)]
  comptimes <- time[index.comptimes]
  
  beta <- model$coef
  
  idxtime=which(time==time)
  otime<-cbind(time,idxtime)
  otime<-otime[!duplicated(otime[,1]),] 
  index.otime<-otime[,2]
  otime<-otime[,1];
  m=length(index.otime)
  
  data.time<-data.frame(time=time)
  
  KM.cens <- summary(survfit(Surv(time,status==0)~1,se.fit=F),times=otime,censored=T)
  G<-KM.cens$surv
  G <- c(1,G[1:m-1])
  data.time.G<-data.frame(time=otime,G)
  data.time.G.sort<-merge(data.time,data.time.G,by="time")
  
  G<-data.time.G.sort$G
  
  
  if(any(is.na(beta))) stop("Over-parametrized model")
  
  if(is.null(variable)==TRUE){
  variable=c(names)}
  
  if(length(variable)!=p) stop("Variables names must have same length than number of variables in model")
  
  myvars <- variable
  myvars.idx <- 1:length(names)
  
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
  
    output <- .C("Wfcovcrr",
                 R=as.integer(R),
                 n=as.integer(n),
                 m=as.integer(m),
                 nd=as.integer(nd),
                 ncom=as.integer(ncom),
                 nc=as.integer(nc),
                 p=as.integer(p),
                 l_data=as.integer(l),
                 G_data=as.double(G),
                 seed=seed,
                 X_data_sort=as.double(X.sort),
                 beta_data=as.double(beta),
                 time_data=as.double(time),
                 index_otime_data=as.integer(index.otime-1),
                 index_dtimes_data=as.integer(index.dtimes-1),
                 index_comptimes_data=as.integer(index.comptimes-1),
                 index_censtimes_data=as.integer(index.censtimes-1),
                 X_data=as.double(X), # nxp
                 index_ox_data=as.integer(index.oX-1),
                 Mt_data=as.double(as.numeric(n)),
                 plotnum=as.integer(plots),
                 type_test_num=as.integer(type.test.num),
                 KS=as.double(numeric(p)),
                 Wsd=as.double(numeric(p*max(l))),
                 cvalues=as.double(numeric(p*R)),
                 Ws=as.double(numeric(p*max(l)*plots)),
                 W=as.double(numeric(p*max(l))),
                 pkg="goftte")
  
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
              R=R, sd=Wsd, type=mytype, model="crr",type.test=type.test,assumption="covariate(s) functional form assumption")
  class(res) <- "scproc"
  res
                }