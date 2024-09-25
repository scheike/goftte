// STL
#include <iostream>
#include <algorithm> // random_shuffle, reverse, sort, ...
#include <cmath>
// SCYTHE
#include "mersenne.h"
#include "rng.h"
#include "distributions.h"
#include "ide.h" 
#include "la.h"
#include "matrix.h" 
#include "stat.h" 
#include "smath.h" 
// R interface
#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
#include "cumres.h"
#include "extra.h"

using namespace scythe;
using namespace std;

Matrix<double> Wscorerate_crr(const unsigned &k,
			      const Matrix<double> &X, 
            const Matrix<double> &qY,
            const Matrix<double> &dNi,
            const Matrix<double> &dNic,
            const Matrix<double> &It,
            const Matrix<double> &Betaiid,
            const Matrix<double> &E,
            const Matrix<double> &dMi,
            const Matrix<double> &dMic,
            const unsigned &type) {
    
  unsigned n=X.rows(); unsigned p=X.cols(); unsigned m=dNi.cols();

  /*Calcul de PSYi variable par variable*/
  
  Matrix<double> PSYi(n,m);
  if (type==1){
    for (unsigned i=0; i<n; i++){
      PSYi(i,0)=(X(i,k)-E(0,k))*dNi(i,0)+qY(0,k)*dNic(i,0);
      for (unsigned j=1; j<m; j++){
      PSYi(i,j)=PSYi(i,j-1)+(X(i,k)-E(j,k))*dNi(i,j)+qY(j,k)*dNic(i,j);
      }
    }
  }
  
  if (type==2){
    for (unsigned i=0; i<n; i++){
      PSYi(i,0)=(X(i,k)-E(0,k))*dMi(i,0)+qY(0,k)*dMic(i,0);
      for (unsigned j=1; j<m; j++){
      PSYi(i,j)=PSYi(i,j-1)+(X(i,k)-E(j,k))*dMi(i,j)+qY(j,k)*dMic(i,j);
      }
    }
  }
  
  /*Calcul du score*/
    
  Matrix<double> Wscorerate(n,m); 

  for (unsigned i=0; i<n; i++) {
     Matrix<double> Betaiidj = Betaiid(_,i);
    for (unsigned j=0; j<m; j++) {
      Matrix<double> Itd = It(j,_); Itd.resize(p,p);
      Wscorerate(i,j) = (PSYi(i,j)) - (Itd*Betaiidj)[k];
    }
  } 
  
  return(Wscorerate);
}

extern "C" {
  void crrscoreW(const int *R, // Number of realizations
		 const unsigned *n, // Number of individuals
     const unsigned *m,
		 const unsigned *nd, // Number of interest events
     const unsigned *ncom, // Number of competitive events
     const unsigned *nc, // Number of censory times
		 const unsigned *p, // Number of  parameters
     const double *G_data,
     const unsigned long *seed,
		 const double *beta_data,  // nxp, parameter vector
		 const double *time_data, // 
     const unsigned *index_otime_data,
		 const unsigned *index_dtimes_data, // interest events times
     const unsigned *index_comptimes_data, // competitive events times
     const unsigned *index_censtimes_data, // censory times
		 const double *X_data,  // nxp, Design matrix
		 const double *Mt_data, // Martingale residuals
		 const unsigned *plotnum,
     const unsigned *type_test_num,
		 double *KS,
		 double *CvM,
     double *AS,
		 double *Wsd,
		 double *cvalues,
		 double *Ws,
		 double *W
		 ) {

   Matrix<double, Col> X(*n, *p, X_data);
   Matrix<double> G(*n, 1, G_data);
   const unsigned type=*type_test_num;
   Matrix<double> beta(*p, 1, beta_data);
   Matrix<double> time(*n, 1, time_data);
   Matrix<unsigned> index_comptimes(*ncom, 1, index_comptimes_data); 
   Matrix<unsigned> index_dtimes(*nd, 1, index_dtimes_data); 
   Matrix<unsigned> index_censtimes(*nc, 1, index_censtimes_data); 
   Matrix<unsigned> index_otime(*m, 1, index_otime_data);
   Matrix<double> comptimes = chrows(time, index_comptimes);
   Matrix<double> dtimes = chrows(time, index_dtimes);
   Matrix<double> censtimes = chrows(time, index_censtimes);
   Matrix<unsigned> index_times = seqa(*n-1,-1,*n);
   Matrix<double> xbeta = ProdMat(X,beta);
   Matrix<double> RR = exp(xbeta);
    
    /*Calcul de WYi=Wi*Yi*/
    
  Matrix<double> Gi =chrows(G, index_comptimes);
  Matrix<double> WYiprep(*n,*n);
  Matrix<double> Yiprep(*n,*n);
  Matrix<double> dNiprep(*n,*n);
  Matrix<double> dNicprep(*n,*n);
    
  for (unsigned j=0; j<*n; j++){
    for (unsigned i=0; i<*nd; i++){
      WYiprep(j,index_dtimes[i]) =  (time[j] <= dtimes[i]); 
      dNiprep(j,index_dtimes[i]) =  (time[j] == dtimes[i]);
      Yiprep(j,index_dtimes[i]) =  (time[j] <= dtimes[i]);
    }
    for (unsigned i=0; i<*nc; i++){
      WYiprep(j,index_censtimes[i]) =  (time[j] <= censtimes[i]);
      dNicprep(j,index_censtimes[i])= (time[j] == censtimes[i]);
      Yiprep(j,index_censtimes[i]) =  (time[j] <= censtimes[i]);
    }
    for (unsigned i=0; i<*ncom; i++){
      WYiprep(j,index_comptimes[i]) =  (time[j] <= comptimes[i]); 
      Yiprep(j,index_comptimes[i]) =  (time[j] <= comptimes[i]);
    }
    for (unsigned i=0; i<*ncom; i++){
      if (time[j] > comptimes[i]) {WYiprep(j,index_comptimes[i])=G[j]/Gi[i];
      }
    }
  }
  
  Matrix<double> WYiprept=chrows(WYiprep,index_otime);
  Matrix<double> WYi=t(WYiprept);
  
  /*Calcul de Yi*/
  
  Matrix<double> Yiprept=chrows(Yiprep,index_otime);
  Matrix<double> Yi=t(Yiprept);
  
  /*Calcul de Y*/
  
   Matrix<double>Y=SumMat(Yi);
  
  /* Calcul de dNi*/
  
   Matrix<double> dNiprept=chrows(dNiprep,index_otime);
   Matrix<double> dNi=t(dNiprept);
  
  /*Calcul de Ni*/
  
   Matrix<double> Ni=SumMat(t(dNi));
  
  /*Calcul de dN*/
  
   Matrix<double> dN=SumMat(dNi);
   
   /* Calcul de dNic*/
  
  Matrix<double> dNicprept=chrows(dNicprep,index_otime);
  Matrix<double> dNic=t(dNicprept);
  
  /*Calcul de dNc*/
  
   Matrix<double> dNc=SumMat(dNic);
  
  /*Calcul de S_0*/
    
   Matrix<double> Si=multCol(WYi,RR);
   Matrix<double> S_0=SumMat(Si);
   
   /*Calcul de S_1*/
   
   Matrix<double> S_1=ProdMat(t(Si),X);
  
   /*Calcul de S_2*/
  
  unsigned p2 = (*p)*(*p);
  Matrix<double>X_2(*n,p2,false);
  for (unsigned i=0; i<*n; i++) { 
      X_2(i,_) = crossprod(X(i,_));
  }
  Matrix<double>S_2=ProdMat(t(Si),X_2);
  
 /*Calcul de dMi*/
  
  Matrix<double> lambda=multCol(t(dN),t(1/S_0));
  Matrix<double> lambdaSi=t(multCol(t(Si),lambda));
  Matrix<double> dMi=dNi-lambdaSi;
  
  /*Calcul de dMic*/
  
  Matrix<double> lambdac=multCol(t(dNc),t(1/Y));
  Matrix<double> lambdacYi=t(multCol(t(Yi),lambdac));
  Matrix<double> dMic=dNic-lambdacYi;
 
 /*Calcul de Mi*/

  Matrix<double> Mi=SumMat(t(dMi));
  Matrix<double> Mit=t(cumsum(t(dMi)));
  
   /*Calcul de It et I-1*/
 
  Matrix<double> E= multCol(S_1,t(1/S_0));
  Matrix<double> otime=chrows(time, index_otime);   
  
  Matrix<double> E_2(*m,p2,false);
  for (unsigned i=0; i<*m; i++) { 
    E_2(i,_) = crossprod(E(i,_));
  }
 
  Matrix<double> Itprep = multCol(S_2, t(1/S_0)); 
  Itprep = Itprep-E_2; // Martinussen & Scheike p. 184
  Matrix<double>It=multCol(Itprep,t(dN));
  It=cumsum(It);  
  Matrix<double> Itau = It(*m-1,_); Itau.resize(*p,*p);
  Matrix<double> Iinv=inv(Itau);
  
  /*Calcul de q*/
  
  Matrix<double> Ius(*n,*m);
    for (unsigned j=0; j<*m; j++){
      for (unsigned i=0; i<*n; i++){
      Ius(i,j)= (otime[j] > time[i]);
    }
  }
  
  Matrix<double> qiprep(*n,*m);
  Matrix<double> q(*m,*p);
  for (unsigned c=0; c<*p; c++){
    for (unsigned i=0; i<*n; i++){
    qiprep(i,*m-1)=(X(i,c)-E(*m-1,c))*dMi(i,*m-1)*Ius(i,*m-1);
      for (unsigned j=1; j<*m; j++){
          qiprep(i,*m-1-j)=(qiprep(i,*m-j)+(X(i,c)-E(*m-1-j,c))*dMi(i,*m-1-j))*Ius(i,*m-1-j);
      }
    }
    q(_,c)=-SumMat(qiprep);
  }
 
  /*Calcul de Betaiid*/
  
  Matrix<double> Siinf(*n,*p);
  Matrix<double> qY=multCol(q,t(1/Y));

  if (type==1){
    Matrix<double> Siinf_prep1=multCol(X,t(Ni));
    Matrix<double> Siinf_prep2=ProdMat(dNi,E); 
    Matrix<double> Siinf_prep4=ProdMat(dNic,qY);
    Siinf=Siinf_prep1-Siinf_prep2+Siinf_prep4;
  }

  if (type==2){
    Matrix<double> Siinf_prep1=multCol(X,t(Mi));
    Matrix<double> Siinf_prep2=ProdMat(dMi,E);
    Matrix<double> Siinf_prep4=ProdMat(dMic,qY);
    Siinf=Siinf_prep1-Siinf_prep2+Siinf_prep4;
  }

  Matrix<double> Betaiid =ProdMat(Iinv,t(Siinf));
  
    Matrix<double> Ws_prep((*plotnum)*(*m),*p);
    Matrix<double> Wsd_prep(*m,*p);
    Matrix<double> W_prep(*m,*p);
    Matrix<double> cvalues_prep(*R,*p);
    mersenne myrng; myrng.initialize((unsigned long) *seed);
    
     for (unsigned k=0; k<*p; k++){
       
    /*Calcul du score résiduel de Schoenfeld*/

    Matrix<double> Ut=ProdMat(t(X(_,k)),Mit);

    /*Calcul du score attendu WW*/
            
    Matrix<double> WW = Wscorerate_crr(k,X,qY,dNi,dNic,It,Betaiid,E,dMi,dMic,type);
  
    Matrix<double> sdW = apply(WW,2,ss2);
    
    unsigned M1=0;
       for (unsigned l=0; l<*m; l++){
       if (sdW[l]!=0)
        {M1=M1+1;}
        }
  
    double KSobs = KolmogorovSmirnov(Ut);
    double CvMobs = CramerVonMises(Ut,It,k);
    double ASobs = AndersonDarling(Ut,It,k);
    unsigned KScount=0; 
    unsigned CvMcount=0;
    unsigned AScount=0;
    Matrix<double> Res(min((double)*plotnum,(double)*R),*m);
    
    Matrix<double>Gi_(*R,*n,false);
    for (int i=0; i<*R; i++) {
      for (unsigned j=0; j<*n; j++) {
      Gi_(i,j)=myrng.rnorm(0, 1);
      }
    }
      Matrix<double>ScoreSimprep=ProdMat(Gi_,WW);

  for (int i=0; i<*R; i++) {
      Matrix<double> ScoreSim=ScoreSimprep(i,_);
      Matrix<double>ScoreSim_bis(M1,1);
      Matrix<double>sdW_bis(M1,1);
      
      unsigned M2=0;
       for (unsigned l=0; l<*m; l++){
       if (sdW[l]!=0)
        {sdW_bis[M2]=sdW[l];
        ScoreSim_bis[M2]=ScoreSim[l];
        M2=M2+1;}
        }
        
      cvalues_prep(i,k) = max(fabs(ScoreSim_bis/sdW_bis));
      double KShat = KolmogorovSmirnov(ScoreSim);
      double CvMhat = CramerVonMises(ScoreSim,It,k);
      double AShat = AndersonDarling(ScoreSim,It,k);
      if (KShat>KSobs) KScount++;
      if (CvMhat>CvMobs) CvMcount++;
      if (AShat>ASobs) AScount++;
      if ((unsigned)i< *plotnum) { Res(i,_) = ScoreSim; }
    }
    KS[k] = (double)KScount/(double)(*R);
    CvM[k] = (double)CvMcount/(double)(*R);
    AS[k] = (double)AScount/(double)(*R);

    for (unsigned s=0; s< Res.cols(); s++) {  
      for (unsigned r=0; r<Res.rows(); r++) {
  unsigned pos = r*(Res.cols())+s;
	Ws_prep(pos,k) = Res(r,s);
      }
    } 
    for (unsigned i=0; i< *m; i++) { 
      Wsd_prep(i,k) = sdW[i];
    } 
    for (unsigned r=0; r< *m; r++) 
      W_prep(r,k) = Ut[r];
		 }

    unsigned m3=0;
    for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<*m;j++){
        W[m3]=W_prep(j,i);
        m3=m3+1;
      }
    }
    m3=0;
     for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<*R;j++){
        cvalues[m3]=cvalues_prep(j,i);
        m3=m3+1;
      }
    }
    m3=0;
     for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<*m;j++){
        Wsd[m3]=Wsd_prep(j,i);
        m3=m3+1;
      }
    }
    
    m3=0;
    unsigned mplots=(*m)*(*plotnum);
     for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<mplots;j++){
        Ws[m3]=Ws_prep(j,i);
        m3=m3+1;
      }
    }
	} 
} // extern "C"
