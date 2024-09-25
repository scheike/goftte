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

Matrix<double> Whatfcov_cox(const unsigned &k,
            const Matrix<unsigned> &l,
            const Matrix<double> &Iz,
            const Matrix<double> &Itz,
            const Matrix<double> &Betaiid,
            const Matrix<double> &g,
            const Matrix<double> &dNi,
            const Matrix<double> &dMi,
            const Matrix<double> &Mi,
            const Matrix<double> &Ni,
            const unsigned &type){
    
  unsigned n=dNi.rows();
  
/*Calcul de Wz*/

  Matrix<double>Wz(n,l[k]);
  if (type==1){
    Matrix<double> Wz_prep1=multCol(Iz,t(Ni));
    Matrix<double> Wz_prep2=ProdMat(dNi,g);
    Wz=Wz_prep1-Wz_prep2;
  }

  if (type==2){
    Matrix<double> Wz_prep1=multCol(Iz,t(Mi));
    Matrix<double> Wz_prep2=ProdMat(dMi,g);
    Wz=Wz_prep1-Wz_prep2;
  }
  /*Calcul du score*/
    
  Matrix<double> Wscorerate(n,l[k]); 

       
  for (unsigned i=0; i<n; i++) {
    Matrix<double> Betaiidj = Betaiid(_,i);
    for (unsigned j=0; j<l[k]; j++) {
      Matrix<double> Itzj = Itz(j,_);
      Wscorerate(i,j) = Wz(i,j) - (Itzj*Betaiidj)[0];
    }
  }
  return(Wscorerate);
}

extern "C" {
  void Wfcov(const int *R, // Number of realizations
  	 const unsigned *n, // Number of individuals
     const unsigned *m,
		 const unsigned *nd, // Number of interest events
     const unsigned *nc, // Number of censory times
		 const unsigned *p, // Number of  parameters
     const unsigned *l_data,
		 const double *beta_data,  // nxp, parameter vector
		 const double *time_data, // 
     const unsigned *index_otime_data,
		 const unsigned *index_dtimes_data, // interest events times
     const unsigned *index_censtimes_data, // censory times
		 const double *X_data,  // nxp, Design matrix
     const unsigned long *seed,
     const unsigned *index_ox_data,
     const double *X_data_sort,
		 const double *Mt_data, // Martingale residuals
		 const unsigned *plotnum,
     const unsigned *type_test_num,
		 double *KS,
		 double *Wsd,
		 double *cvalues,
		 double *Ws,
		 double *W
    ) {
                 
   Matrix<double, Col> X(*n, *p, X_data);
   Matrix<double, Col> X_sort(*n, *p, X_data_sort);  
   const unsigned type=*type_test_num;
   Matrix<double> beta(*p, 1, beta_data);
   Matrix<double> time(*n, 1, time_data);
   Matrix<unsigned> index_dtimes(*nd, 1, index_dtimes_data); 
   Matrix<unsigned> index_censtimes(*nc, 1, index_censtimes_data); 
   Matrix<unsigned> index_otime(*m, 1, index_otime_data);
   Matrix<unsigned> index_ox_prep(*p,*n, index_ox_data);
   Matrix<unsigned> l(*p, 1, l_data);
   Matrix<double> dtimes = chrows(time, index_dtimes);
   Matrix<double> censtimes = chrows(time, index_censtimes);
   Matrix<unsigned> index_times = seqa(*n-1,-1,*n);
   Matrix<double> xbeta = ProdMat(X,beta);
   Matrix<double> RR = exp(xbeta);
    
    /*Calcul de Yi*/
    
   Matrix<double> Yiprep(*n,*n);
   Matrix<double> dNiprep(*n,*n);
    
   for (unsigned j=0; j<*n; j++){
      for (unsigned i=0; i<*nd; i++){
        dNiprep(j,index_dtimes[i]) =  (time[j] == dtimes[i]);
        Yiprep(j,index_dtimes[i]) =  (time[j] <= dtimes[i]);
      }
      for (unsigned i=0; i<*nc; i++){
      Yiprep(j,index_censtimes[i]) =  (time[j] <= censtimes[i]);
      }
   }
  
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
  
  /*Calcul de S_0*/
    
   Matrix<double> Si=multCol(Yi,RR);
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
 
 /*Calcul de Mi*/

  Matrix<double> Mi=SumMat(t(dMi));
  
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
    
  /*Calcul de Betaiid*/
  
  Matrix<double> Siinf(*n,*p);

if (type==1){
  Matrix<double> Siinf_prep1=multCol(X,t(Ni));
  Matrix<double> Siinf_prep2=ProdMat(dNi,E);
  Siinf=Siinf_prep1-Siinf_prep2;
}

if (type==2){
  Matrix<double> Siinf_prep1=multCol(X,t(Mi));
  Matrix<double> Siinf_prep2=ProdMat(dMi,E);
  Siinf=Siinf_prep1-Siinf_prep2;
}

  Matrix<double> Betaiid =ProdMat(Iinv,t(Siinf));
   
    /*Calcul de Iz*/
 
  unsigned lmax=max(l);
  Matrix<double> Ws_prep((*plotnum)*(lmax),*p);
  Matrix<double> Wsd_prep(lmax,*p);
  Matrix<double> W_prep(lmax,*p);
  Matrix<double> cvalues_prep(*R,*p);
  mersenne myrng; myrng.initialize(*seed);

  for (unsigned k=0; k<*p; k++){
    Matrix<unsigned> index_ox(l[k],1);
    for (unsigned i=0; i<l[k]; i++){
      index_ox[i]=index_ox_prep(k,i);
    }
  
    Matrix<double> Iz_prep(*n,*n);
    for (unsigned i=0; i<*n; i++){
      for (unsigned j=0; j<*n; j++){
      Iz_prep(i,j) =  (X(i,k) <= X_sort(j,k));
      }
    }
    
    Matrix<double>Izprep=t(Iz_prep);
    Matrix<double> Izprept(l[k],*n);
    Izprept=chrows(Izprep,index_ox);
    Matrix<double> Iz(*n,l[k]);
    Iz=t(Izprept);
  
    /*Calcul de g(B,t,z)*/
  
    Matrix<double>  g_prep(*m,*n);
    Matrix<double>  g(*m,l[k]);
  
    for (unsigned i=0; i<*n; i++){
      for (unsigned j=0; j<*m; j++){
      g_prep(j,i)=Si(i,j)/S_0[j];
      }
    }

    g=ProdMat(g_prep,Iz);
    
    /*Calcul de Itz*/
    
    Matrix<double>Itz(l[k],*p);
    for (unsigned c=0; c<l[k]; c++){
     Matrix<double> Izc=Iz(_,c);
     Matrix<double> ZIz=multCol(X,Izc);
     Matrix<double> Sz2=ProdMat(t(Si),ZIz);
     Matrix<double> Sz2lambda=ProdMat(t(Sz2),lambda);
    
     Matrix<double> Sz1=ProdMat(t(Si),Izc);
     Matrix<double> Sz1lambdaprep=multCol(E,Sz1);
     Matrix<double> Sz1lambda=ProdMat(t(Sz1lambdaprep),lambda);
     Matrix<double> Itzc=Sz2lambda-Sz1lambda;
     Itz(c,_)=t(Itzc);
   }
    
    /*Calcul du score résiduel de Schoenfeld*/
    
    Matrix<double> Ut=ProdMat(Mi,Iz);
    
    /* Calcul de Wz*/

    Matrix<double> WW=Whatfcov_cox(k,l,Iz,Itz,Betaiid,g,dNi,dMi,Mi,Ni,type);
    Matrix<double> sdW = apply(WW,2,ss2);
    Matrix<double> Score=Ut;
    
    unsigned M1=0;
    for (unsigned i=0; i<l[k]; i++){
      if (sdW[i]!=0)
        {M1=M1+1;}
    }
        
    double KSobs = KolmogorovSmirnov(Score);
    unsigned KScount=0; 
    Matrix<double> Res(min((double)*plotnum,(double)*R),l[k]);
  
    Matrix<double>Gi_(*R,*n,false);
    for (int i=0; i<*R; i++) {
      for (unsigned j=0; j<*n; j++) {
      Gi_(i,j)=myrng.rnorm(0, 1);
      }
    }
  
    Matrix<double>ScoreSimprep=ProdMat(Gi_,WW);

    for (int i=0; i<*R; i++) {
      Matrix<double> ScoreSim = ScoreSimprep(i,_);
      Matrix<double>ScoreSim_bis(M1,1);
      Matrix<double>sdW_bis(M1,1);
      
      unsigned M2=0;
       for (unsigned c=0; c<l[k]; c++){
       if (sdW[c]!=0)
        {sdW_bis[M2]=sdW[c];
        ScoreSim_bis[M2]=ScoreSim[c];
        M2=M2+1;}
        }
        
      cvalues_prep(i,k) = max(fabs(ScoreSim_bis/sdW_bis));
      double KShat = KolmogorovSmirnov(ScoreSim);
      if (KShat>KSobs) KScount++;
      if ((unsigned)i< *plotnum) { Res(i,_) = ScoreSim; }
    }
    
    KS[k] = (double)KScount/(double)(*R);
                 
    for (unsigned s=0; s< lmax; s++) {  
      for (unsigned r=0; r<Res.rows(); r++) {
        unsigned pos = r*lmax+s;
        if (s<Res.cols()){
        Ws_prep(pos,k) = Res(r,s);}
        else {Ws_prep(pos,k) =0;}
      }
    }
    
    for (unsigned i=0; i< l[k]; i++) { 
      Wsd_prep(i,k) = sdW[i];
    } 
    
    for (unsigned r=0; r< l[k]; r++) 
      W_prep(r,k) = Score[r];
		 }

    unsigned m3=0;
    for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<lmax;j++){
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
      for (unsigned j=0;j<lmax;j++){
        Wsd[m3]=Wsd_prep(j,i);
        m3=m3+1;
      }
    }
    
    m3=0;
    unsigned mplots=lmax*(*plotnum);
     for (unsigned i=0;i<*p;i++){
      for (unsigned j=0;j<mplots;j++){
        Ws[m3]=Ws_prep(j,i);
        m3=m3+1;
      }
    }
    
  }  
} // extern "C"
