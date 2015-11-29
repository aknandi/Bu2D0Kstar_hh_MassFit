#include <iostream>
#include <cmath>
#include "TRandom3.h"
#include "TDecompChol.h"
#include "CorrGauss.h"

using namespace std;

CorrGauss::CorrGauss(int seed):
  randgen(0),
  V(0),
  C(0)
{randseed = seed;}

CorrGauss::CorrGauss(TMatrixD &covar, int seed):
  randgen(0),
  C(0)
{
  V = new TMatrixD(covar);
  randseed=seed;
}

CorrGauss::CorrGauss(Int_t nrows, Int_t ncols, const double* data, int seed):
  randgen(0),
  C(0)
{
  randseed=seed;
  V = new TMatrixD(nrows,ncols,data);
}

CorrGauss::~CorrGauss()
{

  if(randgen){
    delete randgen;
  }
  randgen=0;
  if(V){
    delete V;
  }
  V=0;
  if(C){
    delete C;
  }
  C=0;

}


void CorrGauss::SetCovarMatrix(TMatrixD &covar)
{
  //set a covariance matrix
  if(V){
    delete V;
  }
  V = new TMatrixD(covar);
  if(C){
    delete C;
  }
  C=0;
}

void CorrGauss::SetCovarMatrix(Int_t nrows, Int_t ncols, const double* data)
{
  //set a new covariance matrix
  if(V){
    delete V;
  }
  V = new TMatrixD(nrows,ncols,data);
  if(C){
    delete C;
  }
  C=0;
}

void CorrGauss::DoCholesky(int method)
{
  //finds a matrix C such that C=C*C_T
  //The matrix C can then be used to transform a vector of 
  //independent gaussian distruibuted random numbers into a
  //vector of correlated random numbers.

  //two methods have been implemented,default method==1 uses root functions
  //to decompose the matrix.  Note that the matrix returned by root
  //is the transpose of the one that we want

  //method==2 uses code taken from:
  //Monte Carolo Theory and Practice, 
  //F. James, Rep. Prog. Phys. Vol 43, 1980, p. 1182 
  //(except upper limit on sum in last eq. should be j-1, not i-1)

  //method==3 does both and compares the two.

  //root class is probably better protected against strange cases, but 
  //I don't know how speedy it is, method==2 is provided first to check against
  //previous work, but also in case root method proves to be too slow.


  if(C){
    delete C;
  }
  C=0;

  if(V==0){
    cerr<<"You have not input the covariance matrix"<<endl;
    cerr<<"You must provide a covariance matrix before attempting to decompose it"<<endl;
    return;
  }
  if(V->GetNrows()!=V->GetNcols()){
    cerr<<"Nrows!=NCols"<<endl;
    cerr<<"Covariance matrix must be square"<<endl;
    return;
  }

  int NPAR = V->GetNrows();
  if(NPAR<1){
    cerr<<"Nrows==0"<<endl;
    cerr<<"Covariance matrix must have some rows"<<endl;
    return;
  }

  if(method>3||method<1) method=1;

  if(method==1){//use root to do it
    TDecompChol *tdc=new TDecompChol(*V);
    bool worked=tdc->Decompose();
    if(!worked){
      cerr<<"Root can not Decompose this covariance matrix"<<endl;
      return;
    }
    TMatrixD *tcheck = new TMatrixD(tdc->GetU());
    C=new TMatrixD(tcheck->T());
    if(tcheck){
      delete tcheck;
    }
    tcheck=0;
    if(tdc){
      delete tdc;
    }
    tdc=0;
  }
  else if(method==2||method==3){ //use my old code
    C=new TMatrixD(V->GetNrows(),V->GetNcols());
    for(int i=0;i<NPAR;i++){
      double denom=(*V)(0,0);
      if(denom<=0){
	cerr<<"Can not do Cholesky since v[0][0]="<<denom<<endl;
	return;
      }
      (*C)(i,0)=(*V)(i,0)/sqrt(denom);
      for(int j=1;j<NPAR;j++){
	if(j<i){
	  double s=0.;
	  for(int k=0;k<=j-1;k++){
	    s+=(*C)(i,k)*(*C)(j,k);
	  }
	  if((*C)(j,j)!=0){
	    (*C)(i,j)=((*V)(i,j)-s)/(*C)(j,j);
	  }
	  else{
	    (*C)(i,j)=0.;
	  }
	}      
	if(i==j){
	  double s=0.;
	  for(int k=0;k<=i-1;k++){
	    s+=(*C)(i,k)*(*C)(i,k);
	  }
	  if((*V)(i,i)-s>=0){
	    (*C)(j,j)=sqrt((*V)(i,i)-s);
	  }
	  else{
	    (*C)(j,j)=0;
	  }
	}
	if(j>i){
	  (*C)(i,j)=0.;
	}
      }
    }
    //double check, c*c_T=V
    for(int i=0;i<NPAR;i++){
      for(int j=0;j<NPAR;j++){
	double s=0.;
	for(int k=0;k<NPAR;k++){
	  s+=(*C)(i,k)*(*C)(j,k);
	}
	if(fabs(s-(*V)(i,j))>0.000001){
	  cout<<"Does not check! "<<s<<" v["<<i<<"]["<<j<<"]="<<(*V)(i,j)<<endl;
	}
      }
    }
  }

  if(method==3){//do both, then compare
    TDecompChol *tdc=new TDecompChol(*V);
    bool worked=tdc->Decompose();
    if(!worked){
      cerr<<"Root can not Decompose this covariance matrix"<<endl;
      return;
    }
    TMatrixD *tcheck = new TMatrixD(tdc->GetU());
    TMatrixD *check=new TMatrixD(tcheck->T());
    
    if(tdc){
      delete tdc;
    }
    tdc=0;    
    bool match=true;
    for(int i=0;i<C->GetNrows();i++){
      for(int j=0;j<C->GetNrows();j++){
	if(fabs((*C)(i,j)-(*check)(i,j))>1e-11){
	  match=false;
	  cerr<<"Matrices are not the same!"<<endl;
	  cerr<<"C["<<i<<","<<j<<"]="<<(*C)(i,j)
	      <<" check["<<i<<","<<j<<"]="<<(*check)(i,j)
	      <<"Diff "<<(*C)(i,j)-(*check)(i,j)<<endl;
	}	
      }
    }
    if(match){
      cout<<"The methods match"<<endl;
    }
    if(tcheck){
      delete tcheck;
    }
    tcheck=0;
    if(check){
      delete check;
    }
    check=0;

  }
  
}

TMatrixD *CorrGauss::GetCovarMat() const
{
  //returns the covariance matrix

  if(V){
    return V;
  }
  cerr<<"No covariance matrix"<<endl;
  return 0;
}

double CorrGauss::GetCovarElement(int rindx, int cindx) const 
{ 
  //returns the (rindx,cindx) element of the covariance matrix

  if(V){
    if(rindx>=V->GetNrows()||cindx>=V->GetNcols()){
      cerr<<"Wrong indices, asking for ("<<rindx<<","<<cindx<<")"<<endl;
      cerr<<"From a "<<V->GetNrows()<<"x"<<V->GetNcols()<<" matrix"<<endl;
      return -999.;
    }
    return (*V)(rindx,cindx);
  }
  cerr<<"No covariance matrix"<<endl;
  return -999.;
 

}

TMatrixD *CorrGauss::GetDecompMat()
{
  //returns the decomposition matrix
  //if it doesn't exist, trys to find it using DoCholesky
  //if is still can't find it, it gives up.

  if(C){
    return C;
  }
  DoCholesky();
  if(C){
    return C;
  }
  cerr<<"No decomposed matrix"<<endl;
  return 0;
}


double CorrGauss::GetDecompElement(int rindx, int cindx)
{
  //returns the (rindx,cindx) element of the decomposition matrix
  //if it doesn't exist, trys to find it using DoCholesky
  //if is still can't find it, it gives up.


  if(!C){
    DoCholesky();
  }
  if(C){
    if(rindx>=C->GetNrows()||cindx>=C->GetNcols()){
      cerr<<"Wrong indices, asking for ("<<rindx<<","<<cindx<<")"<<endl;
      cerr<<"From a "<<C->GetNrows()<<"x"<<C->GetNcols()<<" matrix"<<endl;
      return -999.;
    }
    return (*C)(rindx,cindx);
  }
  cerr<<"No covariance matrix"<<endl;
  return -999.;
}

void CorrGauss::GetRandomParSet(vector<double> &pars,bool correlated)
{
  //fills the vector pars with gausian distributed random numbers
  //if correlated==true, the elements of pars will be correlated
  //according to the covariance matrix provided.
  //if correlated==false, the elements will not be correlated
  //the gaussians are centered around 0, user should add the nominal value of 
  //the parameter to shift the mean outside of this code.

  pars.clear();
  vector<double> uncorr;
  if(!randgen){
    randgen = new TRandom3(randseed);
  }
  for(int i=0;i<C->GetNrows();i++){
    //pick uncorrelated random numbers
    if(correlated){
      uncorr.push_back(randgen->Gaus(0,1));
    }
    else{
      pars.push_back(randgen->Gaus(0,sqrt((*V)(i,i))));
    }
  }

  if(!correlated){
    return;
  }
  
  //  cout << "SR do I get here? " << endl ;
  

    //should do this with matrix multiplication, but I don't 
  //know how fast root is
  //transform the vector to get correlated random numbers
  for(int j=0;j<C->GetNrows();j++){
      double row=0.;
      for(int k=0;k<C->GetNcols();k++){
        //           cout << "SR do I get here? 3 k" << uncorr[k] << "but it is " <<(*C)(j,k) << endl ;
      row+=(*C)(j,k)*uncorr[k];
      //      cout << "SR do I get here? 4 row" << row << endl ;
    }
      //    cout << "SR do I get here? 5" << endl ;
    //	cout<<"row "<<row<<" beam_syst "<<beam_syst[j]<<endl;
    pars.push_back(row);      
  }
  //  cout << "SR do I get here? 6" << endl ;
}
