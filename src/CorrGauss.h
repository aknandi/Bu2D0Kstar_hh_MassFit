#include <vector>
#include "TMatrixD.h"

class TRandom3;

class CorrGauss
{

 public:
  CorrGauss(int seed);
  CorrGauss(TMatrixD &covar, int seed);
  CorrGauss(Int_t nrows, Int_t ncols, const double* data, int seed);
  ~CorrGauss();
  void SetCovarMatrix(TMatrixD &covar);
  void SetCovarMatrix(Int_t nrows, Int_t ncols, const double* data);

  void DoCholesky(int method=1);

  TMatrixD *GetCovarMat() const;
  double GetCovarElement(int rindx, int cindx) const;

  TMatrixD *GetDecompMat();
  double GetDecompElement(int rindx, int cindx);

  void GetRandomParSet(std::vector<double> &pars,bool correlated=true);
  void SetSeed(int seed){ randseed = seed; }

 private:
  TRandom3 *randgen;
  
  TMatrixD *V; //covariance matrix
  TMatrixD *C; //decomposed matrix, V=C*C_T

  int randseed;
};
