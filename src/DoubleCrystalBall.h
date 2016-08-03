#ifndef DoubleCrystalBall_h
#define DoubleCrystalBall_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class DoubleCrystalBall : public Pdf_Base
{
 public:
  DoubleCrystalBall(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);

  void setMean(RooAbsReal*);
  void setWidth(RooAbsReal*);
  void setWidthRatio(RooAbsReal*);
  void setAlpha(RooAbsReal*);
  void setN(RooAbsReal*);
  void setFrac(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _sigma1, _alpha, _n;
  double _frac;
  double _sigmaRatio;

};
#endif 
