#ifndef DoubleJohnson_h
#define DoubleJohnson_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class DoubleJohnson : public Pdf_Base
{
 public:
  DoubleJohnson(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);

  void setMean(RooAbsReal*);
  void setWidth(RooAbsReal*);
  void setWidthRatio(RooAbsReal*);
  void setGamma(RooAbsReal*);
  void setDelta(RooAbsReal*);
  void setFrac(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _sigma1, _gamma, _delta;
  double _frac;
  double _sigmaRatio;

};
#endif 
