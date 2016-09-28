#ifndef DoubleGaussian_h
#define DoubleGaussian_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class DoubleGaussian : public Pdf_Base
{
 public:
  DoubleGaussian(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);

  void setMean(RooAbsReal*);
  void setWidth(RooAbsReal*);
  void setWidthRatio(RooAbsReal*);
  void setFrac(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _sigma1;
  double _frac;
  double _sigmaRatio;

};
#endif 
