#ifndef myGaussian_h
#define myGaussian_h

#include "RooRealVar.h"
#include "Pdf_Base.h"

class myGaussian : public Pdf_Base
{
 public:
  myGaussian(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
  void setMean(RooAbsReal*);
  void setWidth(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _width;

};


#endif
