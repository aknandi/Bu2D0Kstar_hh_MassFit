#ifndef myCruijff_h
#define myCruijff_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class myCruijff : public Pdf_Base
{
 public:
  myCruijff(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
  void setMean(RooAbsReal*);
  void setSigmaL(RooAbsReal*);
  void setSigmaR(RooAbsReal*);
  void setAlphaL(RooAbsReal*);
  void setAlphaR(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _sigmaL, _sigmaR, _alphaL, _alphaR;

};
#endif 
