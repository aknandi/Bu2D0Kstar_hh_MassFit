#ifndef DoubleExponential_h
#define DoubleExponential_h

#include "RooRealVar.h"
#include "Pdf_Base.h"

class DoubleExponential: public Pdf_Base
{
 public:
  DoubleExponential(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  double _coef1;
  double _coef2;
  double _frac1;

};
#endif 
