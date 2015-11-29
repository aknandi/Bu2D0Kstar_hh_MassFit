#ifndef Exponential_h
#define Exponential_h

#include "RooRealVar.h"
#include "Pdf_Base.h"

class Exponential: public Pdf_Base
{
 public:
  Exponential(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  double _coef;

};
#endif 
