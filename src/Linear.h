#ifndef Linear_h
#define Linear_h

#include "RooRealVar.h"
#include "Pdf_Base.h"

class Linear: public Pdf_Base
{
 public:
  Linear(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  double _coef;

};
#endif 
