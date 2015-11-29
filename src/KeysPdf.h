#ifndef KeysPdf_h
#define KeysPdf_h

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class KeysPdf : public Pdf_Base
{
 public:
  KeysPdf(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  RooAbsReal *v_cutoff;
  RooAbsPdf* pdf_low;
		
};
#endif 
