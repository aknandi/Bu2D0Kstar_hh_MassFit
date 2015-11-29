#ifndef Lambda_h
#define Lambda_h

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"
#include "RooAddPdf.h"
#include "KeysPdf.h"

class Lambda : public Pdf_Base
{
 public:
  Lambda(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  RooAbsPdf* pdf_pK;
  RooAbsPdf* pdf_pPi;
  double _fracpK;
};
#endif 
