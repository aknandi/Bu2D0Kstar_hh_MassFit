#ifndef PartRecoDstKst_h
#define PartRecoDstKst_h

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"
#include "RooAddPdf.h"

class PartRecoDstKst : public Pdf_Base
{
 public:
  PartRecoDstKst(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  RooAbsPdf* pdf_g_001;
  RooAbsPdf* pdf_g_010;
	RooAbsPdf* pdf_pi_001;
  RooAbsPdf* pdf_pi_010;

  RooRealVar* var_G_001;
  RooRealVar* var_G_010;
  RooAbsPdf* pdf_001;
  RooAbsPdf* pdf_010;
  double _frac010;
};
#endif 
