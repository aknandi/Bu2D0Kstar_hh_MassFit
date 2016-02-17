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
  void setFraction(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  RooAbsPdf* pdf_g_101;
  RooAbsPdf* pdf_g_010;
	RooAbsPdf* pdf_pi_101;
  RooAbsPdf* pdf_pi_010;

  RooRealVar* var_G_101;
  RooRealVar* var_G_010;
  RooAbsPdf* pdf_101;
  RooAbsPdf* pdf_010;
  double _frac010;
};
#endif 
