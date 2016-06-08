#ifndef PartRecoDstKst_h
#define PartRecoDstKst_h

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"
#include "RooAddPdf.h"

class PartRecoDstKst : public Pdf_Base
{
 public:
  PartRecoDstKst(RooRealVar*, std::string,std::string,std::string,std::string,std::string);
  void setCoef010(RooAbsReal*);
  void setCoef101(RooAbsReal*);

  RooAbsPdf* getPdf();
 private:
  RooAbsPdf* pdf_bu_gamma_101;
  RooAbsPdf* pdf_bu_gamma_010;
  RooAbsPdf* pdf_bu_pi_101;
  RooAbsPdf* pdf_bu_pi_010;

  RooAbsPdf* pdf_bd_pi_101;
  RooAbsPdf* pdf_bd_pi_010;

  double _coef010;
  double _coef101;

  RooRealVar* frac_pi_101;
  RooRealVar* frac_pi_010;
  RooRealVar* frac_gamma_101;
  RooRealVar* frac_gamma_010;
  RooRealVar* frac_bd_101;
  RooRealVar* frac_bd_010;

  RooAbsPdf* pdf_101;
  RooAbsPdf* pdf_010;
};
#endif 
