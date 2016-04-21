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
  void setFractionBu010(RooAbsReal*);
  void setFractionBd010(RooAbsReal*);
  void setFractionBd(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  RooAbsPdf* pdf_bu_g_101;
  RooAbsPdf* pdf_bu_g_010;
  RooAbsPdf* pdf_bu_pi_101;
  RooAbsPdf* pdf_bu_pi_010;

  RooAbsPdf* pdf_bd_pi_101;
  RooAbsPdf* pdf_bd_pi_010;

  RooRealVar* var_G_101;
  RooRealVar* var_G_010;
  RooAbsPdf* pdf_bu_101;
  RooAbsPdf* pdf_bu_010;
  RooAbsPdf* pdf_bd;
  RooAbsPdf* pdf_bu;
  double _frac010_bu;
  double _frac010_bd;
  double _frac_bd;
};
#endif 
