#ifndef DoubleCrystalBall_h
#define DoubleCrystalBall_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class DoubleCrystalBall : public Pdf_Base
{
 public:
  DoubleCrystalBall(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
		
  RooAbsPdf* getPdf();
 private:
  double _mean1, _width1, _alpha1, _n1;
  double _frac1;
  double _mean2, _width2, _alpha2, _n2;

};
#endif 
