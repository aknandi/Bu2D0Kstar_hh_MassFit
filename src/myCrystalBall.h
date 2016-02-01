#ifndef myCrystalBall_h
#define myCrystalBall_h

#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "Pdf_Base.h"

class myCrystalBall : public Pdf_Base
{
 public:
  myCrystalBall(RooRealVar*, std::string,std::string,std::string,std::string,std::string,std::string);
  void setMean(RooAbsReal*);
  void setWidth(RooAbsReal*);
  void setAlpha(RooAbsReal*);
  void setN(RooAbsReal*);
  RooAbsPdf* getPdf();
 private:
  double _mean, _width, _alpha, _n;

};
#endif 
