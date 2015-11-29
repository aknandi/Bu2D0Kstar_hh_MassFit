#ifndef Pdf_Base_h
#define Pdf_Base_h

#include "RooAbsReal.h"
#include "RooRealVar.h"
using namespace std;

class Pdf_Base {
 public:
  Pdf_Base(){};
  virtual RooAbsPdf* getPdf()=0;
  void setRelation(std::string, RooAbsReal*);
 protected:
  RooRealVar *_mB;
  std::map<std::string,RooAbsReal*> _intVars;
  std::string _name;
};
#endif 
