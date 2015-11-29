#include "RooAbsReal.h"
#include <string>

#include "Pdf_Base.h"

void Pdf_Base::setRelation(std::string varName, RooAbsReal* newVar)
{
  if (_intVars.find(varName)==_intVars.end()) {std::cout<<"Pdf_Base::setRelation requested key " << varName << " not found. Exiting!"<<std::endl;exit(2);}
  // Clear old memory allocation for variable
  delete _intVars[varName];
  // Assign new
  _intVars[varName]=newVar;
}
