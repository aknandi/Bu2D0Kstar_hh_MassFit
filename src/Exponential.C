#include "RooRealVar.h"
#include "RooExponential.h"

#include "Exponential.h"
#include "Settings.h"

Exponential::Exponential(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="Exponential_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("Exponential_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters
  _coef=mySettings.getD(m+"_"+p+"_"+t+"_combs_slope");
  std::cout << "Exponential" << std::endl;
  std::cout << " slope: " << _coef << std::endl;

  // Create any RooRealVars you'll need
  _intVars["slope"]=0;
}

void Exponential::setSlope(RooAbsReal* newSlope)
{
	setRelation("slope",newSlope);
}

RooAbsPdf* Exponential::getPdf()
{
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["slope"]==0) _intVars["slope"]=new RooRealVar(Form("%s_Var_Coef",_name.c_str()),"",_coef);

  RooAbsPdf* returnPdf = new RooExponential(_name.c_str(),"", *_mB, *_intVars["slope"]);
  return returnPdf;
}
