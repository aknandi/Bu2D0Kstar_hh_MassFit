#include "RooRealVar.h"
#include "RooChebychev.h"

#include "Linear.h"
#include "Settings.h"

Linear::Linear(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="CombinatoricShape_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("CombinatoricShape_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters
  _coef=mySettings.getD(m+"_"+p+"_"+t+"_combs_slope");

  // Create any RooRealVars you'll need
  _intVars["slope"]=0;
}

RooAbsPdf* Linear::getPdf()
{
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["slope"]==0) _intVars["slope"]=new RooRealVar(Form("%s_Var_Coef",_name.c_str()),"",_coef);

  RooAbsPdf* returnPdf = new RooChebychev(_name.c_str(),"", *_mB, *_intVars["slope"]);
  return returnPdf;
}
