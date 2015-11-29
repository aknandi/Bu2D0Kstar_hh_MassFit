#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "DoubleExponential.h"
#include "Settings.h"

DoubleExponential::DoubleExponential(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="CombinatoricShape_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("DoubleExponential_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters
  _coef1=mySettings.getD(m+"_"+t+"_dblexp_coef1");
  _coef2=mySettings.getD(m+"_"+t+"_dblexp_coef2");
  _frac1=mySettings.getD(m+"_"+t+"_dblexp_frac1");

  // Create any RooRealVars you'll need
  _intVars["coef1"]=0;
  _intVars["coef2"]=0;
  _intVars["frac1"]=0;
}

RooAbsPdf* DoubleExponential::getPdf()
{
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["coef1"]==0) _intVars["coef1"]=new RooRealVar(Form("%s_Var_Coef1",_name.c_str()),"",_coef1);
  if(_intVars["coef2"]==0) _intVars["coef2"]=new RooRealVar(Form("%s_Var_Coef2",_name.c_str()),"",_coef2);
  if(_intVars["frac1"]==0) _intVars["frac1"]=new RooRealVar(Form("%s_Var_Frac1",_name.c_str()),"",_frac1);

  RooAbsPdf* exp_pdf1 = new RooExponential(Form("%s_exp_pdf1",_name.c_str()),"", *_mB, *_intVars["coef1"]);
  RooAbsPdf* exp_pdf2 = new RooExponential(Form("%s_exp_pdf2",_name.c_str()),"", *_mB, *_intVars["coef2"]);

  RooAddPdf* returnPdf= new RooAddPdf(_name.c_str(),"",*exp_pdf1,*exp_pdf2,*_intVars["frac1"]);

  return returnPdf;
}
