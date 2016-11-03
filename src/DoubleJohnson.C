#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooJohnsonSU.h"

#include "Settings.h"
#include "DoubleJohnson.h"

DoubleJohnson::DoubleJohnson(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="DoubleJohnson_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("DoubleJohnson_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters (shared for all modes)
  /*
   _n			=	mySettings.getD(p+"_n_"+t);
  _alpha = mySettings.getD(p+"_alpha_"+t);
  _frac  = mySettings.getD(p+"_frac");
  _sigmaRatio = mySettings.getD(p+"_width_ratio_"+t);
  */

  // Create any RooRealVars you'll need
  _intVars["mean"]		= 0;
  _intVars["sigma1"]	= 0;
  _intVars["sigma2"] = 0;
  _intVars["sigma_ratio"]	= 0;
  _intVars["gamma"] 	= 0;
  _intVars["delta"]	 		= 0;
  _intVars["frac"]		= 0;

}

void DoubleJohnson::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void DoubleJohnson::setWidth(RooAbsReal* newSigma)
{
	setRelation("sigma1",newSigma);
}

void DoubleJohnson::setWidthRatio(RooAbsReal* newSigmaRatio)
{
	setRelation("sigma_ratio",newSigmaRatio);
}

void DoubleJohnson::setGamma(RooAbsReal* newGamma)
{
	setRelation("gamma",newGamma);
}

void DoubleJohnson::setDelta(RooAbsReal* newDelta)
{
	setRelation("delta",newDelta);
}

void DoubleJohnson::setFrac(RooAbsReal* newFrac)
{
	setRelation("frac",newFrac);
}

RooAbsPdf* DoubleJohnson::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) 	_intVars["mean"] 		= new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["sigma1"]==0) 	_intVars["sigma1"] 		= new RooRealVar(Form("%s_Var_Sigma1",_name.c_str()),"",_sigma1);
  if(_intVars["sigma_ratio"]==0) _intVars["sigma_ratio"] = new RooRealVar(Form("%s_Var_Sigma_ratio",_name.c_str()),"",_sigmaRatio);
  if(_intVars["gamma"]==0)	_intVars["gamma"] 		= new RooRealVar(Form("%s_Var_Gamma",_name.c_str()),"",_gamma);
  if(_intVars["delta"]==0)			_intVars["delta"]		 		= new RooRealVar(Form("%s_Var_Delta",_name.c_str()),"",_delta);
  if(_intVars["frac"]==0)	_intVars["frac"]		 	= new RooRealVar(Form("%s_Var_frac",_name.c_str()),"",_frac);

  _intVars["sigma2"] = new RooFormulaVar(Form("%s_Var_Sigma2",_name.c_str()),"@0*@1",RooArgList(*_intVars["sigma1"],*_intVars["sigma_ratio"]));

  RooJohnsonSU* john1 = new RooJohnsonSU(Form("%s_j1",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma1"],*_intVars["gamma"],*_intVars["delta"]);
  RooJohnsonSU* john2 = new RooJohnsonSU(Form("%s_j2",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma2"],*_intVars["gamma"],*_intVars["delta"]);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " CB1 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma1"]->getVal() << ", " << _intVars["delta"]->getVal() << ", " << _intVars["gamma"]->getVal()
              << std::endl;
    std::cout << " CB2 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma2"]->getVal() << ", " << _intVars["delta"]->getVal() << ", " << _intVars["gamma"]->getVal()
              << std::endl;
    std::cout << "  frac: " << _intVars["frac"]->getVal() << endl;
  }
  RooAddPdf*  john  = new RooAddPdf(_name.c_str(),"",*john1,*john2,*_intVars["frac"]);
  return john;
}
