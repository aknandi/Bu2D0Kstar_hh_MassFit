#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"

#include "Settings.h"
#include "DoubleGaussian.h"

DoubleGaussian::DoubleGaussian(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="DoubleGaussian_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("DoubleGaussian_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
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
  _intVars["frac"]		= 0;

}

void DoubleGaussian::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void DoubleGaussian::setWidth(RooAbsReal* newSigma)
{
	setRelation("sigma1",newSigma);
}

void DoubleGaussian::setWidthRatio(RooAbsReal* newSigmaRatio)
{
	setRelation("sigma_ratio",newSigmaRatio);
}

void DoubleGaussian::setFrac(RooAbsReal* newFrac)
{
	setRelation("frac",newFrac);
}

RooAbsPdf* DoubleGaussian::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) 	_intVars["mean"] 		= new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["sigma1"]==0) 	_intVars["sigma1"] 		= new RooRealVar(Form("%s_Var_Sigma1",_name.c_str()),"",_sigma1);
  if(_intVars["sigma_ratio"]==0) _intVars["sigma_ratio"] = new RooRealVar(Form("%s_Var_Sigma_ratio",_name.c_str()),"",_sigmaRatio);
  if(_intVars["frac"]==0)	_intVars["frac"]		 	= new RooRealVar(Form("%s_Var_frac",_name.c_str()),"",_frac);

  _intVars["sigma2"] = new RooFormulaVar(Form("%s_Var_Sigma2",_name.c_str()),"@0*@1",RooArgList(*_intVars["sigma1"],*_intVars["sigma_ratio"]));

  RooGaussian* g1 = new RooGaussian(Form("%s_g1",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma1"]);
  RooGaussian* g2 = new RooGaussian(Form("%s_g2",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma2"]);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " G1 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma1"]->getVal() << std::endl;
    std::cout << " G2 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma2"]->getVal() << std::endl;
    std::cout << "  frac: " << _intVars["frac"]->getVal() << endl;
  }
  RooAddPdf*  g  = new RooAddPdf(_name.c_str(),"",*g1,*g2,*_intVars["frac"]);
  return g;
}
