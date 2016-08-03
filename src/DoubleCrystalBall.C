#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"

#include "Settings.h"
#include "DoubleCrystalBall.h"

DoubleCrystalBall::DoubleCrystalBall(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="DoubleCrystalBall_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("DoubleCrystalBall_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
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
  _intVars["alpha"] 	= 0;
  _intVars["n"]	 		= 0;
  _intVars["frac"]		= 0;

}

void DoubleCrystalBall::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void DoubleCrystalBall::setWidth(RooAbsReal* newSigma)
{
	setRelation("sigma1",newSigma);
}

void DoubleCrystalBall::setWidthRatio(RooAbsReal* newSigmaRatio)
{
	setRelation("sigma_ratio",newSigmaRatio);
}

void DoubleCrystalBall::setAlpha(RooAbsReal* newAlpha)
{
	setRelation("alpha",newAlpha);
}

void DoubleCrystalBall::setN(RooAbsReal* newN)
{
	setRelation("n",newN);
}

void DoubleCrystalBall::setFrac(RooAbsReal* newFrac)
{
	setRelation("frac",newFrac);
}

RooAbsPdf* DoubleCrystalBall::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) 	_intVars["mean"] 		= new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["sigma1"]==0) 	_intVars["sigma1"] 		= new RooRealVar(Form("%s_Var_Sigma1",_name.c_str()),"",_sigma1);
  if(_intVars["sigma_ratio"]==0) _intVars["sigma_ratio"] = new RooRealVar(Form("%s_Var_Sigma_ratio",_name.c_str()),"",_sigmaRatio);
  if(_intVars["alpha"]==0)	_intVars["alpha"] 		= new RooRealVar(Form("%s_Var_Alpha",_name.c_str()),"",_alpha);
  if(_intVars["n"]==0)			_intVars["n"]		 		= new RooRealVar(Form("%s_Var_n",_name.c_str()),"",_n);
  if(_intVars["frac"]==0)	_intVars["frac"]		 	= new RooRealVar(Form("%s_Var_frac",_name.c_str()),"",_frac);

  _intVars["sigma2"] = new RooFormulaVar(Form("%s_Var_Sigma2",_name.c_str()),"@0*@1",RooArgList(*_intVars["sigma1"],*_intVars["sigma_ratio"]));

  RooCBShape* cb1 = new RooCBShape(Form("%s_cb1",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma1"],*_intVars["alpha"],*_intVars["n"]);
  RooCBShape* cb2 = new RooCBShape(Form("%s_cb2",_name.c_str()),"",*_mB,*_intVars["mean"],*_intVars["sigma2"],*_intVars["alpha"],*_intVars["n"]);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " CB1 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma1"]->getVal() << ", " << _intVars["n"]->getVal() << ", " << _intVars["alpha"]->getVal()
              << std::endl;
    std::cout << " CB2 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["sigma2"]->getVal() << ", " << _intVars["n"]->getVal() << ", " << _intVars["alpha"]->getVal()
              << std::endl;
    std::cout << "  frac: " << _intVars["frac"]->getVal() << endl;
  }
  RooAddPdf*  cb  = new RooAddPdf(_name.c_str(),"",*cb1,*cb2,*_intVars["frac"]);
  return cb;
}
