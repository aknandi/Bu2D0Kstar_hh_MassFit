#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"

#include "Settings.h"
#include "myCrystalBall.h"

myCrystalBall::myCrystalBall(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="myCrystalBall_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("myCrystalBall_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters (shared for all modes)
/*
  _mean	=	mySettings.getD(p+"_mean");
  _width	= mySettings.getD(p+"_width");
  _alpha	=	mySettings.getD(p+"_alpha");
  _n			=	mySettings.getD(p+"_n");
*/

  // Create any RooRealVars you'll need
  _intVars["mean"]		= 0;
  _intVars["width"]	= 0;
  _intVars["alpha"] 	= 0;
  _intVars["n"]	 		= 0;

}

void myCrystalBall::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void myCrystalBall::setWidth(RooAbsReal* newWidth)
{
	setRelation("width",newWidth);
}

void myCrystalBall::setAlpha(RooAbsReal* newAlpha)
{
	setRelation("alpha",newAlpha);
}

void myCrystalBall::setN(RooAbsReal* newN)
{
	setRelation("n",newN);
}

RooAbsPdf* myCrystalBall::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) 	_intVars["mean1"] 		= new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["width"]==0)	_intVars["width"]		= new RooRealVar(Form("%s_Var_Width",_name.c_str()),"",_width);
  if(_intVars["alpha"]==0)	_intVars["alpha"] 		= new RooRealVar(Form("%s_Var_Alpha",_name.c_str()),"",_alpha);
  if(_intVars["n"]==0)			_intVars["n"]		 		= new RooRealVar(Form("%s_Var_n",_name.c_str()),"",_n);


  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " CB1 (mean,width,n,alpha): " << _intVars["mean"]->getVal() << ", " << _intVars["width"]->getVal() << ", " << _intVars["n"]->getVal() << ", " << _intVars["alpha"]->getVal()
              << std::endl;
  }

  RooCBShape* cb = new RooCBShape(_name.c_str(),"",*_mB,*_intVars["mean"],*_intVars["width"],*_intVars["alpha"],*_intVars["n"]);

  return cb;
}
