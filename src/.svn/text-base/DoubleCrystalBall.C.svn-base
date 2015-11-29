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
  _mean1	=	mySettings.getD(p+"_mean");
  _width1	= mySettings.getD(p+"_width");
   _alpha1	=	mySettings.getD(p+"_alpha1");
  _n1			=	mySettings.getD(p+"_n1");
  _frac1  = mySettings.getD(p+"_frac1");
  _mean2  = mySettings.getD(p+p+"_mean2");
  _width2 = mySettings.getD(p+"_width2");
  _alpha2	=	mySettings.getD(p+"_alpha2");
  _n2			=	mySettings.getD(p+"_n2");
  */

  // Create any RooRealVars you'll need
  _intVars["mean1"]		= 0;
  _intVars["mean2"]		= 0;
  _intVars["width1"]	= 0;
  _intVars["width2"]	= 0;
  _intVars["alpha1"] 	= 0;
  _intVars["n1"]	 		= 0;
  _intVars["frac1"]		= 0;
  _intVars["alpha2"] 	= 0;
  _intVars["n2"]	 		= 0;
}

RooAbsPdf* DoubleCrystalBall::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean1"]==0) 	_intVars["mean1"] 		= new RooRealVar(Form("%s_Var_Mean1",_name.c_str()),"",_mean1);
  if(_intVars["mean2"]==0) 	_intVars["mean2"] 		= new RooRealVar(Form("%s_Var_Mean2",_name.c_str()),"",_mean2);
  if(_intVars["width2"]==0)	_intVars["width2"]		= new RooRealVar(Form("%s_Var_Width2",_name.c_str()),"",_width2);
  if(_intVars["alpha1"]==0)	_intVars["alpha1"] 		= new RooRealVar(Form("%s_Var_Alpha1",_name.c_str()),"",_alpha1);
  if(_intVars["n1"]==0)			_intVars["n1"]		 		= new RooRealVar(Form("%s_Var_n1",_name.c_str()),"",_n1);
  if(_intVars["frac1"]==0)	_intVars["frac1"]		 	= new RooRealVar(Form("%s_Var_frac1",_name.c_str()),"",_frac1);
  if(_intVars["alpha2"]==0)	_intVars["alpha2"] 		= new RooRealVar(Form("%s_Var_Alpha2",_name.c_str()),"",_alpha2);
  if(_intVars["n2"]==0)			_intVars["n2"]		 		= new RooRealVar(Form("%s_Var_n2",_name.c_str()),"",_n2);

  RooCBShape* cb1 = new RooCBShape(Form("%s_cb1",_name.c_str()),"",*_mB,*_intVars["mean1"],*_intVars["width1"],*_intVars["alpha1"],*_intVars["n1"]);
  RooCBShape* cb2 = new RooCBShape(Form("%s_cb2",_name.c_str()),"",*_mB,*_intVars["mean2"],*_intVars["width2"],*_intVars["alpha2"],*_intVars["n2"]);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " CB1 (mean,width,n,alpha): " << _intVars["mean1"]->getVal() << ", " << _intVars["width1"]->getVal() << ", " << _intVars["n1"]->getVal() << ", " << _intVars["alpha1"]->getVal()
              << std::endl;
    std::cout << " CB2 (mean,width,n,alpha): " << _intVars["mean2"]->getVal() << ", " << _intVars["width2"]->getVal() << ", " << _intVars["n2"]->getVal() << ", " << _intVars["alpha2"]->getVal()
              << std::endl;
    std::cout << "  frac1: " << _intVars["frac1"]->getVal() << endl;
  }
  RooAddPdf*  cb  = new RooAddPdf(_name.c_str(),"",*cb1,*cb2,*_intVars["frac1"]);
  return cb;
}
