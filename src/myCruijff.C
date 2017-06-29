#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCruijff.h"
#include "TRandom.h"

#include "Settings.h"
#include "myCruijff.h"

myCruijff::myCruijff(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="myCruijff_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("myCruijff_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  Settings genConfs = Settings("GenSettings");
  genConfs.readPairStringsToMap("Settings/GeneralSettings.txt");

  // Read in the relevant parameters (shared for all modes)
  bool lambdasystematic;
  if(genConfs.get("lckst")=="1") {
    lambdasystematic = true;
  }
  else {
    lambdasystematic = false;
  }

  _mean	  = 5280  + (lambdasystematic?(gRandom->Gaus(0,18)):0.);
  _sigmaL = 221   + (lambdasystematic?(gRandom->Gaus(0,26)):0.);
  _sigmaR = 96    + (lambdasystematic?(gRandom->Gaus(0,16)):0.);
  _alphaL = -0.19 + (lambdasystematic?(gRandom->Gaus(0,0.19)):0.);
  _alphaR = -0.04 + (lambdasystematic?(gRandom->Gaus(0,0.06)):0.);

  // Create any RooRealVars you'll need
  _intVars["mean"]		= 0;
  _intVars["sigmaL"]	= 0;
  _intVars["sigmaR"]	= 0;
  _intVars["alphaL"] 	= 0;
  _intVars["alphaR"] 	= 0;

}

void myCruijff::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void myCruijff::setSigmaL(RooAbsReal* newSigmaL)
{
	setRelation("sigmaL",newSigmaL);
}

void myCruijff::setSigmaR(RooAbsReal* newSigmaR)
{
	setRelation("sigmaR",newSigmaR);
}

void myCruijff::setAlphaL(RooAbsReal* newAlphaL)
{
	setRelation("alphaL",newAlphaL);
}

void myCruijff::setAlphaR(RooAbsReal* newAlphaR)
{
	setRelation("alphaR",newAlphaR);
}

RooAbsPdf* myCruijff::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) 	_intVars["mean"] 		= new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["sigmaL"]==0)	_intVars["sigmaL"]		= new RooRealVar(Form("%s_Var_SigmaL",_name.c_str()),"",_sigmaL);
  if(_intVars["sigmaR"]==0)	_intVars["sigmaR"]		= new RooRealVar(Form("%s_Var_SigmaR",_name.c_str()),"",_sigmaR);
  if(_intVars["alphaL"]==0)	_intVars["alphaL"]		= new RooRealVar(Form("%s_Var_AlphaL",_name.c_str()),"",_alphaL);
  if(_intVars["alphaR"]==0)	_intVars["alphaR"]		= new RooRealVar(Form("%s_Var_AlphaR",_name.c_str()),"",_alphaR);


  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " Cruijff (mean,sigmaL,sigmaR,alphaL,alphaR): " << _intVars["mean"]->getVal() << ", " << _intVars["sigmaL"]->getVal() << ", " << _intVars["sigmaR"]->getVal() << ", " << _intVars["alphaL"]->getVal() << ", " << _intVars["alphaR"]->getVal()
              << std::endl;
  }

  RooCruijff* cr = new RooCruijff(_name.c_str(),"",*_mB,*_intVars["mean"],*_intVars["sigmaL"],*_intVars["sigmaR"],*_intVars["alphaL"],*_intVars["alphaR"]);

  return cr;
}
