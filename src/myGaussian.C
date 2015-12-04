#include "RooRealVar.h"
#include "RooGaussian.h"

#include "myGaussian.h"
#include "Settings.h"

myGaussian::myGaussian(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="myGaussian_"+m+"_"+p+"_"+c+"_"+t+"_"+a;

  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("myGaussian_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the relevant parameters
  _mean=mySettings.getD(m+"_"+p+"_"+t+"_gaussian_mean");
  _width=mySettings.getD(m+"_"+p+"_"+t+"_gaussian_width");
  std::cout << "myGaussian" << std::endl;
  std::cout << " mean: " << _mean <<  " width: " << _width << std::endl;

  // Create any RooRealVars you'll need
  _intVars["mean"]=0;
  _intVars["width"]=0;
}


void myGaussian::setMean(RooAbsReal* newMean)
{
	setRelation("mean",newMean);
}

void myGaussian::setWidth(RooAbsReal* newWidth)
{
	setRelation("width",newWidth);
}


RooAbsPdf* myGaussian::getPdf()
{
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["mean"]==0) _intVars["mean"]=new RooRealVar(Form("%s_Var_Mean",_name.c_str()),"",_mean);
  if(_intVars["width"]==0) _intVars["width"]=new RooRealVar(Form("%s_Var_Width",_name.c_str()),"",_width);

  RooAbsPdf* returnPdf = new RooGaussian(_name.c_str(),"", *_mB, *_intVars["mean"], *_intVars["width"]);
  return returnPdf;
}
