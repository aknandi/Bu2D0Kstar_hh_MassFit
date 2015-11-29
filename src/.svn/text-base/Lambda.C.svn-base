#include "RooRealVar.h"
#include "RooAbsPdf.h"

#include "Settings.h"
#include "Lambda.h"

Lambda::Lambda(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="Lambda_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("LambdaSettings_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  Settings genConfs = Settings("GenSettings");
  genConfs.readPairStringsToMap("Settings/GeneralSettings.txt");
  string fitlow = genConfs.get("fit_limit_low");

  // Read in the relevant parameters 
  // Fix the ratio between the p-K, p-Pi misID from MC
  _fracpK = mySettings.getD("fracpK_"+t);

  // Create any RooRealVars you'll need
  _intVars["fracpK"] = 0;

  // ------------------------------------
  // Read in the KEYS pdfs
  KeysPdf* keys_pK  = new KeysPdf(_mB, m,p+"_pK",c,t,a,fileName);
  KeysPdf* keys_pPi = new KeysPdf(_mB, m,p+"_pPi",c,t,a,fileName);

  pdf_pK  = keys_pK->getPdf();
  pdf_pPi = keys_pPi->getPdf();

  // ------------------------------------
  if(!pdf_pK || !pdf_pPi ) std::cout << "Problem obtaining Lambda Pdfs!" << std::endl;

}

RooAbsPdf* Lambda::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["fracpK"]==0)	_intVars["fracpK"]		 	= new RooRealVar(Form("%s_Var_frac1",_name.c_str()),"",_fracpK);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " Lambda (fracpK): " << _intVars["fracpK"]->getVal() 
              << std::endl;
  }

  // Float the ratio between helamp 010 and 100/001
  RooAddPdf *pdf = new RooAddPdf(_name.c_str(),"",*pdf_pK,*pdf_pPi,*_intVars["fracpK"]);
  return pdf;
}
