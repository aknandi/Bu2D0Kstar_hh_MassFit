#include "RooRealVar.h"
#include "TFile.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooEffProd.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"

#include "Settings.h"
#include "KeysPdf.h"
KeysPdf::KeysPdf(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="KeysPdf_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("LowMassSettings_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  // Read in the KEYS low mass pdf (old)
  std::string	_keysPdfFile=mySettings.get("lowmass_fileName_"+t);
  TFile* file = TFile::Open(_keysPdfFile.c_str());
  RooWorkspace *ws = (RooWorkspace*) file->Get((mySettings.get("ws_name_"+t)).c_str());
  if(!ws) {
    cout << " RooWorkspace in " << _keysPdfFile << " NOT FOUND!!!" << endl;
    exit(1);
  }

  std::string _pdf_name = mySettings.get(p+"_"+t); // p contains the option for pdf name

  //std::cout << "Getting " << _pdf_name << std::endl;

  RooAbsPdf* lowmassPdf_preClean = (RooAbsPdf*) ws->pdf(_pdf_name.c_str());
  if(!lowmassPdf_preClean) {
    cout << _pdf_name << " NOT FOUND!!!" << endl;
    exit(1);
  }
  //v_cutoff=new RooFormulaVar(Form("%s_cutoff",_name.c_str()),Form("(@0<%f||@0>%f)?0:1",_mB->getMin(),mySettings.getD(m+"_"+p+"_lowmass_cutoff") ),RooArgSet(*_mB));
  v_cutoff=new RooFormulaVar(Form("%s_cutoff",_name.c_str()),Form("(@0<%f||@0>%f)?0:1",_mB->getMin(),_mB->getMax() ),RooArgSet(*_mB));

  pdf_low = new RooEffProd(_name.c_str(),"",*lowmassPdf_preClean,*v_cutoff);
  //pdf_low = lowmassPdf_preClean;
  pdf_low->SetName(_name.c_str());
  file->Close();

}

RooAbsPdf* KeysPdf::getPdf()
{
  return pdf_low;
}
