#include "RooRealVar.h"
#include "RooAbsPdf.h"

#include "Settings.h"
#include "PartRecoDstKst.h"
#include "PartRecoShapes.h"
PartRecoDstKst::PartRecoDstKst(RooRealVar* pmB, std::string m,std::string p, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="PartRecoDstKst_"+m+"_"+p+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("PartRecoDstKstSettings_"+m+"_"+p+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  Settings genConfs = Settings("GenSettings");
  genConfs.readPairStringsToMap("Settings/GeneralSettings.txt");
  string fitlow = genConfs.get("fit_limit_low");

  // Read in the relevant parameters 
  //_frac010 = mySettings.getD("frac010_"+p+"_"+fitlow);

  // Create any RooRealVars you'll need
  _intVars["frac010"] = 0;

  // Read in PARAMETRIC low mass pdf
  PartRecoShapes* prs = new PartRecoShapes(_mB, true, t);
  pdf_g_001  = prs->pdf_DstKst_D0gamma_101[p];
  pdf_g_010  = prs->pdf_DstKst_D0gamma_010[p];
  pdf_pi_001 = prs->pdf_DstKst_D0pi0_101[p];
  pdf_pi_010 = prs->pdf_DstKst_D0pi0_010[p];
  // ------------------------------------
  
  if(!pdf_g_001 || !pdf_g_010 || !pdf_pi_001 || !pdf_pi_010) {
    std::cout << "Problem obtaining " << p << " Low Mass Pdfs!" << std::endl;
    exit(1);
  }

  double G_001 = mySettings.getD("gamma_frac101_"+t);
  double G_010 = mySettings.getD("gamma_frac010_"+t);

  std::cout << "PartRecoDstKst " << p << " " << t << std::endl;
  std::cout << " G_001 " << G_001 << "\n" << " G_010 " << G_010 << std::endl;
  var_G_001 = new RooRealVar((p+"_"+t+"_gamma_frac001").c_str(),"",G_001);
  var_G_010 = new RooRealVar((p+"_"+t+"_gamma_frac010").c_str(),"",G_010);
  //var_G_001->setConstant(kTRUE);
  //var_G_010->setConstant(kTRUE);

  pdf_001 = new RooAddPdf((_name+"_001").c_str(),"",*pdf_g_001,*pdf_pi_001,*var_G_001);
  pdf_010 = new RooAddPdf((_name+"_010").c_str(),"",*pdf_g_010,*pdf_pi_010,*var_G_010);

}

void PartRecoDstKst::setFraction(RooAbsReal* newFraction)
{
	setRelation("frac010",newFraction);
}

RooAbsPdf* PartRecoDstKst::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["frac010"]==0)	_intVars["frac010"]		 	= new RooRealVar(Form("%s_Var_frac010",_name.c_str()),"",_frac010);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " PartRecoDstKst (frac010): " << _intVars["frac010"]->getVal() 
              << std::endl;
  }

  // Float the ratio between helamp 010 and 100/001
  RooAddPdf *pdf_partreco = new RooAddPdf(_name.c_str(),"",*pdf_010,*pdf_001,*_intVars["frac010"]);
  return pdf_partreco;
}
