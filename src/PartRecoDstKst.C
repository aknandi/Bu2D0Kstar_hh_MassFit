#include "RooRealVar.h"
#include "RooAbsPdf.h"

#include "Settings.h"
#include "PartRecoDstKst.h"
#include "PartRecoShapes.h"
PartRecoDstKst::PartRecoDstKst(RooRealVar* pmB, std::string m, std::string c, std::string t, std::string a, std::string fileName)
{
  _mB=pmB;
  _name="PartRecoDstKst_"+m+"_"+c+"_"+t+"_"+a;
	
  // Make a settings object to find the relevant parameters
  Settings mySettings = Settings("PartRecoDstKstSettings_"+m+"_"+c+"_"+t+"_"+a);
  mySettings.readPairStringsToMap(fileName);

  Settings genConfs = Settings("GenSettings");
  genConfs.readPairStringsToMap("Settings/GeneralSettings.txt");
  string fitlow = genConfs.get("fit_limit_low");

  // Read in the relevant parameters 
  //_frac010 = mySettings.getD("frac010_"+p+"_"+fitlow);

  // Create any RooRealVars you'll need
  _intVars["frac010_bu"] = 0;
  _intVars["frac010_bd"] = 0;
  _intVars["frac_bd"] = 0;

  // Read in PARAMETRIC low mass pdf
  PartRecoShapes* prs = new PartRecoShapes(_mB, true, t);
  pdf_bu_g_101  = prs->pdf_Bu_DstKst_D0gamma_101;
  pdf_bu_g_010  = prs->pdf_Bu_DstKst_D0gamma_010;
  pdf_bu_pi_101 = prs->pdf_Bu_DstKst_D0pi0_101;
  pdf_bu_pi_010 = prs->pdf_Bu_DstKst_D0pi0_010;

  pdf_bd_pi_101 = prs->pdf_Bd_DstKst_D0pi0_101;
  pdf_bd_pi_010 = prs->pdf_Bd_DstKst_D0pi0_010;
  // ------------------------------------

  if(!pdf_bu_g_101 || !pdf_bu_g_010 || !pdf_bu_pi_101 || !pdf_bu_pi_010 || !pdf_bd_pi_101 || !pdf_bd_pi_010) {
    std::cout << "Problem obtaining Low Mass Pdfs!" << std::endl;
    exit(1);
  }

  // Bu shapes
  double G_101 = mySettings.getD("gamma_frac101_"+t);
  double G_010 = mySettings.getD("gamma_frac010_"+t);

  std::cout << "PartRecoDstKst " << t << std::endl;
  std::cout << " G_101 " << G_101 << "\n" << " G_010 " << G_010 << std::endl;
  var_G_101 = new RooRealVar(("bu_"+t+"_gamma_frac101").c_str(),"",G_101);
  var_G_010 = new RooRealVar(("bu_"+t+"_gamma_frac010").c_str(),"",G_010);
  //var_G_001->setConstant(kTRUE);
  //var_G_010->setConstant(kTRUE);

  pdf_bu_101 = new RooAddPdf((_name+"bu_101").c_str(),"",*pdf_bu_g_101,*pdf_bu_pi_101,*var_G_101);
  pdf_bu_010 = new RooAddPdf((_name+"bu_010").c_str(),"",*pdf_bu_g_010,*pdf_bu_pi_010,*var_G_010);


}

void PartRecoDstKst::setFractionBu010(RooAbsReal* newFractionBu010)
{
	setRelation("frac010_bu",newFractionBu010);
}

void PartRecoDstKst::setFractionBd010(RooAbsReal* newFractionBd010)
{
	setRelation("frac010_bd",newFractionBd010);
}

void PartRecoDstKst::setFractionBd(RooAbsReal* newFractionBd)
{
	setRelation("frac_bd",newFractionBd);
}

RooAbsPdf* PartRecoDstKst::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["frac010_bu"]==0)	_intVars["frac010_bu"]		 	= new RooRealVar(Form("%s_Var_frac010_bu",_name.c_str()),"",_frac010_bu);
  if(_intVars["frac010_bd"]==0)	_intVars["frac010_bd"]		 	= new RooRealVar(Form("%s_Var_frac010_bd",_name.c_str()),"",_frac010_bd);
  if(_intVars["frac_bd"]==0)		_intVars["frac_bd"]		 	= new RooRealVar(Form("%s_Var_frac_bd",_name.c_str()),"",_frac_bd);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " PartRecoDstKst (frac010_bu): " << _intVars["frac010_bu"]->getVal()
			  << " PartRecoDstKst (frac010_bd): " << _intVars["frac010_bd"]->getVal()
			  << " PartRecoDstKst (frac_bd): " << _intVars["frac_bd"]->getVal()
              << std::endl;
  }

  // Float the ratio between helamp 010 and 101 in Bu
  pdf_bu = new RooAddPdf((_name+"bu").c_str(),"",*pdf_bu_010,*pdf_bu_101,*_intVars["frac010_bu"]);
  pdf_bd = new RooAddPdf((_name+"bd").c_str(),"",*pdf_bd_pi_010,*pdf_bd_pi_101,*_intVars["frac010_bd"]);

  // Combine Bd and Bu and float ratio
  RooAddPdf *pdf_partreco = new RooAddPdf(_name.c_str(),"",*pdf_bd,*pdf_bu,*_intVars["frac_bd"]);
  return pdf_partreco;
}
