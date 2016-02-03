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

  // ------------------------------------
  // Read in the KEYS low mass pdf -- comment this out if using parametric shapes
  /*
  KeysPdf* keys_g_001 = new KeysPdf(_mB, m,p+"_dgam_001",c,t,a,fileName);
  KeysPdf* keys_g_010 = new KeysPdf(_mB, m,p+"_dgam_010",c,t,a,fileName);
  KeysPdf* keys_pi_001 = new KeysPdf(_mB, m,p+"_dpi0_001",c,t,a,fileName);
  KeysPdf* keys_pi_010 = new KeysPdf(_mB, m,p+"_dpi0_010",c,t,a,fileName);

  pdf_g_001 = keys_g_001->getPdf();
  pdf_g_010 = keys_g_010->getPdf();
  pdf_pi_001 = keys_pi_001->getPdf();
  pdf_pi_010 = keys_pi_010->getPdf();
  */
  // ------------------------------------
  // Read in PARAMETRIC low mass pdf -- comment this out if using keys pdfs
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

  // Fix the ratio between the gamma and pi0 comments from MC, separately for helamps
  const double g_001 = mySettings.getD("BR_Dg")*mySettings.getD("effGen_Dg_001")*mySettings.getD("effSel_Dg_001_"+t)*mySettings.getD("integfrac_"+p+"_g001_"+fitlow);
  const double g_010 = mySettings.getD("BR_Dg")*mySettings.getD("effGen_Dg_010")*mySettings.getD("effSel_Dg_010_"+t)*mySettings.getD("integfrac_"+p+"_g010_"+fitlow);
  const double p_001 = mySettings.getD("BR_Dpi")*mySettings.getD("effGen_Dpi_001")*mySettings.getD("effSel_Dpi_001_"+t)*mySettings.getD("integfrac_"+p+"_pi001_"+fitlow);
  const double p_010 = mySettings.getD("BR_Dpi")*mySettings.getD("effGen_Dpi_010")*mySettings.getD("effSel_Dpi_010_"+t)*mySettings.getD("integfrac_"+p+"_pi010_"+fitlow);

  //std::cout << "g_001 " << g_001 << std::endl;
  //std::cout << "p_001 " << p_001 << std::endl;
  //std::cout << "g_010 " << g_010 << std::endl;
  //std::cout << "p_010 " << p_010 << std::endl;

  double G_001 = g_001 / (g_001+p_001);
  double G_010 = g_010 / (g_010+p_010);

  // Hard code Alexis' numbers (4900 MeV)
  //if(t=="LL"){ G_001=0.37434; G_010=0.38006; }
  //if(t=="DD"){ G_001=0.37633; G_010=0.38107; }

  std::cout << "PartRecoDstKst " << p << " " << t << std::endl;
  std::cout << " G_001 " << G_001 << "\n" << " G_010 " << G_010 << std::endl;
  var_G_001 = new RooRealVar((p+"_"+t+"_gamma_frac001").c_str(),"",G_001);
  var_G_010 = new RooRealVar((p+"_"+t+"_gamma_frac010").c_str(),"",G_010);
  //var_G_001->setConstant(kTRUE);
  //var_G_010->setConstant(kTRUE);

  pdf_001 = new RooAddPdf((_name+"_001").c_str(),"",*pdf_g_001,*pdf_pi_001,*var_G_001);
  pdf_010 = new RooAddPdf((_name+"_010").c_str(),"",*pdf_g_010,*pdf_pi_010,*var_G_010);

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
