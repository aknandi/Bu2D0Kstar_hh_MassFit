#include "RooRealVar.h"
#include "RooAbsPdf.h"

#include "Settings.h"
#include "PartRecoDstKst.h"
#include "PartRecoShapes.h"
PartRecoDstKst::PartRecoDstKst(RooRealVar* pmB, std::string m, std::string c, std::string t, std::string a, std::string fileName, bool pdfgen)
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
  _intVars["coef010"] = 0;
  _intVars["coef101"] = 0;

  // Read in PARAMETRIC low mass pdf
  PartRecoShapes* prs = new PartRecoShapes(_mB, true, t, pdfgen);
  pdf_bu_gamma_101  = prs->pdf_Bu_DstKst_D0gamma_101;
  pdf_bu_gamma_010  = prs->pdf_Bu_DstKst_D0gamma_010;
  pdf_bu_pi_101 = prs->pdf_Bu_DstKst_D0pi0_101;
  pdf_bu_pi_010 = prs->pdf_Bu_DstKst_D0pi0_010;

  pdf_bd_pi_101 = prs->pdf_Bd_DstKst_D0pi0_101;
  pdf_bd_pi_010 = prs->pdf_Bd_DstKst_D0pi0_010;
  // ------------------------------------

  if(!pdf_bu_gamma_101 || !pdf_bu_gamma_010 || !pdf_bu_pi_101 || !pdf_bu_pi_010 || !pdf_bd_pi_101 || !pdf_bd_pi_010) {
    std::cout << "Problem obtaining Low Mass Pdfs!" << std::endl;
    exit(1);
  }

  // Bu shapes
  double gamma_101 = mySettings.getD("gamma_frac101_"+t);
  double gamma_010 = mySettings.getD("gamma_frac010_"+t);
  double bd_101 = mySettings.getD("bd_frac101_"+t);
  double bd_010 = mySettings.getD("bd_frac010_"+t);

  // Calculate the fraction on each 010 and 101 pdf such that they normalise to 1
  // We know the ratio of gamma/pi and Bdpi/Bupi
  double totalFraction101 = 1.0 + gamma_101 + bd_101;
  double pi_101 = 1/totalFraction101;
  gamma_101 /= totalFraction101;
  bd_101 /= totalFraction101;
  double totalFraction010 = 1.0 + gamma_010 + bd_010;
  double pi_010 = 1/totalFraction010;
  gamma_010 /= totalFraction010;
  bd_010 /= totalFraction010;

  std::cout << "PartRecoDstKst " << t << std::endl;
  std::cout << " pi_101 " << pi_101 << "\n" << " pi_010 " << pi_010 << std::endl;
  std::cout << " gamma_101 " << gamma_101 << "\n" << " gamma_010 " << gamma_010 << std::endl;
  std::cout << " bd_101 " << bd_101 << "\n" << " bd_010 " << bd_010 << std::endl;

  frac_pi_101 = new RooRealVar(("bu_"+t+"_pi_frac101").c_str(),"",pi_101);
  frac_pi_010 = new RooRealVar(("bu_"+t+"_pi_frac010").c_str(),"",pi_010);
  frac_gamma_101 = new RooRealVar(("bu_"+t+"_gamma_frac101").c_str(),"",gamma_101);
  frac_gamma_010 = new RooRealVar(("bu_"+t+"_gamma_frac010").c_str(),"",gamma_010);
  frac_bd_101 = new RooRealVar((t+"_bd_frac101").c_str(),"",bd_101);
  frac_bd_010 = new RooRealVar((t+"_bd_frac010").c_str(),"",bd_010);

  RooArgSet pdfList101;
  RooArgSet pdfList010;
  RooArgSet yieldFractions101;
  RooArgSet yieldFractions010;

  pdfList010.add(*pdf_bu_pi_010);
  pdfList010.add(*pdf_bu_gamma_010);
  pdfList010.add(*pdf_bd_pi_010);
  yieldFractions010.add(*frac_pi_010);
  yieldFractions010.add(*frac_gamma_010);
  yieldFractions010.add(*frac_bd_010);

  pdfList101.add(*pdf_bu_pi_101);
  pdfList101.add(*pdf_bu_gamma_101);
  pdfList101.add(*pdf_bd_pi_101);
  yieldFractions101.add(*frac_pi_101);
  yieldFractions101.add(*frac_gamma_101);
  yieldFractions101.add(*frac_bd_101);

  pdf_010 = new RooAddPdf((_name+"010").c_str(),"",pdfList010,yieldFractions010);
  pdf_101 = new RooAddPdf((_name+"101").c_str(),"",pdfList101,yieldFractions101);

}

void PartRecoDstKst::setCoef010(RooAbsReal* newCoef010)
{
	setRelation("coef010",newCoef010);
}

void PartRecoDstKst::setCoef101(RooAbsReal* newCoef101)
{
	setRelation("coef101",newCoef101);
}


RooAbsPdf* PartRecoDstKst::getPdf()
{
  bool dbThis(false);
  // If any RooFit variables still have 0 pointers (ie have not been set from outside (e.g. sharing parameters (like mean) across multiple shapes) then
  // assume that the values should be those passed in the constructor
  if(_intVars["coef010"]==0)	_intVars["coef010"]		 	= new RooRealVar(Form("%s_Var_coef010",_name.c_str()),"",_coef010);
  if(_intVars["coef101"]==0)	_intVars["coef101"]		 	= new RooRealVar(Form("%s_Var_coef101",_name.c_str()),"",_coef101);

  if(dbThis) {
    std::cout << "For name: " << _name << endl;
    std::cout << " PartRecoDstKst (frac010): " << _intVars["frac010"]->getVal()
              << std::endl;
  }

  RooArgSet pdflist;
  RooArgSet coefficients;

  pdflist.add(*pdf_010);
  pdflist.add(*pdf_101);
  coefficients.add(*_intVars["coef010"]);
  coefficients.add(*_intVars["coef101"]);

  // Combine Bd and Bu and float ratio
  RooAddPdf *pdf_partreco = new RooAddPdf(_name.c_str(),"",pdflist,coefficients);
  return pdf_partreco;
}
