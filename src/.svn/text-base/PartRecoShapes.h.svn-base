#ifndef PartRecoShapes_h
#define PartRecoShapes_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "Base.h"

class RooRealVar;
class RooGaussian;

class PartRecoShapes : public Base {
public :
  PartRecoShapes(RooRealVar*,bool=false,std::string=""); 
  ~PartRecoShapes(){}

  std::map<std::string,RooAbsPdf*> pdf_Dst0Kst0_D0pi0_010;
  std::map<std::string,RooAbsPdf*> pdf_Dst0Kst0_D0pi0_001;
  std::map<std::string,RooAbsPdf*> pdf_Dst0Kst0_D0gamma_010;
  std::map<std::string,RooAbsPdf*> pdf_Dst0Kst0_D0gamma_001;

  void fix();
protected :
  void buildShapes(RooRealVar*);

private:
  bool fixed;

  std::string t; // FC

  double lowmass_MASScut;

  RooRealVar* a_Bd_Dst0Kst0_D0pi0;
  RooRealVar* a_Bd_Dst0Kst0_D0gamma;
  RooRealVar* b_Bd_Dst0Kst0_D0pi0;
  RooRealVar* b_Bd_Dst0Kst0_D0gamma;
  RooRealVar* a_Bs_Dst0Kst0_D0pi0;
  RooRealVar* a_Bs_Dst0Kst0_D0gamma;
  RooRealVar* b_Bs_Dst0Kst0_D0pi0;
  RooRealVar* b_Bs_Dst0Kst0_D0gamma;

  RooRealVar* csi_HORNS;
  RooRealVar* csi_LITTLEHORNS;
  RooRealVar* csi_HILL;
  RooRealVar* csi_HILLg;

  RooRealVar* f_MC;
  RooRealVar* ratio_sigma_MC;
  RooRealVar* alpha;
  RooRealVar* n;

  RooRealVar* low_shift;
  RooRealVar* low_shiftg;
  RooRealVar* low_sigma;

};


#endif
