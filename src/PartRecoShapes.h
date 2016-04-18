#ifndef PartRecoShapes_h
#define PartRecoShapes_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "Base.h"

class RooRealVar;
class RooGaussian;

// todo Does this need to inherit from Base?
class PartRecoShapes : public Base {
public :
  PartRecoShapes(RooRealVar*,bool=false,std::string=""); 
  ~PartRecoShapes(){}

  std::map<std::string,RooAbsPdf*> pdf_DstKst_D0pi0_010;
  std::map<std::string,RooAbsPdf*> pdf_DstKst_D0pi0_101;
  std::map<std::string,RooAbsPdf*> pdf_DstKst_D0gamma_010;
  std::map<std::string,RooAbsPdf*> pdf_DstKst_D0gamma_101;

  void fix();
protected :
  void buildShapes(RooRealVar*);

private:
  bool fixed;

  std::string t; // FC

  double lowmass_MASScut;

  RooRealVar* a_Bu_DstKst_D0pi0;
  RooRealVar* a_Bu_DstKst_D0gamma;
  RooRealVar* b_Bu_DstKst_D0pi0;
  RooRealVar* b_Bu_DstKst_D0gamma;

  RooRealVar* a_Bd_DstKst_D0pi0;
  RooRealVar* b_Bd_DstKst_D0pi0;

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
