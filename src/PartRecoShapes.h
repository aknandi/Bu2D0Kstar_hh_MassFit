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
  PartRecoShapes(RooRealVar*,bool=false,std::string="",bool=true);
  ~PartRecoShapes(){}

  RooAbsPdf* pdf_Bu_DstKst_D0pi0_010;
  RooAbsPdf* pdf_Bu_DstKst_D0pi0_101;
  RooAbsPdf* pdf_Bu_DstKst_D0gamma_010;
  RooAbsPdf* pdf_Bu_DstKst_D0gamma_101;

  RooAbsPdf* pdf_Bd_DstKst_D0pi0_010;
  RooAbsPdf* pdf_Bd_DstKst_D0pi0_101;

  void fix();
  void deleteVariables();
protected :
  void buildShapes(RooRealVar*);

private:
  bool fixed;
  bool partrecosystematic;

  std::string t; // FC

  double lowmass_MASScut;

  RooRealVar* a_Bu_DstKst_D0pi0;
  RooRealVar* a_Bu_DstKst_D0gamma;
  RooRealVar* b_Bu_DstKst_D0pi0;
  RooRealVar* b_Bu_DstKst_D0gamma;
  RooRealVar* a_Bd_DstKst_D0pi0;
  RooRealVar* b_Bd_DstKst_D0pi0;

  RooRealVar* csi_bu_p010;
  RooRealVar* csi_bu_g010;
  RooRealVar* csi_bu_p101;
  RooRealVar* csi_bu_g101;
  RooRealVar* csi_bd_p010;
  RooRealVar* csi_bd_p101;

  RooRealVar* shift_bu_p010;
  RooRealVar* shift_bu_g010;
  RooRealVar* shift_bu_p101;
  RooRealVar* shift_bu_g101;
  RooRealVar* shift_bd_p010;
  RooRealVar* shift_bd_p101;

  RooRealVar* sigma_bu_p010;
  RooRealVar* sigma_bu_g010;
  RooRealVar* sigma_bu_p101;
  RooRealVar* sigma_bu_g101;
  RooRealVar* sigma_bd_p010;
  RooRealVar* sigma_bd_p101;

  RooRealVar* ratio_sigma_bu_p010;
  RooRealVar* ratio_sigma_bu_g010;
  RooRealVar* ratio_sigma_bu_p101;
  RooRealVar* ratio_sigma_bu_g101;
  RooRealVar* ratio_sigma_bd_p010;
  RooRealVar* ratio_sigma_bd_p101;

  RooRealVar* fraction_sigma_bu_p010;
  RooRealVar* fraction_sigma_bu_g010;
  RooRealVar* fraction_sigma_bu_p101;
  RooRealVar* fraction_sigma_bu_g101;
  RooRealVar* fraction_sigma_bd_p010;
  RooRealVar* fraction_sigma_bd_p101;

  RooRealVar* shiftg_bu_p010;
  RooRealVar* shiftg_bu_g010;
  RooRealVar* shiftg_bu_p101;
  RooRealVar* shiftg_bu_g101;
  RooRealVar* shiftg_bd_p010;
  RooRealVar* shiftg_bd_p101;

};


#endif
