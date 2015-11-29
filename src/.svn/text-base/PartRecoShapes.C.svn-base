#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"

#include "PartRecoShapes.h"
#include "RooHORNSdini.h"
#include "RooHILLdini.h"
#include "RooLITTLEHORNSdini.h"

PartRecoShapes::PartRecoShapes(RooRealVar *x,bool f, std::string _t)
: Base()
, fixed(f)
, t(_t) // Ks Type 
, lowmass_MASScut(x->getMin())
{ 
  buildShapes(x);

//  modeToUseFor[d2kspipi] = d2kspipi;
//  modeToUseFor[d2kskk]   = d2kspipi;
  
}

void PartRecoShapes::buildShapes(RooRealVar *x)
{
  /////////////////////
  // Endpoints are allowed to vary within range +-TOLERANCE
  double TOLERANCE = 5;
  /////////////////////

  // Kinematic endpoints /////////////////
  double smearshift = 3.1; // due to MC-data difference and fitting to smeared MC
  double Min_Bs_Dst0Kst0_D0gamma=4999.396+smearshift;
  double Max_Bs_Dst0Kst0_D0gamma=5313.523+smearshift;
  double Min_Bs_Dst0Kst0_D0pi0  =5104.125+smearshift;
  double Max_Bs_Dst0Kst0_D0pi0  =5202.589+smearshift;

  double Min_Bd_Dst0Kst0_D0pi0  =5020.901+smearshift;
  double Max_Bd_Dst0Kst0_D0pi0  =5117.039+smearshift;
  double Min_Bd_Dst0Kst0_D0gamma=4918.666+smearshift;
  double Max_Bd_Dst0Kst0_D0gamma=5225.373+smearshift;

  // Kinematic endpoints /////////////////
  a_Bs_Dst0Kst0_D0pi0  = new RooRealVar("low_a_Bs_Dst0Kst0_D0pi0"  ,"a_Bs_Dst0Kst0_D0pi0"  , Min_Bs_Dst0Kst0_D0pi0  , Min_Bs_Dst0Kst0_D0pi0  -TOLERANCE, Min_Bs_Dst0Kst0_D0pi0  +TOLERANCE);
  b_Bs_Dst0Kst0_D0pi0  = new RooRealVar("low_b_Bs_Dst0Kst0_D0pi0"  ,"b_Bs_Dst0Kst0_D0pi0"  , Max_Bs_Dst0Kst0_D0pi0  , Max_Bs_Dst0Kst0_D0pi0  -TOLERANCE, Max_Bs_Dst0Kst0_D0pi0  +TOLERANCE);
  a_Bs_Dst0Kst0_D0gamma= new RooRealVar("low_a_Bs_Dst0Kst0_D0gamma","a_Bs_Dst0Kst0_D0gamma", Min_Bs_Dst0Kst0_D0gamma, Min_Bs_Dst0Kst0_D0gamma-TOLERANCE, Min_Bs_Dst0Kst0_D0gamma+TOLERANCE);
  b_Bs_Dst0Kst0_D0gamma= new RooRealVar("low_b_Bs_Dst0Kst0_D0gamma","b_Bs_Dst0Kst0_D0gamma", Max_Bs_Dst0Kst0_D0gamma, Max_Bs_Dst0Kst0_D0gamma-TOLERANCE, Max_Bs_Dst0Kst0_D0gamma+TOLERANCE);

  a_Bd_Dst0Kst0_D0pi0  = new RooRealVar("low_a_Bd_Dst0Kst0_D0pi0"  ,"a_Bd_Dst0Kst0_D0pi0"  , Min_Bd_Dst0Kst0_D0pi0  , Min_Bd_Dst0Kst0_D0pi0  -TOLERANCE, Min_Bd_Dst0Kst0_D0pi0  +TOLERANCE);
  b_Bd_Dst0Kst0_D0pi0  = new RooRealVar("low_b_Bd_Dst0Kst0_D0pi0"  ,"b_Bd_Dst0Kst0_D0pi0"  , Max_Bd_Dst0Kst0_D0pi0  , Max_Bd_Dst0Kst0_D0pi0  -TOLERANCE, Max_Bd_Dst0Kst0_D0pi0  +TOLERANCE);
  a_Bd_Dst0Kst0_D0gamma= new RooRealVar("low_a_Bd_Dst0Kst0_D0gamma","a_Bd_Dst0Kst0_D0gamma", Min_Bd_Dst0Kst0_D0gamma, Min_Bd_Dst0Kst0_D0gamma-TOLERANCE, Min_Bd_Dst0Kst0_D0gamma+TOLERANCE);
  b_Bd_Dst0Kst0_D0gamma= new RooRealVar("low_b_Bd_Dst0Kst0_D0gamma","b_Bd_Dst0Kst0_D0gamma", Max_Bd_Dst0Kst0_D0gamma, Max_Bd_Dst0Kst0_D0gamma-TOLERANCE, Max_Bd_Dst0Kst0_D0gamma+TOLERANCE);

  // Efficiency effects: fixed value based on MC
  // Initial values from individual fits
  csi_HORNS = new RooRealVar("low_csi_HORNS" ,"",             0.827, 0., 1.0); 
  csi_LITTLEHORNS = new RooRealVar("low_csi_LITTLEHORNS" ,"", 0.569, 0., 1.0); 
  csi_HILL = new RooRealVar("low_csi_HILL" ,"",               0.821, 0., 1.0); 
  csi_HILLg = new RooRealVar("low_csi_HILLg" ,"",             0.554, 0., 1.0); 

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Global shift and width parameters
  low_shift       = new RooRealVar("low_shift"         ,""  ,   0.); 
  low_sigma       = new RooRealVar("low_sigma"         ,""  ,  13.93,    5., 25.); // Base gaussian width

  // Double Gaussian parameters
  low_shiftg      = new RooRealVar("low_shiftg"        ,""  , -35.10, -100.,  0.); // shift of secondary Gaussian
  f_MC            = new RooRealVar("low_f_MC"          ,""  ,   0.928,   0.,  1.); // fraction of secondary Gaussian
  ratio_sigma_MC  = new RooRealVar("low_ratio_sigma_MC",""  ,   6.26,    0., 10.); // low_sigma * ratio = width of secondary Gaussian
 
  //-----------------------------------------------------------------------------------------------------------------------------------
  //The next block of fixing fixes to the Kπ result and can be toggled.
  //For use with low-statistics modes
  if(fixed) fix();

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bd    
  pdf_Dst0Kst0_D0pi0_010["bd"]   = new RooHORNSdini("Bd_Dst0Kst0_D0pi0_010"         ,"",*x,*a_Bd_Dst0Kst0_D0pi0,  *b_Bd_Dst0Kst0_D0pi0,  *csi_HORNS,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0gamma_010["bd"] = new RooHILLdini("bd_Dst0Kst0_D0gamma_010"        ,"",*x,*a_Bd_Dst0Kst0_D0gamma,*b_Bd_Dst0Kst0_D0gamma,*csi_HILLg,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0pi0_001["bd"]   = new RooHILLdini("Bd_Dst0Kst0_D0pi0_001"          ,"",*x,*a_Bd_Dst0Kst0_D0pi0,  *b_Bd_Dst0Kst0_D0pi0,  *csi_HILL,       *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0gamma_001["bd"] = new RooLITTLEHORNSdini("Bd_Dst0Kst0_D0gamma_001" ,"",*x,*a_Bd_Dst0Kst0_D0gamma,*b_Bd_Dst0Kst0_D0gamma,*csi_LITTLEHORNS,*low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bs 
  pdf_Dst0Kst0_D0pi0_010["bs"]   = new RooHORNSdini("Bs_Dst0Kst0_D0pi0_010"         ,"",*x,*a_Bs_Dst0Kst0_D0pi0,  *b_Bs_Dst0Kst0_D0pi0,  *csi_HORNS,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0gamma_010["bs"] = new RooHILLdini("Bs_Dst0Kst0_D0gamma_010"        ,"",*x,*a_Bs_Dst0Kst0_D0gamma,*b_Bs_Dst0Kst0_D0gamma,*csi_HILLg,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0pi0_001["bs"]   = new RooHILLdini("Bs_Dst0Kst0_D0pi0_001"          ,"",*x,*a_Bs_Dst0Kst0_D0pi0,  *b_Bs_Dst0Kst0_D0pi0,  *csi_HILL,       *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_Dst0Kst0_D0gamma_001["bs"] = new RooLITTLEHORNSdini("Bs_Dst0Kst0_D0gamma_001" ,"",*x,*a_Bs_Dst0Kst0_D0gamma,*b_Bs_Dst0Kst0_D0gamma,*csi_LITTLEHORNS,*low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);

  //-----------------------------------------------------------------------------------------------------------------------------------


}

void PartRecoShapes::fix()
{
  a_Bs_Dst0Kst0_D0pi0->setConstant();
  b_Bs_Dst0Kst0_D0pi0->setConstant();
  a_Bs_Dst0Kst0_D0gamma->setConstant();
  b_Bs_Dst0Kst0_D0gamma->setConstant();

  a_Bd_Dst0Kst0_D0pi0->setConstant();
  b_Bd_Dst0Kst0_D0pi0->setConstant();
  a_Bd_Dst0Kst0_D0gamma->setConstant();
  b_Bd_Dst0Kst0_D0gamma->setConstant();

  csi_HORNS->setConstant();
  csi_LITTLEHORNS->setConstant();
  csi_HILL->setConstant();
  csi_HILLg->setConstant();

  low_shift->setConstant();
  low_sigma->setConstant();
  low_shiftg->setConstant();
  ratio_sigma_MC->setConstant(); 
  f_MC->setConstant(); 

}

