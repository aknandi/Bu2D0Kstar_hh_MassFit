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
	// Problem with so many parameters
	// Which ones are the same. Different for LL and DD
	/////////////////////
  // Endpoints are allowed to vary within range +-TOLERANCE
  double TOLERANCE = 0;
  /////////////////////

  // Kinematic endpoints /////////////////
  double smearshift = 0; // due to MC-data difference and fitting to smeared MC
  double Min_Bu_DstKst_D0gamma=4920+smearshift;
  double Max_Bu_DstKst_D0gamma=5220+smearshift;
  double Min_Bu_DstKst_D0pi0  =5020+smearshift;
  double Max_Bu_DstKst_D0pi0  =5120+smearshift;

  double Min_Bd_DstKst_D0pi0  =5020.901+smearshift;
  double Max_Bd_DstKst_D0pi0  =5117.039+smearshift;

  // Kinematic endpoints /////////////////
  a_Bu_DstKst_D0pi0  = new RooRealVar("low_a_Bu_DstKst_D0pi0"  ,"a_Bu_DstKst_D0pi0"  , Min_Bu_DstKst_D0pi0  , Min_Bu_DstKst_D0pi0  -TOLERANCE, Min_Bu_DstKst_D0pi0  +TOLERANCE);
  b_Bu_DstKst_D0pi0  = new RooRealVar("low_b_Bu_DstKst_D0pi0"  ,"b_Bu_DstKst_D0pi0"  , Max_Bu_DstKst_D0pi0  , Max_Bu_DstKst_D0pi0  -TOLERANCE, Max_Bu_DstKst_D0pi0  +TOLERANCE);
  a_Bu_DstKst_D0gamma= new RooRealVar("low_a_Bu_DstKst_D0gamma","a_Bu_DstKst_D0gamma", Min_Bu_DstKst_D0gamma, Min_Bu_DstKst_D0gamma-TOLERANCE, Min_Bu_DstKst_D0gamma+TOLERANCE);
  b_Bu_DstKst_D0gamma= new RooRealVar("low_b_Bu_DstKst_D0gamma","b_Bu_DstKst_D0gamma", Max_Bu_DstKst_D0gamma, Max_Bu_DstKst_D0gamma-TOLERANCE, Max_Bu_DstKst_D0gamma+TOLERANCE);

  a_Bd_DstKst_D0pi0  = new RooRealVar("low_a_Bd_Dst0Kst0_D0pi0"  ,"a_Bd_DstKst_D0pi0"  , Min_Bd_DstKst_D0pi0  , Min_Bd_DstKst_D0pi0  -TOLERANCE, Min_Bd_DstKst_D0pi0  +TOLERANCE);
  b_Bd_DstKst_D0pi0  = new RooRealVar("low_b_Bd_Dst0Kst0_D0pi0"  ,"b_Bd_DstKst_D0pi0"  , Max_Bd_DstKst_D0pi0  , Max_Bd_DstKst_D0pi0  -TOLERANCE, Max_Bd_DstKst_D0pi0  +TOLERANCE);

  // Efficiency effects: fixed value based on MC
  // Initial values from individual fits
  csi_HORNS = new RooRealVar("low_csi_HORNS" ,"",             0.8);
  csi_LITTLEHORNS = new RooRealVar("low_csi_LITTLEHORNS" ,"", 0.4);
  csi_HILL = new RooRealVar("low_csi_HILL" ,"",               0.8);
  csi_HILLg = new RooRealVar("low_csi_HILLg" ,"",             0.6);

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Global shift and width parameters
  low_shift       = new RooRealVar("low_shift"         ,""  ,   0.); 
  low_sigma       = new RooRealVar("low_sigma"         ,""  ,  12.); // Base gaussian width

  // Double Gaussian parameters
  low_shiftg      = new RooRealVar("low_shiftg"        ,""  , -50.0); // shift of secondary Gaussian
  f_MC            = new RooRealVar("low_f_MC"          ,""  ,   0.93); // fraction of secondary Gaussian
  ratio_sigma_MC  = new RooRealVar("low_ratio_sigma_MC",""  ,   5.0); // low_sigma * ratio = width of secondary Gaussian
 
  //-----------------------------------------------------------------------------------------------------------------------------------
  //The next block of fixing fixes to the Kπ result and can be toggled.
  //For use with low-statistics modes
  if(fixed) fix();

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bu
  pdf_DstKst_D0pi0_010["bu"]   = new RooHORNSdini("Bu_DstKst_D0pi0_010"         ,"",*x,*a_Bu_DstKst_D0pi0,  *b_Bu_DstKst_D0pi0,  *csi_HORNS,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_DstKst_D0gamma_010["bu"] = new RooHILLdini("Bu_DstKst_D0gamma_010"        ,"",*x,*a_Bu_DstKst_D0gamma,*b_Bu_DstKst_D0gamma,*csi_HILLg,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_DstKst_D0pi0_101["bu"]   = new RooHILLdini("Bu_DstKst_D0pi0_101"          ,"",*x,*a_Bu_DstKst_D0pi0,  *b_Bu_DstKst_D0pi0,  *csi_HILL,       *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_DstKst_D0gamma_101["bu"] = new RooLITTLEHORNSdini("Bu_DstKst_D0gamma_101" ,"",*x,*a_Bu_DstKst_D0gamma,*b_Bu_DstKst_D0gamma,*csi_LITTLEHORNS,*low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bd
  pdf_DstKst_D0pi0_010["bd"]   = new RooHORNSdini("Bd_DstKst_D0pi0_010"         ,"",*x,*a_Bd_DstKst_D0pi0,  *b_Bd_DstKst_D0pi0,  *csi_HORNS,      *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);
  pdf_DstKst_D0pi0_101["bd"]   = new RooHILLdini("Bd_DstKst_D0pi0_101"          ,"",*x,*a_Bd_DstKst_D0pi0,  *b_Bd_DstKst_D0pi0,  *csi_HILL,       *low_shift,*low_sigma,*ratio_sigma_MC,*f_MC,*low_shiftg);

  //-----------------------------------------------------------------------------------------------------------------------------------


}

void PartRecoShapes::fix()
{
  a_Bu_DstKst_D0pi0->setConstant();
  b_Bu_DstKst_D0pi0->setConstant();
  a_Bu_DstKst_D0gamma->setConstant();
  b_Bu_DstKst_D0gamma->setConstant();

  a_Bd_DstKst_D0pi0->setConstant();
  b_Bd_DstKst_D0pi0->setConstant();

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

