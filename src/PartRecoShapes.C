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
  double TOLERANCE = 0;
  /////////////////////

  // Kinematic endpoints /////////////////
  double smearshift = 0; // due to MC-data difference and fitting to smeared MC
  double Min_Bu_DstKst_D0gamma=4930+smearshift;
  double Max_Bu_DstKst_D0gamma=5235+smearshift;
  double Min_Bu_DstKst_D0pi0  =5020+smearshift;
  double Max_Bu_DstKst_D0pi0  =5120+smearshift;

  double Min_Bd_DstKst_D0pi0  =5020+smearshift;
  double Max_Bd_DstKst_D0pi0  =5110+smearshift;

  // Kinematic endpoints /////////////////
  a_Bu_DstKst_D0pi0  = new RooRealVar("low_a_Bu_DstKst_D0pi0"  ,"a_Bu_DstKst_D0pi0"  , Min_Bu_DstKst_D0pi0  , Min_Bu_DstKst_D0pi0  -TOLERANCE, Min_Bu_DstKst_D0pi0  +TOLERANCE);
  b_Bu_DstKst_D0pi0  = new RooRealVar("low_b_Bu_DstKst_D0pi0"  ,"b_Bu_DstKst_D0pi0"  , Max_Bu_DstKst_D0pi0  , Max_Bu_DstKst_D0pi0  -TOLERANCE, Max_Bu_DstKst_D0pi0  +TOLERANCE);
  a_Bu_DstKst_D0gamma= new RooRealVar("low_a_Bu_DstKst_D0gamma","a_Bu_DstKst_D0gamma", Min_Bu_DstKst_D0gamma, Min_Bu_DstKst_D0gamma-TOLERANCE, Min_Bu_DstKst_D0gamma+TOLERANCE);
  b_Bu_DstKst_D0gamma= new RooRealVar("low_b_Bu_DstKst_D0gamma","b_Bu_DstKst_D0gamma", Max_Bu_DstKst_D0gamma, Max_Bu_DstKst_D0gamma-TOLERANCE, Max_Bu_DstKst_D0gamma+TOLERANCE);

  a_Bd_DstKst_D0pi0  = new RooRealVar("low_a_Bd_Dst0Kst0_D0pi0"  ,"a_Bd_DstKst_D0pi0"  , Min_Bd_DstKst_D0pi0  , Min_Bd_DstKst_D0pi0  -TOLERANCE, Min_Bd_DstKst_D0pi0  +TOLERANCE);
  b_Bd_DstKst_D0pi0  = new RooRealVar("low_b_Bd_Dst0Kst0_D0pi0"  ,"b_Bd_DstKst_D0pi0"  , Max_Bd_DstKst_D0pi0  , Max_Bd_DstKst_D0pi0  -TOLERANCE, Max_Bd_DstKst_D0pi0  +TOLERANCE);

  // Setting up a LOT of parameter values that are fixed from MC
  double a_p010, b_p010, csi_p010, shift_p010, sigma_p010, ratio_sigma_p010, fraction_sigma_p010, shiftg_p010;
  double a_p101, b_p101, csi_p101, shift_p101, sigma_p101, ratio_sigma_p101, fraction_sigma_p101, shiftg_p101;
  double a_g010, b_g010, csi_g010, shift_g010, sigma_g010, ratio_sigma_g010, fraction_sigma_g010, shiftg_g010;
  double a_g101, b_g101, csi_g101, shift_g101, sigma_g101, ratio_sigma_g101, fraction_sigma_g101, shiftg_g101;

  double a_p010_bd, b_p010_bd, csi_p010_bd, shift_p010_bd, sigma_p010_bd, ratio_sigma_p010_bd, fraction_sigma_p010_bd, shiftg_p010_bd;
  double a_p101_bd, b_p101_bd, csi_p101_bd, shift_p101_bd, sigma_p101_bd, ratio_sigma_p101_bd, fraction_sigma_p101_bd, shiftg_p101_bd;

  if(t=="LL") {
    a_p010 = 5025; b_p010 = 5121; csi_p010 = 0.80; shift_p010 = -1.1; sigma_p010 = 11.3; ratio_sigma_p010 = 7.1; fraction_sigma_p010 = 0.96; shiftg_p010 = -55;
    a_p101 = 5012; b_p101 = 5110; csi_p101 = 0.43; shift_p101 = 8.4;  sigma_p101 = 11.7; ratio_sigma_p101 = 5.0; fraction_sigma_p101 = 0.91; shiftg_p101 = -90.0;
    a_g010 = 4908; b_g010 = 5217; csi_g010 = 0.85; shift_g010 = 18.4; sigma_g010 = 5.0;  ratio_sigma_g010 = 7.1; fraction_sigma_g010 = 0.19; shiftg_g010 = -95.8;
    a_g101 = 4930; b_g101 = 5236; csi_g101 = 0.63; shift_g101 = -8.0; sigma_g101 = 12.8; ratio_sigma_g101 = 5.6; fraction_sigma_g101 = 0.88; shiftg_g101 = -53;

    a_p010_bd = 5022; b_p010_bd = 5110; csi_p010_bd = 1.00; shift_p010_bd = 0.05; sigma_p010_bd = 12.4; ratio_sigma_p010_bd = 6.1; fraction_sigma_p010_bd = 0.96; shiftg_p010_bd = -89;
    a_p101_bd = 5017; b_p101_bd = 5104; csi_p101_bd = 0.86; shift_p101_bd = 5.1;  sigma_p101_bd = 14.7; ratio_sigma_p101_bd = 4.8; fraction_sigma_p101_bd = 0.97; shiftg_p101_bd = -97;
  }
  else if(t=="DD") {
    a_p010 = 5025; b_p010 = 5120; csi_p010 = 0.81;   shift_p010 = -1.2; sigma_p010 = 11.8; ratio_sigma_p010 = 7.3; fraction_sigma_p010 = 0.97; shiftg_p010 = -60;
    a_p101 = 5005; b_p101 = 5113; csi_p101 = 0.0015; shift_p101 = -2.7; sigma_p101 = 14.9; ratio_sigma_p101 = 5.1; fraction_sigma_p101 = 0.99; shiftg_p101 = -69;
    a_g010 = 4905; b_g010 = 5206; csi_g010 = 0.52;   shift_g010 = 19.8; sigma_g010 = 12.4; ratio_sigma_g010 = 3.5; fraction_sigma_g010 = 0.97; shiftg_g010 = -100;
    a_g101 = 4929; b_g101 = 5234; csi_g101 = 0.52;   shift_g101 = -7.8; sigma_g101 = 14.4; ratio_sigma_g101 = 5.0; fraction_sigma_g101 = 0.88; shiftg_g101 = -64;

    a_p010_bd = 5022; b_p010_bd = 5111; csi_p010_bd = 1.00; shift_p010_bd = 0.12; sigma_p010_bd = 11.8; ratio_sigma_p010_bd = 7.6; fraction_sigma_p010_bd = 0.96; shiftg_p010_bd = -51;
    a_p101_bd = 5014; b_p101_bd = 5103; csi_p101_bd = 0.45; shift_p101_bd = 5.0;  sigma_p101_bd = 11.6; ratio_sigma_p101_bd = 5.8; fraction_sigma_p101_bd = 0.92; shiftg_p101_bd = -87;
  }
  else if(t=="mix") {
    a_p010 = 5026; b_p010 = 5120; csi_p010 = 0.81;   shift_p010 = -1.6; sigma_p010 = 11.7; ratio_sigma_p010 = 7.3; fraction_sigma_p010 = 0.97; shiftg_p010 = -58;
    a_p101 = 5019; b_p101 = 5113; csi_p101 = 0.47;   shift_p101 = 1.7;  sigma_p101 = 12.5; ratio_sigma_p101 = 5.1; fraction_sigma_p101 = 0.92; shiftg_p101 = -88;
    a_g010 = 4905; b_g010 = 5211; csi_g010 = 0.63;   shift_g010 = 19.8; sigma_g010 = 8.41; ratio_sigma_g010 = 4.1; fraction_sigma_g010 = 0.97; shiftg_g010 = -100;
    a_g101 = 4930; b_g101 = 5236; csi_g101 = 0.58;   shift_g101 = -8.4; sigma_g101 = 13.6; ratio_sigma_g101 = 5.3; fraction_sigma_g101 = 0.88; shiftg_g101 = -58;

    a_p010_bd = 5022; b_p010_bd = 5111; csi_p010_bd = 1.00; shift_p010_bd = 0.12; sigma_p010_bd = 11.8; ratio_sigma_p010_bd = 7.6; fraction_sigma_p010_bd = 0.96; shiftg_p010_bd = -51;
    a_p101_bd = 5014; b_p101_bd = 5103; csi_p101_bd = 0.45; shift_p101_bd = 5.0;  sigma_p101_bd = 11.6; ratio_sigma_p101_bd = 5.8; fraction_sigma_p101_bd = 0.92; shiftg_p101_bd = -87;
  }

  csi_bu_p010 = new RooRealVar("csi_bu_pi010" ,"", 		csi_p010);
  csi_bu_g010 = new RooRealVar("csi_bu_gamma010" ,"", 	csi_g010);
  csi_bu_p101 = new RooRealVar("csi_bu_pi101" ,"", 		csi_p101);
  csi_bu_g101 = new RooRealVar("csi_bu_gamma101" ,"", 	csi_g101);
  csi_bd_p010 = new RooRealVar("csi_bd_pi010" ,"",		csi_p010_bd);
  csi_bd_p101 = new RooRealVar("csi_bd_pi101" ,"",		csi_p101_bd);

  shift_bu_p010 = new RooRealVar("shift_bu_pi010" ,"", 		shift_p010);
  shift_bu_g010 = new RooRealVar("shift_bu_gamma010" ,"", 	shift_g010);
  shift_bu_p101 = new RooRealVar("shift_bu_pi101" ,"", 		shift_p101);
  shift_bu_g101 = new RooRealVar("shift_bu_gamma101" ,"", 	shift_g101);
  shift_bd_p010 = new RooRealVar("shift_bd_pi010" ,"",		shift_p010_bd);
  shift_bd_p101 = new RooRealVar("shift_bd_pi101" ,"",		shift_p101_bd);

  sigma_bu_p010 = new RooRealVar("sigma_bu_pi010" ,"", 		sigma_p010);
  sigma_bu_g010 = new RooRealVar("sigma_bu_gamma010" ,"", 	sigma_g010);
  sigma_bu_p101 = new RooRealVar("sigma_bu_pi101" ,"", 		sigma_p101);
  sigma_bu_g101 = new RooRealVar("sigma_bu_gamma101" ,"", 	sigma_g101);
  sigma_bd_p010 = new RooRealVar("sigma_bd_pi010" ,"",		sigma_p010_bd);
  sigma_bd_p101 = new RooRealVar("sigma_bd_pi101" ,"",		sigma_p101_bd);

  ratio_sigma_bu_p010 = new RooRealVar("ratio_sigma_bu_pi010" ,"", 		ratio_sigma_p010);
  ratio_sigma_bu_g010 = new RooRealVar("ratio_sigma_bu_gamma010" ,"", 	ratio_sigma_g010);
  ratio_sigma_bu_p101 = new RooRealVar("ratio_sigma_bu_pi101" ,"", 		ratio_sigma_p101);
  ratio_sigma_bu_g101 = new RooRealVar("ratio_sigma_bu_gamma101" ,"", 	ratio_sigma_g101);
  ratio_sigma_bd_p010 = new RooRealVar("ratio_sigma_bd_pi010" ,"",		ratio_sigma_p010_bd);
  ratio_sigma_bd_p101 = new RooRealVar("ratio_sigma_bd_pi101" ,"",		ratio_sigma_p101_bd);

  fraction_sigma_bu_p010 = new RooRealVar("fraction_sigma_bu_pi010" ,"", 		fraction_sigma_p010);
  fraction_sigma_bu_g010 = new RooRealVar("fraction_sigma_bu_gamma010" ,"", 	fraction_sigma_g010);
  fraction_sigma_bu_p101 = new RooRealVar("fraction_sigma_bu_pi101" ,"", 		fraction_sigma_p101);
  fraction_sigma_bu_g101 = new RooRealVar("fraction_sigma_bu_gamma101" ,"", 	fraction_sigma_g101);
  fraction_sigma_bd_p010 = new RooRealVar("fraction_sigma_bd_pi010" ,"",		fraction_sigma_p010_bd);
  fraction_sigma_bd_p101 = new RooRealVar("fraction_sigma_bd_pi101" ,"",		fraction_sigma_p101_bd);

  shiftg_bu_p010 = new RooRealVar("shiftg_bu_pi010" ,"", 		shiftg_p010);
  shiftg_bu_g010 = new RooRealVar("shiftg_bu_gamma010" ,"", 	shiftg_g010);
  shiftg_bu_p101 = new RooRealVar("shiftg_bu_pi101" ,"", 		shiftg_p101);
  shiftg_bu_g101 = new RooRealVar("shiftg_bu_gamma101" ,"", 	shiftg_g101);
  shiftg_bd_p010 = new RooRealVar("shiftg_bd_pi010" ,"",		shiftg_p010_bd);
  shiftg_bd_p101 = new RooRealVar("shiftg_bd_pi101" ,"",		shiftg_p101_bd);

  //-----------------------------------------------------------------------------------------------------------------------------------
  //The next block of fixing fixes to the Kπ result and can be toggled.
  //For use with low-statistics modes
  if(fixed) fix();

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bu
  pdf_Bu_DstKst_D0pi0_010   = new RooHORNSdini("Bu_DstKst_D0pi0_010"         ,"",*x,*a_Bu_DstKst_D0pi0,  *b_Bu_DstKst_D0pi0,  *csi_bu_p010,*shift_bu_p010,*sigma_bu_p010,*ratio_sigma_bu_p010,*fraction_sigma_bu_p010,*shiftg_bu_p010);
  pdf_Bu_DstKst_D0gamma_010 = new RooHILLdini("Bu_DstKst_D0gamma_010"        ,"",*x,*a_Bu_DstKst_D0gamma,*b_Bu_DstKst_D0gamma,*csi_bu_g010,*shift_bu_g010,*sigma_bu_g010,*ratio_sigma_bu_g010,*fraction_sigma_bu_g010,*shiftg_bu_g010);
  pdf_Bu_DstKst_D0pi0_101   = new RooHILLdini("Bu_DstKst_D0pi0_101"          ,"",*x,*a_Bu_DstKst_D0pi0,  *b_Bu_DstKst_D0pi0,  *csi_bu_p101,*shift_bu_p101,*sigma_bu_p101,*ratio_sigma_bu_p101,*fraction_sigma_bu_p101,*shiftg_bu_p101);
  pdf_Bu_DstKst_D0gamma_101 = new RooLITTLEHORNSdini("Bu_DstKst_D0gamma_101" ,"",*x,*a_Bu_DstKst_D0gamma,*b_Bu_DstKst_D0gamma,*csi_bu_g101,*shift_bu_g101,*sigma_bu_g101,*ratio_sigma_bu_g101,*fraction_sigma_bu_g101,*shiftg_bu_g101);

  //-----------------------------------------------------------------------------------------------------------------------------------
  // Low Mass Shapes : Bd
  pdf_Bd_DstKst_D0pi0_010   = new RooHORNSdini("Bd_DstKst_D0pi0_010"         ,"",*x,*a_Bd_DstKst_D0pi0,  *b_Bd_DstKst_D0pi0,  *csi_bd_p010,*shift_bd_p010,*sigma_bd_p010,*ratio_sigma_bd_p010,*fraction_sigma_bd_p010,*shiftg_bd_p010);
  pdf_Bd_DstKst_D0pi0_101   = new RooHILLdini("Bd_DstKst_D0pi0_101"          ,"",*x,*a_Bd_DstKst_D0pi0,  *b_Bd_DstKst_D0pi0,  *csi_bd_p101,*shift_bd_p101,*sigma_bd_p101,*ratio_sigma_bd_p101,*fraction_sigma_bd_p101,*shiftg_bd_p101);

  //-----------------------------------------------------------------------------------------------------------------------------------

  // Delete RooRealVars
  //deleteVariables();

}

void PartRecoShapes::fix()
{
  a_Bu_DstKst_D0pi0->setConstant();
  b_Bu_DstKst_D0pi0->setConstant();
  a_Bu_DstKst_D0gamma->setConstant();
  b_Bu_DstKst_D0gamma->setConstant();
  a_Bd_DstKst_D0pi0->setConstant();
  b_Bd_DstKst_D0pi0->setConstant();

  csi_bu_p010->setConstant();
  csi_bu_g010->setConstant();
  csi_bu_p101->setConstant();
  csi_bu_g101->setConstant();
  csi_bd_p010->setConstant();
  csi_bd_p101->setConstant();

  shift_bu_p010->setConstant();
  shift_bu_g010->setConstant();
  shift_bu_p101->setConstant();
  shift_bu_g101->setConstant();
  shift_bd_p010->setConstant();
  shift_bd_p101->setConstant();

  sigma_bu_p010->setConstant();
  sigma_bu_g010->setConstant();
  sigma_bu_p101->setConstant();
  sigma_bu_g101->setConstant();
  sigma_bd_p010->setConstant();
  sigma_bd_p101->setConstant();

  ratio_sigma_bu_p010->setConstant();
  ratio_sigma_bu_g010->setConstant();
  ratio_sigma_bu_p101->setConstant();
  ratio_sigma_bu_g101->setConstant();
  ratio_sigma_bd_p010->setConstant();
  ratio_sigma_bd_p101->setConstant();

  fraction_sigma_bu_p010->setConstant();
  fraction_sigma_bu_g010->setConstant();
  fraction_sigma_bu_p101->setConstant();
  fraction_sigma_bu_g101->setConstant();
  fraction_sigma_bd_p010->setConstant();
  fraction_sigma_bd_p101->setConstant();

  shiftg_bu_p010->setConstant();
  shiftg_bu_g010->setConstant();
  shiftg_bu_p101->setConstant();
  shiftg_bu_g101->setConstant();
  shiftg_bd_p010->setConstant();
  shiftg_bd_p101->setConstant();

}

void PartRecoShapes::deleteVariables() {

	delete a_Bu_DstKst_D0pi0;
	delete b_Bu_DstKst_D0pi0;
	delete a_Bu_DstKst_D0gamma;
	delete b_Bu_DstKst_D0gamma;
	delete a_Bd_DstKst_D0pi0;
	delete b_Bd_DstKst_D0pi0;

	delete csi_bu_p010;
	delete csi_bu_g010;
	delete csi_bu_p101;
	delete csi_bu_g101;
	delete csi_bd_p010;
	delete csi_bd_p101;

	delete shift_bu_p010;
	delete shift_bu_g010;
	delete shift_bu_p101;
	delete shift_bu_g101;
	delete shift_bd_p010;
	delete shift_bd_p101;

	delete sigma_bu_p010;
	delete sigma_bu_g010;
	delete sigma_bu_p101;
	delete sigma_bu_g101;
	delete sigma_bd_p010;
	delete sigma_bd_p101;

	delete ratio_sigma_bu_p010;
	delete ratio_sigma_bu_g010;
	delete ratio_sigma_bu_p101;
	delete ratio_sigma_bu_g101;
	delete ratio_sigma_bd_p010;
	delete ratio_sigma_bd_p101;

	delete fraction_sigma_bu_p010;
	delete fraction_sigma_bu_g010;
	delete fraction_sigma_bu_p101;
	delete fraction_sigma_bu_g101;
	delete fraction_sigma_bd_p010;
	delete fraction_sigma_bd_p101;

	delete shiftg_bu_p010;
	delete shiftg_bu_g010;
	delete shiftg_bu_p101;
	delete shiftg_bu_g101;
	delete shiftg_bd_p010;
	delete shiftg_bd_p101;
}

