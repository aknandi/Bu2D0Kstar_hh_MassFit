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
	    // Low mass parameters for Bu shapes
	    a_p010 = 5023.58; b_p010 = 5120.12; csi_p010 = 0.889402; shift_p010 = -1.68224; sigma_p010 = 11.0057; ratio_sigma_p010 = 8.90725; fraction_sigma_p010 = 0.962802; shiftg_p010 = -49.8499;
	    a_p101 = 5007.72; b_p101 = 5115.44; csi_p101 = 0.117698; shift_p101 = -1.90476; sigma_p101 = 12.9643; ratio_sigma_p101 = 4.06787; fraction_sigma_p101 = 0.937404; shiftg_p101 = -47.1479;
	    a_g010 = 4928.7; b_g010 = 5232.39; csi_g010 = 0.771971; shift_g010 = -0.791192; sigma_g010 = 7.39423; ratio_sigma_g010 = -0.0222271; fraction_sigma_g010 = 0.961891; shiftg_g010 = -74.2666;
	    a_g101 = 4924.72; b_g101 = 5234.22; csi_g101 = 0.564192; shift_g101 = -7.23824; sigma_g101 = 11.0883; ratio_sigma_g101 = 5.95837; fraction_sigma_g101 = 0.80563; shiftg_g101 = -38.422;

	    // Low mass parameters for Bd shapes
	    a_p010_bd = 5020.62; b_p010_bd = 5111.38; csi_p010_bd = 1.24536; shift_p010_bd = -0.00181885; sigma_p010_bd = 12.4348; ratio_sigma_p010_bd = 4.91123; fraction_sigma_p010_bd = 0.97683; shiftg_p010_bd = -99.9999;
	    a_p101_bd = 5018.91; b_p101_bd = 5112.84; csi_p101_bd = 0.71047; shift_p101_bd = -3.64037; sigma_p101_bd = 11.7453; ratio_sigma_p101_bd = 5.70878; fraction_sigma_p101_bd = 0.961203; shiftg_p101_bd = -125.209;
  }
  else if(t=="DD") {
	  // Low mass parameters for Bu shapes
	      a_p010 = 5023.4; b_p010 = 5116.4; csi_p010 = 0.822349; shift_p010 = -0.313756; sigma_p010 = 11.6947; ratio_sigma_p010 = 6.18631; fraction_sigma_p010 = 0.971878; shiftg_p010 = -61.749;
	      a_p101 = 5015.11; b_p101 = 5107.34; csi_p101 = 0.789736; shift_p101 = 6.91286; sigma_p101 = 12.3113; ratio_sigma_p101 = 9.14458; fraction_sigma_p101 = 0.965268; shiftg_p101 = -112.507;
	      a_g010 = 4913.49; b_g010 = 5211.71; csi_g010 = 0.387976; shift_g010 = 7.71486; sigma_g010 = 16.7375; ratio_sigma_g010 = 4.51313e-05; fraction_sigma_g010 = 0.972398; shiftg_g010 = -151.836;
	      a_g101 = 4926.5; b_g101 = 5233.76; csi_g101 = 0.547237; shift_g101 = -7.78492; sigma_g101 = 12.2249; ratio_sigma_g101 = 5.5341; fraction_sigma_g101 = 0.919503; shiftg_g101 = -50.1922;

	      // Low mass parameters for Bd shapes
	      a_p010_bd = 5021.15; b_p010_bd = 5109.96; csi_p010_bd = 0.938805; shift_p010_bd = -0.675983; sigma_p010_bd = 12.541; ratio_sigma_p010_bd = 6.33228; fraction_sigma_p010_bd = 0.960769; shiftg_p010_bd = -50.2524;
	      a_p101_bd = 5013.66; b_p101_bd = 5102.23; csi_p101_bd = 0.613511; shift_p101_bd = 5.13505; sigma_p101_bd = 13.1807; ratio_sigma_p101_bd = 4.50561; fraction_sigma_p101_bd = 0.953672; shiftg_p101_bd = -151.892;
  }
  else if(t=="mix") {
	  // Low mass parameters for Bu shapes
	  a_p010 = 5021; b_p010 = 5114; csi_p010 = 0.87;   shift_p010 = 0.2; sigma_p010 = 12.0; ratio_sigma_p010 = 7.2; fraction_sigma_p010 = 0.97; shiftg_p010 = -60;
	  a_p101 = 5013; b_p101 = 5105; csi_p101 = 0.60;   shift_p101 = 5.7; sigma_p101 = 13.6; ratio_sigma_p101 = 6.1; fraction_sigma_p101 = 0.93; shiftg_p101 = -114;
	  a_g010 = 4917; b_g010 = 5218; csi_g010 = 0.46;   shift_g010 = 0.0; sigma_g010 = 15.0; ratio_sigma_g010 = 6.4; fraction_sigma_g010 = 0.96; shiftg_g010 = -100;
	  a_g101 = 4928; b_g101 = 5234; csi_g101 = 0.65;   shift_g101 = -9.0; sigma_g101 = 12.5; ratio_sigma_g101 = 5.9; fraction_sigma_g101 = 0.91; shiftg_g101 = -50;

	  // Low mass parameters for Bd shapes
	  a_p010_bd = 5021; b_p010_bd = 5109; csi_p010_bd = 0.90; shift_p010_bd = -0.9; sigma_p010_bd = 12.7; ratio_sigma_p010_bd = 7.2; fraction_sigma_p010_bd = 0.96; shiftg_p010_bd = -50;
	  a_p101_bd = 5032; b_p101_bd = 5116; csi_p101_bd = 0.93; shift_p101_bd = -10;  sigma_p101_bd = 14.0; ratio_sigma_p101_bd = 6.3; fraction_sigma_p101_bd = 0.97; shiftg_p101_bd = -103;
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

