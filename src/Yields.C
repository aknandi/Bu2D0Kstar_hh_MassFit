#include <vector>
#include <iostream>
#include <map>
#include <math.h>
#include "Yields.h"
//#include "InternalStorage.h"
//#include "CommonTools.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgSet.h"
//#include "RooUnblindUniform.h"
using namespace std;
using namespace RooFit;

class InternalStorage;

Yields::Yields(Settings* genConfs, Settings* fileList, std::vector<std::string> allmodeList, std::vector<std::string>allchargeList, std::vector<std::string>alltrackList, std::vector<std::string>allrunList, std::string unblind)
{
  _genConfs = genConfs;
  _fileList = fileList;
  //stuff so that the bin lists and categories are defined and also the paramter object p.

  modeList=allmodeList;
  chargeList=allchargeList;
  trackList=alltrackList;
  runList=allrunList;

  input = new Settings("SetRatios");
  input->readPairStringsToMap(_fileList->get("PathnameToRatios"));
  input->readPairStringsToMap(_fileList->get("PathnameToYieldCorrections"));
  input->readPairStringsToMap(_fileList->get("PathnameToTotals"));
  input->readPairStringsToMap(_fileList->get("PathnameToFitRatios"));

  //myMaps.Initialize(fileList, modeList,chargeList,trackList,binList);

  genscale = 1.; // scale yields for generation
  limitlow = _genConfs->get("fit_limit_low");
  unblind=_genConfs->get("UNBLIND");

  std::cout << "*****************************************************" << std::endl;
  std::cout << "Initialising yields ... " << std::endl;
  SetOtherBkgs(); // this has to come before setyieldsGenandfit()
  SetDstKstGenandFit();
  SetYieldRatios();
  SetYieldsGenandFit(); // must be last

  std::cout<<" Yields done "<<std::endl;

}

void Yields::SetOtherBkgs()
{

/*  double scalefactor = input->getD("lowmass_factor");

  std::cout << "Scaling D3h backgrounds by factor: " << scalefactor << std::endl;
  // Shared kspipi and kskk 
  for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
    for(std::vector<std::string>::iterator a=binList.begin(); a!=binList.end();a++){
      for(std::vector<std::string>::iterator c=chargeList.begin(); c!=chargeList.end();c++){

        // --- DKpipi ratio --- 
        ratio_bs_dkpipi[*c][*t][*a] = new RooRealVar(Form("ratio_bs_dkpipi_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"", scalefactor*input->getD(Form("Bu_DKpipi_ratio_to_Bs_%s",limitlow.c_str())),0.,5.0);
        // Gaussian constrain
        gausratio_bs_dkpipi[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_dkpipi_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_dkpipi[*c][*t][*a], RooFit::RooConst(scalefactor*input->getD(Form("Bu_DKpipi_ratio_to_Bs_%s",limitlow.c_str()))), RooFit::RooConst(scalefactor*input->getD(Form("Bu_DKpipi_ratio_to_Bs_%s_err",limitlow.c_str()))));

        // --- Dpipipi ratio --- 
        ratio_bs_dpipipi[*c][*t][*a] = new RooRealVar(Form("ratio_bs_dpipipi_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"", scalefactor*input->getD(Form("Bu_Dpipipi_ratio_to_Bs_%s",limitlow.c_str())),0.,1.0);
        // Gaussian constrain
        gausratio_bs_dpipipi[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_dpipipi_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_dpipipi[*c][*t][*a], RooFit::RooConst(scalefactor*input->getD(Form("Bu_Dpipipi_ratio_to_Bs_%s",limitlow.c_str()))), RooFit::RooConst(scalefactor*input->getD(Form("Bu_Dpipipi_ratio_to_Bs_%s_err",limitlow.c_str()))));

        // --- Lambda_b -> Dppi ratio ---
        double lambda_ratio = input->getD(Form("Lb_Dppi_misIDpK_ratio_to_Bs_%s",limitlow.c_str())) * input->getD("Lb_Dppi_factor");
        double lambda_ratio_err = input->getD(Form("Lb_Dppi_misIDpK_ratio_to_Bs_%s_err",limitlow.c_str())) * input->getD("Lb_Dppi_factor");
        ratio_bs_lambda[*c][*t][*a] = new RooRealVar(Form("ratio_bs_lambda_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"",lambda_ratio,0.0,1.0); 
        // Gaussian constrain       
        gausratio_bs_lambda[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_lambda_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_lambda[*c][*t][*a], RooFit::RooConst(lambda_ratio), RooFit::RooConst(lambda_ratio_err));

		}
      }
    }*/
  }


void Yields::SetDstKstGenandFit()
{
 /* input->readPairStringsToMap(_fileList->get("PathnameToRatiosBsDstKst"));
  input->readPairStringsToMap(_fileList->get("PathnameToRatiosBdDstKst"));

  // Shared Kspipi and KsKK
  for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
    for(std::vector<std::string>::iterator a=binList.begin(); a!=binList.end();a++){
      for(std::vector<std::string>::iterator c=chargeList.begin(); c!=chargeList.end();c++){

        double r_dstkst=0.;
        if(_genConfs->get("genToys")=="false"){ 
          r_dstkst = input->getD(Form("Bs_DstKst_ratio_to_Bs_%s",limitlow.c_str())); 
        }
        else {
          // 1. Calculated
          //double r_dstkst = input->getD(Form("Bs_DstKst_ratio_to_Bs_%s",limitlow.c_str())); 
          // 2. Fitted
          r_dstkst = input->getD(Form("ratio_bs_dstkst_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()));
        }

        // --- B_s -> D*K* ratio --- 
        ratio_bs_dstkst[*c][*t][*a] = new RooRealVar(Form("ratio_bs_dstkst_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"",r_dstkst,0.,5.0); 
        // Gaussian constrain
        gausratio_bs_dstkst[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_dstkst_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_dstkst[*c][*t][*a], 
                                                          RooFit::RooConst(r_dstkst),
                                                          RooFit::RooConst(input->getD(Form("Bs_DstKst_ratio_to_Bs_%s_err",limitlow.c_str())))); 
        
        // --- B_s -> D*K* ratio split by Helamp ---

        double ratio_tot = input->getD(Form("Bs_DstKst_ratio_to_Bs_%s",limitlow.c_str()));
        double ratio_tot_err = input->getD(Form("Bs_DstKst_ratio_to_Bs_%s_err",limitlow.c_str()));
        double frac010 = input->getD(Form("frac010_bs_%s",limitlow.c_str()));
        double frac010_err = input->getD(Form("frac010_bs_%s_err",limitlow.c_str()));
        double frac001 = 1-frac010;
        double frac001_err = frac010_err * frac001 / frac010;
        double ratio_010 = ratio_tot * frac010;
        double ratio_001 = ratio_tot * frac001;
        double ratio_010_err = sqrt( pow(ratio_tot_err/ratio_tot,2) + pow(frac010_err/frac010,2) ) * ratio_010;
        double ratio_001_err = sqrt( pow(ratio_tot_err/ratio_tot,2) + pow(frac001_err/frac001,2) ) * ratio_001;

        ratio_bs_dstkst_010[*c][*t][*a] = new RooRealVar(Form("ratio_bs_dstkst_010_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"",
                                                         ratio_010,0,5.0); 
        ratio_bs_dstkst_001[*c][*t][*a] = new RooRealVar(Form("ratio_bs_dstkst_001_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"",
                                                         ratio_001,0,5.0); 

        // Gaussian constrain
        gausratio_bs_dstkst_010[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_dstkst_010_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_dstkst_010[*c][*t][*a], 
                                                              RooFit::RooConst(ratio_010),RooFit::RooConst(ratio_010_err));
        gausratio_bs_dstkst_001[*c][*t][*a] = new RooGaussian(Form("gausratio_bs_dstkst_001_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bs_dstkst_001[*c][*t][*a], 
                                                              RooFit::RooConst(ratio_001),RooFit::RooConst(ratio_001_err)); 


        // --- B_d -> D*K* ratio --- 
        ratio_bd_dstkst[*c][*t][*a] = new RooRealVar(Form("ratio_bd_dstkst_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()),"",
                                                          input->getD(Form("Bd_DstKst_ratio_to_Bs_%s",limitlow.c_str())),0.,5.0); 

        // Gaussian constrain
        gausratio_bd_dstkst[*c][*t][*a] = new RooGaussian(Form("gausratio_bd_dstkst_%s_%s_%s",(*c).c_str(), (*t).c_str(), (*a).c_str()), "", *ratio_bd_dstkst[*c][*t][*a], 
                                                          RooFit::RooConst(input->getD(Form("Bd_DstKst_ratio_to_Bs_%s",limitlow.c_str()))), 
                                                          RooFit::RooConst(input->getD(Form("Bd_DstKst_ratio_to_Bs_%s_err",limitlow.c_str())))); 

      }
    }
  }
*/

}

void::Yields::SetYieldRatios()
{
   	for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end();m++){

	    A[*m] = new RooRealVar(Form("A_%s",(*m).c_str()),"",0.5,0.0,1.0);

	    if (*m != "d2kpi") {
	    	R[*m] = new RooRealVar(Form("R_%s",(*m).c_str()),"",0.5,0.0,1.0);
	    }
	    else {
	        for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
	          for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){

	       		  N_kpi[*t][*a] = new RooRealVar(Form("N_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",0,100000);

	          }
	        }
	    }
	  }

}

void Yields::SetYieldsGenandFit()
{
  for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){
    for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
      for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){
        for(std::vector<std::string>::iterator c=chargeList.begin(); c!=chargeList.end();c++){

        	const char* identifier = Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str());

			// If fit is charge separated fit to the asymmetries
			// Need to write the signal yields in terms of the fit parameters (A, R, Ni)
			if(_genConfs->isChargeSeparated())
			{

				if(*m == "d2kpi") {
					if(*c == "plus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-@1)",RooArgList(*N_kpi[*t][*a],*A[*m]));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-@1)",RooArgList(*N_kpi[*t][*a],*A[*m]));
					}
					else if (*c == "minus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+@1)",RooArgList(*N_kpi[*t][*a],*A[*m]));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+@1)",RooArgList(*N_kpi[*t][*a],*A[*m]));
					}
				}
				else {
					if(*c == "plus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-@1)*@2",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m]));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-@1)*@2",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m]));
					}
					else if (*c == "minus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+@1)*@2",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m]));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+@1)*@2",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m]));
					}
				}

			}

			// If fit is not charge separated fit to the signal yields
			if(!_genConfs->isChargeSeparated())
			{

				double N_bu   = input->getD(Form("N_bu_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
				n_bu_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_gen_%s",identifier),"",N_bu,0,10000);

				n_bu_fit[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_fit_%s",identifier),"",N_bu,-10.,10000.);
			}

			// The rest of the pdf shapes are always fitted for yields
			// --- Gen yields ---
			double N_comb = input->getD(Form("N_comb_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
			n_comb_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_gen_%s",identifier),"",N_comb,0,10000);
			double N_bu_dstkst = input->getD(Form("N_bu_dstkst_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
			n_bu_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_dstkst_gen_%s",identifier),"",N_bu_dstkst,0,10000);

			// --- Fit yields ---
			n_comb[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_%s",identifier),"",N_comb,0,10000.);
			n_bu_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_dstkst_%s",identifier),"",N_bu_dstkst,0.,10000.);


        }
      }
    }
  }
}



