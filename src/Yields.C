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
#include "RooUnblindPrecision.h"
#include "RooUnblindUniform.h"
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
  genscale = 10.; // scale yields for generation
  limitlow = _genConfs->get("fit_limit_low");
  unblind=_genConfs->get("UNBLIND");
  std::string kstmasscut = _genConfs->get("Kstmasscut");
  std::string kshelcut = _genConfs->get("Kshelcut");
  std::string bdtcutLL = _genConfs->get("BdtcutLL");
  std::string bdtcutDD = _genConfs->get("BdtcutDD");
  std::string bdtadsLL = _genConfs->get("BdtadsLL");
  std::string bdtadsDD = _genConfs->get("BdtadsDD");

  std::cout << "*****************************************************" << std::endl;
  std::cout << "Initialising yields ... " << std::endl;
  SetOtherBkgs(); // this has to come before setyieldsGenandfit()
  SetDstKstGenandFit();
  SetYieldRatios(kstmasscut,kshelcut,bdtcutLL,bdtcutDD);
  SetYieldsGenandFit(kstmasscut,kshelcut,bdtcutLL,bdtcutDD,bdtadsLL,bdtadsDD); // must be last

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

void::Yields::SetYieldRatios(std::string kstmasscut,std::string kshelcut,std::string bdtcutLL,std::string bdtcutDD)
{
   	for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end();m++){

   		if(*m != "d2pik") {
   			double AfromFit = input->getD(Form("A_%s",(*m).c_str()));
   			A[*m] = new RooRealVar(Form("A_%s",(*m).c_str()),"",AfromFit,-5.0,5.0);
   		}

	    if (*m == "d2kk" || *m == "d2pipi") {
	    	double RfromFit = input->getD(Form("R_%s",(*m).c_str()));
	    	R[*m] = new RooRealVar(Form("R_%s",(*m).c_str()),"",RfromFit,-5.0,10.0);
	    }
	    else if (*m == "d2pik") {
    		double Rplus_pik = input->getD(Form("Rplus_%s",(*m).c_str()));
    		double Rminus_pik = input->getD(Form("Rminus_%s",(*m).c_str()));
    		A[*m] = new RooRealVar(Form("Rplus_%s",(*m).c_str()),"",Rplus_pik,-5.0,10.0);
    		R[*m] = new RooRealVar(Form("Rminus_%s",(*m).c_str()),"",Rminus_pik,-5.0,10.0);

	    	if(_genConfs->get("UNBLIND")=="false" && _genConfs->get("genToys")=="false") {
	    		Rplus = new RooUnblindUniform(Form("Rplus_%s_unblind",(*m).c_str()),"Rplus unblind","StringToBlindRplus",0.1,*A[*m]);
	    		Rminus = new RooUnblindUniform(Form("Rminus_%s_unblind",(*m).c_str()),"Rminus unblind","StringToBlindRminus",0.1,*R[*m]);
	    	}
	    	else {
	    		Rplus = A[*m];
	    		Rminus = R[*m];
	    	}

	    }
	    else if(*m == "d2kpi") {
	        for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
	          for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){

	        	  //double N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))*genscale;
	        	  double N_kpifromFit;
	        	  if(*t == "LL") N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEffSig_%s_%s",(*t).c_str(),(*m).c_str()))*input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutLL.c_str()))*genscale;
	        	  if(*t == "DD") N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEffSig_%s_%s",(*t).c_str(),(*m).c_str()))*input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutDD.c_str()))*genscale;
	        	  N_kpi[*t][*a] = new RooRealVar(Form("N_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",N_kpifromFit,-1000.,1000000.);

	          }
	        }
	    }
	  }

}

void Yields::SetYieldsGenandFit(std::string kstmasscut,std::string kshelcut,std::string bdtcutLL,std::string bdtcutDD,std::string bdtadsLL,std::string bdtadsDD)
{
	RooRealVar* effCorrection;
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
				else if(*m == "d2kk" || *m == "d2pipi") {
					//double efficiencyCorrection = input->getD("BR_d2kpi")/input->getD(Form("BR_%s",(*m).c_str()));
					double efficiencyCorrection = input->getD(Form("kpitod2kpi_%s",(*t).c_str()))/input->getD(Form("kpito%s_%s",(*m).c_str(),(*t).c_str()));
					effCorrection = new RooRealVar(Form("effCorrection_%s_%s",(*m).c_str(),(*t).c_str()),"",efficiencyCorrection);

					if(*c == "plus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-@1)*@2/@3",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m],*effCorrection));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-@1)*@2/@3",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m],*effCorrection));
					}
					else if (*c == "minus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+@1)*@2/@3",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m],*effCorrection));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+@1)*@2/@3",RooArgList(*N_kpi[*t][*a],*A[*m],*R[*m],*effCorrection));
					}

				}
				else if(*m == "d2pik") {
					// Account for the efficiencies of the double misID veto adn different BDT cut in ADS mode
					double efficiencyVeto = input->getD(Form("effVeto_%s",(*t).c_str()));
					double efficiencyBdt;
					if (*t == "LL") efficiencyBdt = input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsLL.c_str()))/input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutLL.c_str()));
					if (*t == "DD") efficiencyBdt = input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsDD.c_str()))/input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutDD.c_str()));
					RooRealVar* effVeto = new RooRealVar(Form("effVeto_%s",(*t).c_str()),"",efficiencyVeto);
					RooRealVar* effBdt = new RooRealVar(Form("effBdt_%s",(*t).c_str()),"",efficiencyBdt);

					if(*c == "plus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-@1)*@2*@3*@4",RooArgList(*N_kpi[*t][*a],*A["d2kpi"],*effVeto,*effBdt,*Rplus));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-@1)*@2*@3*@4",RooArgList(*N_kpi[*t][*a],*A["d2kpi"],*effVeto,*effBdt,*Rplus));
					}
					else if (*c == "minus") {
						n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+@1)*@2*@3*@4",RooArgList(*N_kpi[*t][*a],*A["d2kpi"],*effVeto,*effBdt,*Rminus));
						n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+@1)*@2*@3*@4",RooArgList(*N_kpi[*t][*a],*A["d2kpi"],*effVeto,*effBdt,*Rminus));
					}
				}
			}


			// If fit is not charge separated fit to the signal yields
			if(!_genConfs->isChargeSeparated())
			{

				double N_bu   = input->getD(Form("N_bu_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
				n_bu_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_gen_%s",identifier),"",N_bu,0.,100000.);

				n_bu_fit[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_fit_%s",identifier),"",N_bu,0.,100000.);
			}

			// The rest of the pdf shapes are always fitted for yields
			// --- Gen yields ---
			//double N_comb = input->getD(Form("N_comb_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
			// ncomb = 0.5 * kpi comb yield * efficiency of kst selection (account for split by charge)
			double N_comb;
			if(*t == "LL") {
				if(*m != "d2pik") N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpi_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutLL.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
				else N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpi_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsLL.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
			}
			if(*t == "DD") {
				if(*m != "d2pik") N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpi_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutDD.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
				else N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpi_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsDD.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
			}
			n_comb_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_gen_%s",identifier),"",N_comb,-1000.,100000.);

			//double N_dstkst = input->getD(Form("N_dstkst_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
			// ndstkst = 0.5 * (dstkst010*eff010 + dstkst101*eff101) * efficiency of kst selection * D mode fraction (account for split by charge)
			double N_dstkst;
			if(*t == "LL") {
				if(*m != "d2pik") N_dstkst = 0.5*(input->getD(Form("N_dstkst010_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff010_%s",(*t).c_str()))*input->getD(Form("eff010_%s_%s",(*t).c_str(),kshelcut.c_str())) + input->getD(Form("N_dstkst101_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff101_%s",(*t).c_str()))*input->getD(Form("eff101_%s_%s",(*t).c_str(),kshelcut.c_str())))*input->getD(Form("effSig_%s_%s_0_%s",(*t).c_str(),kstmasscut.c_str(),bdtcutLL.c_str()))*input->getD(Form("kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
				else N_dstkst = 0.5*(input->getD(Form("N_dstkst010_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff010_%s",(*t).c_str()))*input->getD(Form("eff010_%s_%s",(*t).c_str(),kshelcut.c_str())) + input->getD(Form("N_dstkst101_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff101_%s",(*t).c_str()))*input->getD(Form("eff101_%s_%s",(*t).c_str(),kshelcut.c_str())))*input->getD(Form("effSig_%s_%s_0_%s",(*t).c_str(),kstmasscut.c_str(),bdtadsLL.c_str()))*input->getD(Form("kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
			}
			if(*t == "DD") {
				if(*m != "d2pik") N_dstkst = 0.5*(input->getD(Form("N_dstkst010_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff010_%s",(*t).c_str()))*input->getD(Form("eff010_%s_%s",(*t).c_str(),kshelcut.c_str())) + input->getD(Form("N_dstkst101_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff101_%s",(*t).c_str()))*input->getD(Form("eff101_%s_%s",(*t).c_str(),kshelcut.c_str())))*input->getD(Form("effSig_%s_%s_0_%s",(*t).c_str(),kstmasscut.c_str(),bdtcutDD.c_str()))*input->getD(Form("kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
				else N_dstkst = 0.5*(input->getD(Form("N_dstkst010_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff010_%s",(*t).c_str()))*input->getD(Form("eff010_%s_%s",(*t).c_str(),kshelcut.c_str())) + input->getD(Form("N_dstkst101_d2kpi_%s_%s",(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEff101_%s",(*t).c_str()))*input->getD(Form("eff101_%s_%s",(*t).c_str(),kshelcut.c_str())))*input->getD(Form("effSig_%s_%s_0_%s",(*t).c_str(),kstmasscut.c_str(),bdtadsDD.c_str()))*input->getD(Form("kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
			}
			n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",N_dstkst,-1000.,100000.);

			// --- Fit yields ---
			n_comb[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_%s",identifier),"",N_comb,-1000.,100000.);
			n_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_%s",identifier),"",N_dstkst);//,0.,100000.);


        }
      }
    }
  }
}



