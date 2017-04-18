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

Yields::Yields(Settings* genConfs, Settings* fileList, std::vector<std::string> allmodeList, std::vector<std::string>allchargeList, std::vector<std::string>alltrackList, std::vector<std::string>allrunList, std::string unblind) : Base()
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

   		// Asymmetry observable in all modes except suppressed ones
   		if(*m != "d2pik" && *m != "d2pikpipi") {
   			double AfromFit = input->getD(Form("A_%s",(*m).c_str()));
   			A[*m] = new RooRealVar(Form("A_%s",(*m).c_str()),"",AfromFit,-5.0,5.0);
   		}
   		// Yield ratio in CP modes
	    if (*m == "d2kk" || *m == "d2pipi" || *m == "d2pipipipi") {
	    	double RfromFit = input->getD(Form("R_%s",(*m).c_str()));
	    	R[*m] = new RooRealVar(Form("R_%s",(*m).c_str()),"",RfromFit,-5.0,10.0);
	    }
	    // R+ and R- in ADS mode
	    else if (*m == "d2pik") {
    		double Rplus_pik = input->getD(Form("Rplus_%s",(*m).c_str()));
    		double Rminus_pik = input->getD(Form("Rminus_%s",(*m).c_str()));
    		A[*m] = new RooRealVar(Form("Rplus_%s",(*m).c_str()),"",Rplus_pik,-5.0,10.0);
    		R[*m] = new RooRealVar(Form("Rminus_%s",(*m).c_str()),"",Rminus_pik,-5.0,10.0);


	    	if(_genConfs->get("UNBLIND")=="false" && _genConfs->get("genToys")=="false") {
	    		Rplus = new RooUnblindUniform(Form("Rplus_%s_unblind",(*m).c_str()),"Rplus unblind","StringToBlindRplus",0.01,*A[*m]);
	    		Rminus = new RooUnblindUniform(Form("Rminus_%s_unblind",(*m).c_str()),"Rminus unblind","StringToBlindRminus",0.01,*R[*m]);
	    	}
	    	else {
	    		Rplus = A[*m];
	    		Rminus = R[*m];
	    	}
	    }
	    // R+ and R- in K3pi mode
	    else if (*m == "d2pikpipi") {
	    	double Rplus_pikpipi = input->getD(Form("Rplus_%s",(*m).c_str()));
	    	double Rminus_pikpipi = input->getD(Form("Rminus_%s",(*m).c_str()));
	    	A[*m] = new RooRealVar(Form("Rplus_%s",(*m).c_str()),"",Rplus_pikpipi,-5.0,10.0);
	    	R[*m] = new RooRealVar(Form("Rminus_%s",(*m).c_str()),"",Rminus_pikpipi,-5.0,10.0);
	    }
	    // Yield in 2-body favoured mode for each category (Summed over charge)
	    else if(*m == "d2kpi") {
	        for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
	          for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){

	        	  double N_kpifromFit;
	        	  N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*genscale;
	        	  /*
	        		if(*t == "LL") N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEffSig_%s_%s",(*t).c_str(),(*m).c_str()))*input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutLL.c_str()))*genscale;
	        		if(*t == "DD") N_kpifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("pidEffSig_%s_%s",(*t).c_str(),(*m).c_str()))*input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutDD.c_str()))*genscale;
	        	  */
	        	  N_kpi[*t][*a] = new RooRealVar(Form("N_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",N_kpifromFit,-10.,1000000.);

	          }
	        }
	    }
	    // Yield in 4-body favoured mode for each category (Summed over charge)
	    else if(*m == "d2kpipipi") {
	    	for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
	    		for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){

	    			double N_kpipipifromFit;
	    			//N_kpipipifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*genscale;
	    			if(*t == "LL") N_kpipipifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effSig_%s_%s_%s_%s",(*m).c_str(),(*t).c_str(),bdtcutLL.c_str(),(*a).c_str()))*genscale;
	    			if(*t == "DD") N_kpipipifromFit = input->getD(Form("N_bu_%s_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effSig_%s_%s_%s_%s",(*m).c_str(),(*t).c_str(),bdtcutDD.c_str(),(*a).c_str()))*genscale;
	    			N_kpipipi[*t][*a] = new RooRealVar(Form("N_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",N_kpipipifromFit,-10.,1000000.);

	    		}
	    	}
	    }

	  }
}

void Yields::SetYieldsGenandFit(std::string kstmasscut,std::string kshelcut,std::string bdtcutLL,std::string bdtcutDD,std::string bdtadsLL,std::string bdtadsDD)
{
	RooRealVar* effCorrection;
	RooRealVar* prodAsymmetry;
	RooRealVar* charmlessContributionPlus;
	RooRealVar* charmlessContributionMinus;
	RooRealVar* pidAsymmetry;
	double Akpi = input->getD("A_det_kpi");
	double Api = input->getD("A_det_pi");
	RooRealVar* detAkpi = new RooRealVar("detAsymmetry_kpi","",Akpi + (_genConfs->get("detectionAsymmetry")=="1"?gRandom->Gaus(0,input->getD("A_det_kpi_err")):0.));;
	RooRealVar* detApi = new RooRealVar("detAsymmetry_pi","",Api + (_genConfs->get("detectionAsymmetry")=="1"?gRandom->Gaus(0,input->getD("A_det_pi_err")):0.));

	for(std::vector<std::string>::iterator a=runList.begin(); a!=runList.end();a++){

		double productionAsymmetry = input->getD(Form("A_prod_%s",(*a).c_str()));
		double Apid = input->getD(Form("A_pid_%s",(*a).c_str()));
		// If switch in on in GeneralSettings, run the systematic
		prodAsymmetry = new RooRealVar(Form("prodAsymmetry_%s",(*a).c_str()),"",productionAsymmetry + (_genConfs->get("productionAsymmetry")=="1"?gRandom->Gaus(0,input->getD(Form("A_prod_%s_err",(*a).c_str()))):0.));
		pidAsymmetry = new RooRealVar(Form("pidAsymmetry_%s",(*a).c_str()),"",Apid + (_genConfs->get("pidAsymmetry")=="1"?gRandom->Gaus(0,input->getD(Form("A_pid_%s_err",(*a).c_str()))):0.));

		for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){
			for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
				for(std::vector<std::string>::iterator c=chargeList.begin(); c!=chargeList.end();c++){

					const char* identifier = Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str());

					// If fit is charge separated fit to the asymmetries
					// Need to write the signal yields in terms of the fit parameters (A, R, Ni)
					if(_genConfs->isChargeSeparated())
					{
						if(*m == "d2kpi") {
							if(*c == "plus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4+@5))",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4+@5))",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
							}
							else if (*c == "minus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4+@5))",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4+@5))",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
							}
						}
						else if(*m == "d2kk" || *m == "d2pipi") {
							double efficiencyCorrection = (input->getD("BR_d2kpi") * input->getD(Form("effSel_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_d2kpi_%s_%s",(*a).c_str(),(*t).c_str()))) / (input->getD(Form("BR_%s",(*m).c_str())) * input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())));
							double errmcEfficiencyCorrection = efficiencyCorrection * sqrt(pow(input->getD(Form("effSel_d2kpi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effSel_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effSel_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2));
							double errpidEfficiencyCorrection = efficiencyCorrection * sqrt(pow(0.002/input->getD(Form("effPid_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())),2) + pow(0.002/input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2));
							double errBranchingRatio = efficiencyCorrection * sqrt(pow(input->getD("BR_d2kpi_err")/input->getD("BR_d2kpi"),2) + pow(input->getD(Form("BR_%s_err",(*m).c_str()))/input->getD(Form("BR_%s",(*m).c_str())),2));
							double charmlessEvents = round(gRandom->Gaus(input->getD(Form("charmless_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str())),input->getD(Form("err_charmless_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()))));
							if(charmlessEvents<=0) charmlessEvents = 0;
							double charmlessPlus = round(gRandom->Poisson(0.5*charmlessEvents));
							double charmlessMinus = charmlessEvents - charmlessPlus;

							effCorrection = new RooRealVar(Form("effCorrection_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyCorrection + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errmcEfficiencyCorrection):0.) + (_genConfs->get("pidefficiencies")=="1"?gRandom->Gaus(0,errpidEfficiencyCorrection):0.) + (_genConfs->get("branchingRatios")=="1"?gRandom->Gaus(0,errBranchingRatio):0.));

							if(*c == "plus") {
								charmlessContributionPlus = new RooRealVar(Form("charmlessEvents_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",(_genConfs->get("charmless")=="1"?charmlessPlus:0.));
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4))*(@5/@6) + @7",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection,*charmlessContributionPlus));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4))*(@5/@6) + @7",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection,*charmlessContributionPlus));
							}
							else if (*c == "minus") {
								charmlessContributionMinus = new RooRealVar(Form("charmlessEvents_%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str()),"",(_genConfs->get("charmless")=="1"?charmlessMinus:0.));
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4))*(@5/@6) + @7",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection,*charmlessContributionMinus));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4))*(@5/@6) + @7",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection,*charmlessContributionMinus));
							}
						}
						else if(*m == "d2pik") {
							// Account for the efficiencies of the double misID veto and different BDT cut in ADS mode
							double efficiencyVeto = input->getD(Form("effVeto_%s_%s",(*a).c_str(),(*t).c_str()));
							double errEfficiencyVeto = input->getD(Form("effVeto_%s_%s_err",(*a).c_str(),(*t).c_str()));

							double efficiencyBdt;
							efficiencyBdt = input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str()));
							/*
							if (*t == "LL") efficiencyBdt = input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsLL.c_str()))/input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutLL.c_str()));
							if (*t == "DD") efficiencyBdt = input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtadsDD.c_str()))/input->getD(Form("effSig_%s_%s_%s_%s_25",(*t).c_str(),kstmasscut.c_str(),kshelcut.c_str(),bdtcutDD.c_str()));
							*/

							double errEfficiencyBdt = efficiencyBdt * sqrt(pow(input->getD(Form("effBdt_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effBdt_d2kpi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())),2));
							RooRealVar* effVeto = new RooRealVar(Form("effVeto_%s_%s",(*a).c_str(),(*t).c_str()),"",efficiencyVeto + (_genConfs->get("vetoefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyVeto):0.));
							RooRealVar* effBdt = new RooRealVar(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyBdt + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyBdt):0.));

							if(*c == "plus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1*@2*@3/(2*@4+1)",RooArgList(*n_bu_gen["d2kpi"]["plus"][*t][*a],*effVeto,*effBdt,*Rplus,*detAkpi));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1*@2*@3/(2*@4+1)",RooArgList(*n_bu_fit["d2kpi"]["plus"][*t][*a],*effVeto,*effBdt,*Rplus,*detAkpi));
							}
							else if (*c == "minus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1*@2*@3*(2*@4+1)",RooArgList(*n_bu_gen["d2kpi"]["minus"][*t][*a],*effVeto,*effBdt,*Rminus,*detAkpi));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1*@2*@3*(2*@4+1)",RooArgList(*n_bu_fit["d2kpi"]["minus"][*t][*a],*effVeto,*effBdt,*Rminus,*detAkpi));
							}
						}
						else if(*m == "d2kpipipi") {
							if(*c == "plus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4+@5))",RooArgList(*N_kpipipi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4+@5))",RooArgList(*N_kpipipi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
							}
							else if (*c == "minus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4+@5))",RooArgList(*N_kpipipi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4+@5))",RooArgList(*N_kpipipi[*t][*a],*A[*m],*prodAsymmetry,*detAkpi,*detApi,*pidAsymmetry));
							}
						}
						else if(*m == "d2pipipipi") {
							double efficiencyCorrection = (input->getD("BR_d2kpipipi") * input->getD(Form("effSel_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str()))) / (input->getD(Form("BR_%s",(*m).c_str())) * input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())));
							double errmcEfficiencyCorrection = efficiencyCorrection * sqrt(pow(input->getD(Form("effSel_d2kpipipi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effSel_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effSel_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2));
							double errpidEfficiencyCorrection = efficiencyCorrection * sqrt(pow(0.002/input->getD(Form("effPid_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())),2) + pow(0.002/input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2));
							double errBranchingRatio = efficiencyCorrection * sqrt(pow(input->getD("BR_d2kpipipi_err")/input->getD("BR_d2kpipipi"),2) + pow(input->getD(Form("BR_%s_err",(*m).c_str()))/input->getD(Form("BR_%s",(*m).c_str())),2));

							effCorrection = new RooRealVar(Form("effCorrection_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyCorrection + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errmcEfficiencyCorrection):0.) + (_genConfs->get("pidefficiencies")=="1"?gRandom->Gaus(0,errpidEfficiencyCorrection):0.) + (_genConfs->get("branchingRatios")=="1"?gRandom->Gaus(0,errBranchingRatio):0.));

							if(*c == "plus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4))*(@5/@6)",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1-(@1+@2+@3+@4))*(@5/@6)",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection));
							}
							else if (*c == "minus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4))*(@5/@6)",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"0.5*@0*(1+(@1+@2+@3+@4))*(@5/@6)",RooArgList(*N_kpi[*t][*a],*A[*m],*prodAsymmetry,*detApi,*pidAsymmetry,*R[*m],*effCorrection));
							}
						}
						else if(*m == "d2pikpipi") {
							/*
							// Only to be used if LL and DD bdt cuts are different between favoured and suppressed K3pi
							double efficiencyBdt = input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str()));
							double errEfficiencyBdt = efficiencyBdt * sqrt(pow(input->getD(Form("effBdt_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effBdt_d2kpipipi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())),2));
							RooRealVar* effBdt = new RooRealVar(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyBdt + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyBdt):0.));
							*/
							if(*c == "plus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1/(2*@2+1)",RooArgList(*n_bu_gen["d2kpipipi"]["plus"][*t][*a],*A[*m],*detAkpi));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1/(2*@2+1)",RooArgList(*n_bu_fit["d2kpipipi"]["plus"][*t][*a],*A[*m],*detAkpi));
							}
							else if (*c == "minus") {
								n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1*(2*@2+1)",RooArgList(*n_bu_gen["d2kpipipi"]["minus"][*t][*a],*R[*m],*detAkpi));
								n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1*(2*@2+1)",RooArgList(*n_bu_fit["d2kpipipi"]["minus"][*t][*a],*R[*m],*detAkpi));
							}
						}
					}


					// If fit is not charge separated fit to the signal yields
					if(!_genConfs->isChargeSeparated()) {

						//double N_bu   = input->getD(Form("N_bu_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))*genscale;
						//n_bu_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_gen_%s",identifier),"",N_bu,0.,100000.);
						//n_bu_fit[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_fit_%s",identifier),"",N_bu,0.,100000.);

						if(*m!="d2pik" && *m!="d2pikpipi") {
							double N_bu   = input->getD(Form("N_bu_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))*genscale;
							n_bu_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_gen_%s",identifier),"",N_bu,0.,100000.);
							n_bu_fit[*m][*c][*t][*a] = new RooRealVar(Form("n_bu_fit_%s",identifier),"",N_bu,0.,100000.);
						}
						else if(*m=="d2pik") {
							double efficiencyVeto = input->getD(Form("effVeto_%s_%s",(*a).c_str(),(*t).c_str()));
							double errEfficiencyVeto = input->getD(Form("effVeto_%s_%s_err",(*a).c_str(),(*t).c_str()));
							double efficiencyBdt = input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str()));
							double errEfficiencyBdt = efficiencyBdt * sqrt(pow(input->getD(Form("effBdt_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effBdt_d2kpi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())),2));
							RooRealVar* effVeto = new RooRealVar(Form("effVeto_%s_%s",(*a).c_str(),(*t).c_str()),"",efficiencyVeto + (_genConfs->get("vetoefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyVeto):0.));
							RooRealVar* effBdt = new RooRealVar(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyBdt + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyBdt):0.));

							RooRealVar *Rads = new RooRealVar("Rads","",0,1);
							n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1*@2*@3",RooArgList(*Rads,*n_bu_gen["d2kpi"][*c][*t][*a],*effVeto,*effBdt));
							n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1*@2*@3",RooArgList(*Rads,*n_bu_fit["d2kpi"][*c][*t][*a],*effVeto,*effBdt));
						}
						else if(*m=="d2pikpipi") {
							/* Only to be used if LL and DD bdt cuts are different between favoured and suppressed K3pi
							double efficiencyBdt = input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str()));
							double errEfficiencyBdt = efficiencyBdt * sqrt(pow(input->getD(Form("effBdt_%s_%s_%s_err",(*m).c_str(),(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())),2) + pow(input->getD(Form("effBdt_d2kpipipi_%s_%s_err",(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())),2));
							RooRealVar* effBdt = new RooRealVar(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()),"",efficiencyBdt + (_genConfs->get("mcefficiencies")=="1"?gRandom->Gaus(0,errEfficiencyBdt):0.));
							 */
							RooRealVar *Rads_k3pi = new RooRealVar("Rads_k3pi","",0,1);
							n_bu_gen[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_gen_%s",identifier),"@0*@1",RooArgList(*Rads_k3pi,*n_bu_gen["d2kpipipi"][*c][*t][*a]));
							n_bu_fit[*m][*c][*t][*a] = new RooFormulaVar(Form("n_bu_fit_%s",identifier),"@0*@1",RooArgList(*Rads_k3pi,*n_bu_fit["d2kpipipi"][*c][*t][*a]));
						}

					}

					// The rest of the pdf shapes are always fitted for yields
					// --- Gen yields ---
					//double N_comb = input->getD(Form("N_comb_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
					// ncomb = 0.5 * kpi comb yield * efficiency of kst selection (account for split by charge)
					double N_comb;
					// Used for optimisation of cuts
					if(*t == "LL") {
						if( *m=="d2kpipipi" || *m=="d2pikpipi" || *m=="d2pipipipi") N_comb = 0.5*input->getD(Form("N_comb_d2kpipipi_both_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpipipi_%s_%s",(*t).c_str(),bdtcutLL.c_str(),(*a).c_str()))*input->getD(Form("effComb_kpipipito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
						else N_comb = N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
					}
					if(*t == "DD") {
						if( *m=="d2kpipipi" || *m=="d2pikpipi" || *m=="d2pipipipi") N_comb = 0.5*input->getD(Form("N_comb_d2kpipipi_both_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_%s_d2kpipipi_%s_%s",(*t).c_str(),bdtcutDD.c_str(),(*a).c_str()))*input->getD(Form("effComb_kpipipito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
						else N_comb = N_comb = 0.5*input->getD(Form("N_comb_d2kpi_both_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effComb_kpito%s_%s",(*m).c_str(),(*t).c_str()))*genscale;
					}

					n_comb_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_gen_%s",identifier),"",N_comb,-10.,100000.);

					// double N_dstkst = input->getD(Form("N_dstkst_%s_both_%s",(*m).c_str(),(*t).c_str()))*genscale;
					// ndstkst = 0.5 * (dstkst010*eff010 + dstkst101*eff101) * efficiency of kst selection * D mode fraction (account for split by charge)
					// Total D*K* NOT split by charge- below need to add asymmetry as systematic
					double N_dstkst;
					if( *m=="d2kpipipi" || *m=="d2pikpipi" || *m=="d2pipipipi") {
						double efficiencyCorrection;
						if(*m == "d2kpipipi") efficiencyCorrection = 1.0;
						else if (*m == "d2pipipipi") efficiencyCorrection = (input->getD("BR_d2kpipipi") * input->getD(Form("effSel_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str()))) / (input->getD(Form("BR_%s",(*m).c_str())) * input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())));
						else if (*m == "d2pikpipi") efficiencyCorrection = input->getD(Form("effBdt_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str()))/(input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))*input->getD(Form("effVeto_d2kpipipi_%s_%s",(*a).c_str(),(*t).c_str())));

						if(*t == "LL") N_dstkst = input->getD(Form("N_dstkst_d2kpipipi_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effSig_d2kpipipi_%s_%s_%s",(*t).c_str(),bdtcutLL.c_str(),(*a).c_str()))*(input->getD(Form("R_%s",(*m).c_str()))/efficiencyCorrection)*genscale;
						if(*t == "DD") N_dstkst = input->getD(Form("N_dstkst_d2kpipipi_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*input->getD(Form("effSig_d2kpipipi_%s_%s_%s",(*t).c_str(),bdtcutDD.c_str(),(*a).c_str()))*(input->getD(Form("R_%s",(*m).c_str()))/efficiencyCorrection)*genscale;
					}
					else {
						double efficiencyCorrection;
						if(*m == "d2kpi") efficiencyCorrection = 1.0;
						else if (*m == "d2kk" || *m == "d2pipi") efficiencyCorrection = (input->getD("BR_d2kpi") * input->getD(Form("effSel_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_d2kpi_%s_%s",(*a).c_str(),(*t).c_str()))) / (input->getD(Form("BR_%s",(*m).c_str())) * input->getD(Form("effSel_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())) * input->getD(Form("effPid_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str())));
						else if (*m == "d2pik") efficiencyCorrection = input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str()))/(input->getD(Form("effBdt_%s_%s_%s",(*m).c_str(),(*a).c_str(),(*t).c_str()))*input->getD(Form("effVeto_%s_%s",(*a).c_str(),(*t).c_str())));

						if(*m == "d2pik" && *t == "DD") N_dstkst = input->getD(Form("N_dstkst_d2kpi_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*(input->getD(Form("effBdt_d2pik_%s_%s",(*a).c_str(),(*t).c_str()))/input->getD(Form("effBdt_d2kpi_%s_%s",(*a).c_str(),(*t).c_str())))*(input->getD(Form("R_%s",(*m).c_str()))/efficiencyCorrection)*genscale;
						else N_dstkst = input->getD(Form("N_dstkst_d2kpi_%s_%s_%s",(*a).c_str(),(*t).c_str(),limitlow.c_str()))*(input->getD(Form("R_%s",(*m).c_str()))/efficiencyCorrection)*genscale;
					}

					if(_genConfs->get("partrecoShape") == "0") {
						if(!_genConfs->isChargeSeparated()) {
							n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",N_dstkst,-10.,100000.);
							n_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_%s",identifier),"",N_dstkst);
						}
						else {
							// Assume no asymmetry
							n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",0.5*N_dstkst,-10.,100000.);
							n_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_%s",identifier),"",0.5*N_dstkst);
						}
					}
					else {
						if(!_genConfs->isChargeSeparated()) {
							n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",1.2*N_dstkst,-10.,100000.);
						}
						else {
							// Increase yield by 20% and add 10% asymmetry
							// N+ = 1.2 * 0.5*(1+asym) * N_tot
							double partRecoAsym = 0.1;
							if(*c == "plus") {
								n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",1.2*0.5*(1+partRecoAsym)*N_dstkst,-10.,100000.);
								n_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_%s",identifier),"",0.5*N_dstkst);
							}
							if(*c == "minus") {
								n_dstkst_gen[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_gen_%s",identifier),"",1.2*0.5*(1-partRecoAsym)*N_dstkst,-10.,100000.);
								n_dstkst[*m][*c][*t][*a] = new RooRealVar(Form("n_dstkst_%s",identifier),"",0.5*N_dstkst);
							}
						}
					}
					// --- Fit yields ---
					n_comb[*m][*c][*t][*a] = new RooRealVar(Form("n_comb_%s",identifier),"",N_comb,0.,100000.);

				}
			}
		}
	}
}



