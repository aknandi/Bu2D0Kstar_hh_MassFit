#include "Pdf_Gen.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "DoubleCrystalBall.h"
#include "DoubleJohnson.h"
#include "PartRecoDstKst.h"
#include "myCruijff.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "RooWorkspace.h"

Pdf_Gen::Pdf_Gen(Settings* fileList, RooRealVar* pmB, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackTypeList, std::vector<std::string> runList)
{
  _fileList=fileList;
  _modeList=modeList;
  _chargeList=chargeList;
  _trackTypeList=trackTypeList;
  _runList=runList;

  // Initialise the PDFs
  for(std::vector<std::string>::iterator mode=_modeList.begin();mode!=_modeList.end();mode++){
    for(std::vector<std::string>::iterator charge=_chargeList.begin();charge!=_chargeList.end();charge++){
      for(std::vector<std::string>::iterator trackType=_trackTypeList.begin();trackType!=_trackTypeList.end();trackType++){
        for(std::vector<std::string>::iterator run=_runList.begin();run!=_runList.end();run++){
        	bu[*mode][*charge][*trackType][*run]  = new DoubleCrystalBall(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("gen_signal"));
        	//bu[*mode][*charge][*trackType][*run]  = new DoubleGaussian(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("gen_signal"));
        	//bu[*mode][*charge][*trackType][*run]  = new DoubleJohnson(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("gen_signal"));
        	comb[*mode][*charge][*trackType][*run]   = new Exponential(pmB, *mode,"exp",*charge,*trackType,*run,_fileList->get("gen_combs"));
        	dstkst[*mode][*charge][*trackType][*run]    = new PartRecoDstKst(pmB, *mode,*charge,*trackType,*run,_fileList->get("gen_partreco"),true);
        	if(*mode=="d2kk") lckst[*mode][*charge][*trackType][*run] = new myCruijff(pmB,*mode,"bu",*charge,*trackType,*run,_fileList->get("gen_signal"));
           }
      }
    }
  }
  // Now set relations, which prepares the PDFs to be returned
  setRelations();
}

void Pdf_Gen::setRelations()
{
  // Set up the configuration files
  Settings relConfs("Pdf_Gen::SetRelations");
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));
  relConfs.readPairStringsToMap(_fileList->get("gen_signal"));
  relConfs.readPairStringsToMap(_fileList->get("gen_combs"));
  relConfs.readPairStringsToMap(_fileList->get("gen_partreco"));
  relConfs.readPairStringsToMap(_fileList->get("PathnameToTotals"));
  relConfs.readPairStringsToMap(_fileList->get("PathnameToYieldCorrections"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
  
  //Signal -- Double Crystal Ball
  RooRealVar* bu_mean_kpi = new RooRealVar("bu_mean_kpi","",relConfs.getD("bu_mean_kpi"),
                                       relConfs.getD("bu_mean_kpi_LimL"),relConfs.getD("bu_mean_kpi_LimU") );
  RooRealVar* bu_mean_kpipipi = new RooRealVar("bu_mean_kpipipi","",relConfs.getD("bu_mean_kpipipi"),
                                       relConfs.getD("bu_mean_kpipipi_LimL"),relConfs.getD("bu_mean_kpipipi_LimU") );
  RooRealVar* bu_n_LL = new RooRealVar("bu_n_LL","",relConfs.getD("bu_n_LL"),
                                        relConfs.getD("bu_n_LL_LimL"),relConfs.getD("bu_n_LL_LimU") );
  RooRealVar* bu_n_DD = new RooRealVar("bu_n_DD","",relConfs.getD("bu_n_DD"),
                                        relConfs.getD("bu_n_DD_LimL"),relConfs.getD("bu_n_DD_LimU") );
  RooRealVar* bu_width_kpi = new RooRealVar("bu_width_kpi","",relConfs.getD("bu_width_kpi"),
                                        relConfs.getD("bu_width_kpi_LimL"),relConfs.getD("bu_width_kpi_LimU") );
  RooRealVar* bu_width_kpipipi = new RooRealVar("bu_width_kpipipi","",relConfs.getD("bu_width_kpipipi"),
                                        relConfs.getD("bu_width_kpipipi_LimL"),relConfs.getD("bu_width_kpipipi_LimU") );
  RooRealVar* bu_alpha_kpi_LL = new RooRealVar("bu_alpha_kpi_LL","",relConfs.getD("bu_alpha_kpi_LL"),
                                       relConfs.getD("bu_alpha_kpi_LL_LimL"),relConfs.getD("bu_alpha_kpi_LL_LimU") );
  RooRealVar* bu_alpha_kpipipi_LL = new RooRealVar("bu_alpha_kpipipi_LL","",relConfs.getD("bu_alpha_kpipipi_LL"),
                                       relConfs.getD("bu_alpha_kpipipi_LL_LimL"),relConfs.getD("bu_alpha_kpipipi_LL_LimU") );
  RooRealVar* bu_alpha_kpi_DD = new RooRealVar("bu_alpha_kpi_DD","",relConfs.getD("bu_alpha_kpi_DD"),
                                       relConfs.getD("bu_alpha_kpi_DD_LimL"),relConfs.getD("bu_alpha_kpi_DD_LimU") );
  RooRealVar* bu_alpha_kpipipi_DD = new RooRealVar("bu_alpha_kpipipi_DD","",relConfs.getD("bu_alpha_kpipipi_DD"),
                                       relConfs.getD("bu_alpha_kpipipi_DD_LimL"),relConfs.getD("bu_alpha_kpipipi_DD_LimU") );
  RooRealVar* bu_width_ratio_kpi_LL = new RooRealVar("bu_width_ratio_kpi_LL","",relConfs.getD("bu_width_ratio_kpi_LL"),
                                        relConfs.getD("bu_width_ratio_kpi_LL_LimL"),relConfs.getD("bu_width_ratio_kpi_LL_LimU") );
  RooRealVar* bu_width_ratio_kpi_DD = new RooRealVar("bu_width_ratio_kpi_DD","",relConfs.getD("bu_width_ratio_kpi_DD"),
                                        relConfs.getD("bu_width_ratio_kpi_DD_LimL"),relConfs.getD("bu_width_ratio_kpi_DD_LimU") );
  RooRealVar* bu_width_ratio_kpipipi_LL = new RooRealVar("bu_width_ratio_kpipipi_LL","",relConfs.getD("bu_width_ratio_kpipipi_LL"),
                                        relConfs.getD("bu_width_ratio_kpipipi_LL_LimL"),relConfs.getD("bu_width_ratio_kpipipi_LL_LimU") );
  RooRealVar* bu_width_ratio_kpipipi_DD = new RooRealVar("bu_width_ratio_kpipipi_DD","",relConfs.getD("bu_width_ratio_kpipipi_DD"),
                                        relConfs.getD("bu_width_ratio_kpipipi_DD_LimL"),relConfs.getD("bu_width_ratio_kpipipi_DD_LimU") );
  RooRealVar* bu_frac_kpi_LL = new RooRealVar("bu_frac_kpi_LL","",relConfs.getD("bu_frac_kpi_LL"),
                                        relConfs.getD("bu_frac_kpi_LL_LimL"),relConfs.getD("bu_frac_kpi_LL_LimU") );
  RooRealVar* bu_frac_kpi_DD = new RooRealVar("bu_frac_kpi_DD","",relConfs.getD("bu_frac_kpi_DD"),
                                        relConfs.getD("bu_frac_kpi_DD_LimL"),relConfs.getD("bu_frac_kpi_DD_LimU") );
  RooRealVar* bu_frac_kpipipi_LL = new RooRealVar("bu_frac_kpipipi_LL","",relConfs.getD("bu_frac_kpipipi_LL"),
                                        relConfs.getD("bu_frac_kpipipi_LL_LimL"),relConfs.getD("bu_frac_kpipipi_LL_LimU") );
  RooRealVar* bu_frac_kpipipi_DD = new RooRealVar("bu_frac_kpipipi_DD","",relConfs.getD("bu_frac_kpipipi_DD"),
                                        relConfs.getD("bu_frac_kpipipi_DD_LimL"),relConfs.getD("bu_frac_kpipipi_DD_LimU") );

  // parameters for LL and DD combined
  RooRealVar* bu_n_mix = new RooRealVar("bu_n_mix","",relConfs.getD("bu_n_mix"),
                                        relConfs.getD("bu_n_mix_LimL"),relConfs.getD("bu_n_mix_LimU") );
  RooRealVar* bu_alpha_kpi_mix = new RooRealVar("bu_alpha_kpi_mix","",relConfs.getD("bu_alpha_kpi_mix"),
                                       relConfs.getD("bu_alpha_kpi_mix_LimL"),relConfs.getD("bu_alpha_kpi_mix_LimU") );
  RooRealVar* bu_alpha_kpipipi_mix = new RooRealVar("bu_alpha_kpipipi_mix","",relConfs.getD("bu_alpha_kpipipi_mix"),
                                       relConfs.getD("bu_alpha_kpipipi_mix_LimL"),relConfs.getD("bu_alpha_kpipipi_mix_LimU") );
  RooRealVar* bu_width_ratio_kpi_mix = new RooRealVar("bu_width_ratio_kpi_mix","",relConfs.getD("bu_width_ratio_kpi_mix"),
                                          relConfs.getD("bu_width_ratio_kpi_mix_LimL"),relConfs.getD("bu_width_ratio_kpi_mix_LimU") );
  RooRealVar* bu_width_ratio_kpipipi_mix = new RooRealVar("bu_width_ratio_kpipipi_mix","",relConfs.getD("bu_width_ratio_kpipipi_mix"),
                                          relConfs.getD("bu_width_ratio_kpipipi_mix_LimL"),relConfs.getD("bu_width_ratio_kpipipi_mix_LimU") );
  RooRealVar* bu_frac_kpi_mix = new RooRealVar("bu_frac_kpi_mix","",relConfs.getD("bu_frac_kpi_mix"),
                                         relConfs.getD("bu_frac_kpi_mix_LimL"),relConfs.getD("bu_frac_kpi_mix_LimU") );
  RooRealVar* bu_frac_kpipipi_mix = new RooRealVar("bu_frac_kpipipi_mix","",relConfs.getD("bu_frac_kpipipi_mix"),
                                          relConfs.getD("bu_frac_kpipipi_mix_LimL"),relConfs.getD("bu_frac_kpipipi_mix_LimU") );

  // For Johnson
  RooRealVar* bu_delta_kpi_LL = new RooRealVar("bu_delta_kpi_LL","",relConfs.getD("bu_delta_kpi_LL"));
  RooRealVar* bu_delta_kpipipi_LL = new RooRealVar("bu_delta_kpipipi_LL","",relConfs.getD("bu_delta_kpipipi_LL"));
  RooRealVar* bu_delta_kpi_DD = new RooRealVar("bu_delta_kpi_DD","",relConfs.getD("bu_delta_kpi_DD"));
  RooRealVar* bu_delta_kpipipi_DD = new RooRealVar("bu_delta_kpipipi_DD","",relConfs.getD("bu_delta_kpipipi_DD"));
  RooRealVar* bu_gamma = new RooRealVar("bu_gamma","",relConfs.getD("bu_gamma"));

  RooRealVar* bu_delta_kpi_mix = new RooRealVar("bu_delta_mix","",relConfs.getD("bu_delta_kpi_DD"));
  RooRealVar* bu_delta_kpipipi_mix = new RooRealVar("bu_delta_mix","",relConfs.getD("bu_delta_kpipipi_DD"));


  // Combs- exponential
  RooRealVar *combs_slope_kpi_mix = new RooRealVar("exp_kpi_mix_combs_slope","",relConfs.getD("exp_kpi_mix_combs_slope"),
  		  	  	  	  	  	  	  	  relConfs.getD("exp_kpi_mix_combs_slope_LimL"), relConfs.getD("exp_kpi_mix_combs_slope_LimU") );
  RooRealVar *combs_slope_kpi_LL = new RooRealVar("exp_kpi_LL_combs_slope","",relConfs.getD("exp_kpi_LL_combs_slope"),
  		  	  	  	  	  	  	  	  relConfs.getD("exp_kpi_LL_combs_slope_LimL"), relConfs.getD("exp_kpi_LL_combs_slope_LimU") );
  RooRealVar *combs_slope_kpi_DD = new RooRealVar("exp_kpi_DD_combs_slope","",relConfs.getD("exp_kpi_DD_combs_slope"),
  		  	  	  	  	  	  	  	  relConfs.getD("exp_kpi_DD_combs_slope_LimL"), relConfs.getD("exp_kpi_DD_combs_slope_LimU") );
  RooRealVar *combs_slope_kpipipi_mix = new RooRealVar("exp_kpipipi_mix_combs_slope","",relConfs.getD("exp_kpipipi_mix_combs_slope"),
  		  	  	  	  	  	  	  	  	  relConfs.getD("exp_kpipipi_mix_combs_slope_LimL"), relConfs.getD("exp_kpipipi_mix_combs_slope_LimU") );
  RooRealVar *combs_slope_kpipipi_LL = new RooRealVar("exp_kpipipi_LL_combs_slope","",relConfs.getD("exp_kpipipi_LL_combs_slope"),
  		  	  	  	  	  	  	  	  	  relConfs.getD("exp_kpipipi_LL_combs_slope_LimL"), relConfs.getD("exp_kpipipi_LL_combs_slope_LimU") );
  RooRealVar *combs_slope_kpipipi_DD = new RooRealVar("exp_kpipipi_DD_combs_slope","",relConfs.getD("exp_kpipipi_DD_combs_slope"),
  		  	  	  	  	  	  	  	  	  relConfs.getD("exp_kpipipi_DD_combs_slope_LimL"), relConfs.getD("exp_kpipipi_DD_combs_slope_LimU") );

  //Lb->LcK*
  RooRealVar *lambda_mean = new RooRealVar("lambda_mean","",5269);
  RooRealVar *lambda_sigmaR = new RooRealVar("lambda_sigmaR","",221);
  RooRealVar *lambda_sigmaL = new RooRealVar("lambda_sigmaL","",96);
  RooRealVar *lambda_alphaR = new RooRealVar("lambda_alphaR","",-0.19);
  RooRealVar *lambda_alphaL = new RooRealVar("lambda_alphaL","",-0.04);

  //Get Ks helicity angle selection from general settings
  std::string kshelcut = relConfs.get("Kshelcut");
  std::string lowRange = relConfs.get("fit_limit_low");
  // frac010 = (n010*eff010)/(n101*eff101), different for DD and LL

  double frac010LL = (relConfs.getD(Form("N_dstkst010_d2kpi_LL_%s",lowRange.c_str()))*relConfs.getD(Form("eff010_LL_%s",kshelcut.c_str())))/(relConfs.getD(Form("N_dstkst101_d2kpi_LL_%s",lowRange.c_str()))*relConfs.getD(Form("eff101_LL_%s",kshelcut.c_str())));
  double frac010DD = (relConfs.getD(Form("N_dstkst010_d2kpi_DD_%s",lowRange.c_str()))*relConfs.getD(Form("eff010_DD_%s",kshelcut.c_str())))/(relConfs.getD(Form("N_dstkst101_d2kpi_DD_%s",lowRange.c_str()))*relConfs.getD(Form("eff101_DD_%s",kshelcut.c_str())));

  double coef010LL = frac010LL/(1+frac010LL);
  double coef101LL = 1/(1+frac010LL);
  double coef010DD = frac010DD/(1+frac010DD);
  double coef101DD= 1/(1+frac010DD);

  RooRealVar *coef010_LL = new RooRealVar("coef010_LL","",coef010LL);//, 0.0, 10.0);
  RooRealVar *coef010_DD = new RooRealVar("coef010_DD","",coef010DD);//, 0.0, 10.0);
  RooRealVar *coef101_LL = new RooRealVar("coef101_LL","",coef101LL);//, 0.0, 10.0);
  RooRealVar *coef101_DD = new RooRealVar("coef101_DD","",coef101DD);//, 0.0, 10.0);

  //RooRealVar *frac010_LL = new RooRealVar("frac010_LL","",frac010LL);//, 0.0, 10.0);
  //RooRealVar *frac010_DD = new RooRealVar("frac010_DD","",frac010DD);//, 0.0, 10.0);

  std::cout << std::endl << "PdfGen : floating parameters " << std::endl;
  std::vector<RooRealVar*> *floatParams = new std::vector <RooRealVar*>;
  floatParams->push_back(bu_mean_kpi);
  floatParams->push_back(bu_mean_kpipipi);
  floatParams->push_back(bu_width_kpi);
  floatParams->push_back(bu_width_kpipipi);


  for (Int_t n_p = 0; n_p < (Int_t)floatParams->size(); ++n_p)
    {
      RooRealVar *par = floatParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
    }
 
   //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfGen : fixed parameters " << std::endl;
  fixedParams = new std::vector <RooRealVar*>;

  //fixedParams->push_back(bu_mean);
  //fixedParams->push_back(bu_width);
  fixedParams->push_back(bu_alpha_kpi_LL);
  fixedParams->push_back(bu_alpha_kpipipi_LL);
  fixedParams->push_back(bu_n_LL);
  fixedParams->push_back(bu_alpha_kpi_DD);
  fixedParams->push_back(bu_alpha_kpipipi_DD);
  fixedParams->push_back(bu_n_DD);
  fixedParams->push_back(bu_width_ratio_kpi_LL);
  fixedParams->push_back(bu_width_ratio_kpi_DD);
  fixedParams->push_back(bu_width_ratio_kpipipi_LL);
  fixedParams->push_back(bu_width_ratio_kpipipi_DD);
  fixedParams->push_back(bu_frac_kpi_LL);
  fixedParams->push_back(bu_frac_kpi_DD);
  fixedParams->push_back(bu_frac_kpipipi_LL);
  fixedParams->push_back(bu_frac_kpipipi_DD);
  if(relConfs.get("signalShape")=="1") {
	  fixedParams->push_back(bu_delta_kpi_LL);
	  fixedParams->push_back(bu_delta_kpipipi_LL);
	  fixedParams->push_back(bu_delta_kpi_DD);
	  fixedParams->push_back(bu_delta_kpipipi_DD);
	  fixedParams->push_back(bu_gamma);
  }
  //fixedParams->push_back(bu_alpha_mix);
  //fixedParams->push_back(bu_n_mix);
  //fixedParams->push_back(combs_slope_LL);
  //fixedParams->push_back(combs_slope_DD);
  //fixedParams->push_back(frac010_LL);
  //fixedParams->push_back(frac010_DD);
  fixedParams->push_back(coef010_LL);
  fixedParams->push_back(coef010_DD);
  fixedParams->push_back(coef101_LL);
  fixedParams->push_back(coef101_DD);
  fixedParams->push_back(lambda_mean);
  fixedParams->push_back(lambda_sigmaL);
  fixedParams->push_back(lambda_sigmaR);
  fixedParams->push_back(lambda_alphaL);
  fixedParams->push_back(lambda_alphaR);

  for (Int_t n_p = 0; n_p < (Int_t)fixedParams->size(); ++n_p)
    {
      RooRealVar *par = fixedParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
      par->setConstant(kTRUE);
    }

  std::cout << "... done" << std::endl;
  
  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
//  std::map<std::string,std::map<std::string,RooRealVar*> > dk_signal_width;
	
  // Now loop over the pdfs and set the shared variables which you've just defined
  for(std::vector<std::string>::iterator mode=_modeList.begin();mode!=_modeList.end();mode++){
    for(std::vector<std::string>::iterator charge=_chargeList.begin();charge!=_chargeList.end();charge++){
      for(std::vector<std::string>::iterator trackType=_trackTypeList.begin();trackType!=_trackTypeList.end();trackType++){
        for(std::vector<std::string>::iterator run=_runList.begin();run!=_runList.end();run++){

          //bu[*mode][*charge][*trackType][*run]->setGamma(bu_gamma);

          if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
              bu[*mode][*charge][*trackType][*run]->setMean(bu_mean_kpipipi);
        	  bu[*mode][*charge][*trackType][*run]->setWidth(bu_width_kpipipi);
          }
          else {
              bu[*mode][*charge][*trackType][*run]->setMean(bu_mean_kpi);
        	  bu[*mode][*charge][*trackType][*run]->setWidth(bu_width_kpi);
          }

          if(*trackType=="LL")  {

              if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
            	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_LL); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_LL);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpipipi_LL);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_LL);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_LL);
        	  }
        	  else {
        		  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_LL); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_LL);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpi_LL);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpi_LL);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpi_LL);
        	  }

        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_LL);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_LL);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef010_LL);

          }
          else if(*trackType=="DD") {

              if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
            	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_DD); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_DD);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpipipi_DD);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_DD);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_DD);
        	  }
        	  else {
        		  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_DD); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_DD);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpi_DD);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpi_DD);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpi_DD);
        	  }

        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_DD);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_DD);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef101_DD);
          }
          else if(*trackType=="mix") {

              if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
            	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_mix); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_mix);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpipipi_mix);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_mix);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_mix);
        	  }
        	  else {
        		  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_mix); }
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_mix);
        		  //bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_kpi_mix);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpi_mix);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpi_mix);
        	  }

              bu[*mode][*charge][*trackType][*run]->setN(bu_n_mix);
              dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_DD);
              dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef101_DD);
          }
          if(*mode=="d2kk") {
        	  lckst[*mode][*charge][*trackType][*run]->setMean(lambda_mean);
        	  lckst[*mode][*charge][*trackType][*run]->setSigmaL(lambda_sigmaL);
        	  lckst[*mode][*charge][*trackType][*run]->setSigmaR(lambda_sigmaR);
        	  lckst[*mode][*charge][*trackType][*run]->setAlphaL(lambda_alphaL);
        	  lckst[*mode][*charge][*trackType][*run]->setAlphaR(lambda_alphaR);
          }
        }
      }
    }
  }


  // Finally, create map of RooPdfs to return
  for(std::vector<std::string>::iterator mode=_modeList.begin();mode!=_modeList.end();mode++){
    for(std::vector<std::string>::iterator charge=_chargeList.begin();charge!=_chargeList.end();charge++){
      for(std::vector<std::string>::iterator trackType=_trackTypeList.begin();trackType!=_trackTypeList.end();trackType++){
        for(std::vector<std::string>::iterator run=_runList.begin();run!=_runList.end();run++){
          roopdf_bu[*mode][*charge][*trackType][*run]     = bu[*mode][*charge][*trackType][*run]->getPdf();
          roopdf_comb[*mode][*charge][*trackType][*run]      = comb[*mode][*charge][*trackType][*run]->getPdf();
          roopdf_dstkst[*mode][*charge][*trackType][*run]  = dstkst[*mode][*charge][*trackType][*run]->getPdf();
          if(*mode=="d2kk") roopdf_lckst[*mode][*charge][*trackType][*run] = lckst[*mode][*charge][*trackType][*run]->getPdf();
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Gen" << endl;
}

