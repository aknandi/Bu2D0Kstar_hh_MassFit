#include "Pdf_Gen.h"
//#include "KeysPdf.h"
//#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "DoubleCrystalBall.h"
#include "PartRecoDstKst.h"
//#include "CorrGauss.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"

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
        	comb[*mode][*charge][*trackType][*run]   = new Exponential(pmB, *mode,"exp",*charge,*trackType,*run,_fileList->get("gen_combs"));
        	dstkst[*mode][*charge][*trackType][*run]    = new PartRecoDstKst(pmB, *mode,*charge,*trackType,*run,_fileList->get("gen_partreco"),true);
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
  RooRealVar* bu_mean = new RooRealVar("bu_mean","",relConfs.getD("bu_mean"),
                                       relConfs.getD("bu_mean_LimL"),relConfs.getD("bu_mean_LimU") );
  RooRealVar* bu_width_all = new RooRealVar("bu_width_all","",relConfs.getD("bu_width_all"),
                                        relConfs.getD("bu_width_all_LimL"),relConfs.getD("bu_width_all_LimU") );
  RooRealVar* bu_alpha_all_LL = new RooRealVar("bu_alpha_all_LL","",relConfs.getD("bu_alpha_all_LL"),
                                       relConfs.getD("bu_alpha_all_LL_LimL"),relConfs.getD("bu_alpha_all_LL_LimU") );
  RooRealVar* bu_n_LL = new RooRealVar("bu_n_LL","",relConfs.getD("bu_n_LL"),
                                        relConfs.getD("bu_n_LL_LimL"),relConfs.getD("bu_n_LL_LimU") );
  RooRealVar* bu_alpha_all_DD = new RooRealVar("bu_alpha_all_DD","",relConfs.getD("bu_alpha_all_DD"),
                                       relConfs.getD("bu_alpha_all_DD_LimL"),relConfs.getD("bu_alpha_all_DD_LimU") );
  RooRealVar* bu_n_DD = new RooRealVar("bu_n_DD","",relConfs.getD("bu_n_DD"),
                                        relConfs.getD("bu_n_DD_LimL"),relConfs.getD("bu_n_DD_LimU") );
  RooRealVar* bu_width_ratio_all_LL = new RooRealVar("bu_width_ratio_all_LL","",relConfs.getD("bu_width_ratio_all_LL"),
                                        relConfs.getD("bu_width_ratio_all_LL_LimL"),relConfs.getD("bu_width_ratio_all_LL_LimU") );
  RooRealVar* bu_width_ratio_all_DD = new RooRealVar("bu_width_ratio_all_DD","",relConfs.getD("bu_width_ratio_all_DD"),
                                        relConfs.getD("bu_width_ratio_all_DD_LimL"),relConfs.getD("bu_width_ratio_all_DD_LimU") );
  RooRealVar* bu_frac_all_LL = new RooRealVar("bu_frac_all_LL","",relConfs.getD("bu_frac_all_LL"),
                                        relConfs.getD("bu_frac_all_LL_LimL"),relConfs.getD("bu_frac_all_LL_LimU") );
  RooRealVar* bu_frac_all_DD = new RooRealVar("bu_frac_all_DD","",relConfs.getD("bu_frac_all_DD"),
                                        relConfs.getD("bu_frac_all_DD_LimL"),relConfs.getD("bu_frac_all_DD_LimU") );


  // RooRealVars with either trackType combined or run separated
  RooRealVar* bu_width_run1 = new RooRealVar("bu_width_run1","",relConfs.getD("bu_width_run1"),
                                        relConfs.getD("bu_width_run1_LimL"),relConfs.getD("bu_width_run1_LimU") );
  RooRealVar* bu_width_run2 = new RooRealVar("bu_width_run2","",relConfs.getD("bu_width_run2"),
                                        relConfs.getD("bu_width_run2_LimL"),relConfs.getD("bu_width_run2_LimU") );
  RooRealVar* bu_alpha_run1_LL = new RooRealVar("bu_alpha_run1_LL","",relConfs.getD("bu_alpha_run1_LL"),
                                       relConfs.getD("bu_alpha_run1_LL_LimL"),relConfs.getD("bu_alpha_run1_LL_LimU") );
  RooRealVar* bu_alpha_run2_LL = new RooRealVar("bu_alpha_run2_LL","",relConfs.getD("bu_alpha_run2_LL"),
                                       relConfs.getD("bu_alpha_run2_LL_LimL"),relConfs.getD("bu_alpha_run2_LL_LimU") );
  RooRealVar* bu_alpha_run1_DD = new RooRealVar("bu_alpha_run1_DD","",relConfs.getD("bu_alpha_run1_DD"),
                                       relConfs.getD("bu_alpha_run1_DD_LimL"),relConfs.getD("bu_alpha_run1_DD_LimU") );
  RooRealVar* bu_alpha_run2_DD = new RooRealVar("bu_alpha_run2_DD","",relConfs.getD("bu_alpha_run2_DD"),
                                       relConfs.getD("bu_alpha_run2_DD_LimL"),relConfs.getD("bu_alpha_run2_DD_LimU") );
  RooRealVar* bu_width_ratio_run1_LL = new RooRealVar("bu_width_ratio_run1_LL","",relConfs.getD("bu_width_ratio_run1_LL"),
                                        relConfs.getD("bu_width_ratio_run1_LL_LimL"),relConfs.getD("bu_width_ratio_run1_LL_LimU") );
  RooRealVar* bu_width_ratio_run1_DD = new RooRealVar("bu_width_ratio_run1_DD","",relConfs.getD("bu_width_ratio_run1_DD"),
                                        relConfs.getD("bu_width_ratio_run1_DD_LimL"),relConfs.getD("bu_width_ratio_run1_DD_LimU") );
  RooRealVar* bu_width_ratio_run2_LL = new RooRealVar("bu_width_ratio_run2_LL","",relConfs.getD("bu_width_ratio_run2_LL"),
                                        relConfs.getD("bu_width_ratio_run2_LL_LimL"),relConfs.getD("bu_width_ratio_run2_LL_LimU") );
  RooRealVar* bu_width_ratio_run2_DD = new RooRealVar("bu_width_ratio_run2_DD","",relConfs.getD("bu_width_ratio_run2_DD"),
                                        relConfs.getD("bu_width_ratio_run2_DD_LimL"),relConfs.getD("bu_width_ratio_run2_DD_LimU") );
  RooRealVar* bu_frac_run1_LL = new RooRealVar("bu_frac_run1_LL","",relConfs.getD("bu_frac_run1_LL"),
                                        relConfs.getD("bu_frac_run1_LL_LimL"),relConfs.getD("bu_frac_run1_LL_LimU") );
  RooRealVar* bu_frac_run1_DD = new RooRealVar("bu_frac_run1_DD","",relConfs.getD("bu_frac_run1_DD"),
                                        relConfs.getD("bu_frac_run1_DD_LimL"),relConfs.getD("bu_frac_run1_DD_LimU") );
  RooRealVar* bu_frac_run2_LL = new RooRealVar("bu_frac_run2_LL","",relConfs.getD("bu_frac_run2_LL"),
                                        relConfs.getD("bu_frac_run2_LL_LimL"),relConfs.getD("bu_frac_run2_LL_LimU") );
  RooRealVar* bu_frac_run2_DD = new RooRealVar("bu_frac_run2_DD","",relConfs.getD("bu_frac_run2_DD"),
                                        relConfs.getD("bu_frac_run2_DD_LimL"),relConfs.getD("bu_frac_run2_DD_LimU") );

  RooRealVar* bu_alpha_all_mix = new RooRealVar("bu_alpha_all_mix","",relConfs.getD("bu_alpha_all_mix"),
                                       relConfs.getD("bu_alpha_all_mix_LimL"),relConfs.getD("bu_alpha_all_mix_LimU") );
  RooRealVar* bu_n_mix = new RooRealVar("bu_n_mix","",relConfs.getD("bu_n_mix"),
                                        relConfs.getD("bu_n_mix_LimL"),relConfs.getD("bu_n_mix_LimU") );
  RooRealVar* bu_width_ratio_all_mix = new RooRealVar("bu_width_ratio_all_mix","",relConfs.getD("bu_width_ratio_all_mix"),
                                        relConfs.getD("bu_width_ratio_all_mix_LimL"),relConfs.getD("bu_width_ratio_all_mix_LimU") );

  // Combs- exponential
  RooRealVar *combs_slope_mix = new RooRealVar("exp_combs_slope","",relConfs.getD("exp_mix_combs_slope"),
                                                   relConfs.getD("exp_mix_combs_slope_LimL"), relConfs.getD("exp_mix_combs_slope_LimU") );

  RooRealVar *combs_slope_LL = new RooRealVar("exp_LL_combs_slope","",relConfs.getD("exp_LL_combs_slope"),
                                                 relConfs.getD("exp_LL_combs_slope_LimL"), relConfs.getD("exp_LL_combs_slope_LimU") );
  RooRealVar *combs_slope_DD = new RooRealVar("exp_DD_combs_slope","",relConfs.getD("exp_DD_combs_slope"),
                                                 relConfs.getD("exp_DD_combs_slope_LimL"), relConfs.getD("exp_DD_combs_slope_LimU") );

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
  floatParams->push_back(bu_mean);
  floatParams->push_back(bu_width_all);
  floatParams->push_back(bu_width_run1);
  floatParams->push_back(bu_width_run2);


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

  fixedParams->push_back(bu_alpha_all_LL);
  fixedParams->push_back(bu_alpha_run1_LL);
  fixedParams->push_back(bu_alpha_run2_LL);
  fixedParams->push_back(bu_n_LL);
  fixedParams->push_back(bu_alpha_all_DD);
  fixedParams->push_back(bu_alpha_run1_DD);
  fixedParams->push_back(bu_alpha_run2_DD);
  fixedParams->push_back(bu_n_DD);
  fixedParams->push_back(bu_width_ratio_all_LL);
  fixedParams->push_back(bu_width_ratio_all_DD);
  fixedParams->push_back(bu_width_ratio_run1_LL);
  fixedParams->push_back(bu_width_ratio_run1_DD);
  fixedParams->push_back(bu_width_ratio_run2_LL);
  fixedParams->push_back(bu_width_ratio_run2_DD);
  fixedParams->push_back(bu_frac_all_LL);
  fixedParams->push_back(bu_frac_all_DD);
  fixedParams->push_back(bu_frac_run1_LL);
  fixedParams->push_back(bu_frac_run1_DD);
  fixedParams->push_back(bu_frac_run2_LL);
  fixedParams->push_back(bu_frac_run2_DD);

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

          //Signal
          bu[*mode][*charge][*trackType][*run]->setMean(bu_mean);
          if(*run=="all") {
        	  bu[*mode][*charge][*trackType][*run]->setWidth(bu_width_all);
          }
          else if(*run=="run1") {
        	  bu[*mode][*charge][*trackType][*run]->setWidth(bu_width_run1);
          }
          else if(*run=="run2") {
        	  bu[*mode][*charge][*trackType][*run]->setWidth(bu_width_run2);
          }

          if(*trackType=="LL")  {

        	  if(*run=="all") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_all_LL);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_all_LL);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_all_LL);
        	  }
        	  else if(*run=="run1") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_run1_LL);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_run1_LL);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_run1_LL);
        	  }
        	  else if(*run=="run2") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_run2_LL);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_run2_LL);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_run2_LL);
        	  }

        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_LL);
        	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_LL); }
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_LL);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef010_LL);

          }
          else if(*trackType=="DD") {

        	  if(*run=="all") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_all_DD);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_all_DD);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_all_DD);
        	  }
        	  else if(*run=="run1") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_run1_DD);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_run1_DD);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_run1_DD);
        	  }
        	  else if(*run=="run2") {
        		  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_run2_DD);
        		  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_run2_DD);
        		  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_run2_DD);
        	  }

        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_DD);
        	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_DD); }
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_DD);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef101_DD);
          }
          else if(*trackType=="mix") {

        	  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_all_mix);
        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_mix);
        	  bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_all_mix);
        	  bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_all_LL);

        	  if(relConfs.get("combinatoricShape")=="0") { comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_mix); }
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_DD);
        	  dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef101_DD);
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
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Gen" << endl;
}

