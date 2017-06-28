#include "Pdf_Fit.h"
//#include "KeysPdf.h"
//#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "DoubleCrystalBall.h"
#include "myGaussian.h"
#include "PartRecoDstKst.h"
//#include "RooKeysPdf.h"
#include "myCruijff.h"
//#include "CorrGauss.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"

Pdf_Fit::Pdf_Fit(Settings* fileList, Settings* genConfs, RooRealVar* pmB, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackTypeList, std::vector<std::string> runList, int systematicFactor, std::string MCsimfit)
{
  _fileList=fileList;
  _genConfs=genConfs;
  _MCsimfit=MCsimfit;
  _modeList=modeList;
  _chargeList=chargeList;
  _trackTypeList=trackTypeList;
  _runList=runList;
  //_shiftFactor = systematicFactor;

  /*
  TFile LambdaFile("Settings/PDFShapes/Gen/lckst.root");
  RooWorkspace* workspace = (RooWorkspace*)LambdaFile.Get("workspace");
  RooKeysPdf* lckst_all = (RooKeysPdf*)workspace->pdf("lckst");
  LambdaFile.Close();
*/
  // Initialise the PDFs
  // Note that if you uncomment that KeysPdfs but don't use them they use up unnecessary memory!
  for(std::vector<std::string>::iterator mode=_modeList.begin();mode!=_modeList.end();mode++){
    for(std::vector<std::string>::iterator charge=_chargeList.begin();charge!=_chargeList.end();charge++){
      for(std::vector<std::string>::iterator trackType=_trackTypeList.begin();trackType!=_trackTypeList.end();trackType++){
        for(std::vector<std::string>::iterator run=_runList.begin();run!=_runList.end();run++){
          bu[*mode][*charge][*trackType][*run]  = new DoubleCrystalBall(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("fit_signal"));
          //bu[*mode][*charge][*trackType][*run]  = new myGaussian(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("fit_signal"));
          comb[*mode][*charge][*trackType][*run]   = new Exponential(pmB, *mode,"exp",*charge,*trackType,*run,_fileList->get("fit_combs"));
          dstkst[*mode][*charge][*trackType][*run]    = new PartRecoDstKst(pmB, *mode,*charge,*trackType,*run,_fileList->get("fit_partreco"),false);
          //if(*mode=="d2kk") lckst[*mode][*charge][*trackType][*run] = new RooKeysPdf(*lckst_all,Form("lckst_%s_%s_%s",(*charge).c_str(),(*trackType).c_str(),(*run).c_str()));
          if(*mode=="d2kk") lckst[*mode][*charge][*trackType][*run] = new myCruijff(pmB,*mode,"bu",*charge,*trackType,*run,_fileList->get("fit_signal"));
        }
      }
    }
  }
  
  // Now set relations, which prepares the PDFs to be returned
  setRelations();
}

void Pdf_Fit::setRelations()
{
  // Set up the configuration files
  Settings relConfs("Pdf_Fit::SetRelations");
  relConfs.readPairStringsToMap(_fileList->get("fit_signal"));
  relConfs.readPairStringsToMap(_fileList->get("fit_combs"));
  relConfs.readPairStringsToMap(_fileList->get("fit_partreco"));
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));
  relConfs.readPairStringsToMap(_fileList->get("PathnameToTotals"));
  relConfs.readPairStringsToMap(_fileList->get("PathnameToYieldCorrections"));
  Settings pdfGenConfs("pdfGenConfs");
  pdfGenConfs.readPairStringsToMap(_fileList->get("gen_partreco"));


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
  RooRealVar* bu_width_ratio_kpi_LL = new RooRealVar("bu_width_ratio_kpi_LL","",relConfs.getD("bu_width_ratio_kpi_LL") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_width_ratio_kpi_LL_err")):0.));
  RooRealVar* bu_width_ratio_kpi_DD = new RooRealVar("bu_width_ratio_kpi_DD","",relConfs.getD("bu_width_ratio_kpi_DD") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_width_ratio_kpi_DD_err")):0.));
  RooRealVar* bu_width_ratio_kpipipi_LL = new RooRealVar("bu_width_ratio_kpipipi_LL","",relConfs.getD("bu_width_ratio_kpipipi_LL") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_width_ratio_kpipipi_LL_err")):0.));
  RooRealVar* bu_width_ratio_kpipipi_DD = new RooRealVar("bu_width_ratio_kpipipi_DD","",relConfs.getD("bu_width_ratio_kpipipi_DD") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_width_ratio_kpipipi_DD_err")):0.));
  RooRealVar* bu_frac_kpi_LL = new RooRealVar("bu_frac_kpi_LL","",relConfs.getD("bu_frac_kpi_LL") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_frac_kpi_LL_err")):0.));
  RooRealVar* bu_frac_kpi_DD = new RooRealVar("bu_frac_kpi_DD","",relConfs.getD("bu_frac_kpi_DD") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_frac_kpi_DD_err")):0.));
  RooRealVar* bu_frac_kpipipi_LL = new RooRealVar("bu_frac_kpipipi_LL","",relConfs.getD("bu_frac_kpipipi_LL") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_frac_kpipipi_LL_err")):0.));
  RooRealVar* bu_frac_kpipipi_DD = new RooRealVar("bu_frac_kpipipi_DD","",relConfs.getD("bu_frac_kpipipi_DD") + (_genConfs->get("signalShapeData")=="1"?gRandom->Gaus(0,relConfs.getD("bu_frac_kpipipi_DD_err")):0.));

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
  RooRealVar *lambda_mean = new RooRealVar("lambda_mean","",5269  + (_genConfs->get("lckst")=="1"?(gRandom->Gaus(0,18)):0.));
  RooRealVar *lambda_sigmaR = new RooRealVar("lambda_sigmaR","",221  + (_genConfs->get("lckst")=="1"?(gRandom->Gaus(0,26)):0.));
  RooRealVar *lambda_sigmaL = new RooRealVar("lambda_sigmaL","",96  + (_genConfs->get("lckst")=="1"?(gRandom->Gaus(0,16)):0.));
  RooRealVar *lambda_alphaR = new RooRealVar("lambda_alphaR","",-0.19  + (_genConfs->get("lckst")=="1"?(gRandom->Gaus(0,0.19)):0.));
  RooRealVar *lambda_alphaL = new RooRealVar("lambda_alphaL","",-0.04  + (_genConfs->get("lckst")=="1"?(gRandom->Gaus(0,0.06)):0.));

  //Get Ks helicity angle selection from general settings
  std::string kshelcut = relConfs.get("Kshelcut");
  std::string lowRange = relConfs.get("fit_limit_low");
  // frac010 = (n010*eff010)/(n101*eff101), different for DD and LL
  double frac010LL, frac010DD;

  if(_genConfs->get("genToys")=="true") {
	  frac010LL = (relConfs.getD(Form("N_dstkst010_d2kpi_LL_%s",lowRange.c_str()))*relConfs.getD(Form("eff010_LL_%s",kshelcut.c_str())))/(relConfs.getD(Form("N_dstkst101_d2kpi_LL_%s",lowRange.c_str()))*relConfs.getD(Form("eff101_LL_%s",kshelcut.c_str())));
	  frac010DD = (relConfs.getD(Form("N_dstkst010_d2kpi_DD_%s",lowRange.c_str()))*relConfs.getD(Form("eff010_DD_%s",kshelcut.c_str())))/(relConfs.getD(Form("N_dstkst101_d2kpi_DD_%s",lowRange.c_str()))*relConfs.getD(Form("eff101_DD_%s",kshelcut.c_str())));
  }
  else {
	  frac010LL = relConfs.getD("frac010_LL");
	  frac010DD = relConfs.getD("frac010_DD");
  }

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

  //PartReco
  string limitlow = relConfs.get("fit_limit_low");

  //set certain parameters constant
  //this list has to be maintained manually
  std::cout << std::endl << "PdfFit: Setting parameters constant" << std::endl;
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

  // todo set parameters in low mass shape constant

  for (Int_t n_p = 0; n_p < (Int_t)fixedParams->size(); ++n_p)
    {
      RooRealVar *par = fixedParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
      par->setConstant(kTRUE);
    }

  std::cout << "... done" << std::endl;
  
  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
	
  // Now loop over the pdfs and set the shared variables which you've just defined
  for(std::vector<std::string>::iterator mode=_modeList.begin();mode!=_modeList.end();mode++){
    for(std::vector<std::string>::iterator charge=_chargeList.begin();charge!=_chargeList.end();charge++){
      for(std::vector<std::string>::iterator trackType=_trackTypeList.begin();trackType!=_trackTypeList.end();trackType++){
        for(std::vector<std::string>::iterator run=_runList.begin();run!=_runList.end();run++){

        	//Signal
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
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_LL);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_LL);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_all_LL);
        			bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_LL);
        			bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_LL);
        		}
        		else {
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_LL);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_LL);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_run1_LL);
        			bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpi_LL);
        			bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpi_LL);
        		}

        		bu[*mode][*charge][*trackType][*run]->setN(bu_n_LL);
        		dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_LL);
        		dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef010_LL);

        	}
        	else if(*trackType=="DD") {

        		if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_DD);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_DD);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_all_DD);
        			bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_DD);
        			bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_DD);
        		}
        		else {
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_DD);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_DD);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_run1_DD);
        			bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpi_DD);
        			bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpi_DD);
        		}

        		bu[*mode][*charge][*trackType][*run]->setN(bu_n_DD);
        		dstkst[*mode][*charge][*trackType][*run]->setCoef010(coef010_DD);
        		dstkst[*mode][*charge][*trackType][*run]->setCoef101(coef101_DD);
        	}
        	else if(*trackType=="mix") {

        		if(*mode=="d2kpipipi" || *mode=="d2pikpipi" || *mode=="d2pipipipi") {
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpipipi_mix);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpipipi_mix);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_all_mix);
        			bu[*mode][*charge][*trackType][*run]->setWidthRatio(bu_width_ratio_kpipipi_mix);
        			bu[*mode][*charge][*trackType][*run]->setFrac(bu_frac_kpipipi_mix);
        		}
        		else {
        			comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_kpi_mix);
        			bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_kpi_mix);
        			//bu[*mode][*charge][*trackType][*run]->setDelta(bu_delta_run1_mix);
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
        	std::cout << "Start getting pdfs" << std::endl;
          roopdf_bu[*mode][*charge][*trackType][*run]     = bu[*mode][*charge][*trackType][*run]->getPdf();
          roopdf_comb[*mode][*charge][*trackType][*run]   = comb[*mode][*charge][*trackType][*run]->getPdf();
          roopdf_dstkst[*mode][*charge][*trackType][*run] = dstkst[*mode][*charge][*trackType][*run]->getPdf();
          if(*mode=="d2kk") roopdf_lckst[*mode][*charge][*trackType][*run] = lckst[*mode][*charge][*trackType][*run]->getPdf();
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Fit" << endl;
}

