#include "Pdf_Gen.h"
//#include "KeysPdf.h"
//#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "myCrystalBall.h"
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
        	bu[*mode][*charge][*trackType][*run]  = new myCrystalBall(pmB, *mode,"bu",*charge,*trackType,*run,_fileList->get("gen_signal"));
        	comb[*mode][*charge][*trackType][*run]   = new Exponential(pmB, *mode,"exp",*charge,*trackType,*run,_fileList->get("gen_combs"));
        	dstkst[*mode][*charge][*trackType][*run]    = new PartRecoDstKst(pmB, *mode,*charge,*trackType,*run,_fileList->get("gen_partreco"));
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
  relConfs.readPairStringsToMap(_fileList->get("gen_signal"));
  relConfs.readPairStringsToMap(_fileList->get("gen_combs"));
  relConfs.readPairStringsToMap(_fileList->get("gen_partreco"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
  
  //Signal -- Crystal Ball
  RooRealVar* bu_mean = new RooRealVar("bu_mean","",relConfs.getD("bu_mean"),
                                       relConfs.getD("bu_mean_LimL"),relConfs.getD("bu_mean_LimU") );
  RooRealVar* bu_width = new RooRealVar("bu_width","",relConfs.getD("bu_width"),
                                        relConfs.getD("bu_width_LimL"),relConfs.getD("bu_width_LimU") );
  RooRealVar* bu_alpha_LL = new RooRealVar("bu_alpha_LL","",relConfs.getD("bu_alpha_LL"),
                                       relConfs.getD("bu_alpha_LL_LimL"),relConfs.getD("bu_alpha_LL_LimU") );
  RooRealVar* bu_n_LL = new RooRealVar("bu_n_LL","",relConfs.getD("bu_n_LL"),
                                        relConfs.getD("bu_n_LL_LimL"),relConfs.getD("bu_n_LL_LimU") );
  RooRealVar* bu_alpha_DD = new RooRealVar("bu_alpha_DD","",relConfs.getD("bu_alpha_DD"),
                                       relConfs.getD("bu_alpha_DD_LimL"),relConfs.getD("bu_alpha_DD_LimU") );
  RooRealVar* bu_n_DD = new RooRealVar("bu_n_DD","",relConfs.getD("bu_n_DD"),
                                        relConfs.getD("bu_n_DD_LimL"),relConfs.getD("bu_n_DD_LimU") );

  RooRealVar* bu_alpha_mix = new RooRealVar("bu_alpha_mix","",relConfs.getD("bu_alpha_mix"),
                                       relConfs.getD("bu_alpha_mix_LimL"),relConfs.getD("bu_alpha_mix_LimU") );
  RooRealVar* bu_n_mix = new RooRealVar("bu_n_mix","",relConfs.getD("bu_n_mix"),
                                        relConfs.getD("bu_n_mix_LimL"),relConfs.getD("bu_n_mix_LimU") );


  RooRealVar *combs_slope_mix = new RooRealVar("exp_combs_slope","",relConfs.getD("d2kpi_exp_mix_combs_slope"),
                                                   relConfs.getD("d2kpi_exp_mix_combs_slope_LimL"), relConfs.getD("d2kpi_exp_mix_combs_slope_LimU") );

  RooRealVar *combs_slope_LL = new RooRealVar("d2kpi_exp_LL_combs_slope","",relConfs.getD("d2kpi_exp_LL_combs_slope"),
                                                 relConfs.getD("d2kpi_exp_LL_combs_slope_LimL"), relConfs.getD("d2kpi_exp_LL_combs_slope_LimU") );
  RooRealVar *combs_slope_DD = new RooRealVar("d2kpi_exp_DD_combs_slope","",relConfs.getD("d2kpi_exp_DD_combs_slope"),
                                                 relConfs.getD("d2kpi_exp_DD_combs_slope_LimL"), relConfs.getD("d2kpi_exp_DD_combs_slope_LimU") );

  RooRealVar *bu_frac010 = new RooRealVar("frac010_bu","",relConfs.getD("frac010_bu"), 0.0, 1.0);
  RooRealVar *bd_frac010 = new RooRealVar("frac010_bd","",relConfs.getD("frac010_bd"), 0.0, 1.0);
  RooRealVar *lowmass_fracBd = new RooRealVar("frac_bd","",relConfs.getD("frac_bd"), 0.0, 1.0);

  std::cout << std::endl << "PdfGen : floating parameters " << std::endl;
  std::vector<RooRealVar*> *floatParams = new std::vector <RooRealVar*>;
  floatParams->push_back(bu_mean);
  floatParams->push_back(bu_width);
  floatParams->push_back(bu_frac010);

  for (Int_t n_p = 0; n_p < (Int_t)floatParams->size(); ++n_p)
    {
      RooRealVar *par = floatParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
    }
 
   //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfGen : fixed parameters " << std::endl;
  fixedParams = new std::vector <RooRealVar*>;
  fixedParams->push_back(bu_alpha_LL);
  fixedParams->push_back(bu_n_LL);
  fixedParams->push_back(bu_alpha_DD);
  fixedParams->push_back(bu_n_DD);
  //fixedParams->push_back(bu_alpha_mix);
  //fixedParams->push_back(bu_n_mix);
  fixedParams->push_back(combs_slope_LL);
  fixedParams->push_back(combs_slope_DD);


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
          bu[*mode][*charge][*trackType][*run]->setWidth(bu_width);

          if(*trackType=="LL")  {
        	  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_LL);
        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_LL);
        	  comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_LL);

          }
          else if(*trackType=="DD") {
        	  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_DD);
        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_DD);
        	  comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_DD);
          }
          else if(*trackType=="mix") {
        	  bu[*mode][*charge][*trackType][*run]->setAlpha(bu_alpha_mix);
        	  bu[*mode][*charge][*trackType][*run]->setN(bu_n_mix);
        	  comb[*mode][*charge][*trackType][*run]->setSlope(combs_slope_mix);
          }


          dstkst[*mode][*charge][*trackType][*run]->setFractionBu010(bu_frac010);
          dstkst[*mode][*charge][*trackType][*run]->setFractionBd010(bd_frac010);
          dstkst[*mode][*charge][*trackType][*run]->setFractionBd(lowmass_fracBd);
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

