#include "Pdf_Gen.h"
#include "DoubleCrystalBall.h"
#include "KeysPdf.h"
#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "CorrGauss.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"

Pdf_Gen::Pdf_Gen(Settings* fileList, RooRealVar* pmB, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackTypeList, std::vector<std::string> binList)
{
  _fileList=fileList;
  _modeList=modeList;
  _chargeList=chargeList;
  _trackTypeList=trackTypeList;
  _binList=binList;

  // Initialise the PDFs
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){ 
          bd[*m][*c][*t][*a]  = new DoubleCrystalBall(pmB, *m,"bd",*c,*t,*a,_fileList->get("gen_signal"));
          bs[*m][*c][*t][*a]  = new DoubleCrystalBall(pmB, *m,"bs",*c,*t,*a,_fileList->get("gen_signal"));
          comb[*m][*c][*t][*a]   = new Exponential(pmB, *m,"exp",*c,*t,*a,_fileList->get("gen_combs"));
          //drho[*m][*c][*t][*a]    = new KeysPdf(pmB, *m,"drho",*c,*t,*a,_fileList->get("gen_drho"));
          drho[*m][*c][*t][*a]    = new DoubleCrystalBall(pmB, *m,"drho",*c,*t,*a,_fileList->get("gen_drho"));
          bs_dstkst[*m][*c][*t][*a]    = new PartRecoDstKst(pmB, *m,"bs",*c,*t,*a,_fileList->get("gen_partreco"));
          bd_dstkst[*m][*c][*t][*a]    = new PartRecoDstKst(pmB, *m,"bd",*c,*t,*a,_fileList->get("gen_partreco"));
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
  relConfs.readPairStringsToMap(_fileList->get("gen_partreco"));
  relConfs.readPairStringsToMap(_fileList->get("gen_combs"));
  relConfs.readPairStringsToMap(_fileList->get("gen_drho"));
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
  
  //Signal -- Double Crystal Ball
  RooRealVar* bs_mean = new RooRealVar("bs_mean","",relConfs.getD("bs_mean"), 
                                       relConfs.getD("bs_mean_LimL"),relConfs.getD("bs_mean_LimU") );
  RooRealVar* bs_width = new RooRealVar("bs_width","",relConfs.getD("bs_width"), 
                                        relConfs.getD("bs_width_LimL"),relConfs.getD("bs_width_LimU") );
  RooRealVar* bs_minus = new RooRealVar("bs_minus","",-1);
  RooFormulaVar* bs_width2 = new RooFormulaVar("bs_width2","@0*@1",RooArgList(*bs_width,*bs_minus));
  RooRealVar *bs_alpha1 = new RooRealVar("bs_alpha1","",relConfs.getD("bs_alpha1"),
                                            relConfs.getD("bs_alpha1_LimL"), relConfs.getD("bs_alpha1_LimU") );
  RooRealVar *bs_alpha2 = new RooRealVar("bs_alpha2","",relConfs.getD("bs_alpha2"),
                                            relConfs.getD("bs_alpha2_LimL"), relConfs.getD("bs_alpha2_LimU") );
  RooRealVar *bs_n1 = new RooRealVar("bs_n1","",relConfs.getD("bs_n1"),
                                            relConfs.getD("bs_n1_LimL"), relConfs.getD("bs_n1_LimU") );
  RooRealVar *bs_n2 = new RooRealVar("bs_n2","",relConfs.getD("bs_n2"),
                                            relConfs.getD("bs_n2_LimL"), relConfs.getD("bs_n2_LimU") );
  RooRealVar *bs_frac1 = new RooRealVar("bs_frac1","",relConfs.getD("bs_frac1"),
                                            relConfs.getD("bs_frac1_LimL"), relConfs.getD("bs_frac1_LimU") );

  RooRealVar *bd_bs_shift = new RooRealVar("bd_bs_shift","",relConfs.getD("bd_bs_shift"));
  RooFormulaVar *bd_mean = new RooFormulaVar("bd_mean","@0-@1",RooArgSet(*bs_mean,*bd_bs_shift));

  //Combinatoric
  RooRealVar *combs_slope_mix = new RooRealVar("d2kspipi_exp_mix_combs_slope","",relConfs.getD("d2kspipi_exp_mix_combs_slope"),
                                                 relConfs.getD("d2kspipi_exp_mix_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_mix_combs_slope_LimU") );
  RooRealVar *combs_slope_mix_kskk = new RooRealVar("d2kskk_exp_mix_combs_slope","",relConfs.getD("d2kskk_exp_mix_combs_slope"),
                                                 relConfs.getD("d2kskk_exp_mix_combs_slope_LimL"), relConfs.getD("d2kskk_exp_mix_combs_slope_LimU") );

//  RooRealVar *combs_slope_LL = new RooRealVar("d2kspipi_exp_LL_combs_slope","",relConfs.getD("d2kspipi_exp_LL_combs_slope"),
//                                                 relConfs.getD("d2kspipi_exp_LL_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_LL_combs_slope_LimU") );
//  RooRealVar *combs_slope_DD = new RooRealVar("d2kspipi_exp_DD_combs_slope","",relConfs.getD("d2kspipi_exp_DD_combs_slope"),
//                                                 relConfs.getD("d2kspipi_exp_DD_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_DD_combs_slope_LimU") );

  //PartReco
  string limitlow = relConfs.get("fit_limit_low");
  RooRealVar *bs_frac010 = new RooRealVar(Form("frac010_bs_%s",limitlow.c_str()),"",relConfs.getD(Form("frac010_bs_%s",limitlow.c_str())), 0.0, 1.0);

  // Drho -- Double Crystal Ball
  RooRealVar* drho_mean = new RooRealVar("drho_mean","",relConfs.getD("drho_mean"),
                                         relConfs.getD("drho_mean_LimL"),relConfs.getD("drho_mean_LimU") );
  RooRealVar* drho_mean2 = new RooRealVar("drho_mean2","",relConfs.getD("drho_mean2"),
                                          relConfs.getD("drho_mean2_LimL"),relConfs.getD("drho_mean2_LimU") );
  RooRealVar* drho_width = new RooRealVar("drho_width","",relConfs.getD("drho_width"),
                                          relConfs.getD("drho_width_LimL"),relConfs.getD("drho_width_LimU") );
  RooRealVar* drho_width2 = new RooRealVar("drho_width2","",relConfs.getD("drho_width2"),
                                           relConfs.getD("drho_width2_LimL"),relConfs.getD("drho_width2_LimU") );
  RooRealVar *drho_alpha1 = new RooRealVar("drho_alpha1","",relConfs.getD("drho_alpha1"),
                                           relConfs.getD("drho_alpha1_LimL"), relConfs.getD("drho_alpha1_LimU") );
  RooRealVar *drho_alpha2 = new RooRealVar("drho_alpha2","",relConfs.getD("drho_alpha2"),
                                           relConfs.getD("drho_alpha2_LimL"), relConfs.getD("drho_alpha2_LimU") );
  RooRealVar *drho_n1 = new RooRealVar("drho_n1","",relConfs.getD("drho_n1"),
                                       relConfs.getD("drho_n1_LimL"), relConfs.getD("drho_n1_LimU") );
  RooRealVar *drho_n2 = new RooRealVar("drho_n2","",relConfs.getD("drho_n2"),
                                       relConfs.getD("drho_n2_LimL"), relConfs.getD("drho_n2_LimU") );
  RooRealVar *drho_frac1 = new RooRealVar("drho_frac1","",relConfs.getD("drho_frac1"),
                                          relConfs.getD("drho_frac1_LimL"), relConfs.getD("drho_frac1_LimU") );

  std::cout << std::endl << "PdfGen : floating parameters " << std::endl;
  std::vector<RooRealVar*> *floatParams = new std::vector <RooRealVar*>;
  floatParams->push_back(bs_mean);
  floatParams->push_back(bs_width);
  floatParams->push_back(bs_frac010);
  floatParams->push_back(combs_slope_mix);
  floatParams->push_back(combs_slope_mix_kskk);

  for (Int_t n_p = 0; n_p < (Int_t)floatParams->size(); ++n_p)
    {
      RooRealVar *par = floatParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
    }
 
   //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfGen : fixed parameters " << std::endl;
  fixedParams = new std::vector <RooRealVar*>;
  fixedParams->push_back(bs_alpha1);
  fixedParams->push_back(bs_alpha2);
  fixedParams->push_back(bs_n1);
  fixedParams->push_back(bs_n2);
  fixedParams->push_back(bs_frac1);
  fixedParams->push_back(bd_bs_shift);
  fixedParams->push_back(drho_mean);
  fixedParams->push_back(drho_mean2);
  fixedParams->push_back(drho_width);
  fixedParams->push_back(drho_width2);
  fixedParams->push_back(drho_alpha1);
  fixedParams->push_back(drho_alpha2);
  fixedParams->push_back(drho_n1);
  fixedParams->push_back(drho_n2);
  fixedParams->push_back(drho_frac1);

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
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){

          //Signal
          bd[*m][*c][*t][*a]->setRelation("mean1",bd_mean);
          bd[*m][*c][*t][*a]->setRelation("mean2",bd_mean);
          bd[*m][*c][*t][*a]->setRelation("width1",bs_width);
          bd[*m][*c][*t][*a]->setRelation("width2",bs_width2);
          bd[*m][*c][*t][*a]->setRelation("alpha1",bs_alpha1);
          bd[*m][*c][*t][*a]->setRelation("alpha2",bs_alpha2);
          bd[*m][*c][*t][*a]->setRelation("n1",bs_n1);
          bd[*m][*c][*t][*a]->setRelation("n2",bs_n2);
          bd[*m][*c][*t][*a]->setRelation("frac1",bs_frac1);

          bs[*m][*c][*t][*a]->setRelation("mean1",bs_mean);
          bs[*m][*c][*t][*a]->setRelation("mean2",bs_mean);
          bs[*m][*c][*t][*a]->setRelation("width1",bs_width);
          bs[*m][*c][*t][*a]->setRelation("width2",bs_width2);
          bs[*m][*c][*t][*a]->setRelation("alpha1",bs_alpha1);
          bs[*m][*c][*t][*a]->setRelation("alpha2",bs_alpha2);
          bs[*m][*c][*t][*a]->setRelation("n1",bs_n1);
          bs[*m][*c][*t][*a]->setRelation("n2",bs_n2);
          bs[*m][*c][*t][*a]->setRelation("frac1",bs_frac1);
 
          //Drho (Double Crystal Ball)
          drho[*m][*c][*t][*a]->setRelation("mean1",drho_mean);
          drho[*m][*c][*t][*a]->setRelation("mean2",drho_mean2);
          drho[*m][*c][*t][*a]->setRelation("width1",drho_width);
          drho[*m][*c][*t][*a]->setRelation("width2",drho_width2);
          drho[*m][*c][*t][*a]->setRelation("alpha1",drho_alpha1);
          drho[*m][*c][*t][*a]->setRelation("alpha2",drho_alpha2);
          drho[*m][*c][*t][*a]->setRelation("n1",drho_n1);
          drho[*m][*c][*t][*a]->setRelation("n2",drho_n2);
          drho[*m][*c][*t][*a]->setRelation("frac1",drho_frac1);
 
          //PartReco
          bs_dstkst[*m][*c][*t][*a]->setRelation("frac010",bs_frac010);
          bd_dstkst[*m][*c][*t][*a]->setRelation("frac010",bs_frac010);

          //Comb
          if(*m=="d2kspipi") {
            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix);
            //if(*t=="LL")  comb[*m][*c]["LL"][*a]->setRelation("slope",combs_slope_LL);
            //if(*t=="DD")  comb[*m][*c]["DD"][*a]->setRelation("slope",combs_slope_DD);
          }
          else {
            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix_kskk);
          }

        }
      }
    }
  }


  // Finally, create map of RooPdfs to return
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){ 
          roopdf_bs[*m][*c][*t][*a]     = bs[*m][*c][*t][*a]->getPdf();
          roopdf_bd[*m][*c][*t][*a]     = bd[*m][*c][*t][*a]->getPdf();
          roopdf_comb[*m][*c][*t][*a]      = comb[*m][*c][*t][*a]->getPdf();
          roopdf_drho[*m][*c][*t][*a]       = drho[*m][*c][*t][*a]->getPdf();
          roopdf_bs_dstkst[*m][*c][*t][*a]  = bs_dstkst[*m][*c][*t][*a]->getPdf();
          roopdf_bd_dstkst[*m][*c][*t][*a]  = bd_dstkst[*m][*c][*t][*a]->getPdf();

        }
      }
    }
  }
	
  cout << "Finished making Pdf_Gen" << endl;
}

