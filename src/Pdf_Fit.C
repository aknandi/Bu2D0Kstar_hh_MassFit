#include "Pdf_Fit.h"
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
#include "RooGaussian.h"
#include "RooConstVar.h"

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

  // Initialise the PDFs
  // Note that if you uncomment that KeysPdfs but don't use them they use up unnecessary memory!
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){
          bu[*m][*c][*t][*a]  = new myCrystalBall(pmB, *m,"bu",*c,*t,*a,_fileList->get("fit_signal"));
          comb[*m][*c][*t][*a]   = new Exponential(pmB, *m,"exp",*c,*t,*a,_fileList->get("fit_combs"));
          bu_dstkst[*m][*c][*t][*a]    = new PartRecoDstKst(pmB, *m,"bu",*c,*t,*a,_fileList->get("fit_partreco"));
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
//  relConfs.readPairStringsToMap(_fileList->get("fit_drho"));
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));
  Settings pdfGenConfs("pdfGenConfs");
  pdfGenConfs.readPairStringsToMap(_fileList->get("gen_partreco"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)

  //Signal -- Double Crystal Ball
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

  RooRealVar *combs_slope_mix = new RooRealVar("d2kpi_exp_mix_combs_slope","",relConfs.getD("d2kpi_exp_mix_combs_slope"),
                                                   relConfs.getD("d2kpi_exp_mix_combs_slope_LimL"), relConfs.getD("d2kpi_exp_mix_combs_slope_LimU") );

  //  RooRealVar *combs_slope_LL = new RooRealVar("d2kspipi_exp_LL_combs_slope","",relConfs.getD("d2kspipi_exp_LL_combs_slope"),
  //                                                 relConfs.getD("d2kspipi_exp_LL_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_LL_combs_slope_LimU") );
  //  RooRealVar *combs_slope_DD = new RooRealVar("d2kspipi_exp_DD_combs_slope","",relConfs.getD("d2kspipi_exp_DD_combs_slope"),
  //                                                 relConfs.getD("d2kspipi_exp_DD_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_DD_combs_slope_LimU") );

  RooRealVar *bu_frac010 = new RooRealVar("frac010_bu","",relConfs.getD("frac010_bu"), 0.0, 1.0);

  //PartReco
  string limitlow = relConfs.get("fit_limit_low");

  // Ensure that when running toys we generate and fit the same frac010 parameter
/*  double frac010_val=0.;
  if(_genConfs->get("genToys")=="true") {
    frac010_val = pdfGenConfs.getD(Form("frac010_bs_%s",limitlow.c_str()));
  }
  else {
    frac010_val = relConfs.getD(Form("frac010_bs_%s",limitlow.c_str()));
  }
  bs_frac010 = new RooRealVar(Form("frac010_bs_%s",limitlow.c_str()),"",frac010_val, 0.0, 1.0);

  RooRealVar *bd_frac010 = new RooRealVar(Form("frac010_bd_%s",limitlow.c_str()),"",relConfs.getD(Form("frac010_bd_%s",limitlow.c_str())), 0.0, 1.0);
  //RooRealVar *bs_010_frac010 = new RooRealVar(Form("frac010_bs_010_%s",limitlow.c_str()),"",1.0);
  //RooRealVar *bs_001_frac010 = new RooRealVar(Form("frac010_bs_001_%s",limitlow.c_str()),"",0.0);


  // Initialise the PDFs for RooGaussian
  for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
    for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
      for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){

        gaus_frac010_bs[*c][*t][*a] = new RooGaussian(Form("gaus_frac010_bs_%s_%s_%s",(*c).c_str(),(*t).c_str(),(*a).c_str()),"",*bs_frac010,
                                                      RooFit::RooConst(frac010_val),
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bs_%s_err",limitlow.c_str()))));
        gaus_frac010_bd[*c][*t][*a] = new RooGaussian(Form("gaus_frac010_bd_%s_%s_%s",(*c).c_str(),(*t).c_str(),(*a).c_str()),"",*bd_frac010,
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bd_%s",limitlow.c_str()))),
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bd_%s_err",limitlow.c_str()))));
      }
    }
  }*/


  //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfFit: Setting parameters constant" << std::endl;
  fixedParams = new std::vector <RooRealVar*>;

  fixedParams->push_back(bu_alpha_LL);
  fixedParams->push_back(bu_n_LL);
  fixedParams->push_back(bu_alpha_DD);
  fixedParams->push_back(bu_n_DD);

 
  for (Int_t n_p = 0; n_p < (Int_t)fixedParams->size(); ++n_p)
    {
      RooRealVar *par = fixedParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
      par->setConstant(kTRUE);
    }

  std::cout << "... done" << std::endl;
  
  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
	
  // Now loop over the pdfs and set the shared variables which you've just defined
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){

          //Signal
          bu[*m][*c][*t][*a]->setMean(bu_mean);
          bu[*m][*c][*t][*a]->setWidth(bu_width);

          if(*t=="LL")  {
        	  bu[*m][*c][*t][*a]->setAlpha(bu_alpha_LL);
        	  bu[*m][*c][*t][*a]->setN(bu_n_LL);
          }
          else if(*t=="DD") {
        	  bu[*m][*c][*t][*a]->setAlpha(bu_alpha_DD);
        	  bu[*m][*c][*t][*a]->setN(bu_n_DD);
          }
          else if(*t=="mix") {
        	  bu[*m][*c][*t][*a]->setAlpha(bu_alpha_mix);
        	  bu[*m][*c][*t][*a]->setN(bu_n_mix);
          }

          comb[*m][*c][*t][*a]->setSlope(combs_slope_mix);

          bu_dstkst[*m][*c][*t][*a]->setFraction(bu_frac010);
        }
      }
    }
  }


  // Finally, create map of RooPdfs to return
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){
          roopdf_bu[*m][*c][*t][*a]     = bu[*m][*c][*t][*a]->getPdf();
          roopdf_comb[*m][*c][*t][*a]      = comb[*m][*c][*t][*a]->getPdf();
          roopdf_bu_dstkst[*m][*c][*t][*a]  = bu_dstkst[*m][*c][*t][*a]->getPdf();
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Fit" << endl;
}

