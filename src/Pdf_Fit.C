#include "Pdf_Fit.h"
//#include "DoubleCrystalBall.h"
//#include "KeysPdf.h"
//#include "PartRecoDstKst.h"
#include "Exponential.h"
//#include "DoubleExponential.h"
//#include "Linear.h"
#include "RooFormulaVar.h"
#include "myGaussian.h"
//#include "CorrGauss.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"
//#include "Lambda.h"
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
          bu[*m][*c][*t][*a]  = new myGaussian(pmB, *m,"bu",*c,*t,*a,_fileList->get("fit_signal"));

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
  //relConfs.readPairStringsToMap(_fileList->get("fit_combs"));
//  relConfs.readPairStringsToMap(_fileList->get("fit_drho"));
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));
  //Settings pdfGenConfs("pdfGenConfs");
  //pdfGenConfs.readPairStringsToMap(_fileList->get("gen_partreco"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)

  //Signal -- Double Crystal Ball
  RooRealVar* bu_mean = new RooRealVar("bu_mean","",relConfs.getD("bu_mean"),
                                              relConfs.getD("bu_mean_LimL"),relConfs.getD("bu_mean_LimU") );

  RooRealVar* bu_width = new RooRealVar("bu_width","",relConfs.getD("bu_width"),
                                              relConfs.getD("bu_width_LimL"),relConfs.getD("bu_width_LimU") );



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
//  fixedParams->push_back(bs_alpha1);
//  fixedParams->push_back(bs_alpha2);
//  fixedParams->push_back(bs_n1);
//  fixedParams->push_back(bs_n2);
//  fixedParams->push_back(bs_frac1);
//  fixedParams->push_back(bd_bs_shift);
//  fixedParams->push_back(drho_mean);
//  fixedParams->push_back(drho_mean2);
//  fixedParams->push_back(drho_width);
//  fixedParams->push_back(drho_width2);
//  fixedParams->push_back(drho_alpha1);
//  fixedParams->push_back(drho_alpha2);
//  fixedParams->push_back(drho_n1);
//  fixedParams->push_back(drho_n2);
//  fixedParams->push_back(drho_frac1);
  //fixedParams->push_back(bs_frac010);
  //fixedParams->push_back(bs_010_frac010);
  //fixedParams->push_back(bs_001_frac010);
  //fixedParams->push_back(bd_frac010);
 
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
          bu[*m][*c][*t][*a]->setRelation("mean",bu_mean);
          bu[*m][*c][*t][*a]->setRelation("width",bu_width);

//          //Combinatoric - floating Exponential
//          if(*m=="d2kspipi") {
//            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix);
//            //if(*t=="LL")  comb[*m][*c]["LL"][*a]->setRelation("slope",combs_slope_LL);
//            //if(*t=="DD")  comb[*m][*c]["DD"][*a]->setRelation("slope",combs_slope_DD);
//          }
//          else
//            comb[*m][*c][*t][*a]->setRelation("slope",combs_slope_mix_kskk);
//
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
          //roopdf_comb[*m][*c][*t][*a]      = comb[*m][*c][*t][*a]->getPdf();
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Fit" << endl;
}

