#include "Pdf_Gen.h"
//#include "DoubleCrystalBall.h"
//#include "KeysPdf.h"
//#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "RooFormulaVar.h"
#include "myGaussian.h"
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
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){
          bu[*m][*c][*t][*a]  = new myGaussian(pmB, *m,"bu",*c,*t,*a,_fileList->get("gen_signal"));

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
//  relConfs.readPairStringsToMap(_fileList->get("gen_combs"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
  
  //Signal -- Double Crystal Ball
  RooRealVar* bu_mean = new RooRealVar("bu_mean","",relConfs.getD("bu_mean"),
                                       relConfs.getD("bu_mean_LimL"),relConfs.getD("bu_mean_LimU") );
  RooRealVar* bu_width = new RooRealVar("bu_width","",relConfs.getD("bu_width"),
                                        relConfs.getD("bu_width_LimL"),relConfs.getD("bu_width_LimU") );


  std::cout << std::endl << "PdfGen : floating parameters " << std::endl;
  std::vector<RooRealVar*> *floatParams = new std::vector <RooRealVar*>;
  floatParams->push_back(bu_mean);
  floatParams->push_back(bu_width);/*
  floatParams->push_back(bs_frac010);
  floatParams->push_back(combs_slope_mix);
  floatParams->push_back(combs_slope_mix_kskk);*/

  for (Int_t n_p = 0; n_p < (Int_t)floatParams->size(); ++n_p)
    {
      RooRealVar *par = floatParams->at(n_p);
      std::cout << " " << par->GetName() << " " << par->getVal() << std::endl;
    }
 
   //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfGen : fixed parameters " << std::endl;
  fixedParams = new std::vector <RooRealVar*>;
/*
  fixedParams->push_back(bs_alpha1);
  fixedParams->push_back(bs_alpha2);
  fixedParams->push_back(bs_n1);
  fixedParams->push_back(bs_n2);
  fixedParams->push_back(bs_frac1);*/

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
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){

          //Signal
          bu[*m][*c][*t][*a]->setRelation("mean",bu_mean);
          bu[*m][*c][*t][*a]->setRelation("width",bu_width);

/*          //Comb
          if(*m=="d2kspipi") {
            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix);
            //if(*t=="LL")  comb[*m][*c]["LL"][*a]->setRelation("slope",combs_slope_LL);
            //if(*t=="DD")  comb[*m][*c]["DD"][*a]->setRelation("slope",combs_slope_DD);
          }
          else {
            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix_kskk);
          }*/

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
//          roopdf_comb[*m][*c][*t][*a]      = comb[*m][*c][*t][*a]->getPdf();

        }
      }
    }
  }
	
  cout << "Finished making Pdf_Gen" << endl;
}

