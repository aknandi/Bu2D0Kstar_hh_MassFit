#include "Pdf_Fit.h"
#include "DoubleCrystalBall.h"
#include "KeysPdf.h"
#include "PartRecoDstKst.h"
#include "Exponential.h"
#include "DoubleExponential.h"
#include "Linear.h"
#include "RooFormulaVar.h"
#include "CorrGauss.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "Lambda.h"
#include "RooGaussian.h"
#include "RooConstVar.h"

Pdf_Fit::Pdf_Fit(Settings* fileList, Settings* genConfs, RooRealVar* pmB, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackTypeList, std::vector<std::string> binList, int systematicFactor, std::string MCsimfit)
{
  _fileList=fileList;
  _genConfs=genConfs;
  _MCsimfit=MCsimfit;
  _modeList=modeList;
  _chargeList=chargeList;
  _trackTypeList=trackTypeList;
  _binList=binList;
  _shiftFactor = systematicFactor;

  // Initialise the PDFs
  // Note that if you uncomment that KeysPdfs but don't use them they use up unnecessary memory!
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackTypeList.begin();t!=_trackTypeList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){ 
          bd[*m][*c][*t][*a]  = new DoubleCrystalBall(pmB, *m,"bd",*c,*t,*a,_fileList->get("fit_signal"));
          bs[*m][*c][*t][*a]  = new DoubleCrystalBall(pmB, *m,"bs",*c,*t,*a,_fileList->get("fit_signal"));
          comb[*m][*c][*t][*a]   = new Exponential(pmB, *m,"exp",*c,*t,*a,_fileList->get("fit_combs")); 
          //comb[*m][*c][*t][*a]   = new Linear(pmB, *m,"combs",*c,*t,*a,_fileList->get("fit_combs")); 
          //comb[*m][*c][*t][*a]   = new DoubleExponential(pmB, *m,"combs",*c,*t,*a,_fileList->get("fit_combs")); 
          //drho[*m][*c][*t][*a]   = new KeysPdf(pmB, *m,"drho",*c,*t,*a,_fileList->get("fit_drho"));
          drho[*m][*c][*t][*a]   = new DoubleCrystalBall(pmB, *m,"drho",*c,*t,*a,_fileList->get("fit_drho"));
          bs_dstkst[*m][*c][*t][*a]    = new PartRecoDstKst(pmB, *m,"bs",*c,*t,*a,_fileList->get("fit_partreco"));
          //bs_dstkst_010[*m][*c][*t][*a]= new PartRecoDstKst(pmB, *m,"bs_010",*c,*t,*a,_fileList->get("fit_partreco_010"));
          //bs_dstkst_001[*m][*c][*t][*a]= new PartRecoDstKst(pmB, *m,"bs_001",*c,*t,*a,_fileList->get("fit_partreco_001"));
          bd_dstkst[*m][*c][*t][*a]    = new PartRecoDstKst(pmB, *m,"bd",*c,*t,*a,_fileList->get("fit_partreco"));
          //dkpipi[*m][*c][*t][*a]   = new KeysPdf(pmB, *m,"dkpipi",*c,*t,*a,_fileList->get("fit_d3h"));
          //dpipipi[*m][*c][*t][*a]   = new KeysPdf(pmB, *m,"dpipipi",*c,*t,*a,_fileList->get("fit_d3h"));
          //lambda[*m][*c][*t][*a]   = new Lambda(pmB, *m,"lambda",*c,*t,*a,_fileList->get("fit_Lb"));
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
  relConfs.readPairStringsToMap(_fileList->get("fit_partreco"));
  relConfs.readPairStringsToMap(_fileList->get("fit_combs"));
  relConfs.readPairStringsToMap(_fileList->get("fit_drho"));
  relConfs.readPairStringsToMap(_fileList->get("gensettings"));
  Settings pdfGenConfs("pdfGenConfs");
  pdfGenConfs.readPairStringsToMap(_fileList->get("gen_partreco"));


  // Need to have separate sets of related parameters for different modes (d2kspipi) or track types (LL/DD)
  
  // Get correlations matrix
  vector<double>ranShifts;

  //NB - this is out of date, but is not used here - it's used in the binned fit instead
  if(_shiftFactor==1)
    {
      std::cout << "Warning - this is out of date" << std::endl;
      assert(0);
      cout << "**********************************************************************" << endl;
      cout << "* Getting correlation matrix and doing Cholesky decomposition        *" << endl;
      TFile*	f_summedMassFitResult = new TFile(_fileList->get("fit_d2kspipi_shapesFitResults").c_str(),"read");
      TMatrixD* corr = (TMatrixD*)f_summedMassFitResult->Get("CorrelationMatrix");
      CorrGauss *cg = new CorrGauss(*corr, (int)gRandom->GetSeed());
      cg->DoCholesky();
      cg->GetRandomParSet(ranShifts,true);
      TVectorD *fitErrors = (TVectorD*)f_summedMassFitResult->Get("FitParameterErrors");
      f_summedMassFitResult->Close();
      cout << "**********************************************************************" << endl;
      for (std::vector<double>::iterator it_shifts=ranShifts.begin();it_shifts!=ranShifts.end();it_shifts++){
        *it_shifts *= (*fitErrors)( (int)(it_shifts - ranShifts.begin()) );
      }
    }	
  else
    { // If not shifting the shapes, then set 100 shifts to be zero so that the following doesn't add anything
      int i=0;
      do
        {
          ranShifts.push_back(0);
          i++;
        }
      while (i<100);
    }

  //Signal -- Double Crystal Ball
  RooRealVar* bs_mean = new RooRealVar("bs_mean","",relConfs.getD("bs_mean") + ((_shiftFactor==1)?ranShifts.at(0):0),
                                              relConfs.getD("bs_mean_LimL"),relConfs.getD("bs_mean_LimU") );

  RooRealVar* bs_width = new RooRealVar("bs_width","",relConfs.getD("bs_width") + ((_shiftFactor==1)?ranShifts.at(1):0),
                                              relConfs.getD("bs_width_LimL"),relConfs.getD("bs_width_LimU") );

  RooRealVar* bs_minus = new RooRealVar("bs_minus","",-1);

  RooFormulaVar* bs_width2 = new RooFormulaVar("bs_width2","@0*@1",RooArgList(*bs_width,*bs_minus));

  RooRealVar *bs_alpha1 = new RooRealVar("bs_alpha1","",relConfs.getD("bs_alpha1")+((_shiftFactor==1)?ranShifts.at(2):0),
                                            relConfs.getD("bs_alpha1_LimL"), relConfs.getD("bs_alpha1_LimU") );
  RooRealVar *bs_alpha2 = new RooRealVar("bs_alpha2","",relConfs.getD("bs_alpha2")+((_shiftFactor==1)?ranShifts.at(3):0),
                                            relConfs.getD("bs_alpha2_LimL"), relConfs.getD("bs_alpha2_LimU") );
  RooRealVar *bs_n1 = new RooRealVar("bs_n1","",relConfs.getD("bs_n1")+((_shiftFactor==1)?ranShifts.at(4):0),
                                            relConfs.getD("bs_n1_LimL"), relConfs.getD("bs_n1_LimU") );
  RooRealVar *bs_n2 = new RooRealVar("bs_n2","",relConfs.getD("bs_n2")+((_shiftFactor==1)?ranShifts.at(5):0),
                                            relConfs.getD("bs_n2_LimL"), relConfs.getD("bs_n2_LimU") );
  RooRealVar *bs_frac1 = new RooRealVar("bs_frac1","",relConfs.getD("bs_frac1")+((_shiftFactor==1)?ranShifts.at(6):0),
                                            relConfs.getD("bs_frac1_LimL"), relConfs.getD("bs_frac1_LimU") );

  RooRealVar *bd_bs_shift = new RooRealVar("bd_bs_shift","",relConfs.getD("bd_bs_shift"));
  RooFormulaVar *bd_mean = new RooFormulaVar("bd_mean","@0-@1",RooArgSet(*bs_mean,*bd_bs_shift));

  //PartReco
  string limitlow = relConfs.get("fit_limit_low");

  // Ensure that when running toys we generate and fit the same frac010 parameter
  double frac010_val=0.;
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
      for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){ 

        gaus_frac010_bs[*c][*t][*a] = new RooGaussian(Form("gaus_frac010_bs_%s_%s_%s",(*c).c_str(),(*t).c_str(),(*a).c_str()),"",*bs_frac010,
                                                      RooFit::RooConst(frac010_val),
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bs_%s_err",limitlow.c_str()))));
        gaus_frac010_bd[*c][*t][*a] = new RooGaussian(Form("gaus_frac010_bd_%s_%s_%s",(*c).c_str(),(*t).c_str(),(*a).c_str()),"",*bd_frac010,
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bd_%s",limitlow.c_str()))),
                                                      RooFit::RooConst(relConfs.getD(Form("frac010_bd_%s_err",limitlow.c_str()))));
      }
    }
  }

  // --- Combinatorial ---
  
  // Split LL, DD
  //RooRealVar *combs_slope_LL = new RooRealVar("d2kspipi_exp_LL_combs_slope","",relConfs.getD("d2kspipi_exp_LL_combs_slope")+((_shiftFactor==1)?ranShifts.at(8):0),
  //                                               relConfs.getD("d2kspipi_exp_LL_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_LL_combs_slope_LimU") );
  //RooRealVar *combs_slope_DD = new RooRealVar("d2kspipi_exp_DD_combs_slope","",relConfs.getD("d2kspipi_exp_DD_combs_slope")+((_shiftFactor==1)?ranShifts.at(9):0),
  //                                               relConfs.getD("d2kspipi_exp_DD_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_DD_combs_slope_LimU") );
  
  RooRealVar *combs_slope_mix = new RooRealVar("d2kspipi_exp_mix_combs_slope","",relConfs.getD("d2kspipi_exp_mix_combs_slope")+((_shiftFactor==1)?ranShifts.at(8):0),
                                                 relConfs.getD("d2kspipi_exp_mix_combs_slope_LimL"), relConfs.getD("d2kspipi_exp_mix_combs_slope_LimU") );
  RooRealVar *combs_slope_mix_kskk = new RooRealVar("d2kskk_exp_mix_combs_slope","",relConfs.getD("d2kskk_exp_mix_combs_slope")+((_shiftFactor==1)?ranShifts.at(8):0),
                                                 relConfs.getD("d2kskk_exp_mix_combs_slope_LimL"), relConfs.getD("d2kskk_exp_mix_combs_slope_LimU") );

  // Drho -- Double Crystal Ball
  RooRealVar* drho_mean = new RooRealVar("drho_mean","",relConfs.getD("drho_mean") + ((_shiftFactor==1)?ranShifts.at(0):0),
                                              relConfs.getD("drho_mean_LimL"),relConfs.getD("drho_mean_LimU") );
  RooRealVar* drho_mean2 = new RooRealVar("drho_mean2","",relConfs.getD("drho_mean2") + ((_shiftFactor==1)?ranShifts.at(0):0),
                                              relConfs.getD("drho_mean2_LimL"),relConfs.getD("drho_mean2_LimU") );
  RooRealVar* drho_width = new RooRealVar("drho_width","",relConfs.getD("drho_width") + ((_shiftFactor==1)?ranShifts.at(1):0),
                                              relConfs.getD("drho_width_LimL"),relConfs.getD("drho_width_LimU") );
  RooRealVar* drho_width2 = new RooRealVar("drho_width2","",relConfs.getD("drho_width2") + ((_shiftFactor==1)?ranShifts.at(1):0),
                                              relConfs.getD("drho_width2_LimL"),relConfs.getD("drho_width2_LimU") );
  RooRealVar *drho_alpha1 = new RooRealVar("drho_alpha1","",relConfs.getD("drho_alpha1")+((_shiftFactor==1)?ranShifts.at(2):0),
                                            relConfs.getD("drho_alpha1_LimL"), relConfs.getD("drho_alpha1_LimU") );
  RooRealVar *drho_alpha2 = new RooRealVar("drho_alpha2","",relConfs.getD("drho_alpha2")+((_shiftFactor==1)?ranShifts.at(3):0),
                                            relConfs.getD("drho_alpha2_LimL"), relConfs.getD("drho_alpha2_LimU") );
  RooRealVar *drho_n1 = new RooRealVar("drho_n1","",relConfs.getD("drho_n1")+((_shiftFactor==1)?ranShifts.at(4):0),
                                            relConfs.getD("drho_n1_LimL"), relConfs.getD("drho_n1_LimU") );
  RooRealVar *drho_n2 = new RooRealVar("drho_n2","",relConfs.getD("drho_n2")+((_shiftFactor==1)?ranShifts.at(5):0),
                                            relConfs.getD("drho_n2_LimL"), relConfs.getD("drho_n2_LimU") );
  RooRealVar *drho_frac1 = new RooRealVar("drho_frac1","",relConfs.getD("drho_frac1")+((_shiftFactor==1)?ranShifts.at(6):0),
                                            relConfs.getD("drho_frac1_LimL"), relConfs.getD("drho_frac1_LimU") );

  //set certain parameters constant
  //this list has to be maintained manually
  
  std::cout << std::endl << "PdfFit: Setting parameters constant" << std::endl;
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
          //bs_dstkst_010[*m][*c][*t][*a]->setRelation("frac010",bs_010_frac010);
          //bs_dstkst_001[*m][*c][*t][*a]->setRelation("frac010",bs_001_frac010);

          bd_dstkst[*m][*c][*t][*a]->setRelation("frac010",bd_frac010);
          //bd_dstkst[*m][*c][*t][*a]->setRelation("frac010",bs_frac010);

          /*
          if(*t=="LL"){
            bs_dstkst[*m][*c]["LL"][*a]->setRelation("frac010",bs_LL_frac010);
            bd_dstkst[*m][*c]["LL"][*a]->setRelation("frac010",bd_LL_frac010);
          }
          if(*t=="DD"){
            bs_dstkst[*m][*c]["DD"][*a]->setRelation("frac010",bs_DD_frac010);
            bd_dstkst[*m][*c]["DD"][*a]->setRelation("frac010",bd_DD_frac010);
          }
          */

          //Combinatoric - floating Exponential
          if(*m=="d2kspipi") {
            if(*t=="mix") comb[*m][*c]["mix"][*a]->setRelation("slope",combs_slope_mix);
            //if(*t=="LL")  comb[*m][*c]["LL"][*a]->setRelation("slope",combs_slope_LL);
            //if(*t=="DD")  comb[*m][*c]["DD"][*a]->setRelation("slope",combs_slope_DD);
          }
          else
            comb[*m][*c][*t][*a]->setRelation("slope",combs_slope_mix_kskk);

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
          //roopdf_bs_dstkst_010[*m][*c][*t][*a]  = bs_dstkst_010[*m][*c][*t][*a]->getPdf();
          //roopdf_bs_dstkst_001[*m][*c][*t][*a]  = bs_dstkst_001[*m][*c][*t][*a]->getPdf();
          roopdf_bd_dstkst[*m][*c][*t][*a]  = bd_dstkst[*m][*c][*t][*a]->getPdf();
          //roopdf_dkpipi[*m][*c][*t][*a]     = dkpipi[*m][*c][*t][*a]->getPdf();
          //roopdf_dpipipi[*m][*c][*t][*a]    = dpipipi[*m][*c][*t][*a]->getPdf();
          //roopdf_lambda[*m][*c][*t][*a]     = lambda[*m][*c][*t][*a]->getPdf();
        }
      }
    }
  }
	
  cout << "Finished making Pdf_Fit" << endl;
}

