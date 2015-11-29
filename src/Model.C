#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include "TFile.h"

#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooAbsPdf.h"
#include "RooMsgService.h"
#include "Model.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"

Model::Model(Settings* genConfs, RooRealVar* pmB, RooCategory* cat, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackList, std::vector<std::string> binList)
  : Base()
  , _modeList(modeList)
  , _chargeList(chargeList)
  , _trackList(trackList)
  , _binList(binList)									 
{  
  _genConfs=genConfs;
  mB=pmB;
  _cat=cat;
	
  // Set up yields -- contains the yields for gen and fit pdfs
  Settings* yields_fileList = new Settings("Yields settings");
  yields_fileList->readPairStringsToMap("Inputs/ControlFiles.txt");
  std::string unblind=_genConfs->get("UNBLIND");
  if(_genConfs->get("genToys")=="true") unblind="true"; // if generating toys, don't want to blind
  yields = new Yields(genConfs, yields_fileList,_modeList,_chargeList,_trackList,_binList,unblind);

  // Set up shapes
  Settings* pdf_fit_settings = new Settings("pdf_fit_settings");
  pdf_fit_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_FitPdf.txt");
  pdf_fit_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_GenPdf.txt");
  fitPdf = Pdf_Fit(pdf_fit_settings,_genConfs,mB,_modeList,_chargeList,_trackList,_binList,_genConfs->getI("inputTwoStageFitErrors"),_genConfs->get("MCsimfit"));
}

RooSimultaneous* Model::getGenPdf()
{
  // Create simultaneous PDF
  sim = new RooSimultaneous("model_gen","Simultaneous GENERATING model",*_cat);    

  // Set up shapes
  Settings* pdf_gen_settings = new Settings("pdf_gen_settings");
  pdf_gen_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_GenPdf.txt");
  genPdf = Pdf_Gen(pdf_gen_settings,mB,_modeList,_chargeList,_trackList,_binList);

  //Add to master PDF
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){
          std::string tag=(*m)+underscore+(*c)+underscore+(*t)+underscore+(*a);

          RooArgSet pdflist;
          RooArgSet nevents;

          // Bd
          pdflist.add(*(genPdf.roopdf_bd[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bd_gen[*m][*c][*t][*a]));
          // Bs
          pdflist.add(*(genPdf.roopdf_bs[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bs_gen[*m][*c][*t][*a]));
          if(_genConfs->get("MCsimfit")!="true") {
            // Combinatorics
            pdflist.add(*(genPdf.roopdf_comb[*m][*c][*t][*a]));
            nevents.add(*(yields->n_comb[*m][*c][*t][*a]));
            // Part Reco
            if(_genConfs->get("bd_dstkst")=="true"){
              pdflist.add(*(genPdf.roopdf_bd_dstkst[*m][*c][*t][*a]));
              nevents.add(*(yields->n_bd_dstkst_gen[*m][*c][*t][*a]));
            }
            if(_genConfs->get("bs_dstkst")=="true"){
              pdflist.add(*(genPdf.roopdf_bs_dstkst[*m][*c][*t][*a]));
              nevents.add(*(yields->n_bs_dstkst_gen[*m][*c][*t][*a]));
            }
            // Drho
            if(_genConfs->get("bd_drho")=="true"){
              pdflist.add(*(genPdf.roopdf_drho[*m][*c][*t][*a]));
              nevents.add(*(yields->n_drho_gen[*m][*c][*t][*a]));
            }
          }
         

          // --- Print out generated yields ---
          cout << "Generating yields ..." << endl;
          TIter nevtPar(nevents.createIterator());
          RooRealVar * par;
          while ((par = (RooRealVar *)nevtPar())) {
            cout << " " << par->GetName() << " " << par->getVal() << endl;
          }

          RooAddPdf* pdf = new RooAddPdf(Form("GENpdf_%s",tag.c_str()) ,"",pdflist,nevents);
          pdf->Print();
          sim->addPdf(*pdf,Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str()));

        }
      }
    }
  }

  std::cout<<" Have set generating PDF in Model.C "<<std::endl;
  return sim;
}

RooSimultaneous* Model::getFitPdf()
{
  // Create simultaneous PDF
  sim = new RooSimultaneous("model_fit","Simultaneous FITTING model",*_cat);    

         
  std::cout << "creating model" << std::endl;
  
  //Add to master PDF
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){

          std::string tag=(*m)+underscore+(*c)+underscore+(*t)+underscore+(*a);
          RooArgSet pdflist;
          RooArgSet nevents;

          // Bd
          pdflist.add(*(fitPdf.roopdf_bd[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bd_fit[*m][*c][*t][*a]));
          // Bs
          pdflist.add(*(fitPdf.roopdf_bs[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bs_fit[*m][*c][*t][*a]));
          if(_genConfs->get("MCsimfit")!="true")
            {
              // Combinatorics
              pdflist.add(*(fitPdf.roopdf_comb[*m][*c][*t][*a]));
              nevents.add(*(yields->n_comb[*m][*c][*t][*a]));
              // Part Reco
              if(_genConfs->get("bd_dstkst")=="true"){
                pdflist.add(*(fitPdf.roopdf_bd_dstkst[*m][*c][*t][*a]));
                nevents.add(*(yields->n_bd_dstkst[*m][*c][*t][*a]));
              }
              if(_genConfs->get("bs_dstkst")=="true"){
                pdflist.add(*(fitPdf.roopdf_bs_dstkst[*m][*c][*t][*a]));
                nevents.add(*(yields->n_bs_dstkst[*m][*c][*t][*a]));
                //pdflist.add(*(fitPdf.roopdf_bs_dstkst_010[*m][*c][*t][*a]));
                //nevents.add(*(yields->n_bs_dstkst_010[*m][*c][*t][*a]));
                //pdflist.add(*(fitPdf.roopdf_bs_dstkst_001[*m][*c][*t][*a]));
                //nevents.add(*(yields->n_bs_dstkst_001[*m][*c][*t][*a]));
              }
              // Drho
              if(_genConfs->get("bd_drho")=="true"){
                pdflist.add(*(fitPdf.roopdf_drho[*m][*c][*t][*a]));
                nevents.add(*(yields->n_drho[*m][*c][*t][*a]));
              }
              // Lambda
              if(_genConfs->get("lb_dppi")=="true"){ 
                pdflist.add(*(fitPdf.roopdf_lambda[*m][*c][*t][*a]));
                nevents.add(*(yields->n_lambda[*m][*c][*t][*a]));
              }
              // D3h
              if(_genConfs->get("bu_dkpipi")=="true"){ 
                pdflist.add(*(fitPdf.roopdf_dkpipi[*m][*c][*t][*a]));
                nevents.add(*(yields->n_dkpipi[*m][*c][*t][*a]));
              }
              if(_genConfs->get("bu_dpipipi")=="true"){ 
                pdflist.add(*(fitPdf.roopdf_dpipipi[*m][*c][*t][*a]));
                nevents.add(*(yields->n_dpipipi[*m][*c][*t][*a]));
              }
            }
          
          // --- No Gaussian Constraints --- 
          //RooAddPdf* pdf = new RooAddPdf(Form("FITpdf_%s",tag.c_str()) ,"",pdflist,nevents);

          // --- The Gaussian Constraints --- 
          RooAddPdf* addpdf = new RooAddPdf(Form("FITpdf_%s",tag.c_str()) ,"",pdflist,nevents);
          RooArgSet constpdfset(*addpdf);

          if(_genConfs->get("genToys")=="true" && _genConfs->get("fixedGC")=="false")
          {
            int seed = _genConfs->getD("startSeed");
            TRandom3* rand = new TRandom3();
            rand->SetSeed(seed*2);

            // Set up sampling Gaussians to fit toys
            double f010_mean(0.), f010_err(0.);
            double r_dstkst_mean(0.), r_dstkst_err(0.);
            double r_drho_mean(0.), r_drho_err(0.);

            if(_genConfs->getI("fit_limit_low")==5200) 
            {
              f010_mean = 0.75;
              f010_err = 0.13;
              r_dstkst_mean = 0.42;
              r_dstkst_err = 0.14;
              r_drho_mean = 0.035;
              r_drho_err = 0.010;
            }
            if(_genConfs->getI("fit_limit_low")==5160) 
            {
              f010_mean = 0.77;
              f010_err = 0.11;
              r_dstkst_mean = 0.76;
              r_dstkst_err = 0.29;
              r_drho_mean = 0.035;
              r_drho_err = 0.010;
            }

            // frac010
            double f010_mean_fit = rand->Gaus(f010_mean,f010_err);
            gaus_f010 = new RooGaussian("gaus_f010","",*(fitPdf.bs_frac010),RooFit::RooConst(f010_mean_fit), RooFit::RooConst(f010_err));
            constpdfset.add(*gaus_f010);
            // r(D*K*)
            double r_dstkst_mean_fit = rand->Gaus(r_dstkst_mean,r_dstkst_err);
            gaus_r_dstkst = new RooGaussian("gaus_r_dstkst","",*(yields->ratio_bs_dstkst[*c][*t][*a]), RooFit::RooConst(r_dstkst_mean_fit), RooFit::RooConst(r_dstkst_err));
            constpdfset.add(*gaus_r_dstkst);
            // r(Drho) -- don't fluctuate?
            //double r_drho_mean_fit = rand->Gaus(r_drho_mean,r_drho_err);
            //RooGaussian* gaus_r_drho = new RooGaussian("gaus_r_drho","",*(yields->ratio_bs_drho[*c][*t][*a]), RooFit::RooConst(r_drho_mean_fit), RooFit::RooConst(r_drho_err));
            //if(_genConfs->get("bd_drho")=="true") constpdfset.add(*gaus_r_drho);
            constpdfset.add(*(yields->gausratio_bs_drho[*c][*t][*a]));

          }
          else
          {
            // Add gaussian constraint for yields
            if(_genConfs->get("bd_drho")=="true"){
              constpdfset.add(*(yields->gausratio_bs_drho[*c][*t][*a]));
            }
            if(_genConfs->get("bs_dstkst")=="true"){
              constpdfset.add(*(yields->gausratio_bs_dstkst[*c][*t][*a]));
              constpdfset.add(*(fitPdf.gaus_frac010_bs[*c][*t][*a]));
//            // Split helamp yields
//            //constpdfset.add(*(yields->gausratio_bs_dstkst_010[*c][*t][*a]));
//            //constpdfset.add(*(yields->gausratio_bs_dstkst_001[*c][*t][*a]));
            }
            if(_genConfs->get("bd_dstkst")=="true") constpdfset.add(*(yields->gausratio_bd_dstkst[*c][*t][*a]));
            if(_genConfs->get("bd_dstkst")=="true") constpdfset.add(*(fitPdf.gaus_frac010_bd[*c][*t][*a]));
            if(_genConfs->get("lb_dppi")=="true") constpdfset.add(*(yields->gausratio_bs_lambda[*c][*t][*a]));
            if(_genConfs->get("bu_dkpipi")=="true") constpdfset.add(*(yields->gausratio_bs_dkpipi[*c][*t][*a]));
            if(_genConfs->get("bu_dpipipi")=="true") constpdfset.add(*(yields->gausratio_bs_dpipipi[*c][*t][*a]));
          }

          constpdfset.Print("v");
          RooAbsPdf* pdf = new RooProdPdf(Form("CONSTFITpdf_%s",tag.c_str()),"",constpdfset);

          // --- Add to simultaneous pdf ---
          sim->addPdf(*pdf,Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str()));

        }
      }
    }
  }

  std::cout<<" Have set fit PDF in Model.C "<<std::endl;  
  return sim;
}


void Model::printYieldsAndPurities(string b, double integ_limit_low, double integ_limit_high)
{

  mB->setRange("Bsigbox",integ_limit_low, integ_limit_high);

  // --- for manual integration of RooKeysPdf ---
  RooDataSet *drho_integEvents = fitPdf.roopdf_drho[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  RooDataSet *bd_dstkst_integEvents = fitPdf.roopdf_bd_dstkst[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  RooDataSet *bs_dstkst_integEvents = fitPdf.roopdf_bs_dstkst[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  //RooDataSet *bs_dstkst_010_integEvents = fitPdf.roopdf_bs_dstkst_010[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  //RooDataSet *bs_dstkst_001_integEvents = fitPdf.roopdf_bs_dstkst_001[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  RooDataSet *lambda_integEvents = 0;
  if(_genConfs->get("lb_dppi")=="true") lambda_integEvents = fitPdf.roopdf_lambda[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));

  std::string integRange = Form("%s>%f && %s<%f",mB->GetName(), integ_limit_low, mB->GetName(), integ_limit_high);
  std::string fitRange = Form("%s>%f && %s<%f",mB->GetName(), _genConfs->getD("fit_limit_low"), mB->GetName(),_genConfs->getD("fit_limit_high"));

  // --- output yields to a text file --- 
  string ofilename="";
  if(b=="Bs")   ofilename="output/GenTotals_bs_"+_genConfs->get("fit_limit_low")+".txt";
  if(b=="Bd")   ofilename="output/GenTotals_reduced_"+_genConfs->get("fit_limit_low")+".txt";
  if(b=="full") ofilename="output/GenTotals_"+_genConfs->get("fit_limit_low")+".txt";
  if(b=="binnedfit") ofilename="output/GenTotals_binnedfit_"+_genConfs->get("fit_limit_low")+".txt";

  ofstream GenTotals(ofilename);
  GenTotals << "* 2013 (2011 + 2012)" << std::endl;
  GenTotals << "*" << std::endl;
  GenTotals << "DEBUG_runInHighStatsMode 0" << std::endl;

  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){
          //////////////////////////////////////////////////////////////////////
          // Integrals in signal window
          //////////////////////////////////////////////////////////////////////
          double integral_bd  = fitPdf.roopdf_bd[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          double integral_bs  = fitPdf.roopdf_bs[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          double integral_comb   = fitPdf.roopdf_comb[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          // integration of RooKeysPdf needs to be done by hand 
          double integral_drho    = drho_integEvents->sumEntries(integRange.c_str())/drho_integEvents->sumEntries(fitRange.c_str());
          double integral_bd_dstkst    = bd_dstkst_integEvents->sumEntries(integRange.c_str())/bd_dstkst_integEvents->sumEntries(fitRange.c_str());
          double integral_bs_dstkst    = bs_dstkst_integEvents->sumEntries(integRange.c_str())/bs_dstkst_integEvents->sumEntries(fitRange.c_str());
          double integral_lambda = 0;
          if(_genConfs->get("lb_dppi")=="true") integral_lambda    = lambda_integEvents->sumEntries(integRange.c_str())/lambda_integEvents->sumEntries(fitRange.c_str());

          //////////////////////////////////////////////////////////////////////
          // Integrated yields
          //////////////////////////////////////////////////////////////////////
          //std::cout << "Normalised integrals" << std::endl;
          double integyield_bd     = integral_bd * yields->n_bd_fit[*m][*c][*t][*a]->getVal();
          double integyield_bs     = integral_bs * yields->n_bs_fit[*m][*c][*t][*a]->getVal();
          double integyield_comb   = integral_comb * yields->n_comb[*m][*c][*t][*a]->getVal();
          double integyield_bd_dstkst  = integral_bd_dstkst * yields->n_bd_dstkst[*m][*c][*t][*a]->getVal();
          double integyield_bs_dstkst   = integral_bs_dstkst * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_dstkst[*c][*t][*a]->getVal();
          double integyield_drho   = integral_drho * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_drho[*c][*t][*a]->getVal();
          double integyield_lambda   = 0;
          if(_genConfs->get("lb_dppi")=="true") integyield_lambda = integral_lambda * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_lambda[*c][*t][*a]->getVal();

          // Errors
          double integyield_bd_err     = integral_bd * yields->n_bd_fit[*m][*c][*t][*a]->getError();
          double integyield_bs_err     = integral_bs * yields->n_bs_fit[*m][*c][*t][*a]->getError();
          double integyield_comb_err   = integral_comb * yields->n_comb[*m][*c][*t][*a]->getError();
          //double integyield_bd_dstkst_err = integral_bd_dstkst * yields->n_bd_dstkst[*m][*c][*t][*a]->getError();
          double integyield_bd_dstkst_err   = integyield_bd_dstkst * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bd_dstkst[*c][*t][*a]->getError()/yields->ratio_bd_dstkst[*c][*t][*a]->getVal(),2) ); 
          double integyield_bs_dstkst_err   = integyield_bs_dstkst * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_dstkst[*c][*t][*a]->getError()/yields->ratio_bs_dstkst[*c][*t][*a]->getVal(),2) ); 
          double integyield_drho_err   = integyield_drho * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_drho[*c][*t][*a]->getError()/yields->ratio_bs_drho[*c][*t][*a]->getVal(),2) ); 
          double integyield_lambda_err = 0;
          if(_genConfs->get("lb_dppi")=="true") integyield_lambda_err = integyield_lambda * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_lambda[*c][*t][*a]->getError()/yields->ratio_bs_lambda[*c][*t][*a]->getVal(),2) ); 

          // --- If Bs->D*K* yields split by helamp ---
          //double integral_bs_dstkst_010    = bs_dstkst_010_integEvents->sumEntries(integRange.c_str())/bs_dstkst_010_integEvents->sumEntries(fitRange.c_str());
          //double integral_bs_dstkst_001    = bs_dstkst_001_integEvents->sumEntries(integRange.c_str())/bs_dstkst_001_integEvents->sumEntries(fitRange.c_str());
          //double integyield_bs_dstkst_010   = integral_bs_dstkst_010 * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_dstkst_010[*c][*t][*a]->getVal();
          //double integyield_bs_dstkst_001   = integral_bs_dstkst_001 * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_dstkst_001[*c][*t][*a]->getVal();
          //double integyield_bs_dstkst_010_err   = integyield_bs_dstkst_010 * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_dstkst_010[*c][*t][*a]->getError()/yields->ratio_bs_dstkst_010[*c][*t][*a]->getVal(),2) ); 
          //double integyield_bs_dstkst_001_err   = integyield_bs_dstkst_001 * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_dstkst_001[*c][*t][*a]->getError()/yields->ratio_bs_dstkst_001[*c][*t][*a]->getVal(),2) ); 
          //integyield_bs_dstkst = integyield_bs_dstkst_010 + integyield_bs_dstkst_001;
          //integyield_bs_dstkst_err = sqrt( pow(integyield_bs_dstkst_010_err,2) + pow(integyield_bs_dstkst_001_err,2));

          // Set to zero any yields of pdfs that are not in fit
          if(_genConfs->get("bd_dstkst")!="true") { integyield_bd_dstkst=0; integyield_bd_dstkst_err=0; }
          if(_genConfs->get("lb_dppi")!="true") { integyield_lambda=0; integyield_lambda_err=0; }

          //////////////////////////////////////////////////////////////////////
          // Print out yields in mass window
          //////////////////////////////////////////////////////////////////////
          cout.setf(ios::fixed);
          cout.precision(2);

          cout << "//////////////////////////////////////////////////////////////////////\n// Integrated yields in " << b << " window (" << integ_limit_low << " - " << integ_limit_high << ")\n/////////////////////////////////////////////////////////" << endl;
          cout<<"In bin : "<<*m<<", "<<*c<<", "<<*t<<", "<<*a<<endl;
          cout << "Bd signal: " << integyield_bd << " +/- " << integyield_bd_err << std::endl;
          cout << "Bs:        " << integyield_bs << " +/- " << integyield_bs_err << std::endl;
          cout << "Combs:     " << integyield_comb << " +/- " << integyield_comb_err << std::endl;
          cout << "D0rho0:    " << integyield_drho << " +/- " << integyield_drho_err << std::endl;
          cout << "Bs D*K*:   " << integyield_bs_dstkst << " +/- " << integyield_bs_dstkst_err << std::endl;
          cout << "Bd D*K*:   " << integyield_bd_dstkst << " +/- " << integyield_bd_dstkst_err << std::endl;
          cout << "Lambda:    " << integyield_lambda << " +/- " << integyield_lambda_err << std::endl;
          cout << "/////////////////////////////////////////////////////////////////" << endl;
          double total = integyield_bd + integyield_bs + integyield_comb + integyield_drho + integyield_bd_dstkst + integyield_bs_dstkst;
          double nS = (b == "Bd" ? integyield_bd : integyield_bs);
          double nSerr = (b == "Bd" ? integyield_bd_err : integyield_bs_err);
          double nB = total - nS;
          double nBerr =sqrt( pow(integyield_bs,2) + pow(integyield_comb,2) + pow(integyield_drho,2) + pow(integyield_bd_dstkst,2) + pow(integyield_bs_dstkst,2) );
          double purity = nS/total;
          double purity_err = sqrt( (nS*nS*nBerr*nBerr + nB*nB*nSerr*nSerr)/pow(nS+nB,4));
          plotNums[*m][*c][*t][*a]["purity_val"]=purity;
          plotNums[*m][*c][*t][*a]["purity_err"]=purity_err;
          if(b!="full") cout << "PURITY: " << purity << " +- " << purity_err << endl;

          //////////////////////////////////////////////////////////////////////
          // Latex style
          //////////////////////////////////////////////////////////////////////
          cout << "\\begin{tabular}{l r c l}" << endl;
          cout << "\\hline" << endl;
          cout << b << " window & & & \\\\" << endl;
          cout << "$B^0 \\to D^{0} K^{\\ast 0}$ & $" << integyield_bd << "$ & $\\pm$ & $" << integyield_bd_err << "$ \\\\ " << endl;
          cout << "$B_s^0 \\to D^{0} K^{\\ast 0}$ & $" << integyield_bs << "$ & $\\pm$ & $" << integyield_bs_err << "$ \\\\ " << endl;
          cout << "$\\mathrm{Combinatoric}$ & $" << integyield_comb << "$ & $\\pm$ & $" << integyield_comb_err << "$ \\\\ " << endl;
          cout << "$B^0 \\to D^{0} \\rho^{0}$ & $" << integyield_drho << "$ & $\\pm$ & $" << integyield_drho_err << "$ \\\\ " << endl;
          cout << "$B_s^0 \\to D^{\\ast 0} K^{\\ast 0}$ & $" << integyield_bs_dstkst << "$ & $\\pm$ & $" << integyield_bs_dstkst_err << "$ \\\\ " << endl;
          cout << "$B^0 \\to D^{\\ast 0} K^{\\ast 0}$ & $" << integyield_bd_dstkst << "$ & $\\pm$ & $" << integyield_bd_dstkst_err << "$ \\\\ " << endl;
          cout << "\\hline"  << endl;
          cout << "\\end{tabular}"  << endl;

          //////////////////////////////////////////////////////////////////////
          // Output to text files
          //////////////////////////////////////////////////////////////////////
          GenTotals << "N_bd_" << *m << "_both_" << *t << " " << integyield_bd << std::endl;
          GenTotals << "N_bd_" << *m << "_plus_" << *t << " " << integyield_bd/2.0 << std::endl;
          GenTotals << "N_bd_" << *m << "_minus_" << *t << " " << integyield_bd/2.0 << std::endl;
          GenTotals << "N_bs_" << *m << "_both_" << *t << " " << integyield_bs << std::endl;
          GenTotals << "N_bs_" << *m << "_plus_" << *t << " " << integyield_bs/2.0 << std::endl;
          GenTotals << "N_bs_" << *m << "_minus_" << *t << " " << integyield_bs/2.0 << std::endl;
          GenTotals << "N_comb_" << *m << "_both_" << *t << " " << integyield_comb << std::endl;
          GenTotals << "N_comb_" << *m << "_plus_" << *t << " " << integyield_comb/2.0 << std::endl;
          GenTotals << "N_comb_" << *m << "_minus_" << *t << " " << integyield_comb/2.0 << std::endl;
          GenTotals << "N_drho_" << *m << "_both_" << *t << " " << integyield_drho << std::endl;
          GenTotals << "N_drho_" << *m << "_plus_" << *t << " " << integyield_drho/2.0 << std::endl;
          GenTotals << "N_drho_" << *m << "_minus_" << *t << " " << integyield_drho/2.0 << std::endl;
          GenTotals << "N_bs_dstkst_" << *m << "_both_" << *t << " " << integyield_bs_dstkst << std::endl;
          GenTotals << "N_bs_dstkst_" << *m << "_plus_" << *t << " " << integyield_bs_dstkst/2.0 << std::endl;
          GenTotals << "N_bs_dstkst_" << *m << "_minus_" << *t << " " << integyield_bs_dstkst/2.0 << std::endl;
          GenTotals << "N_bd_dstkst_" << *m << "_both_" << *t << " " << integyield_bd_dstkst << std::endl;
          GenTotals << "N_bd_dstkst_" << *m << "_plus_" << *t << " " << integyield_bd_dstkst/2.0 << std::endl;
          GenTotals << "N_bd_dstkst_" << *m << "_minus_" << *t << " " << integyield_bd_dstkst/2.0 << std::endl;

        }
      }
    }
  }

  GenTotals.close();

}
/*
void Model::printYields()
{
   

  // need to integrate RooKeysPdf manually
  RooDataSet *drho_integEvents = fitPdf.roopdf_drho[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  RooDataSet *bd_dstkst_integEvents = fitPdf.roopdf_bd_dstkst[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  RooDataSet *bs_dstkst_integEvents = fitPdf.roopdf_bs_dstkst[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_binList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));
  std::string integRange = Form("%s>%f && %s<%f",mB->GetName(), _genConfs->getD("integ_limit_low_Bd"), mB->GetName(),_genConfs->getD("integ_limit_high_Bd"));
  std::string fitRange = Form("%s>%f && %s<%f",mB->GetName(), _genConfs->getD("fit_limit_low"), mB->GetName(),_genConfs->getD("fit_limit_high"));

  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_binList.begin();a!=_binList.end();a++){
          //////////////////////////////////////////////////////////////////////
          // Write out yields to text files
          //////////////////////////////////////////////////////////////////////
          GenTotals_full << "Nsignal_" << *m << "_fail_plus_" << *t << " " << integyield_full_dpi_sig/2.0 << std::endl;
          GenTotals_full << "Nsignal_" << *m << "_fail_minus_" << *t << " " << integyield_full_dpi_sig/2.0 << std::endl;
          GenTotals_full << "Nsignal_" << *m << "_pass_plus_" << *t << " " << integyield_full_dk_sig/2.0 << std::endl;
          GenTotals_full << "Nsignal_" << *m << "_pass_minus_" << *t << " " << integyield_full_dk_sig/2.0 << std::endl;
          
          GenTotals_full << "NLow_" << *m << "_fail_plus_" << *t << " " << integyield_full_dpi_lowmass/2.0 << std::endl;
          GenTotals_full << "NLow_" << *m << "_fail_minus_" << *t << " " << integyield_full_dpi_lowmass/2.0 << std::endl;
          GenTotals_full << "NLow_" << *m << "_pass_plus_" << *t << " " << integyield_full_dk_lowmass/2.0 << std::endl;
          GenTotals_full << "NLow_" << *m << "_pass_minus_" << *t << " " << integyield_full_dk_lowmass/2.0 << std::endl;
          
          GenTotals_full << "NComb_" << *m << "_fail_plus_" << *t << " " << integyield_full_dpi_comb/2.0 << std::endl;
          GenTotals_full << "NComb_" << *m << "_fail_minus_" << *t << " " << integyield_full_dpi_comb/2.0 << std::endl;
          GenTotals_full << "NComb_" << *m << "_pass_plus_" << *t << " " << integyield_full_dk_comb/2.0 << std::endl;
          GenTotals_full << "NComb_" << *m << "_pass_minus_" << *t << " " << integyield_full_dk_comb/2.0 << std::endl;

          GenTotals_reduced << "Nsignal_" << *m << "_fail_plus_" << *t << " " << integyield_reduced_dpi_sig/2.0 << std::endl;
          GenTotals_reduced << "Nsignal_" << *m << "_fail_minus_" << *t << " " << integyield_reduced_dpi_sig/2.0 << std::endl;
          GenTotals_reduced << "Nsignal_" << *m << "_pass_plus_" << *t << " " << integyield_reduced_dk_sig/2.0 << std::endl;
          GenTotals_reduced << "Nsignal_" << *m << "_pass_minus_" << *t << " " << integyield_reduced_dk_sig/2.0 << std::endl;
          
          GenTotals_reduced << "NLow_" << *m << "_fail_plus_" << *t << " " << integyield_reduced_dpi_lowmass/2.0 << std::endl;
          GenTotals_reduced << "NLow_" << *m << "_fail_minus_" << *t << " " << integyield_reduced_dpi_lowmass/2.0 << std::endl;
          GenTotals_reduced << "NLow_" << *m << "_pass_plus_" << *t << " " << integyield_reduced_dk_lowmass/2.0 << std::endl;
          GenTotals_reduced << "NLow_" << *m << "_pass_minus_" << *t << " " << integyield_reduced_dk_lowmass/2.0 << std::endl;
          
          GenTotals_reduced << "NComb_" << *m << "_fail_plus_" << *t << " " << integyield_reduced_dpi_comb/2.0 << std::endl;
          GenTotals_reduced << "NComb_" << *m << "_fail_minus_" << *t << " " << integyield_reduced_dpi_comb/2.0 << std::endl;
          GenTotals_reduced << "NComb_" << *m << "_pass_plus_" << *t << " " << integyield_reduced_dk_comb/2.0 << std::endl;
          GenTotals_reduced << "NComb_" << *m << "_pass_minus_" << *t << " " << integyield_reduced_dk_comb/2.0 << std::endl;
        }
      }
    }
  }

}
*/

