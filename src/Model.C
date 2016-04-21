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

Model::Model(Settings* genConfs, RooRealVar* pmB, RooCategory* cat, std::vector<std::string> modeList, std::vector<std::string> chargeList, std::vector<std::string> trackList, std::vector<std::string> runList)
  : Base()
  , _modeList(modeList)
  , _chargeList(chargeList)
  , _trackList(trackList)
  , _runList(runList)
{  
  _genConfs=genConfs;
  mB=pmB;
  _cat=cat;
	
  // Set up yields -- contains the yields for gen and fit pdfs
  Settings* yields_fileList = new Settings("Yields settings");
  yields_fileList->readPairStringsToMap("Settings/Yields/ControlFiles.txt");
  std::string unblind=_genConfs->get("UNBLIND");
  if(_genConfs->get("genToys")=="true") unblind="true"; // if generating toys, don't want to blind
  yields = new Yields(genConfs, yields_fileList,_modeList,_chargeList,_trackList,_runList,unblind);

  // Set up shapes
  Settings* pdf_fit_settings = new Settings("pdf_fit_settings");
  pdf_fit_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_FitPdf.txt");
  pdf_fit_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_GenPdf.txt");
  fitPdf = Pdf_Fit(pdf_fit_settings,_genConfs,mB,_modeList,_chargeList,_trackList,_runList,_genConfs->getI("inputTwoStageFitErrors"),_genConfs->get("MCsimfit"));
}

RooSimultaneous* Model::getGenPdf()
{
  // Create simultaneous PDF
  sim = new RooSimultaneous("model_gen","Simultaneous GENERATING model",*_cat);    

  // Set up shapes
  Settings* pdf_gen_settings = new Settings("pdf_gen_settings");
  pdf_gen_settings->readPairStringsToMap("Settings/PDFShapes/ControlFiles_GenPdf.txt");
  genPdf = Pdf_Gen(pdf_gen_settings,mB,_modeList,_chargeList,_trackList,_runList);

  //Add to master PDF
  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){
          std::string tag=(*m)+underscore+(*c)+underscore+(*t)+underscore+(*a);

          RooArgSet pdflist;
          RooArgSet nevents;

          // Bu
          pdflist.add(*(genPdf.roopdf_bu[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bu_gen[*m][*c][*t][*a]));

          if(_genConfs->get("MCsimfit")!="true") {
            // Combinatorics
            pdflist.add(*(genPdf.roopdf_comb[*m][*c][*t][*a]));
            nevents.add(*(yields->n_comb[*m][*c][*t][*a]));
            // DstKst
            pdflist.add(*(genPdf.roopdf_dstkst[*m][*c][*t][*a]));
            nevents.add(*(yields->n_dstkst_gen[*m][*c][*t][*a]));
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
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){

          std::string tag=(*m)+underscore+(*c)+underscore+(*t)+underscore+(*a);
          RooArgSet pdflist;
          RooArgSet nevents;

          // Bu
          pdflist.add(*(fitPdf.roopdf_bu[*m][*c][*t][*a]));
          nevents.add(*(yields->n_bu_fit[*m][*c][*t][*a]));

          if(_genConfs->get("MCsimfit")!="true")
            {
              // Combinatorics
              pdflist.add(*(fitPdf.roopdf_comb[*m][*c][*t][*a]));
              nevents.add(*(yields->n_comb[*m][*c][*t][*a]));
              // DstKst
              pdflist.add(*(fitPdf.roopdf_dstkst[*m][*c][*t][*a]));
              nevents.add(*(yields->n_dstkst[*m][*c][*t][*a]));
            }
          
          // --- No Gaussian Constraints --- 
          RooAddPdf* pdf = new RooAddPdf(Form("FITpdf_%s",tag.c_str()) ,"",pdflist,nevents);

          // --- The Gaussian Constraints --- 
          RooAddPdf* addpdf = new RooAddPdf(Form("FITpdf_%s",tag.c_str()) ,"",pdflist,nevents);
          RooArgSet constpdfset(*addpdf);

          if(_genConfs->get("genToys")=="true")
          {
            int seed = _genConfs->getD("startSeed");
            TRandom3* rand = new TRandom3();
            rand->SetSeed(seed*2);

/*            // Set up sampling Gaussians to fit toys
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
            }*/

            /*// frac010
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
*/
          }
          else
          {
/*            // Add gaussian constraint for yields
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
            if(_genConfs->get("bu_dpipipi")=="true") constpdfset.add(*(yields->gausratio_bs_dpipipi[*c][*t][*a]));*/
          }

          constpdfset.Print("v");
          //RooAbsPdf* pdf = new RooProdPdf(Form("CONSTFITpdf_%s",tag.c_str()),"",constpdfset);

          // --- Add to simultaneous pdf ---
          sim->addPdf(*pdf,Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str()));

        }
      }
    }
  }

  std::cout<<" Have set fit PDF in Model.C "<<std::endl;  
  return sim;
}


void Model::printYieldsAndPurities(string b, double integ_limit_low, double integ_limit_high, RooFitResult* result)
{

  mB->setRange("Bsigbox",integ_limit_low, integ_limit_high);

  // --- for manual integration of RooKeysPdf ---
  //RooDataSet *bu_dstkst_integEvents = fitPdf.roopdf_bu_dstkst[_modeList.at(0)][_chargeList.at(0)][_trackList.at(0)][_runList.at(0)]->generate(*mB, 100000, RooFit::Verbose(kFALSE));


  std::string integRange = Form("%s>%f && %s<%f",mB->GetName(), integ_limit_low, mB->GetName(), integ_limit_high);
  std::string fitRange = Form("%s>%f && %s<%f",mB->GetName(), _genConfs->getD("fit_limit_low"), mB->GetName(),_genConfs->getD("fit_limit_high"));

  // --- output yields to a text file --- 
  string ofilename="";
  if(b=="Bu")   ofilename="output/GenTotals_bu_"+_genConfs->get("fit_limit_low")+".txt";
  if(b=="full") ofilename="output/GenTotals_"+_genConfs->get("fit_limit_low")+".txt";
  //if(b=="binnedfit") ofilename="output/GenTotals_binnedfit_"+_genConfs->get("fit_limit_low")+".txt";

  ofstream GenTotals(ofilename);
  GenTotals << "* 2013 (2011 + 2012)" << std::endl;
  GenTotals << "*" << std::endl;
  GenTotals << "DEBUG_runInHighStatsMode 0" << std::endl;


  if(_genConfs->isChargeSeparated()) {
	  double A, A_error;
	  double R, R_error;
	  double N_kpi, N_kpi_error;

	  for(std::vector<std::string>::iterator m=_modeList.begin(); m!=_modeList.end();m++){

		  A = yields->A[*m]->getVal();
		  A_error = yields->A[*m]->getError();

		  if (*m != "d2kpi") {
			  R = yields->R[*m]->getVal();
			  R_error = yields->R[*m]->getError();
		  }
		  else {
			  for(std::vector<std::string>::iterator t=_trackList.begin(); t!=_trackList.end(); t++){
				  for(std::vector<std::string>::iterator a=_runList.begin(); a!=_runList.end();a++){

					  N_kpi = yields->N_kpi[*t][*a]->getVal();
					  N_kpi_error = yields->N_kpi[*t][*a]->getError();
			          GenTotals << "N_bu_" << *m << "_" << *a << "_" << *t << " " << N_kpi << std::endl;
				  }
			  }
		  }
          GenTotals << "A_" << *m << " " << A << std::endl;
          GenTotals << "R_" << *m << " " << R << std::endl;

	  }

  }


  for(std::vector<std::string>::iterator m=_modeList.begin();m!=_modeList.end();m++){
    for(std::vector<std::string>::iterator c=_chargeList.begin();c!=_chargeList.end();c++){
      for(std::vector<std::string>::iterator t=_trackList.begin();t!=_trackList.end();t++){
        for(std::vector<std::string>::iterator a=_runList.begin();a!=_runList.end();a++){
          //////////////////////////////////////////////////////////////////////
          // Integrals in signal window
          //////////////////////////////////////////////////////////////////////
          double integral_bu  = fitPdf.roopdf_bu[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          double integral_comb   = fitPdf.roopdf_comb[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          double integral_dstkst   = fitPdf.roopdf_dstkst[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          // integration of RooKeysPdf needs to be done by hand 
          //double integral_bd_dstkst    = bu_dstkst_integEvents->sumEntries(integRange.c_str())/bd_dstkst_integEvents->sumEntries(fitRange.c_str());
          //////////////////////////////////////////////////////////////////////
          // Integrated yields
          //////////////////////////////////////////////////////////////////////
          //std::cout << "Normalised integrals" << std::endl;
          double integyield_bu;
          double integyield_bu_err;
          if(_genConfs->isChargeSeparated()) {
        	  RooFormulaVar* n_bu_fit_asRooFormulaVar = static_cast<RooFormulaVar*>(yields->n_bu_fit[*m][*c][*t][*a]);
        	  integyield_bu = integral_bu * n_bu_fit_asRooFormulaVar->getVal();
        	  integyield_bu_err = integral_bu * n_bu_fit_asRooFormulaVar->getPropagatedError(*result);
          }
          if(!_genConfs->isChargeSeparated()) {
              RooRealVar* n_bu_fit_asRooRealVar = static_cast<RooRealVar*>(yields->n_bu_fit[*m][*c][*t][*a]);
              integyield_bu = integral_bu * n_bu_fit_asRooRealVar->getVal();
              integyield_bu_err = integral_bu * n_bu_fit_asRooRealVar->getError();
          }
          //double integyield_bu     = integral_bu * yields->n_bu_fit[*m][*c][*t][*a]->getVal();
          double integyield_comb   = integral_comb * yields->n_comb[*m][*c][*t][*a]->getVal();
          double integyield_dstkst  = integral_dstkst * yields->n_dstkst[*m][*c][*t][*a]->getVal();

          // Errors
          //double integyield_bu_err     = integral_bu * yields->n_bu_fit[*m][*c][*t][*a]->getError();
          double integyield_comb_err   = integral_comb * yields->n_comb[*m][*c][*t][*a]->getError();
          double integyield_dstkst_err = integral_dstkst * yields->n_dstkst[*m][*c][*t][*a]->getError();

          // --- If Bs->D*K* yields split by helamp ---
          //double integral_bs_dstkst_010    = bs_dstkst_010_integEvents->sumEntries(integRange.c_str())/bs_dstkst_010_integEvents->sumEntries(fitRange.c_str());
          //double integral_bs_dstkst_001    = bs_dstkst_001_integEvents->sumEntries(integRange.c_str())/bs_dstkst_001_integEvents->sumEntries(fitRange.c_str());
          //double integyield_bs_dstkst_010   = integral_bs_dstkst_010 * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_dstkst_010[*c][*t][*a]->getVal();
          //double integyield_bs_dstkst_001   = integral_bs_dstkst_001 * yields->n_bs_fit[*m][*c][*t][*a]->getVal() * yields->ratio_bs_dstkst_001[*c][*t][*a]->getVal();
          //double integyield_bs_dstkst_010_err   = integyield_bs_dstkst_010 * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_dstkst_010[*c][*t][*a]->getError()/yields->ratio_bs_dstkst_010[*c][*t][*a]->getVal(),2) ); 
          //double integyield_bs_dstkst_001_err   = integyield_bs_dstkst_001 * sqrt( pow(yields->n_bs_fit[*m][*c][*t][*a]->getError()/yields->n_bs_fit[*m][*c][*t][*a]->getVal(),2) + pow(yields->ratio_bs_dstkst_001[*c][*t][*a]->getError()/yields->ratio_bs_dstkst_001[*c][*t][*a]->getVal(),2) ); 
          //integyield_bs_dstkst = integyield_bs_dstkst_010 + integyield_bs_dstkst_001;
          //integyield_bs_dstkst_err = sqrt( pow(integyield_bs_dstkst_010_err,2) + pow(integyield_bs_dstkst_001_err,2));

         //////////////////////////////////////////////////////////////////////
          // Print out yields in mass window
          //////////////////////////////////////////////////////////////////////
          cout.setf(ios::fixed);
          cout.precision(2);


          cout << "//////////////////////////////////////////////////////////////////////\n// Integrated yields in " << b << " window (" << integ_limit_low << " - " << integ_limit_high << ")\n/////////////////////////////////////////////////////////" << endl;
          cout<<"In bin : "<<*m<<", "<<*c<<", "<<*t<<", "<<*a<<endl;
          cout << "Bu signal: " << integyield_bu << " +/- " << integyield_bu_err << std::endl;
          cout << "Combs:     " << integyield_comb << " +/- " << integyield_comb_err << std::endl;
          cout << "Bu D*K*:   " << integyield_dstkst << " +/- " << integyield_dstkst_err << std::endl;
          cout << "/////////////////////////////////////////////////////////////////" << endl;
          double total = integyield_bu + integyield_comb + integyield_dstkst;
          double nS = integyield_bu;
          double nSerr = integyield_bu_err;
          double nB = total - nS;
          double nBerr =sqrt( pow(integyield_comb,2) + pow(integyield_dstkst,2) );
          double purity = nS/total;
          double purity_err = sqrt( (nS*nS*nBerr*nBerr + nB*nB*nSerr*nSerr)/pow(nS+nB,4));
          // Not sure if the yield give should be in the B mass region or the total yield of signal peak
          // nS and nSerr- inside or outside the if statement?
          if(b!="full") {
          plotNums[*m][*c][*t][*a]["val"]= nS;
          plotNums[*m][*c][*t][*a]["err"]= nSerr;
          plotNums[*m][*c][*t][*a]["purity_val"]=purity;
          plotNums[*m][*c][*t][*a]["purity_err"]=purity_err;
          cout << "PURITY: " << purity << " +- " << purity_err << endl;
          }

          //////////////////////////////////////////////////////////////////////
          // Latex style
          //////////////////////////////////////////////////////////////////////
          cout << "\\begin{tabular}{l r c l}" << endl;
          cout << "\\hline" << endl;
          cout << b << " window & & & \\\\" << endl;
          cout << "$B^+ \\to D^{0} K^{\\ast +}$ & $" << integyield_bu << "$ & $\\pm$ & $" << integyield_bu_err << "$ \\\\ " << endl;
          cout << "$\\mathrm{Combinatoric}$ & $" << integyield_comb << "$ & $\\pm$ & $" << integyield_comb_err << "$ \\\\ " << endl;
          cout << "$B^+ \\to D^{\\ast 0} K^{\\ast +}$ & $" << integyield_dstkst << "$ & $\\pm$ & $" << integyield_dstkst_err << "$ \\\\ " << endl;
          cout << "\\hline"  << endl;
          cout << "\\end{tabular}"  << endl;

          //////////////////////////////////////////////////////////////////////
          // Output to text files
          //////////////////////////////////////////////////////////////////////
          GenTotals << "N_bu_" << *m << "_both_" << *t << " " << integyield_bu << std::endl;
          GenTotals << "N_bu_" << *m << "_plus_" << *t << " " << integyield_bu/2.0 << std::endl;
          GenTotals << "N_bu_" << *m << "_minus_" << *t << " " << integyield_bu/2.0 << std::endl;
          GenTotals << "N_comb_" << *m << "_both_" << *t << " " << integyield_comb << std::endl;
          GenTotals << "N_comb_" << *m << "_plus_" << *t << " " << integyield_comb/2.0 << std::endl;
          GenTotals << "N_comb_" << *m << "_minus_" << *t << " " << integyield_comb/2.0 << std::endl;
          GenTotals << "N_dstkst_" << *m << "_both_" << *t << " " << integyield_dstkst << std::endl;
          GenTotals << "N_dstkst_" << *m << "_plus_" << *t << " " << integyield_dstkst/2.0 << std::endl;
          GenTotals << "N_dstkst_" << *m << "_minus_" << *t << " " << integyield_dstkst/2.0 << std::endl;


        }
      }
    }
  }
  GenTotals.close();
}
