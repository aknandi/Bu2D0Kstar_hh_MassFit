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
            nevents.add(*(yields->n_comb_gen[*m][*c][*t][*a]));
            // DstKst
            pdflist.add(*(genPdf.roopdf_dstkst[*m][*c][*t][*a]));
            nevents.add(*(yields->n_dstkst_gen[*m][*c][*t][*a]));
            // LcKst
            if(*m=="d2kk") {
                pdflist.add(*(genPdf.roopdf_lckst[*m][*c][*t][*a]));
                nevents.add(*(yields->n_lckst_gen[*m][*c][*t][*a]));
            }
          }


          // --- Print out generated yields ---
          cout << "Generating yields ..." << endl;
          TIter nevtPar(nevents.createIterator());
          RooRealVar * par;
          while ((par = (RooRealVar *)nevtPar())) {
            cout << " " << par->GetName() << " " << par->getVal() << endl;
          }

          TIter pdfPar(pdflist.createIterator());
          RooRealVar * par2;
          while ((par2 = (RooRealVar *)pdfPar())) {
            cout << " " << par2->GetName() << " " << par2->getVal() << endl;
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
        	  if(*m=="d2kk") {
        		  pdflist.add(*(fitPdf.roopdf_lckst[*m][*c][*t][*a]));
        		  nevents.add(*(yields->n_lckst[*m][*c][*t][*a]));
        	  }
          }

          
          // --- No Gaussian Constraints --- 
          RooAddPdf* pdf = new RooAddPdf(Form("FITpdf_%s",tag.c_str()) ,"",pdflist,nevents);

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
  mB->setRange("combtest",5400, integ_limit_high);
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
	  double N_kpipipi, N_kpipipi_error;

	  for(std::vector<std::string>::iterator m=_modeList.begin(); m!=_modeList.end();m++){

		  A = yields->A[*m]->getVal();
		  A_error = yields->A[*m]->getError();

		  if(*m == "d2kpi"){
			  for(std::vector<std::string>::iterator t=_trackList.begin(); t!=_trackList.end(); t++){
				  for(std::vector<std::string>::iterator a=_runList.begin(); a!=_runList.end();a++){

					  N_kpi = yields->N_kpi[*t][*a]->getVal();
					  N_kpi_error = yields->N_kpi[*t][*a]->getError();
					  GenTotals << "N_bu_" << *m << "_" << *a << "_" << *t << " " << N_kpi << std::endl;

				  }
			  }
		  }
		  else if(*m == "d2kpipipi"){
			  for(std::vector<std::string>::iterator t=_trackList.begin(); t!=_trackList.end(); t++){
				  for(std::vector<std::string>::iterator a=_runList.begin(); a!=_runList.end();a++){

					  N_kpipipi = yields->N_kpipipi[*t][*a]->getVal();
					  N_kpipipi_error = yields->N_kpipipi[*t][*a]->getError();
					  GenTotals << "N_bu_" << *m << "_" << *a << "_" << *t << " " << N_kpipipi << std::endl;

				  }
			  }
		  }
		  else {
			  R = yields->R[*m]->getVal();
			  R_error = yields->R[*m]->getError();
		  }

          GenTotals << "A_" << *m << " " << A << std::endl;
          if (*m != "d2kpi" && *m != "d2kpipipi") {
        	  GenTotals << "R_" << *m << " " << R << std::endl;
          }

	  }

  }

  double adsSignal = 0, adsBackground = 0;
  double erradsSignalsq = 0, erradsBackgroundsq = 0;
  double erradsSignal = 0, erradsBackground = 0;

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
          double integral_lckst;
          if(*m=="d2kk") integral_lckst  = fitPdf.roopdf_lckst[*m][*c][*t][*a]->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("Bsigbox"))->getVal();
          // integration of RooKeysPdf needs to be done by hand 
          //double integral_bd_dstkst    = bu_dstkst_integEvents->sumEntries(integRange.c_str())/bd_dstkst_integEvents->sumEntries(fitRange.c_str());
          //////////////////////////////////////////////////////////////////////
          // Integrated yields
          //////////////////////////////////////////////////////////////////////
          //std::cout << "Normalised integrals" << std::endl;
          double integyield_bu;
          double integyield_bu_err;
          if(_genConfs->isChargeSeparated()) {
        	  RooFormulaVar* n_bu_fit_asRooFormulaVar;
        	  //if(*m == "d2pik" && _genConfs->get("UNBLIND")=="false") n_bu_fit_asRooFormulaVar = static_cast<RooFormulaVar*>(yields->n_bu_fit[*m][*c][*t][*a]);
        	  //else
        	  n_bu_fit_asRooFormulaVar = static_cast<RooFormulaVar*>(yields->n_bu_fit[*m][*c][*t][*a]);
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

          double integyield_lckst, integyield_lckst_err;
          if(*m=="d2kk") {
        	  integyield_lckst = integral_lckst * static_cast<RooFormulaVar*>(yields->n_lckst[*m][*c][*t][*a])->getVal();
        	  integyield_lckst_err = integral_lckst * static_cast<RooFormulaVar*>(yields->n_lckst[*m][*c][*t][*a])->getPropagatedError(*result);
          }

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
          if(*m=="d2kk") cout << "LcKst:     " << integyield_lckst  << " +/- " << integyield_lckst_err  << std::endl;
          cout << "/////////////////////////////////////////////////////////////////" << endl;
          double total = integyield_bu + integyield_comb + integyield_dstkst + ((*m=="d2kk")?integyield_lckst:0);
          double nS = integyield_bu;
          double nSerr = integyield_bu_err;
          double nB = total - nS;
          double nBerr =sqrt( pow(integyield_comb_err,2) + pow(integyield_dstkst_err,2) + ((*m=="d2kk")?pow(integyield_lckst_err,2):0) );
          double purity = nS/total;
          double purity_err = sqrt( (nS*nS*nBerr*nBerr + nB*nB*nSerr*nSerr)/pow(nS+nB,4));
          double significance = nS/sqrt(total);
          double significance_err = sqrt(nS*nS*(pow(nS+nB+nB,2)*nSerr*nSerr + nS*nS*nBerr*nBerr)/(4*pow(nS+nB,4)));
          // Not sure if the yield give should be in the B mass region or the total yield of signal peak
          // nS and nSerr- inside or outside the if statement?
          if(b!="full") {
        	  plotNums[*m][*c][*t][*a]["val"]= nS;
        	  plotNums[*m][*c][*t][*a]["err"]= nSerr;
        	  plotNums[*m][*c][*t][*a]["purity_val"]=purity;
        	  plotNums[*m][*c][*t][*a]["purity_err"]=purity_err;
        	  cout << "PURITY: " << purity << " +- " << purity_err << endl;
        	  cout << "SIGNIFICANCE: " << significance << " +- " << significance_err <<endl;
        	  if(*m=="d2pik") {
                  adsSignal += nS;
        	  	  erradsSignalsq += pow(integyield_bu_err,2);
                  adsBackground += nB;
                  erradsBackgroundsq += pow(integyield_comb_err,2) + pow(integyield_dstkst_err,2) + ((*m=="d2kk")?pow(integyield_lckst_err,2):0);
        	  }
          }

        if(b!="full") {
        	erradsSignal = sqrt(erradsSignalsq);
        	erradsBackground = sqrt(erradsBackgroundsq);
        	totalSignificance = adsSignal/sqrt(adsSignal+adsBackground);
        	errSignificance = sqrt(adsSignal*adsSignal*(pow(adsSignal+adsBackground+adsBackground,2)*erradsSignal*erradsSignal + adsSignal*adsSignal*erradsBackground*erradsBackground)/(4*pow(adsSignal+adsBackground,4)));
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
