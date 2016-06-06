#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <map>
#include <math.h>
#include <time.h>

#include "TSystem.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TIterator.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TKey.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TTree.h"

#include "RooHist.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooErrorVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooMCStudy.h"
#include "Roo1DTable.h"

#include "Fitting.h"
#include "CommonTools.h"
#include "TH2F.h"
#include "TROOT.h"
#include <boost/algorithm/string.hpp>

// Trace memory issues
#include "RooTrace.h"

Fitting::Fitting(TApplication* app, Settings* genConfs)
  : Base()
  , inputlist("contents of Final ntuple")
  , fulllist("contents, including categories")
  , reducedlist("reduced content")
  , mB("Bu_D0constKS0constPVconst_M","m(B)",genConfs->getD("fit_limit_low"),genConfs->getD("fit_limit_high"),"MeV/c^{2}")
  , assignedCharge("KstK_Q","KstK charge",-1,1.0)
  //, BDTGresponse("BDT","BDT",0.7,1)
  , mode("mode","D^{0} decay mode")
  , run("run","Data set")
  , charge("charge","bachelor charge")
  , track("track","Ks daughters' track type")
  , _genConfs(genConfs)
  , batchMode()
  , drawProjections()
  , doFit()
  , readToys()
  , readData()
  , genToys()
  , modeList()
  , chargeList()
  , trackList()
{
  // Memory problem tracing
  //RooTrace::verbose(kTRUE);
  //RooTrace::active(kTRUE);
  
  // Read true/false parameters from GeneralSettings file and convert the string to lower case
  batchMode = _genConfs->get("batchMode");boost::to_lower(batchMode);
  drawProjections = _genConfs->get("drawProjections");boost::to_lower(drawProjections);
  doFit = _genConfs->get("doFit");boost::to_lower(doFit);
  readToys = _genConfs->get("readToys");boost::to_lower(readToys);
  readData = _genConfs->get("readData");boost::to_lower(readData);
  genToys = _genConfs->get("genToys");boost::to_lower(genToys);

  //drawing options
  gStyle->SetErrorX(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetNdivisions(205,"XYZ");
  gStyle->SetStatFont(132);
  gStyle->SetStatFontSize(0.08);
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTitleSize(0.077,"XYZ");
  gStyle->SetLabelSize(0.08,"XYZ");
  gStyle->SetTitleOffset(0.83,"X");
  gStyle->SetTitleOffset(1.05,"Y");

  //'LHCb preliminary' box
  lhcbpreliminary = new TPaveText(0.7,0.72,0.83,0.9,"TR NDC");
  lhcbpreliminary->SetBorderSize(0); lhcbpreliminary->SetFillStyle(0);
  lhcbpreliminary->SetTextFont(132); lhcbpreliminary->SetTextSize(0.075); 
  //lhcbpreliminary->AddText("LHCb 2013");
  //lhcbpreliminary->AddText("LHCb preliminary");
  lhcbpreliminary->AddText("#scale[0.5]{#int }L d#it{t} = 3.0 fb^{-1}");

  // Setup the file for fit projections output
  saveOutputForPlottingMacro = new TFile("output/saveOutputForPlottingMacro.root","RECREATE");
  // Set up the random generator reproducibly	
  RooRandom::randomGenerator()->SetSeed(_genConfs->getI("startSeed"));
  gRandom=RooRandom::randomGenerator();
	
  // Create model object and super category from which can be retrieved the fit/gen pdfs
  DefineRooCategories();
  cat = new RooSuperCategory("cat","mode/charge/track/run",RooArgSet(mode,charge,track,run));
  //model = new Model(_genConfs,&mB,cat,modeList,chargeList,trackList,runList);
  model = new Model(_genConfs,&mB,catNew,modeList,chargeList,trackList,runList);

  if(readToys=="true")
    {
      DisplayToys();
    }
  else
    {
      if(genToys=="true")
        {
          if(_genConfs->getI("nToys")==1)
            {
              OrderToys(1);
              PrintDataSet(false);
              RunFullFit(true);
            }
          else
            {
              RunManyToys();
              return;
            }
        }
      else if(readData=="true")
        {
          cout << "**********************************************************************" << endl;
          cout << "Ready to read in data" << endl;
          LoadDataSet();
          PrintDataSet(false);
          RunFullFit(true);
        }
      else
        {
          cout << "I'm not sure what you want to do. Check your settings and the logic flow in Fitting.C" << endl;
        }
    }
  if(batchMode=="false")
    {
      std::cout<<"Starting: app->Run()"<<std::endl;
      app->Run(true);
    }
  saveOutputForPlottingMacro->Close();
}

Fitting::~Fitting(){}

void Fitting::DefineRooCategories()
{
  // Define categories. In each case push back into a list (to loop over) and define a new 'type' in that RooCategory.
  
  // Define mode category
  modeList.push_back(d2kpi);
  modeList.push_back(d2kk);
  modeList.push_back(d2pipi);
  //modeList.push_back(d2pik);
  for (std::vector<std::string>::iterator m = modeList.begin(); m != modeList.end(); m++)
    {
      mode.defineType((*m).c_str());
    }

  // Define charge categories - could be minus, plus or both
  if(_genConfs->get("chargeSeparated")=="true") {
	  chargeList.push_back(minus);
	  chargeList.push_back(plus);
  }
  else {
	  chargeList.push_back(both);
  }
  for (std::vector<std::string>::iterator c = chargeList.begin(); c != chargeList.end(); c++)
    {
      charge.defineType((*c).c_str());
    }

 
  // Define Ks track type categories - could be LL, DD or mix
  if(_genConfs->get("trackSeparated")=="true") {
	  trackList.push_back(LL);
	  trackList.push_back(DD);
  }
  else {
	  trackList.push_back(mix);
  }
  for (std::vector<std::string>::iterator t = trackList.begin(); t != trackList.end(); t++)
    {
      track.defineType((*t).c_str());
    }

  // Define run categories - could be run1, run2 or all
  if(_genConfs->get("runSeparated")=="true") {
	  runList.push_back(run1);
	  runList.push_back(run2);
  }
  else {
	  // Do not have combined Ntuples yet. For one run category, just use run 1 (change all->run1)
	  runList.push_back(all);
  }

  for (std::vector<std::string>::iterator a = runList.begin(); a != runList.end(); a++)
    {
      run.defineType((*a).c_str());
    }


  catNew = new RooCategory("catNew", "catNew");
  std::string catNewLabel;

  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++) {
    for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++) {
      for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++) {
        for(std::vector<std::string>::iterator b=runList.begin();b!=runList.end();b++) {
           catNewLabel = *m + "_" + *c + "_" + *t + "_" + *b;
           catNew->defineType(catNewLabel.c_str());
        }
      }
    }
  }

  reducedlist.add(mB);
  //reducedlist.add(BDTGresponse);
  //RooRealVars created above
  reducedlist.add(mode);
  reducedlist.add(run);
  reducedlist.add(charge);
  reducedlist.add(track);
  reducedlist.add(*catNew);

  if(_genConfs->get("MCsimfit")!="true")
    {
      //inputlist.add(assignedCharge);
    }
  fulllist.add(reducedlist);  
  fulllist.add(inputlist);
  inputlist.add(mB);
  //inputlist.add(BDTGresponse);

  std::cout << "reducedlist:" << std::endl;
  reducedlist.Print();
  std::cout << "inputlist:" << std::endl;
  inputlist.Print();
  std::cout << "fulllist:" << std::endl;
  fulllist.Print();
  
}

int Fitting::LoadDataSet()
{
  cout << "Reached LoadDataSet" << endl;
  // Get paths to data
  Settings dataSettings("Data settings");
  dataSettings.readPairStringsToMap(_genConfs->get("dataSetLists"));
  // create the full dataset
  std::cout << "Creating RooDataSet" << std::endl;
  data = new RooDataSet("data","Data",fulllist);

  // Loop over the categories into which the files are split, and send on to FinalDataSet() 
  // which is where the types and flags are assigned.
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++) {
    for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++) {
    	for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++) {
    		for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++) {

				std::string fullPathAndName = dataSettings.get("pathToData_"+(*m)+"_"+(*t)+"_"+(*c)+"_"+(*a));

				cout << "File path name is " << fullPathAndName << endl;
				TFile* tfile = TFile::Open(fullPathAndName.c_str());
				RooDataSet *dataset=0;

				// IS this what you want it to do? What about when it is not charge separated? Will it still be fine?
				// Make RooDataSet the same if it is correct
				RooDataSet *ds = (RooDataSet*)tfile->FindObjectAny("DS");
					  if(ds)
					  {
						// RooDataSet
						  cout << "RooDataSet" <<endl;
						tfile->cd();
						dataset = FinalDataSet(*m, *t, *c, *a, (RooDataSet*)tfile->FindObjectAny("DS"));
					  }
					  else
					  {
						// TTree
						  cout << "TTree" << endl;
						tfile->cd();
						dataset = FinalDataSet(*m, *t, *c, *a, (TTree*)tfile->Get("DecayTree"));
					  }

					  data->append(*dataset);

					  tfile->Close();
    		}
    	}
    }
  }
  
  cout << "DataSetSize: " << data->numEntries() << endl;
  std::cout << "Printing dataset: " << std::endl;
  data->Print("v");
  std::cout << std::endl;
  return 0;
}

//version passed a TTree
RooDataSet* Fitting::FinalDataSet(const std::string s_mode, const std::string s_track, const std::string s_charge, const std::string s_run, TTree* tree) //, TH2F &_h_diagnose)
{
  if(!tree){ std::cout << "\n THE TREE IS A ZERO POINTER IN "<< s_mode <<" "<< s_track << " "<< s_charge << " " <<  s_run << std::endl; return 0; }

  TString exclusionString;
  std::string masscut = "Bu_D0constKS0constPVconst_M > " + _genConfs->get("fit_limit_low") + " && Bu_D0constKS0constPVconst_M < " + _genConfs->get("fit_limit_high");
/*  std::string charmlesscut = "D0_FDsignificance > " + _genConfs->get("D0fdcut");
  std::string kstmasscut = "abs(Kst_M - 891.66) < " + _genConfs->get("Kstmasscut");
  std::string kshelcut = "abs(Ks_cosHelicityAngleA) > " + _genConfs->get("Kshelcut");
  std::string bdtcut = "BDTG > " + _genConfs->get("Bdtcut");*/
  exclusionString = masscut;// + " && " + charmlesscut + " && " + kstmasscut + " && " + kshelcut + " && " + bdtcut;

  std::cout << "Exclusion string: " << exclusionString << std::endl;
  
  gROOT->cd();		// Copy the TTree into memory otherwise new tree will try to write to the old file
  TTree* reducedtree = (TTree*) tree->CopyTree(exclusionString);
  double bm(0);
  double cla(0);
  reducedtree->SetBranchAddress("Bu_D0constKS0constPVconst_M",&bm);
  reducedtree->SetBranchAddress("BDTG",&cla);

  /*
  TTree* newtree = new TTree("TTT","");
  double bm_new(0);
  double bdt(0);
  newtree->Branch("Bu_D0constKS0constPVconst_M",&bm_new,"Bu_D0constKS0constPVconst_M/D");
  newtree->Branch("BDTG",&bdt,"BDTG/D");
  for(Long64_t n=0; n<reducedtree->GetEntries(); ++n) {
    reducedtree->GetEntry(n);
    bm_new = bm;
    bdt = cla;
    newtree->Fill();
  }
  */

  RooDataSet *input = new RooDataSet(Form("%s_%s_%s_%s",s_mode.c_str(),s_track.c_str(),s_charge.c_str(),s_run.c_str()),
                                     Form("Final %s %s %s %s",s_mode.c_str(),s_track.c_str(),s_charge.c_str(),s_run.c_str()),
                                     (TTree*)reducedtree,inputlist); //,exclusionString);
  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(mode,run,track,charge,*catNew));
  
  // Now loop over the entries in the dataset and make sure all the types are assigned 
  int numEntries = input->numEntries();
  std::cout << "Number of entries in input file: " << tree->GetEntries() << std::endl;
  std::cout << "Number of entries after exclusion: " << numEntries << std::endl;

  mode.setLabel(s_mode.c_str());
  track.setLabel(s_track.c_str());
  charge.setLabel(s_charge.c_str());
  run.setLabel(s_run.c_str());
  RooArgSet datadetails(mode, run, track, charge, *catNew);

  std::string catNewLabel;
  for(int i = 0; i < input->numEntries(); i++)
    {
      catNewLabel = s_mode + "_" + s_charge + "_" + s_track + "_" + s_run;
      catNew->setLabel(catNewLabel.c_str());
      extra->add(datadetails);
    }
  
  input->merge(extra);
  return input;
}

//version passed a RooDataSet
RooDataSet* Fitting::FinalDataSet(const std::string s_mode, const std::string s_track, const std::string s_charge, const std::string s_run, RooDataSet* DS)
{
  if (!DS)
    {
      std::cout << "\n THE DATASET IS A ZERO POINTER IN " << s_mode << " " << s_track << " "<< s_charge << " " <<  s_run << std::endl;
      return 0;
    }

  TString exclusionString;

  std::string masscut = "Bu_D0constKS0constPVconst_M > " + _genConfs->get("fit_limit_low") + " && Bu_D0constKS0constPVconst_M < " + _genConfs->get("fit_limit_high");
  std::string charmlesscut = "D0_FDsignificance > " + _genConfs->get("D0fdcut");
  std::string kstmasscut = "abs(Kst_M - 891.66) < " + _genConfs->get("Kstmasscut");
  std::string kshelcut = "abs(Ks_cosHelicityAngleA) > " + _genConfs->get("Kshelcut");
  std::string bdtcut = "BDTG > " + _genConfs->get("Bdtcut");
  exclusionString = masscut;// + " && " + charmlesscut + " && " + kstmasscut + " && " + kshelcut + " && " + bdtcut;

  std::cout << "Exclusion string: " << exclusionString << std::endl;

  RooDataSet *input = (RooDataSet*)DS->reduce(exclusionString);
  int numEntries = input->numEntries();
  
  std::cout << "Number of entries in input file: " << DS->numEntries() << std::endl;
  std::cout << "Number of entries after exclusion: " << numEntries << std::endl;
  
  // Now loop over the entries in the dataset and make sure all the types are assigned 
  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(mode,run,track,charge,*catNew));
  mode.setLabel(s_mode.c_str());
  track.setLabel(s_track.c_str());
  charge.setLabel(s_charge.c_str());
  run.setLabel(s_run.c_str());
  RooArgSet datadetails(mode, run, track, charge, *catNew);
  
  std::string catNewLabel;
  for(int i = 0; i < input->numEntries(); i++)
    {
      catNewLabel = s_mode + "_" + s_charge + "_" + s_track + "_" + s_run;
      catNew->setLabel(catNewLabel.c_str());
      extra->add(datadetails);
    }
  
  input->merge(extra);
  return input;
}


void Fitting::PrintDataSet(bool verbose)
{
  if(verbose)
    {
      for(int i=0;i<data->numEntries();i++)
        {
          const RooArgSet *row = data->get(i);
          RooRealVar*  v1 = (RooRealVar*)  row->find(mB.GetName());
          RooCategory* c0 = (RooCategory*) row->find("run");
          RooCategory* c1 = (RooCategory*) row->find("track");
          RooCategory* c2 = (RooCategory*) row->find("charge");
          RooCategory* c3 = (RooCategory*) row->find("mode");
     
          std::cout <<"run "<<c0->getLabel()<<", "<<", "<<c1->getLabel()<<": "<<c2->getLabel()<<" "
                    <<"\t mB="<<v1->getVal()<<"\t "<<c3->getLabel()<<"   "<<std::endl;
        }
    }
  //data->table(*cat)->Print("v");
  data->table(*catNew)->Print("v");

  if(genToys=="true"){
    std::cout<<"   NUMBER OF GENERATED TOY EVENTS = "<<data->sumEntries()<<std::endl;
  }else{
    std::cout<<"                            TOTAL = "<<data->sumEntries()<<std::endl;
  }
}

void Fitting::RunFullFit(bool draw=true)
{
  std::cout << std::endl << "In RunFullFit" << std::endl;
  //mB.setBins(65);
  RooSimultaneous* sim = model->getFitPdf();
  std::cout << "Simultaneous PDF:" << std::endl;
  sim->Print();
  std::cout << std::endl;

  std::cout << "Fixed parameters" << std::endl;
  std::vector <RooRealVar*> *fixedParams = model->GetFixedParameters();
  //create new RooRealVars
  //this is necessary because saving the default RooRealVars somehow pulls along the PDF dependencies and RooFit complains
  RooArgList fixedVars;
  for (Int_t n_p = 0; n_p < (Int_t)fixedParams->size(); ++n_p)
    {
      std::string fixedVarName = std::string(fixedParams->at(n_p)->GetName());
      Double_t fixedVarVal = fixedParams->at(n_p)->getVal();
      RooRealVar *fixedVar = new RooRealVar(fixedVarName.c_str(), fixedVarName.c_str(), fixedVarVal);
      fixedVar->setConstant(kTRUE);
      std::cout << " " << fixedVarName << " " << fixedVarVal << std::endl;
      fixedVars.add(*fixedVar);
    }
  
  RooFitResult *result = 0;  
  Int_t nFloatParams = 0;
  if (doFit == "true")
    {
      std::cout << "Start fitTo" << std::endl;

      //std::vector< RooRealVar* > listlowcomb = model->getMyLowComb();

      std::cout << "About to actually run 'fitTo' " << std::endl;
      //result = sim->fitTo(*data, RooFit::Save(), RooFit::Extended());
      //RooFitResult* result = sim->fitTo(*dataBinned, RooFit::Save(), RooFit::Extended());
      //, RooFit::Hesse(false), RooFit::InitialHesse(true));//,RooFit::Strategy(0),RooFit::PrintLevel(1),RooFit::Warnings(true));
      //result = sim->fitTo(*data, RooFit::Save(), RooFit::Extended(), RooFit::Minos(true));
      result = sim->fitTo(*data, RooFit::Save(), RooFit::Extended(), RooFit::Minos(true), RooFit::Strategy(2));


      std::cout << "End fitTo" << std::endl;


      // Get mean in order to calculate yields and purities in Bu region
      RooArgList allPars = result->floatParsFinal();
      int listindex = allPars.index("bu_mean");
      double mean = ((RooRealVar*)allPars.at(listindex))->getVal();
      double integ_limit_low  = mean - _genConfs->getD("integ_range_low");
      double integ_limit_high = mean + _genConfs->getD("integ_range_high");
      model->printYieldsAndPurities("Bu",integ_limit_low,integ_limit_high,result);
      integ_limit_low = _genConfs->getD("fit_limit_low");
      integ_limit_high = _genConfs->getD("fit_limit_high");
      model->printYieldsAndPurities("full",integ_limit_low,integ_limit_high,result);
      //model->printYieldsAndPurities("binnedfit",5200,5800);


      std::cout <<"\n-------- FULL FIT, ALL FLOATING ---------"<< std::endl;
      std::cout <<"\n-------- minNLL = " << double(result->minNll()) << std::endl;
      result->Print("v");

      if (_genConfs->get("makeCorrelationMatrix")=="true")
        {
          std::cout << "Floating parameters" << std::endl;
          RooArgList varNames = result->floatParsFinal();
          TIter floatedFinalPars(result->floatParsFinal().createIterator());
          RooRealVar* par;
          std::vector<std::string> v_varNames;
          while ((par=(RooRealVar*)floatedFinalPars()))
            {
              std::cout << "found variable " << par->GetName() << " in rooarglist at entry " << varNames.index(par->GetName())
                        << " with value " << par->getVal() << " +/- " << par->getError() << std::endl;
              v_varNames.push_back(par->GetName());
            }
          std::cout << "Correlation and covariance matrices" << std::endl;
          std::cout << std::endl << "Correlation matrix (symmetric)" << std::endl << std::endl;
          result->correlationMatrix().Print();
          std::cout << std::endl << "Covariance matrix (symmetric)" << std::endl << std::endl;
          result->covarianceMatrix().Print();
          
          const int varNum = varNames.getSize();
          TVectorD varValues(varNum);
          TVectorD varErrors(varNum);
          TMatrixD correlationMatrix(varNum,varNum);
          TMatrixD covarianceMatrix(varNum,varNum);
      
          //extract correlation and covariance matrices directly
          TMatrixDSym corrmat = result->correlationMatrix();
          TMatrixDSym covmat = result->covarianceMatrix();
          //transfer from TMatrixDSym into TMatrixD for compatibility
          int nrows = corrmat.GetNrows(); 
          for (int n_r = 0; n_r < nrows; ++n_r)
            {
              varValues(n_r) = ((RooRealVar*)result->floatParsFinal().at(n_r))->getVal();
              varErrors(n_r) = ((RooRealVar*)result->floatParsFinal().at(n_r))->getError();
              for (int n_c = 0; n_c < nrows; ++n_c)
                {
                  correlationMatrix(n_r, n_c) = corrmat(n_r, n_c);
                  covarianceMatrix(n_r, n_c) = covmat(n_r, n_c);
                }
            }
          std::cout << std::endl << "Correlation matrix" << std::endl << std::endl;
          correlationMatrix.Print();
          std::cout << std::endl << "Covariance matrix" << std::endl << std::endl;
          covarianceMatrix.Print();
          std::cout << std::endl;

          // Save it
          TFile* f_fitOutput = new TFile("output/persistFitOutput.root","RECREATE");
          f_fitOutput->cd();
          correlationMatrix.Write("CorrelationMatrix");
          covarianceMatrix.Write("CovarianceMatrix");
          varValues.Write("FitParameterValues");
          varErrors.Write("FitParameterErrors");
          varNames.Write("ListVariablesFloated");
          fixedVars.Write("ListVariablesFixed");
          f_fitOutput->Close();
        }

    } // end of doFit = true

  if (!draw) return;

  //===========================//
  //=== Drawing starts here ===//
  //===========================//

  std::cout<<" going to draw projections "<<std::endl;

  bool drawpulls = false;

  //create canvases
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++)
    {
      std::map< std::string, std::map<std::string, TCanvas* > > v_canvas;
      std::map< std::string, std::map<std::string, TCanvas* > > v_canRes;
      std::map< std::string, std::map<std::string, TCanvas* > > v_canvaslog;
      for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++) {
        for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++) {

        	const char* identifier = Form("%s_%s_%s",(*m).c_str(),(*t).c_str(),(*a).c_str());

        	// create canvas for normal data and fit
        	TCanvas* canvas = new TCanvas(Form("canvas_%s",identifier),Form("%s",identifier),30,30,500,chargeList.size()*250);
        	if (drawpulls) canvas->Divide(1,chargeList.size()*2);
        	else           canvas->Divide(1,chargeList.size());
        	v_canvas[*t][*a] = canvas;

        	// create canvas for residual distributions- these are actually not plotted at the moment
        	// Are they really needed if they can be chosen to be plotted on the fits canvas?
        	TCanvas* canRes = 0;
        	if(doFit=="true"){
				canRes = new TCanvas(Form("canres_%s",identifier),Form("%s",identifier),30,30,500,chargeList.size()*250);
				canRes->Divide(1,chargeList.size());
				v_canRes[*t][*a] = canRes;
        	}

        	// create canvas for log plot of data and fit
        	TCanvas* canvaslog = new TCanvas(Form("canvaslog_%s",identifier),Form("%s",identifier),30,30,500,chargeList.size()*250);
        	canvaslog->Divide(1,chargeList.size());
        	v_canvaslog[*t][*a] = canvaslog;
        }
      }


      std::cout<<" canvases made "<<std::endl;

      RooHist* residualHist = 0;
      float maxH = 0.;
      std::map<std::string,std::map<std::string,std::map<std::string,RooPlot*> > > plot;

      //draw each fit
      for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
        for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++){
          for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++){

            int numbins = 40;
            plot[*c][*t][*a] = mB.frame(RooFit::Bins(numbins));
            std::string tag = Form("%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str());
            catNew->setLabel(tag.c_str());

            //////////////////////////////////////////////////////
            // Plot the separated, standard LL,DD plots
            //////////////////////////////////////////////////////
            data->plotOn(plot[*c][*t][*a],
                         RooFit::Cut(Form("catNew==catNew::%s",tag.c_str())),
                         RooFit::MarkerStyle(6));
				
            std::cout<<" data to plot: "<<*m<<" "<<*c<<" "<<*t<<" "<<*a<<std::endl;

            if(drawProjections=="true"){
              //plot total PDF
              sim->plotOn( plot[*c][*t][*a],RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*data),
                             RooFit::LineWidth(3) );
              if(v_canRes[*t][*a]) {
            	  cout << "Make residuals: " <<*m<<" "<<*c<<" "<<*t<<" "<<*a<<std::endl;
                residualHist = plot[*c][*t][*a]->pullHist(); residualHist->GetYaxis()->SetNdivisions(515);residualHist->setYAxisLimits(-5,5);
              }
              //Bu - CB
              std::cout<<" plotting Bu "<<std::endl;
              sim->plotOn( plot[*c][*t][*a],RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*data),
                           RooFit::Components(Form("myCrystalBall_%s_bu_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str())),
                           RooFit::LineStyle(kSolid),RooFit::LineColor(kRed), RooFit::LineWidth(3), RooFit::Name("sig") );

              if(_genConfs->get("MCsimfit")!="true") {
              // Combinatoric
              std::cout <<" plotting combinatoric "<<std::endl;
              sim->plotOn( plot[*c][*t][*a],RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*data),
                           RooFit::Components(Form("Exponential_%s_exp_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str())),
                           RooFit::LineStyle(kDotted), RooFit::LineColor(kMagenta), RooFit::LineWidth(3), RooFit::Name("comb"));

              //Bu -> D*K* - regexp to pick up also version split by helamp
              std::cout<<" plotting Bu -> D*K*  "<<std::endl;
              sim->plotOn( plot[*c][*t][*a],RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*data),
                           RooFit::Components(Form("PartRecoDstKst_%s_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str())),
                           //RooFit::Components(Form("PartRecoDstKst_%s_bu_001_%s_%s_%s, PartRecoDstKst_%s_bu_010_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str(),
                           //                                                                                               (*m).c_str(),(*c).c_str(),(*t).c_str(),(*a).c_str())
                           //),
                           RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack), RooFit::LineWidth(3), RooFit::Name("partreco"));
              }
            }

            //plot total PDF again
            sim->plotOn( plot[*c][*t][*a],RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*data),
            		RooFit::LineWidth(3) );
            //std::cout << "RooPlot chi^2: " << plot[*c][*t][*a]->chiSquare() << std::endl; //wrong - need number of floating parameters

            if (result)
            {
            	Double_t chisquareperDoF = plot[*c][*t][*a]->chiSquare(nFloatParams);
            	Int_t nBins = plot[*c][*t][*a]->GetNbinsX();
            	Int_t nDoF = nBins - nFloatParams;
            	Double_t chisquare = chisquareperDoF * (Double_t)nDoF;
            	std::cout << "RooPlot chi^2/DoF: " << chisquareperDoF << std::endl;
            	std::cout << "Number of bins: " << nBins << std::endl;
            	std::cout << "Number of degrees of freedom: " << nDoF << std::endl;
            	std::cout << "chi^2: " << chisquare << std::endl;
            	std::cout << "Probability: " << TMath::Prob(chisquare, nDoF) << std::endl;
            }

            int ipad = 0;
            if(*c==both)
            {
            	ipad = 1;
            	v_canvas[*t][*a]->cd(ipad);
            	if(_genConfs->get("setLogScale")=="true") gPad->SetLogy();
            	gPad->SetTicks(1, 1);//upper and right-hand ticks
            	//resize pad if drawing pulls (makes pads too small if not drawing pulls)
            	//if (drawpulls) gPad->SetPad(Form("fit_%s_%s",t->c_str(),a->c_str()),"",0,0.0,1,1,0.3,0,-1);

            }
            if(*c==plus)
            {
            	ipad = 1;
            	v_canvas[*t][*a]->cd(ipad);
            }
            if(*c==minus)
            {
            	// Need to make this different whether pulls are plotted or not
            	ipad = 2;
            	if(drawpulls) ipad = 3;
            	v_canvas[*t][*a]->cd(ipad);
            }
            plot[*c][*t][*a]->SetMinimum(0.1);
            //if(*t=="LL") plot_combLLDD[*c][*a]->SetMinimum(0.1);
            plot[*c][*t][*a]->Draw();
            // Add yields and purities to the plots
            // Cannot do yet as yields are of type RooAbsArg so cannot call getVal()
            char num[100];
            sprintf(num,"%.1f",model->plotNums[*m][*c][*t][*a]["val"]); TString ts_yield(num);
            sprintf(num,"%.1f",model->plotNums[*m][*c][*t][*a]["err"]); TString ts_yielderr(num);
            sprintf(num,"%.1f",model->plotNums[*m][*c][*t][*a]["purity_val"]*100); TString ts_purity(num);
            sprintf(num,"%.1f",model->plotNums[*m][*c][*t][*a]["purity_err"]*100); TString ts_purityerr(num);
            //write signal yield and purity on plots (not done by default)
            // TLatex *latex_yield = new TLatex(5350,plot[*c][*t][*a]->GetMaximum()*0.70,"N_{Sig} = "+ts_yield+" #pm "+ts_yielderr);
            // if(_genConfs->get("setLogScale")!="true") latex_yield->Draw();
            // TLatex *latex_purity= new TLatex(5350,plot[*c][*t][*a]->GetMaximum()*0.60,"Purity = "+ts_purity+" #pm "+ts_purityerr+" %");
            // if(_genConfs->get("setLogScale")!="true") latex_purity->Draw();
            plot[*c][*t][*a]->SetTitle("");
            plot[*c][*t][*a]->SetXTitle("m(DK^{*}) [MeV/#it{c}^{2}]");
            double binwidth = (_genConfs->getD("fit_limit_high") - _genConfs->getD("fit_limit_low")) / numbins;
            plot[*c][*t][*a]->SetYTitle(Form("Candidates / (%.1f MeV/#it{c}^{2})",binwidth));
            //if(*t=="LL") plot_combLLDD[*c][*a]->SetTitle("");
            plot[*c][*t][*a]->GetXaxis()->SetTitleOffset(1.05);
            plot[*c][*t][*a]->GetYaxis()->SetTitleOffset(1.3);

            //draw legend
            Double_t legtop = 0.75;
            TLegend *leg = new TLegend(0.65, 0.35, 0.85, legtop);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(132);
            leg->AddEntry((TObject*)0,"","");
            leg->AddEntry(plot[*c][*t][*a]->findObject("sig"),"B_{u} #rightarrow D^{0}K*","l");
            if(_genConfs->get("MCsimfit")!="true") {
            leg->AddEntry((TObject*)0,"","");
            leg->AddEntry(plot[*c][*t][*a]->findObject("partreco"),"B #rightarrow D^{*}K^{*}","l");
            leg->AddEntry((TObject*)0,"","");
            leg->AddEntry(plot[*c][*t][*a]->findObject("comb"),"Combinatorial","l");
            }

            //do not draw legend if a log plot or pulls are drawn
            //if(_genConfs->get("setLogScale")!="true" && !drawpulls) leg->Draw();
            if(_genConfs->get("setLogScale")!="true") leg->Draw();

            //do not draw title
            /*
                TPaveLabel *pav = new TPaveLabel(0.13,0.9,0.97,0.99,title[*m][*c][*t][*a].c_str(),"NDC");
                pav->SetBorderSize(0); pav->SetFillStyle(0);
                pav->SetTextFont(132); pav->SetTextSize(0.8);
                pav->Draw();
             */
            //if(genToys=="false") lhcbpreliminary->Draw();
            if(plot[*c][*t][*a]->GetMaximum()>maxH)
            {
            	maxH=plot[*c][*t][*a]->GetMaximum();
            }

            // Draw pulls
            if(residualHist && drawpulls)
            {
            	std::cout<<" plotting pulls: "<<*m<<" "<<*c<<" "<<*t<<" "<<*a<<std::endl;
            	v_canvas[*t][*a]->cd(ipad+1);
            	//gPad->SetPad(Form("fit_%s_%s",t->c_str(),a->c_str()),"",0.,0.0,1.0,0.3,0,0,-1);
            	RooPlot* frame = mB.frame(RooFit::Title("Residual Distribution"));
            	frame->GetYaxis()->SetNdivisions(515);//,kTrue);
            	frame->SetXTitle("");//#it{m}(DK^{#pm}#pi^{#mp}) (MeV/c^{2})");
            	residualHist->SetFillColor(4);
            	residualHist->SetLineColor(4);
            	frame->addPlotable(residualHist,"BEX0");
            	frame->Draw();
            	frame->SetTitle("");
            }

            // Save data and pdf to draw projections in external macro
            saveOutputForPlottingMacro->cd();
            //if(*t=="LL") plot_combLLDD[*c][*a]->Write(Form("%s_%s_comb_%s",m->c_str(),c->c_str(),a->c_str()));
            plot[*c][*t][*a]->Write(Form("%s_%s_%s_%s",m->c_str(),c->c_str(),t->c_str(),a->c_str()));

            //=== log plots ===//

            ipad = 0;
            if(*c==both)
            {
            	ipad = 1;
            	v_canvaslog[*t][*a]->cd(ipad);
            	gPad->SetLogy();
            	gPad->SetTicks(1, 1);//upper and right-hand ticks
            }
            if(*c==plus)
            {
            	ipad = 1;
            	v_canvaslog[*t][*a]->cd(ipad);
            	gPad->SetLogy();
            }
            if(*c==minus)
            {
            	ipad = 2;
            	v_canvaslog[*t][*a]->cd(ipad);
            	gPad->SetLogy();
            }
            plot[*c][*t][*a]->Draw();
          }
        }
      }

      //fonts and y-axis minimum
      for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
        for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++){
          for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++){
            plot[*c][*t][*a]->SetTitleFont(132,"X"); plot[*c][*t][*a]->SetLabelFont(132,"X");
            plot[*c][*t][*a]->SetTitleFont(132,"Y"); plot[*c][*t][*a]->SetLabelFont(132,"Y");
            if(_genConfs->get("setLogScale")=="true") plot[*c][*t][*a]->SetMinimum(0.100001);
            else                                      plot[*c][*t][*a]->SetMinimum(0.000001);
          }
        }
      }
    
      //save plots
      for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++){
        for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++){
          if(_genConfs->get("setLogScale")=="true")
            {
              v_canvas[*t][*a]->Print(Form("figs/fits/%s_log.pdf",v_canvas[*t][*a]->GetName()));
              v_canvas[*t][*a]->Print(Form("figs/fits/%s_log.C",v_canvas[*t][*a]->GetName()));
              if(v_canRes[*t][*a]) v_canRes[*t][*a]->Print(Form("figs/residuals/%s_log.pdf",v_canRes[*t][*a]->GetName()));
              if(v_canRes[*t][*a]) v_canRes[*t][*a]->Print(Form("figs/residuals/%s_log.C",v_canRes[*t][*a]->GetName()));
            }
          else
            {
              v_canvas[*t][*a]->Print(Form("figs/fits/%s.pdf",v_canvas[*t][*a]->GetName()));
              //v_canvas[*t][*a]->Print(Form("figs/fits/%s.eps",v_canvas[*t][*a]->GetName()));
              v_canvas[*t][*a]->Print(Form("figs/fits/%s.C",v_canvas[*t][*a]->GetName()));
              if(v_canRes[*t][*a]) v_canRes[*t][*a]->Print(Form("figs/residuals/%s.pdf",v_canRes[*t][*a]->GetName()));
              if(v_canRes[*t][*a]) v_canRes[*t][*a]->Print(Form("figs/residuals/%s.C",v_canRes[*t][*a]->GetName()));
            }
        }
      }
      
      //=== log plots ===//

      //fonts and y-axis minimum
      for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
        for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++){
          for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++){
            //plot[*c][*t][*a]->SetMinimum(0.100001);
            plot[*c][*t][*a]->SetMinimum(0.001);
          }
        }
      }
      for(std::vector<std::string>::iterator t=trackList.begin();t!=trackList.end();t++){
        for(std::vector<std::string>::iterator a=runList.begin();a!=runList.end();a++){
          v_canvaslog[*t][*a]->Print(Form("figs/fits/%s_log.pdf",v_canvaslog[*t][*a]->GetName()));
          v_canvaslog[*t][*a]->Print(Form("figs/fits/%s_log.C",v_canvaslog[*t][*a]->GetName()));
        }
      }
    }
  
}


void Fitting::RunManyFits()
{
  unsigned int n = 10;
  std::string filebase = "./fit_";
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){filebase+=(*m)+underscore;}
  //  model->R_ADS_b2dk_d2pik_b->setConstant();
  //model->A_ADS_b2dk_d2pik_b->setConstant();
  for(unsigned int i=0;i<n;i++){
    int secs=time(NULL);
   
    //change stuff around here for systematics i guess
    //model->R_ADS_b2dk_d2pik_b->setVal(gRandom.Uniform( 0.,0.03));
    //model->A_ADS_b2dk_d2pik_b->setVal(gRandom.Uniform(-1.,1.));
    std::string autofile = filebase+std::string(Form("%i.root",secs));  
    TFile f(autofile.c_str(),"RECREATE");
    RunFullFit(false);
    f.Close();
  }
}

void Fitting::RunManyToys()
{
  NewOrderToys(_genConfs->getI("nToys"));
}


/********* TOYS **********/


void Fitting::OrderToys(int n)
{
  RooAbsPdf* genPdf = model->getGenPdf();
  RooAbsPdf* fitPdf = model->getFitPdf();
  bool DB=false; // Switch on debug mode, more verbose output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  RooMCStudy* mcstudy = new RooMCStudy(*genPdf,reducedlist,RooFit::FitModel(*fitPdf),RooFit::FitOptions(RooFit::Save(true),RooFit::Extended(true),RooFit::NumCPU(_genConfs->getI("numCPUsToUse")),RooFit::PrintLevel(DB?1:-1)));

  int nEvtsPerSample=(int)genPdf->expectedEvents(*catNew);
  if(1==n&&!DB){
    mcstudy->generate(n,nEvtsPerSample,true);
    data=(RooDataSet*)mcstudy->genData(n-1);
  }else{
    std::cout<<"Generating "<<n<<" toys of "<<nEvtsPerSample<<" events."<<std::endl;
    mcstudy->generateAndFit(n,nEvtsPerSample);
    std::cout<<"Finished generating toys"<<std::endl;
    if(DB){return;}
    
    gSystem->Exec((std::string("mkdir -p ")+_genConfs->get("toyLocation")).c_str());
    std::string toyfile=_genConfs->get("toyLocation")+"/toy_";
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){toyfile+=(*m)+underscore;}
    toyfile+=std::string(Form("%i.root",_genConfs->getI("startSeed")));  

    std::string textfile = _genConfs->get("toyLocation")+"/toy_";
    textfile+=std::string(Form("%i.txt",_genConfs->getI("startSeed")));  

    ofstream resfile(textfile.c_str());

    TFile f(toyfile.c_str(),"RECREATE");
    for(int i=0;i<n;i++){
      mcstudy->fitResult(i)->Write(Form("toy%i",i));
      resfile<<i<<'\n';

      TIter finalPar(mcstudy->fitResult(i)->floatParsFinal().createIterator());
      RooRealVar * par;
      while ((par = (RooRealVar *)finalPar())) {
        double trueVal=0;
        TIter initPar(mcstudy->fitResult(i)->floatParsInit().createIterator());
        RooRealVar *initpar;
        while((initpar=(RooRealVar *)initPar())){
          if(0==strcmp(par->GetName(),initpar->GetName())){ 
            trueVal = initpar->getVal();
          }
        }
	
        resfile << par->GetName() << ' '
                << par->getVal() << ' '
                << par->getError() << ' '
                <<mcstudy->fitResult(i)->covQual()<< ' '
                << _genConfs->getI("startSeed") << ' '
                << i << ' '
                << trueVal << '\n';
      }
    }
    mcstudy->fitParDataSet().Write("toyPars");
    gDirectory->ls();
    resfile.close();	
    f.Close();
		
  }
}

void Fitting::NewOrderToys(int n)
{
  // Only use RooMCStudy when only running a few toys locally
  // For many toys submitted in parallel RooMCStudy won't work
  bool roomcstudy = false;

  gSystem->Exec((std::string("mkdir -p ")+_genConfs->get("toyLocation")).c_str());
  std::string textfile = _genConfs->get("toyLocation")+"/toy_";
  textfile+=std::string(Form("%i.txt",_genConfs->getI("startSeed")));  

  ofstream resfile(textfile.c_str());

  RooAbsPdf* genPdf = model->getGenPdf();
  RooAbsPdf* fitPdf = model->getFitPdf();
  RooAbsPdf* sim_template = (RooAbsPdf*)fitPdf->Clone("Copy of the initial state of the simultaneous pdf");
  RooArgSet* fittedPars = sim_template->getParameters(RooArgSet(mB,mode,charge,track,run));
  RooArgSet* frozen = (RooArgSet*) fittedPars->snapshot(kTRUE); 

  //from here i should save the fitPdf state  

  bool DB=false;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  RooMCStudy* mcstudy = new RooMCStudy(*genPdf,reducedlist,
                                       RooFit::FitModel(*fitPdf),
                                       RooFit::Extended(true),
                                       RooFit::FitOptions(RooFit::Save(true),RooFit::Extended(true),RooFit::NumCPU(_genConfs->getI("numCPUsToUse")),RooFit::PrintLevel(DB?1:-1)));

  int nEvtsPerSample=(int)genPdf->expectedEvents(*catNew);
  
  std::cout<<"Generating "<<n<<" toys of "<<nEvtsPerSample<<" events."<<std::endl;

  // Cannot use RooMCStudy if submitting many jobs to speed up toys
  if(roomcstudy) {
  // ========= Using RooMCStudy ============= 
  mcstudy->generateAndFit(n,nEvtsPerSample,true);
  fittedPars->Print("v");

  TFile* f_dsoutput_fitpars = new TFile((_genConfs->get("toyLocation")+"/toyFitPars.root").c_str(),"recreate");
  RooDataSet ds_fitpars = mcstudy->fitParDataSet();
  //RooDataSet* ds_genpars = (RooDataSet*) mcstudy->genParDataSet();
  //ds_genpars->Print("v");
  //ds_genpars->Write();
  ds_fitpars.Write();
  f_dsoutput_fitpars->Close();

  // Loop over all variables and use RooMCStudy's plotting tool
  TIterator* fittedParIter = fittedPars->createIterator();
  RooRealVar* fitpar;
  std::vector <RooRealVar*> *fixedParams = model->GetFixedParameters();
  while((fitpar = (RooRealVar*)fittedParIter->Next()))
  {
    TString vname = fitpar->GetName();
    if(vname.Contains("catNew")) continue;
    bool isfixed=0;
    // Do not want to work out pulls for fixed params
    // All parameters relating to partreco are fixed but not included in fixedParams
    if(vname.Contains("gamma_frac")) isfixed=1;
    if(vname.Contains("Bu_DstKst")) isfixed=1;
    if(vname.Contains("HORNS")) isfixed=1;
    if(vname.Contains("HILL")) isfixed=1;
    if(vname.Contains("low")) isfixed=1;
    if(vname.Contains("MC")) isfixed=1;
    // Set parameters in fixedParams to fixed
    for(int n=0; n<(int)fixedParams->size(); ++n){
      TString myvar = fixedParams->at(n)->GetName();
      if(vname==myvar) isfixed=1;
    }
    if (!isfixed) {
      TCanvas* cvar = new TCanvas("cvar","",1200,400);
      cvar->Divide(3);
      cvar->cd(1);
      RooPlot* parframe = mcstudy->plotParam(*fitpar);
      parframe->Draw();
      cvar->cd(2);
      RooPlot* errframe = mcstudy->plotError(*fitpar);
      errframe->Draw();
      cvar->cd(3);
      // Plotting the pulls will only work if you have the exact same parameters (and param names) in the genPdf
      // Note that Contains("n_") includes tail parameter of cb (does not matter as it is fixed anyway
      if (!(vname.Contains("n_")) && !(vname.Contains("ratio_"))){
      RooPlot* pullframe = mcstudy->plotPull(*fitpar,-5,5,30,kTRUE);//RooFit::Bins(30)); //,-5,5,30); //,kTRUE);
      pullframe->Draw();
      }
      cvar->SaveAs("TOYS/"+vname+".pdf");
      delete cvar;
    }
  }
  // =========== End using RooMCStudy =============
  }

  if(!roomcstudy) {
  mcstudy->generate(n,nEvtsPerSample,true);
  if(1==n&&!DB){
    data=(RooDataSet*)mcstudy->genData(n-1);
  }
  else{
    RooAbsPdf * toyFitPdf = 0; 
    RooFitResult* result = 0;

    for(int i=0; i<n; i++){
      std::cout<<" I'm going to generate the " << i << "th toy data set and fit pdf!" << std::endl;

      RooAbsData* toyData = (RooAbsData*)mcstudy->genData(i);
      toyFitPdf = model->getFitPdf();  // back in the for loop because GC's implemented here

      // Reset fit params to defaults (not the last fit result)
      RooArgSet* initialPars= toyFitPdf->getParameters(RooArgSet(mB,mode,charge,track,run));
      TIterator* initialIter= initialPars->createIterator();
      TObject* initialObj;  TObject* fittedObj;
      while((initialObj=initialIter->Next())) {
        RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
        if(!initialRrv) continue;
        TIterator* fittedIter = frozen->createIterator();
        while((fittedObj=fittedIter->Next())) {
          RooRealVar* fittedRrv = dynamic_cast<RooRealVar*>(fittedObj);
          if(!fittedRrv) continue;
          if(std::string(initialRrv->GetName())==std::string(fittedRrv->GetName())){
            if(i<5) std::cout<<" "<<initialRrv->GetName()<<" <--- "<<initialRrv->getVal()<<" (default: "<<fittedRrv->getVal()<<")"<<std::endl;
            initialRrv->setVal(fittedRrv->getVal());
          }
        }
        delete fittedIter;
      }
      delete initialPars;
      delete initialObj;
      delete fittedObj;
      delete initialIter;

      toyData->table(*catNew)->Print("v");

      std::cout<<" Sneha - Moving on to Toy "<<i<<std::endl;
      
      result = toyFitPdf->fitTo(*toyData,RooFit::Save(),RooFit::Extended(true), RooFit::PrintLevel(DB?1:-1));
      std::cout<<" the rooo fit status is "<<result->status()<<std::endl;
      // Fit with HESSE
      //result = toyFitPdf->fitTo(*toyData,RooFit::Save(),RooFit::Extended(true), RooFit::Hesse(true),RooFit::PrintLevel(DB?1:-1));
      //if(result->covQual()<3 || result->status()!=0) {
      //	result = toyFitPdf->fitTo(*toyData,RooFit::Save(),RooFit::Extended(true), RooFit::Hesse(true),RooFit::PrintLevel(DB?1:-1));
      //}
      // Fit with MINOS
      result = toyFitPdf->fitTo(*toyData,RooFit::Save(),RooFit::Extended(true), RooFit::Minos(true),RooFit::PrintLevel(-1), RooFit::Strategy(2));
      if(result->covQual()<3 || result->status()!=0) {
      	result = toyFitPdf->fitTo(*toyData,RooFit::Save(),RooFit::Extended(true), RooFit::Minos(true),RooFit::PrintLevel(-1), RooFit::Strategy(2));
      }
 
      std::cout<<" the rooo fit status is "<<result->status()<<std::endl;
      std::cout<<" the EDM is " << result->edm() << std::endl;
      
      //Writing the results file
      resfile<<i<<'\n';
      result->floatParsFinal().Print() ;
      TIter initPar(result->floatParsInit().createIterator());
      TIter finalPar(result->floatParsFinal().createIterator());
      RooRealVar * par;
      while ((par = (RooRealVar *)finalPar())) {
        //std::cout << par->GetName() << " " << par->getVal() << std::endl;
      	double trueVal=0;
	      RooRealVar *initpar;
	      while((initpar=(RooRealVar *)initPar())){
  	      if(0==strcmp(par->GetName(),initpar->GetName())){ 
            trueVal = initpar->getVal();  // doesn't always work because Pdf_Gen might be different from Pdf_Fit
            //std::cout << initpar->GetName() << " " << trueVal << std::endl;
          }
        }
        resfile << par->GetName() << ' '
                << par->getVal() << ' '
                << par->getError() << ' ' // HESSE
                //<< par->getAsymErrorLo() << ' '
                //<< par->getAsymErrorHi() << ' '
                << result->covQual()<< ' '
                << _genConfs->getI("startSeed") << ' '
                << i << ' '
                << result->status() << std::endl;
                //<< result->edm() << std::endl;
        delete initpar;
      }
      delete par;
/*
      if(i==0) {
        std::cout << "Now making a plot" << std::endl;
        catNew->setLabel(Form("%s_%s_%s_%s","d2kpi","both","mix","all"));

        RooPlot* frame = mB.frame(RooFit::Bins(40));
        toyData->plotOn(frame,RooFit::Cut(Form("catNew==catNew::%s_%s_%s_%s",
                                               "d2kpi","both","mix","all")));
                                               //"d2kpi","both","DD","all")));
        toyFitPdf->plotOn(frame,RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*toyData));
        TCanvas* ctoy = new TCanvas("ctoy","",1000,700);
        frame->Draw();
        string filename = _genConfs->get("toyLocation");
        filename += Form("/toy_%i_mix.pdf",_genConfs->getI("startSeed"));
        ctoy->SaveAs(filename.c_str());
        delete ctoy;
        delete frame;
      }

      if(i==0) {
        std::cout << "Now making a plot" << std::endl;
        catNew->setLabel(Form("%s_%s_%s_%s","d2kpi","both","LL","all"));

        RooPlot* frame = mB.frame(RooFit::Bins(40));
        toyData->plotOn(frame,RooFit::Cut(Form("catNew==catNew::%s_%s_%s_%s",
                                               "d2kpi","both","LL","all")));
        toyFitPdf->plotOn(frame,RooFit::Slice(RooArgSet(*catNew)), RooFit::ProjWData(RooArgSet(*catNew),*toyData));
        TCanvas* ctoy = new TCanvas("ctoy","",1000,700);
        frame->Draw();
        string filename = _genConfs->get("toyLocation");
        filename += Form("/toy_%i_LL.pdf",_genConfs->getI("startSeed"));
        ctoy->SaveAs(filename.c_str());
        delete ctoy;
        delete frame;
      }
*/
      std::cout<<" now going to delete some stuff" << std::endl;

      //now to delete all the things created.
      result = 0;
      delete par;
      delete toyData;
      delete toyFitPdf;

      std::cout << "now ready to start next toy loop" << std::endl;
    }
    resfile.close();	
    //gDirectory->ls();
  }
  }
  delete mcstudy;
}



void Fitting::DisplayToys()
{
  RooSimultaneous* genPdf=model->getGenPdf();
  RooSimultaneous* fitPdf=model->getFitPdf();
 	
  // First load the toys into the RooMCStudy and find some general success statistics	
  RooMCStudy* mcstudy = new RooMCStudy(*genPdf,reducedlist,RooFit::FitModel(*fitPdf),RooFit::FitOptions("r")); 
  int nFits=1;
  int nBad=0;
  int nFPD=0;
  std::string toyfolder("toySeed_");
  std::string toyfile("toy_");
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){toyfile+=(*m)+underscore;}
  std::cout <<" Looking for files like: " << _genConfs->get("toyLocation")+slash+toyfolder+"*/"+toyfile <<"*"<< std::endl;
  std::vector<std::string> foldername;
  CommonTools::getdir(_genConfs->get("toyLocation").c_str(),foldername);
  for(unsigned int i=0;i<foldername.size();i++){
    if(foldername[i].find(toyfolder) == std::string::npos) continue;
    cout << "Got to element " << i << " in foldername" << endl;
    cout << "Element is: " << foldername[i] << endl;
    // Extract seed from foldername
    int extractedSeed = atoi(foldername[i].substr(8).c_str());
    cout << "I think the seed is " <<  extractedSeed << endl;
    TFile f((_genConfs->get("toyLocation")+slash+foldername[i]+slash+toyfile+Form("%i.root",extractedSeed)).c_str(),"READ");
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    while((key=(TKey*)nextkey())){
      if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
      RooFitResult *r = (RooFitResult*)key->ReadObj();
      if(r->covQual()<2){nBad++;continue;}
      if(r->covQual()<3){nFPD++;continue;}
      std::cout <<"\r   Collecting "<<nFits<<" toys, of which "<<nBad<<" don't converge and "<<100*float(nFPD)/nFits<<"% are forced positive definite."<<std::flush;
      mcstudy->addFitResult(*r);
      nFits++;
    }
  }
  std::cout<<"\n"<<std::endl;
  if(1==nFits){
    std::cout<<"SOMETHING WRONG. NO TOYS LOADED"<<std::endl;
    return;
  }
  gStyle->SetPalette(1) ; //gStyle->SetOptStat(0) ;

  // Now that you've loaded the toys into the RooMCStudy, draw the pulls
  TCanvas* c_all = new TCanvas("canvas_all_pull","Pulls",0,0,800,800) ; c_all->Divide(4,4) ;
  
  int nv=0;
  RooRealVar* var=0;
  RooArgSet* vars = fitPdf->getVariables();
  TIterator* it = vars->createIterator();
  // Create a file to store the pull plots
  TFile* f_pull = new TFile((_genConfs->get("toyLocation")+slash+Form("F_variableValErrPullDists.root")).c_str(),"RECREATE");
  while((var = (RooRealVar*)it->Next())) {
    std::cout<< std::endl<< var->GetName()<<std::flush;
    if(var->isConstant()) continue;
    if(var->InheritsFrom("RooAbsCategory")) continue;
    // For now I want to look at all pull distributions, so can uncomment this later
    if(0==std::string(var->GetName()).find(mB.GetName())) continue;
    // If the variable is not interesting (so not one of the below) then skip
    if((0!=std::string(var->GetName()).find("n_comb_d2kpi_fail_plus_DD_binp"))
       //if((0!=std::string(var->GetName()).find("n_")) and 
       //(0!=std::string(var->GetName()).find("comb_coef")) and
       //(0!=std::string(var->GetName()).find("xplus")) and
       //(0!=std::string(var->GetName()).find("yplus")) and
       //(0!=std::string(var->GetName()).find("xminus")) and
       //(0!=std::string(var->GetName()).find("yminus"))
       //(0!=std::string(var->GetName()).find("effPID")) and 
       //(0!=std::string(var->GetName()).find("alphaL_dpi_")) and 
       //(0!=std::string(var->GetName()).find("alphaR_dpi_")) and 
       //(0!=std::string(var->GetName()).find("MC_"))
       )
      {std::cout << "About to continue" << std::endl; continue;}
    std::cout<<"  ------->  Plotting PULL " << std::endl;
    TDirectory *currDir = f_pull->mkdir(var->GetName());
    currDir->cd();
    RooPlot* frame1 = mcstudy->plotParam(*var,RooFit::Bins(20));frame1->Write();
    cout << "Got past plotParam" << endl;
    RooPlot* frame2 = mcstudy->plotError(*var,RooFit::Bins(20));frame2->Write();
    cout << "Got past plotError" << endl;
    RooPlot* frame3 = mcstudy->plotPull( *var,RooFit::Bins(20),RooFit::Range(-5,5),RooFit::FitGauss(true));
    cout << "Got past plotPull" << endl;
    frame3->SetTitle(var->GetTitle());frame3->Write();
    nv++;
    TCanvas* c = new TCanvas(Form("canv_%s",var->GetName()),var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;
    c->cd(1) ; gPad->SetLeftMargin(0.15); frame1->GetXaxis()->SetNdivisions(5); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15); frame2->GetXaxis()->SetNdivisions(5); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw() ;
    c->cd(3) ; gPad->SetLeftMargin(0.15); frame3->GetXaxis()->SetNdivisions(5); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->getAttText()->SetTextSize(0.001); frame3->Draw() ;
    gStyle->SetTitleFontSize(0.1);
    c_all->cd(nv); gPad->SetLeftMargin(0.15); frame3->Draw();
    ((TPaveText*)gPad->FindObject("pullGauss_paramBox"))->SetTextSize(0.07);
    gPad->Update();
  }
  f_pull->Close();

  std::cout<< std::endl;
  return;
}

