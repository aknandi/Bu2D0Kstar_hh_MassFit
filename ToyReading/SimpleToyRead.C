#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include "SimpleToyRead.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>

SimpleToyRead::SimpleToyRead(){

  std::cout<<"Bet i get my plots faster than tha plotOn function. And in a more useful format than a CANVAS! "<<std::endl;
}


void SimpleToyRead::MakeNtupleFromTextFile(std::string name_var, float trueval){

  std::cout << "In MakeNtupleFromTextFile" << std::endl;
  std::string textfile = "results/"+name_var+".txt";
  std::string rootfile = "results/"+name_var+"_ntp.root";
  std::ifstream inputfile;
  inputfile.open(textfile.c_str(),std::ifstream::in);
  //if (fp <= 0){
  //  printf("Sorry, cannot open file \n");
  //  printf("%s",textfile.c_str());
  //  return;
 // }

  std::string tempvar;
  float  var, evar, ivar, seed, run, pull, cov, junk, evarPlus, evarMinus;
  //Int_t ncols=0;
  Int_t nlines = 0;
  TFile* f = new TFile(rootfile.c_str(),"RECREATE");

  TTree* tree = new TTree(name_var.c_str(),"");
  tree->Branch("value",&var,"value/F");
  tree->Branch("error",&evar,"error/F");
  //tree->Branch("errorHigh",&evarHigh,"errorHigh/F");
  //tree->Branch("errorLow",&evarLow,"errorLow/F");
  tree->Branch("truevalue",&trueval,"truevalue/F");
  tree->Branch("pull",&pull,"pull/F");
  tree->Branch("run",&run,"run/F");
  tree->Branch("cov",&cov,"cov/F");

  //while (inputfile) {
  while (inputfile) {
    inputfile >> tempvar >> evar >> evarMinus >> evarPlus >> cov >> seed >> run >> ivar;
    if(tempvar!="nan") {
      var = atof(tempvar.c_str());
    if (nlines < 3) std::cout << Form("var=%5f, error=%5f, cov=%5f, seed=%5f, run=%5f, init=%5f\n", var,evar,cov,seed, run, ivar);
    //if(var>trueval) evarHigh = evar;
    //if(var<trueval) evarLow = evar;

    if(var>trueval) pull = (var-trueval)/evarMinus;
    else pull = (trueval-var)/evarPlus;

    //pull = (var-trueval)/evar;

    tree->Fill();
    //if(nlines==1214) nlines++;
    }
    nlines++;

  }

  std::cout << Form(" found %d points\n",nlines) << std::endl;

  f->cd();
  tree->Write();
  f->Close();

  inputfile.close();
}

double SimpleToyRead::MakeSomePlotsFromRootFile(std::string name_var, float trueval, std::string method){

  std::cout << "In MakeSomePlots" << std::endl;
  std::string rootfile = "results/"+name_var+"_ntp.root";
  TFile f(rootfile.c_str(),"READ");
  TTree *mytree = (TTree*)f.Get(name_var.c_str());

  float val=0; float eval=0; float valp=0; float cov=0; float evalHigh=0; float evalLow=0;

  mytree->SetBranchAddress("value",&val);
  mytree->SetBranchAddress("error",&eval);
  //mytree->SetBranchAddress("errorHigh",&evalHigh);
  //mytree->SetBranchAddress("errorLow",&evalLow);
  mytree->SetBranchAddress("pull",&valp);
  mytree->SetBranchAddress("cov",&cov);

  mytree->GetEntry(1);
  TH1F * var = new TH1F("var",name_var.c_str(),40,val-10*eval,val+10*eval);
  TH1F * err = new TH1F("err",(name_var+"_err").c_str(),40,0,2*eval);
  //TH1F * errHigh = new TH1F("errHigh",(name_var+"_errHigh").c_str(),30,0,2*eval);
  //TH1F * errLow = new TH1F("errLow",(name_var+"_errLow").c_str(),30,0,2*eval);
  TH1F * pull = new TH1F("pull",(name_var+"_pull").c_str(),40,-5.0,5.0);

  //int count =0;

  for(int i=0; i<mytree->GetEntries(); ++i){

    mytree->GetEntry(i);

    //   if(eval<0.02)continue;

    if(cov<3)continue;

    var->Fill(val);
    err->Fill(eval);
    //errHigh->Fill(evalHigh);
    //errLow->Fill(evalLow);
    pull->Fill(valp);
    //if(name_var.compare("adsSignificance")==0) {
    //  if ( val < 1.3 ) count += 1;
    //}
  }

  double systematicError;

  if(method.compare("data")==0) {
    systematicError = var->GetRMS();
  }
  else if(method.compare("toys")==0) {
    systematicError = std::abs(var->GetMean() - trueval);
  }
  else {
    systematicError = 0;
  }

  std::cout << "systematic " << name_var << " " << systematicError << std::endl;
  //std::cout << count << std::endl;

  if(name_var.compare("adsSignificance")==0) {
    var->Fit("gaus","EM");
  }

  if(pull->GetEntries()>0)
    pull->Fit("gaus", "EM");

  TCanvas *c1 = new TCanvas(name_var.c_str(),name_var.c_str(),0,0,1800,500);
  c1->Divide(3,1);
  c1->cd(1);
  var->Draw();
  c1->cd(2);
  err->Draw();
  //c1->cd(3);
  //errHigh->Draw();
  //c1->cd(4);
  //errLow->Draw();
  c1->cd(3);
  pull->Draw("E");

  c1->Print(("results/"+name_var+"_plots.eps").c_str());

  delete c1;

  return systematicError;
}
