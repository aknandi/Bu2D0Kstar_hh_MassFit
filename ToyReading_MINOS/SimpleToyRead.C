#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include "SimpleToyRead.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooDataSet.h"

using namespace RooFit;

SimpleToyRead::SimpleToyRead(){

  std::cout<<"Bet i get my plots faster than tha plotOn function. And in a more useful format than a CANVAS! "<<std::endl;
  outputfile.open("bias.txt",std::ofstream::out);

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

  float  var, evar, edm, seed, run, pull, cov, status, evarlo, evarhi;
  //Int_t ncols=0;
  Int_t nlines = 0;
  TFile* f = new TFile(rootfile.c_str(),"RECREATE");

  TTree* tree = new TTree(name_var.c_str(),"");
  tree->Branch("value",&var,"value/F");
  tree->Branch("error_low",&evarlo,"error_low/F");
  tree->Branch("error_high",&evarhi,"error_high/F");
  tree->Branch("truevalue",&trueval,"truevalue/F");
  tree->Branch("pull",&pull,"pull/F");
  tree->Branch("run",&run,"run/F");
  tree->Branch("cov",&cov,"cov/F");
  tree->Branch("edm",&edm,"edm/F");
  tree->Branch("status",&status,"status/F");


  //while (inputfile) {
  while (inputfile  >> var >> evarlo >> evarhi >> cov >> seed >> run >> status >> edm){
    if (nlines < 3) std::cout << Form("var=%5f, errorlow=%5f, errorhigh=%5f, cov=%5f, seed=%5f, run=%5f, init=%5f\n", var,evarlo, evarhi,cov,seed, run, edm);
    if(cov!=3) std::cout << cov << std::endl;

    if( (var-trueval)<0 )
      pull = (var-trueval)/evarhi;
    else
      pull = (var-trueval)/evarlo*-1;

    tree->Fill();
    nlines++;
  }

  std::cout << Form(" found %d points\n",nlines) << std::endl;

  f->cd();
  tree->Write();
  f->Close();

  inputfile.close();
}

void SimpleToyRead::MakeSomePlotsFromRootFile(std::string name_var, float trueval){

  std::cout << "In MakeSomePlots" << std::endl;
  std::string rootfile = "results/"+name_var+"_ntp.root";
  TFile f(rootfile.c_str(),"READ");
  TTree *mytree = (TTree*)f.Get(name_var.c_str());

  float val=0; float eval=0; float valp=0; float cov=0;
  float evarlo=0; float evarhi=0;
  float edm=0; float status=0;

  mytree->SetBranchAddress("value",&val);
  mytree->SetBranchAddress("error_low",&evarlo);
  mytree->SetBranchAddress("error_high",&evarhi);
  mytree->SetBranchAddress("pull",&valp);
  mytree->SetBranchAddress("cov",&cov);
  mytree->SetBranchAddress("edm",&edm);
  mytree->SetBranchAddress("status",&status);
  mytree->GetEntry(0);
  TH1F * var = new TH1F("var",name_var.c_str(),60,val-10*evarhi,val+10*evarhi);
  var->SetBit(TH1::kCanRebin);
  TH1F * err = new TH1F("err",(name_var+"_err").c_str(),60,2*evarlo,2*evarlo*-1);
  TH1F * pull = new TH1F("pull",(name_var+"_pull").c_str(),40,-5,5);


  for(int i=0; i<mytree->GetEntries(); ++i){

    mytree->GetEntry(i);

    //   if(eval<0.02)continue;
    if(cov!=3)continue;
    
    if (status!=0) continue;

    var->Fill(val);
    err->Fill(evarlo);
    err->Fill(evarhi);
    pull->Fill(valp);
  }

  if(pull->GetEntries()>0)
    pull->Fit("gaus", "QE");

  TCanvas *c1 = new TCanvas(name_var.c_str(),name_var.c_str(),0,0,1800,500);
  c1->Divide(3,1);
  c1->cd(1);
  var->Draw();
  c1->cd(2);
  err->Draw();
  c1->cd(3);
  pull->Draw("E");

  c1->Print(("results/"+name_var+"_plots.eps").c_str());

  delete c1;

  RooRealVar rrv_pull("pull","pull",-5,5);
  RooRealVar rrv_status("status","status",0,6);

  RooDataSet data("data","",mytree,RooArgSet(rrv_pull,rrv_status),"status==0");
  RooRealVar mean("mean","mean",0,-5,5);
  RooRealVar sigma("sigma","sigma",1,0,5);
  RooGaussian gaus("gaus","",rrv_pull,mean,sigma);

  RooFitResult* result =  gaus.fitTo(data,Save(),Hesse(),PrintLevel(-1),Strategy(2));
  std::cout << result->status() << " " << result->edm() << " " << result->covQual() << std::endl;
  if(result->status()!=0 || result->covQual()!=3 || result->edm()>1){
    std::cout << "doing the fit again" << std::endl;
    mean.setVal(-0.1);
    sigma.setVal(0.11);
    result = gaus.fitTo(data,Save(),Hesse(),PrintLevel(-1),Strategy(2));
  }
  result->Print("v");
  RooArgSet parset = result->floatParsFinal();
  TIterator* parIter = parset.createIterator();
  RooRealVar* par;
  outputfile << name_var << " " ;
  while(par = (RooRealVar*) parIter->Next()){
    outputfile << par->getVal() << " " << par->getError() << " ";
  }
  outputfile << std::endl;
  RooPlot* frame = rrv_pull.frame(Bins(50));
  data.plotOn(frame);
  gaus.plotOn(frame);
  gaus.paramOn(frame,&data);
  frame->SetTitle(name_var.c_str());
  TCanvas *c2 = new TCanvas("c","",0,0,600,500);
  frame->Draw();
 
  c2->Print(("results_roofit/"+name_var+"_plots.eps").c_str());
  delete c2;
}
