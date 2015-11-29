#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include "SimpleToyRead.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

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

  float  var, evar, ivar, seed, run, pull, cov, junk;
  //Int_t ncols=0;
  Int_t nlines = 0;
  TFile* f = new TFile(rootfile.c_str(),"RECREATE");

  TTree* tree = new TTree(name_var.c_str(),"");
  tree->Branch("value",&var,"value/F");
  tree->Branch("error",&evar,"error/F");
  tree->Branch("truevalue",&trueval,"truevalue/F");
  tree->Branch("pull",&pull,"pull/F");
  tree->Branch("run",&run,"run/F");
  tree->Branch("cov",&cov,"cov/F");

  //while (inputfile) {
  while (inputfile  >> var >> evar >> cov >> seed >> run >> ivar >> junk){
    if (nlines < 3) std::cout << Form("var=%5f, error=%5f, cov=%5f, seed=%5f, run=%5f, init=%5f\n", var,evar,cov,seed, run, ivar);

    pull = (var-trueval)/evar;

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

  mytree->SetBranchAddress("value",&val);
  mytree->SetBranchAddress("error",&eval);
  mytree->SetBranchAddress("pull",&valp);
  mytree->SetBranchAddress("cov",&cov);

  mytree->GetEntry(0);
  TH1F * var = new TH1F("var",name_var.c_str(),30,val-10*eval,val+10*eval);
  TH1F * err = new TH1F("err",(name_var+"_err").c_str(),30,0,2*eval);
  TH1F * pull = new TH1F("pull",(name_var+"_pull").c_str(),40,-5,5);


  for(int i=0; i<mytree->GetEntries(); ++i){

    mytree->GetEntry(i);

    //   if(eval<0.02)continue;

    if(cov<3)continue;

    var->Fill(val);
    err->Fill(eval);
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


}
