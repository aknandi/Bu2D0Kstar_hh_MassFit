// Include files 
#include "TAxis.h"
#include <cmath>
// local
#include "Correlation.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Correlation
//
// 2015-05-14 : Shu-Faye cheungs
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
Correlation::Correlation( ) {

}
//=============================================================================
// Destructor
//=============================================================================
Correlation::~Correlation() {} 

//=============================================================================

void Correlation::MakeHost(TString var){
  f_host = new TFile("results/"+var+"_ntp.root");
  t_host = (TTree*) f_host->Get(var);
}

void Correlation::MakeFriends(TString var){
  t_host->AddFriend(var,"results/"+var+"_ntp.root");
}

void Correlation::Get2DHist(TString var1, TString var2){
  TCanvas* c = new TCanvas("c","",700,700);
  t_host->Draw(var1+".value:"+var2+".value","","colz");
  TH2F* h2 = (TH2F*) t_host->GetHistogram(); // returns the last drawn histogram
  double correlation = h2->GetCorrelationFactor();
  results[var1][var2] = correlation;
  c->SaveAs("correlations/"+var1+"_"+var2+".pdf");
  delete h2;
  delete c;
}

void Correlation::DrawMatrix(){

  map<TString, map<TString,double> >::iterator it1 = results.begin();
  map<TString, double>::iterator it2 = ((*it1).second).begin();

  int nvars = (int) results.size();
  TH2F* hmatrix = new TH2F("matrix","",nvars,0,nvars,nvars,0,nvars);
  hmatrix->SetStats(0);

  int counter1= 1;
  for(it1 = results.begin(); it1!= results.end(); ++it1){

    TString var1 = (*it1).first;
    var1.ReplaceAll("_both_mix_merge","");
    hmatrix->GetXaxis()->SetBinLabel(counter1,var1);
    counter1++;

    int counter2=1;
    map<TString, double> results2 = (*it1).second;
    map<TString, double>::iterator it2 = results2.begin();

    for(it2 = results2.begin(); it2!=results2.end(); ++it2){
      TString var2 = (*it2).first;
      var2.ReplaceAll("_both_mix_merge","");
      //cout << counter2 << " " << var2 << endl;
      if ( counter1==1) hmatrix->GetYaxis()->SetBinLabel(counter2,var2);
      counter2++;
      double value = abs((*it2).second);
      //cout << var1 << " " << var2 << " " << value << endl;
      hmatrix->Fill(var1,var2,value);
      if(var1!=var2) hmatrix->Fill(var2,var1,value);
    }
  }

  TCanvas* c = new TCanvas("c","",700,700);
  hmatrix->GetXaxis()->SetLabelSize(0.03);
  hmatrix->GetYaxis()->SetLabelSize(0.03);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.2);
  hmatrix->Draw("colz");
  c->SaveAs("results/matrix.pdf");
  c->SaveAs("results/matrix.eps");

}
