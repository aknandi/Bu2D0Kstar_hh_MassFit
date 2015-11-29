#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"

#include "SimpleToyRead.h"
#include "Correlation.h"
//#include "boost/algorithm/string.hpp"
using namespace std;
int main(int argc, char *argv[])
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1011);

	SimpleToyRead * str = new SimpleToyRead();
  Correlation* corr = new Correlation();

  // Get true values
  ifstream infile;
  infile.open("truevals_4900_modAlexis.txt",ifstream::in);
  string var;
  float val;
  vector<TString> varList;
  while (infile >> var >> val) {
    varList.push_back(var.c_str());
    cout << var << " " << val << endl; 
    str->MakeNtupleFromTextFile(var,val);
    cout << "hello back to main " << endl;
    str->MakeSomePlotsFromRootFile(var,val);
  }

  vector<TString>::iterator it = varList.begin();
  cout << "Making host: " << *it << endl;
  corr->MakeHost(*it);
  
  for(it = varList.begin()+1; it!=varList.end(); ++it){
    cout << "Adding friend: " << *it << endl;
    corr->MakeFriends(*it);
  }

  for(int i = 0; i<(int) varList.size(); ++i){
    for(int j = 0; j<(int) varList.size(); ++j){
      if(i<=j) corr->Get2DHist(varList.at(i),varList.at(j));
    }
  }

  corr->DrawMatrix();
}
