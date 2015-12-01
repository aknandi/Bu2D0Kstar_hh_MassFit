#ifndef Fitting_h
#define Fitting_h

#include <string>

#include "TH2.h"
#include "TFile.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsCollection.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "Settings.h"

#include "Base.h"
#include "Model.h"

using namespace std;

class TTree;
class TPaveText;
class TApplication;

class Fitting : public Base
{
 public :
  Fitting(TApplication*,Settings*);
  ~Fitting();
  void OrderToys(int);
  void NewOrderToys(int);
  int  LoadDataSet();
  void RunFullFit(bool);
  void RunManyFits();
  void RunManyToys();
  void DefineRooCategories();
  void DisplayToys();
  RooDataSet *FinalDataSet(const std::string, const std::string, TTree*);
  RooDataSet *FinalDataSet(const std::string, const std::string, RooDataSet* DS);
  void PrintDataSet(bool=false);
  void Fit();
  //void ReInitialize(RooAbsPdf*, RooAbsPdf*);	

 private:
  TFile *saveOutputForPlottingMacro;
  Model *model;
  RooDataSet* data;
  RooDataHist* dataBinned;
  RooWorkspace *w;
  RooArgSet inputlist;
  RooArgSet fulllist;
  RooArgSet reducedlist;
 
  RooRealVar mB;
  RooRealVar bid;
  RooRealVar bach_dll;
  RooRealVar assignedCharge;
//  RooRealVar dalitzm2_plus;
//  RooRealVar dalitzm2_minus;
//  RooRealVar DalitzBinNumber_equal;
//  RooRealVar DalitzBinNumber_optimal;
//  RooRealVar DalitzBinNumber_modopt;
//  RooRealVar DalitzBinNumber;
  RooRealVar BDTGresponse;

  RooSuperCategory* cat;	
  RooCategory* catNew;
  RooCategory mode;
  RooCategory run;
  RooCategory charge;
  RooCategory track;
  
  Settings* _genConfs;
  std::string batchMode;
  std::string drawProjections;
  std::string doFit;
  std::string readToys;
  std::string readData;
  std::string genToys;
  
  std::vector<std::string> modeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> trackList;
  std::vector<std::string> runList;
	
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,string> > > > title;
  std::map<std::string,std::map<std::string,int> > nVeto;
  
  TPaveText *lhcbpreliminary;
};


#endif
