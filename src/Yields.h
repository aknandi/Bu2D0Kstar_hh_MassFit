#ifndef Yields_h
#define Yields_h

#include <string>
#include "TRandom3.h"
#include <vector>
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooSuperCategory.h"
#include "InternalStorage.h"
#include "Base.h"
//#include "CommonTools.h"
#include "Settings.h"

class Yields
{

 public :
  Yields(Settings* genConfs, Settings* fileList, std::vector<std::string> allmodeList, std::vector<std::string>allchargeList, std::vector<std::string>alltrackList, std::vector<std::string>allbinList, std::string unblind); //call in the constructor the other functions required
 
  Settings* _genConfs;
  Settings* _fileList;
  Settings* input;

  void SetYieldsGenandFit();
  void SetOtherBkgs();
  void SetDstKstGenandFit();
  void SetDrhoGenandFit();
  
  double genscale;
  std::string limitlow;

  // Yields can be accessed by Model.
  // Bs 
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bu_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bu_fit;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_comb_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_comb;


  //std::map<std::string, std::map<std::string, std::map<std::string,  RooRealVar*> > >n_totsigs; //no pid & no bins
  //std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > stdBinFracs; //no pid part

//  RooAbsReal *xplus;
//  RooAbsReal *xminus;
//  RooAbsReal *yplus;
//  RooAbsReal *yminus;
//


  std::vector< RooRealVar* > ParametersOfFit;




 private:
  // Stuff that is internal to Yields.
//  RooRealVar *xplus_hiding;
//  RooRealVar *xminus_hiding;
//  RooRealVar *yplus_hiding;
//  RooRealVar *yminus_hiding;
//


//  InternalStorage myMaps;

  std::vector<std::string> modeList;
  std::vector<std::string> pidList;
  std::vector<std::string> trackList;
  std::vector<std::string> chargeList;
  std::vector<std::string> binList;





};

#endif

