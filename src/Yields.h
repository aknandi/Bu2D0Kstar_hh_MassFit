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
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bd_gen; 
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bd_fit; 
  // Bs
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bs_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bs_fit;
  // Combinatorics
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_comb;
  // Part Reco
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_bd_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_bs_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_bs_dstkst_010;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_bs_dstkst_001;

  //std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bd_dstkst;
  //std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bs_dstkst;
 
  // Drho
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_drho;

  // D3h
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_dkpipi;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_dpipipi;
  // Lambda
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > n_lambda;

  // GenPdf
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_comb_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bs_dstkst_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_bd_dstkst_gen;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > n_drho_gen;

  // More interesting parameters
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_drho;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_drho;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_dstkst_010;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_dstkst_010;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_dstkst_001;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_dstkst_001;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bd_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bd_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_lambda;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_lambda;


  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_dkpipi;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_dkpipi;
  std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > ratio_bs_dpipipi;
  std::map<std::string, std::map<std::string, std::map<std::string, RooAbsPdf*> > > gausratio_bs_dpipipi;


  std::map<std::string, std::map<std::string, std::map<std::string,  RooRealVar*> > >n_totsigs; //no pid & no bins
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooFormulaVar*> > > > stdBinFracs; //no pid part

  RooAbsReal *xplus;
  RooAbsReal *xminus;
  RooAbsReal *yplus;
  RooAbsReal *yminus;



  std::vector< RooRealVar* > ParametersOfFit;




 private:
  // Stuff that is internal to Yields.
  RooRealVar *xplus_hiding;
  RooRealVar *xminus_hiding;
  RooRealVar *yplus_hiding;
  RooRealVar *yminus_hiding;



  InternalStorage myMaps;

  std::vector<std::string> modeList;
  std::vector<std::string> pidList;
  std::vector<std::string> trackList;
  std::vector<std::string> chargeList;
  std::vector<std::string> binList;





};

#endif

