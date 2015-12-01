#ifndef InternalStorage_h
#define InternalStorage_h

#include <string>
#include "TRandom3.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooSuperCategory.h"
#include "Base.h"
#include <vector>
#include <map>
#include "Settings.h"


class InternalStorage
{
 public:

  void Initialize(Settings* Input , std::vector<std::string> allmodeList ,  std::vector<std::string> allChargeList, std::vector<std::string> alltrackList, std::vector<std::string> allrunList);

  InternalStorage();

  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > >Eff;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double>  > > >PidGenEff; //same over bins
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, RooRealVar*> > > > PidFitEff; //same over bins

  std::map<std::string, double>  ek;
  std::map<std::string, double>  epi;

  std::map<std::string, double>  ek_KsKK;
  std::map<std::string, double>  epi_KsKK;

  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > TotSig;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > TotLowmass;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > TotComb;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > FracBsInLowMass;

  void SetTotals();
  //  void SetBackgroundModel();
  //  void SetCiSiKi();
  //  void SetCiSiKiFit();
  //  void SetEfficiency();
  //  void SetNormKiModEff();
  //void SetPIDEffs();
  //void SetPIDEffsFit();

 private:

  Settings* _Control;
  std::vector<std::string> modeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> trackList;
  std::vector<std::string> runList;

};
#endif

