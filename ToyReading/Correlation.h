#ifndef CORRELATION_H 
#define CORRELATION_H 1

// Include files
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2F.h"
#include <vector>

#include <string>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

/** @class Correlation Correlation.h ToyReading/Correlation.h
 *  
 *
 *  @author Shu-Faye cheungs
 *  @date   2015-05-14
 */
class Correlation {
public: 
  /// Standard constructor
  Correlation( ); 
  void MakeHost(TString);
  void MakeFriends(TString);
  void Get2DHist(TString,TString);
  void DrawMatrix(TString);
  vector<TString> observables;
  virtual ~Correlation( ); ///< Destructor

  std::map<TString, std::map<TString, double> > results;

protected:

private:
  TFile* f_host;
  TTree* t_host;
};
#endif // TOYREADING_CORRELATION_H
