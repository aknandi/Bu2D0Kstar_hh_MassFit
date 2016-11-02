#ifndef SimpleToyRead_h
#define SimpleToyRead_h

#include <string>
#include "TRandom3.h"
#include "TFile.h"
#include <iostream>
#include "TString.h"
#include <vector>

class SimpleToyRead {

 public :
  SimpleToyRead();
 

  void MakeNtupleFromTextFile(std::string name_var, float trueval);
  double MakeSomePlotsFromRootFile(std::string name_var, float trueval, std::string method);


};

#endif

