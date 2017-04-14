#ifndef BASE_H 
#define BASE_H

// Include files
#include<string>
#include<map>
#include<vector>

class Base
{
 public: 
  /// Standard constructor
  Base(); 
  virtual ~Base( ){}
 protected:
  std::string d2kpi;
  std::string d2kk;
  std::string d2pipi;
  std::string d2pik;
  std::string d2kpipipi;
  std::string d2pipipipi;
  std::string d2pikpipi;
  std::string minus;
  std::string plus;
  std::string both;
  std::string LL;
  std::string DD;
  std::string mix;
  std::string run1;
  std::string run2;
  std::string all;
  std::string merge;
  std::string separate;
  std::string slash;
  std::string underscore;
  std::string null;
	
  double pi;
  double twopi;
  double pibytwo;

  std::vector<std::string> allmodeList;
};

#endif // BASE_H
