// Include files
#include <dirent.h>
#include <errno.h>
#include <iostream>

#include "TColor.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CommonTools.h"

bool CommonTools::exists(std::string file)
{
  FileStat_t x;
  return (not gSystem->GetPathInfo(file.c_str(),x));
}

bool CommonTools::isADirectory(std::string dir)
{
  return (bool)opendir(dir.c_str());
}

void CommonTools::defineColorMap()
{
  int nColors=51;
  Double_t s[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
  Double_t r[5] = {0.99, 0.00, 0.87, 1.00, 0.70};
  Double_t g[5] = {0.99, 0.81, 1.00, 0.20, 0.00};
  Double_t b[5] = {0.99, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(5, s, r, g, b, nColors);
  gStyle->SetNumberContours(nColors);
}

int CommonTools::getdir (std::string dir, std::vector<std::string> &files)
{
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }
  while ((dirp = readdir(dp)) != NULL) {
    std::string file(dirp->d_name);
    files.push_back(file);
  }
  closedir(dp);
  std::cout <<"\n" << std::endl;
  return 0;
}

void CommonTools::split(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = " ")
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
	
  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

