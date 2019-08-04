#ifndef H_HIG_WRITE_DATACARDS_LUMINOSITY
#define H_HIG_WRITE_DATACARDS_LUMINOSITY

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <map>

#include "core/gamma_params.hpp"
#include "core/table_row.hpp"

namespace HigWriteDataCards {
  void GetOptions(int argc, char *argv[]);
  void writeDataCardHeader(std::vector<std::pair<std::string, std::string> > sampleBins, std::ofstream & cardFile);
  void setDataCardObserved(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::string const & dataTag, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardSignalBackground(std::string const & processName, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues);
  void setDataCardSignalStatistics(std::string const & processName, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::vector<std::vector<std::string> > & tableValues);
  // Counts number of "sig" to find tags for equations. Uses "bkg" tag. Relies on "xname_yname_dimensionname" tag 
  int countSubstring(const std::string& str, const std::string& sub);
  void setDataCardBackground(std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, std::vector<std::pair<std::string, std::string> > sampleBins, std::string const & mcTag, std::vector<std::vector<std::string> > & tableValues);
  void writeTableValues(std::vector<std::vector<std::string> > & tableValues, std::ofstream & cardFile);
  void setRow(std::vector<std::string> & row, std::string const & value);
}

#endif
