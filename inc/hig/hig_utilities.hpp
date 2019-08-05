#ifndef H_HIG_UTILITIES
#define H_HIG_UTILITIES

#include <iostream>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <memory>

#include "TString.h"

#include "core/named_func.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/gamma_params.hpp"

namespace HigUtilities {
  class RowInformation {
    public:
    std::vector<std::string> labels;
    // To get cut: tableRows[0].cut_.Name()
    std::vector<TableRow> tableRows;
    // Index used by PlotMaker. Set in makePlots()
    int tableIndex;
  };

  int stringToVectorString(std::string const& inString, std::vector<std::string>& outputVector, std::string const & delimiter);
  int vectorStringToString(std::vector<std::string> const & inVector, std::string &outString, std::string const & delimiter);
  std::string removeSpaces(std::string inString);

  extern const NamedFunc pass_2016;
  extern const NamedFunc weight_2016;

  TString nom2sys_bin(TString ibin, size_t shift_index);
  TString nom2genmet(TString ibin);

  std::string getBaseFolder();
  // Example: "400_1, 700_1, 1000_1" => {{400,1}, {700,1}, {1000,1}}
  void parseMassPoints(std::string mass_points_string, std::vector<std::pair<std::string, std::string> > & mass_points);
  // Example: "2016,2017" => {2016,2017}
  void parseYears(std::string years_string, std::set<int> & years);

  std::string setProcessName(std::string const & model, int const & mGluino, int const & mLSP);
  std::string setProcessNameLong(std::string const & model, int const & mGluino, int const & mLSP);
  void getInfoFromProcessName(std::string const & processName, std::string & model, int & mGluino, int & mLSP);

  // sampleBins = { {label, cut} }
  void setABCDBins(std::map<std::string, std::string> xBins, std::map<std::string, std::string> yBins, std::map<std::string, std::vector<std::pair<std::string, std::string> > > dimensionBins, std::vector<std::pair<std::string, std::string> > & sampleBins);
  void combineDimensionBins(std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins);

  void setMcProcesses(std::set<int> years, std::map<std::string, std::string> samplePaths, NamedFunc filters, std::map<std::string, std::vector<std::shared_ptr<Process> > > & sampleProcesses);
  void setDataProcesses(std::set<int> years, std::map<std::string, std::string> samplePaths, NamedFunc filters, std::map<std::string, std::vector<std::shared_ptr<Process> > > & sampleProcesses);
  void setSignalProcesses(std::vector<std::pair<std::string, std::string> > massPoints, std::set<int> years, std::map<std::string, std::string> samplePaths, NamedFunc filters, std::map<std::string, std::vector<std::shared_ptr<Process> > > & sampleProcesses);

  void addBinCuts(std::vector<std::pair<std::string, std::string> > sampleBins, std::string baseline, NamedFunc weight, std::string tag, RowInformation & cutRows);
  void addBinCuts(std::vector<std::pair<std::string, std::string> > sampleBins, std::string baseline, NamedFunc weight, std::string tag, TString (*replaceFunc)(TString) ,RowInformation & cutRows);
  void addBinCuts(std::vector<std::pair<std::string, std::string> > sampleBins, std::string baseline, NamedFunc weight, std::string tag, TString (*replaceFunc)(TString, size_t), RowInformation & cutRows);
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  void makePlots(std::map<std::string, RowInformation > & cutTable, std::map<std::string, std::vector<std::shared_ptr<Process> > > & sampleProcesses, float luminosity, PlotMaker & pm);

  void addToMapYields(std::string const & label, GammaParams & yield, TableRow & yieldMeta, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields);
  void fillDataYields(PlotMaker & pm, RowInformation & dataRow, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, bool verbose=false);
  void fillMcYields(PlotMaker & pm, float luminosity, RowInformation & mcRow, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, bool verbose=false);
  void fillSignalYieldsProcesses(PlotMaker & pm, float luminosity, std::vector<std::shared_ptr<Process> > & signalProcesses, RowInformation & signalRow, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, bool verbose=false);
  void fillAverageGenMetYields(std::vector<std::shared_ptr<Process> > & sampleProcesses, std::vector<std::pair<std::string, std::string> > sampleBins, std::string const & signalTag, std::string const & signalGenMetTag, std::string const & signalAverageGenMetTag, std::map<std::string, std::pair<GammaParams, TableRow> > & mYields, bool verbose=false);

}
#endif
