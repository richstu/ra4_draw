#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>
#include <dirent.h>

#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/palette.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"
#include "hig/write_datacards_luminosity.hpp"
#include "hig/hig_utilities.hpp"

using namespace std;
namespace 
{
  string outFolder = getenv("PWD");
  string mass_points_string = "127_1";
  string years_string = "2016";
  //float luminosity = 137.;
  float luminosity = 35.9;
  string dimensionFilePath = "";
  bool use_data = true;
}

int main(int argc, char *argv[])
{
  
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  HigWriteDataCards::GetOptions(argc, argv);

  string baseFolder = HigUtilities::getBaseFolder();
  cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_sys_abcd!!"<<endl;
  gSystem->mkdir(outFolder.c_str(), kTRUE);
  // massPoints = { {"1000","1"} }
  vector<pair<string, string> > massPoints;
  HigUtilities::parseMassPoints(mass_points_string, massPoints);

  string baseline = "higd_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0";
  NamedFunc metFilters = HigUtilities::pass_2016;
  // TODO below is wrong. above is right
  //NamedFunc metFilters = Functions::hem_veto && Functions::pass_run2;
  //NamedFunc weight = "weight/w_btag*w_bhig_deep" * Higfuncs::eff_higtrig;
  //NamedFunc weight = "weight/w_btag*w_bhig_deep" * Higfuncs::eff_higtrig;
  NamedFunc weight = HigUtilities::weight_2016 * Higfuncs::eff_higtrig;

  set<int> years;
  HigUtilities::parseYears(years_string, years);

  map<string, string> samplePaths;
  samplePaths["mc_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higtight/";
  samplePaths["data_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  samplePaths["signal_2016"] = baseFolder + "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/";

  // sampleProcesses[mc, data, signal]
  map<string, vector<shared_ptr<Process> > > sampleProcesses;
  HigUtilities::setMcProcesses(years, samplePaths, metFilters&&"stitch_met", sampleProcesses);
  // TODO Below is wrong, above is correct
  //setMcProcesses(years, samplePaths, "pass && stitch_met", sampleProcesses);
  // Higfuncs::trig_hig will only work for 2016
  HigUtilities::setDataProcesses(years, samplePaths, metFilters&&Higfuncs::trig_hig>0., sampleProcesses);
  // TODO Below is wrong, above is correct
  //setDataProcesses(years, samplePaths, Higfuncs::trig_hig>0. && "pass", sampleProcesses);
  HigUtilities::setSignalProcesses(massPoints, years, samplePaths, metFilters, sampleProcesses);
  // TODO Below is wrong, above is correct
  //setSignalProcesses(massPoints, years, samplePaths, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet", sampleProcesses);

  // Set ABCDbins
  // "bkg" and "sig" should be kept for setDataCardBackground()
  map<string, string> xBins;
  xBins["bkg"] = "!(higd_am>100&&higd_am<140)&&higd_dm<40&&higd_am<200";
  xBins["sig"] = "higd_am>100&&higd_am<140&&higd_dm<40";
  map<string, string> yBins; // Shares sideband
  yBins["bkg"] = "nbdt==2&&nbdm==2";
  yBins["sig1"] = "nbdt>=2&&nbdm==3&&nbdl==3";
  yBins["sig2"] = "nbdt>=2&&nbdm>=3&&nbdl>=4";
  map<string, vector<pair<string, string> > > dimensionBins;
  if (dimensionFilePath=="")
  {
    dimensionBins["met"].push_back({"met0", "met>150&&met<=200"});
    dimensionBins["met"].push_back({"met1", "met>200&&met<=300"});
    dimensionBins["met"].push_back({"met2", "met>300&&met<=450"});
    dimensionBins["met"].push_back({"met3", "met>450"});
    //dimensionBins["drmax"].push_back({"drmaxL", "higd_drmax<=1.5"});
    //dimensionBins["drmax"].push_back({"drmaxH", "higd_drmax>1.5"});
  } 
  else 
  {
    HigWriteDataCards::readDimensionFile(dimensionFilePath, dimensionBins);
    for (auto & dimension : dimensionBins)
    {
      cout<<"dimension bins"<<endl;
      cout<<dimension.first<<endl;
      for (auto & entry : dimension.second)
      {
        cout<<" "<<entry.first<<" "<<entry.second<<endl;
      }
    }
  }


  // sampleBins = { {label, cut} }
  vector<pair<string, string> > sampleBins;
  HigUtilities::setABCDBins(xBins, yBins, dimensionBins, sampleBins);
  //sampleBins = {{"test","1"}};

  // Set systematics somehow..
  // Consider weight for nonHH somehow..

  // cuts[mc,data,signal] = RowInformation(labels, tableRows, yields)
  map<string, HigUtilities::RowInformation > cutTable;
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "data", cutTable["data"]);
  // TODO below is wrong above is correct
  //addBinCuts(sampleBins, baseline, Higfuncs::weight_higd * Higfuncs::eff_higtrig, "data", cutTable["data"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signal", cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "signalGenMet", HigUtilities::nom2genmet, cutTable["signal"]);
  HigUtilities::addBinCuts(sampleBins, baseline, weight, "mc", cutTable["mc"]);

  PlotMaker pm;
  // Luminosity used for labeling for table
  // Luminosity used for scaling for hist1d
  HigUtilities::makePlots(cutTable, sampleProcesses, luminosity, pm);

  // mYields[process_tag_sampleBinLabel] = GammaParams, TableRow
  map<string, pair<GammaParams, TableRow> > mYields;
  HigUtilities::fillDataYields(pm, cutTable["data"], mYields);
  // Luminosity used for scaling
  HigUtilities::fillMcYields(pm, luminosity, cutTable["mc"], mYields);
  // Luminosity used for scaling
  HigUtilities::fillSignalYieldsProcesses(pm, luminosity, sampleProcesses["signal"], cutTable["signal"], mYields);
  HigUtilities::fillAverageGenMetYields(sampleProcesses["signal"], sampleBins, "signal", "signalGenMet", "signalAverageGenMet", mYields, true);

  for (auto & process : sampleProcesses["signal"])
  {
    cout<<process->name_<<endl;
    string model;
    int mGluino=0, mLSP=0;
    HigUtilities::getInfoFromProcessName(process->name_, model, mGluino, mLSP);
    TString outPath = outFolder+"/datacard-"+HigUtilities::setProcessNameLong(model, mGluino, mLSP)+"_"+years_string+".txt";
    cout<<"open "<<outPath<<endl;
    ofstream cardFile(outPath);
    HigWriteDataCards::writeDataCardHeader(sampleBins,cardFile);

    vector<vector<string> > tableValues;
    if (use_data) HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "data", tableValues);
    else HigWriteDataCards::setDataCardObserved(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::setDataCardSignalBackground(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
    HigWriteDataCards::setDataCardSignalStatistics(process->name_, "signalAverageGenMet", mYields, sampleBins, tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile);
    tableValues.clear();
    HigWriteDataCards::setDataCardBackground(mYields, sampleBins, "mc", tableValues);
    HigWriteDataCards::writeTableValues(tableValues,cardFile, true);
    tableValues.clear();
    //writeDataCardClosure()
    //writeDataCardUncertainties()
    //writeDataCardParam()
  }

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

namespace HigWriteDataCards{
  void readDimensionFile(std::string const & dimensionFilePath, std::map<std::string, std::vector<std::pair<std::string, std::string> > > & dimensionBins)
  {
    map<string, int> dimensionIndex;
    ifstream dimensionFile(dimensionFilePath.c_str());
    if (dimensionFile.is_open())
    {
      string dimension, dimensionCut;
      while (dimensionFile >> dimension >> dimensionCut)
      {
        dimensionBins[dimension].push_back(std::make_pair(dimension+to_string(dimensionIndex[dimension]), dimensionCut));
        dimensionIndex[dimension]++;
      }
      dimensionFile.close();
    }
  }

  void writeDataCardHeader(vector<pair<string, string> > sampleBins, ofstream & cardFile)
  {
    cardFile<<"imax "<<sampleBins.size()<<"  number of channels\n";
    cardFile<<"jmax 1  number of backgrounds\n";
    cardFile<<"kmax *  number of nuisance parameters\n";
    cardFile<<"shapes * * FAKE\n\n";
  }

  void setDataCardObserved(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & dataTag, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "Observation";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      row[iBin*2+3] = to_string(mYields.at(dataTag+"_"+sampleBins[iBin].first).first.Yield());
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }

  void setDataCardSignalBackground(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    row[0] = "bin";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      row[iBin*2+2] = sampleBins[iBin].first;
      row[iBin*2+3] = sampleBins[iBin].first;
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      row[iBin*2+2] = "sig";
      row[iBin*2+3] = "bkg";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "process";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      row[iBin*2+2] = "0";
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    row[0] = "rate";
    for (unsigned iBin=0;iBin<sampleBins.size();++iBin)
    {
      string label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[iBin].first;
      row[iBin*2+2] = to_string(mYields.at(label).first.Yield());
      row[iBin*2+3] = "1";
    }
    tableValues.push_back(row);
    setRow(row,"");
  
    tableValues.push_back(row);
  }
  
  void setDataCardSignalStatistics(string const & processName, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, vector<vector<string> > & tableValues)
  {
    // title + type + nBins * 2
    vector<string> row(2+2*sampleBins.size());
    string label;
  
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      setRow(row,"-");
      row[0] = "stat_"+sampleBins[rBin].first;
      row[1] = "lnN";
      label = processName + "_"+signalAverageGenMetTag+"_" +sampleBins[rBin].first;
      string uncertainty;
      if (mYields.at(label).first.Yield() == 0) uncertainty = "2.00";
      else uncertainty = to_string(1+mYields.at(label).first.Uncertainty()/mYields.at(label).first.Yield());
      row[2+2*rBin] = uncertainty;
      tableValues.push_back(row);
    }
  
    setRow(row,"");
    tableValues.push_back(row);
  }
  
  // returns count of non-overlapping occurrences of 'sub' in 'str'
  int countSubstring(const std::string& str, const std::string& sub)
  {
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos; offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
  }

  void setDataCardBackground(map<string, pair<GammaParams, TableRow> > & mYields, vector<pair<string, string> > sampleBins, string const & mcTag, vector<vector<string> > & tableValues)
  {
    // title + type + tag + tagType + (yield,equation), equation arguments
    vector<string> row(6);
    for (unsigned rBin=0;rBin<sampleBins.size();++rBin)
    {
      row[0] = "rp_"+sampleBins[rBin].first;
      row[1] = "rateParam";
      row[2] = sampleBins[rBin].first;
      row[3] = "bkg";
      if (countSubstring(sampleBins[rBin].first,"sig") != 2) 
      {
        row[4] = to_string(mYields.at(mcTag+"_"+sampleBins[rBin].first).first.Yield());
        // Infinity defined in https://root.cern.ch/doc/master/RooNumber_8cxx_source.html
        row[5] = "[0,1.0e30]";
      }
      else 
      {
        row[4] = "(@0*@1/@2)";
        vector<string> xydimension;
        HigUtilities::stringToVectorString(sampleBins[rBin].first, xydimension, "_");
        vector<string> aNameVector = xydimension;
        aNameVector[0] = "xbkg";
        aNameVector[1] = "ybkg";
        string aName;
        HigUtilities::vectorStringToString(aNameVector, aName, "_");
        vector<string> bNameVector = xydimension;
        bNameVector[1] = "ybkg";
        string bName;
        HigUtilities::vectorStringToString(bNameVector, bName, "_");
        vector<string> cNameVector = xydimension;
        cNameVector[0] = "xbkg";
        string cName;
        HigUtilities::vectorStringToString(cNameVector, cName, "_");
        //cout<<aName<<" "<<bName<<" "<<cName<<endl;
        row[5] = "rp_"+cName+",rp_"+bName+",rp_"+aName;
      }
      tableValues.push_back(row);
      setRow(row,"");
    }
  
  }
  
  void writeTableValues(vector<vector<string> > & tableValues, ofstream & cardFile, bool alignLeft)
  {
    // Set space values
    vector<unsigned> xSpace(tableValues.back().size());
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (xSpace[xIndex] < tableValues[yIndex][xIndex].size()) xSpace[xIndex] = tableValues[yIndex][xIndex].size();
      }
    }
  
    // Write values
    for (unsigned yIndex = 0; yIndex < tableValues.size(); ++yIndex)
    {
      for (unsigned xIndex = 0; xIndex < tableValues[yIndex].size(); ++xIndex)
      {
        if (alignLeft) cardFile << setw(xSpace[xIndex]) << std::left << tableValues[yIndex][xIndex] << " ";
        else cardFile << setw(xSpace[xIndex]) << tableValues[yIndex][xIndex] << " ";
      }
      cardFile<<endl;
    }
  }
  
  void setRow(vector<string> & row, string const & value)
  {
    for(auto & item : row) item = value;
  }

  void GetOptions(int argc, char *argv[])
  {
    string blah;
    while(true){
      static struct option long_options[] = 
      {
        {"output_folder", required_argument, 0, 'o'},
        {"mass_points", required_argument, 0, 'p'},
        {"years", required_argument, 0, 'y'},
        {"luminosity", required_argument, 0, 'l'},
        {"dimension", required_argument, 0, 'd'},
        {"no_data", no_argument, 0, 'n'},
        {0, 0, 0, 0}
      };
  
      char opt = -1;
      int option_index;
      opt = getopt_long(argc, argv, "o:p:y:l:d:n", long_options, &option_index);
      if( opt == -1) break;
  
      string optname;
      switch(opt)
      {
        case 'o': outFolder = optarg; break;
        case 'p': mass_points_string = optarg; break;
        case 'y': years_string = optarg; break;
        case 'l': luminosity = atof(optarg); break;
        case 'd': 
          dimensionFilePath = optarg; 
          if (!FileExists(dimensionFilePath)) 
          {
            cout<<"[Error] No file caled "<<dimensionFilePath<<endl;
            exit(EXIT_FAILURE);
          }
          break;
        case 'n':
          use_data = false;
          break;
        case 0:
          optname = long_options[option_index].name;
          printf("Bad option! Found option name %s\n", optname.c_str());
          break;
        default: 
          printf("Bad option! getopt_long returned character code 0%o\n", opt); 
          break;
      }
    }
  }

}
