#include <algorithm>
#include <stdlib.h>
#include "hig/hig_utilities.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"

namespace HigUtilities {
  using std::string;
  using std::to_string;
  using std::map;
  using std::vector;
  using std::set;
  using std::pair;
  using std::cout;
  using std::endl;
  using std::shared_ptr;
  using std::max;

  int stringToVectorString(std::string const& inString, std::vector<std::string>& outputVector, std::string const & delimiter)
  {
      if( outputVector.size() != 0) {
      std::cout<<"[Error] StringFunctions::divideString() => outputVector size is not 0."<<std::endl;
      return 1;
    }
  
    if(delimiter.size()==0){
      std::cout<<"[Error] StringFunctions::divideString() => delimiter size is 0."<<std::endl;
      return 1;
    }
  
    size_t start = 0;
    size_t end = 0;
    while((end=inString.find(delimiter,start)) != string::npos) {
      if(start!=end) {
        //std::cout<<"found: "<<removeSpaces(inString.substr(start, end-start))<<std::endl;
        outputVector.push_back(removeSpaces(inString.substr(start, end-start)));
      }
      start = end + delimiter.length();
    }
    if(start!=inString.length()) {
      //std::cout<<"found: "<<removeSpaces(inString.substr(start))<<std::endl;
      outputVector.push_back(removeSpaces(inString.substr(start)));
    }
    return 0;
  }
  
  int vectorStringToString(std::vector<std::string> const & inVector, std::string &outString, std::string const & delimiter){
    if (inVector.size()==0) {
      outString = "";
    } else {
      outString = "";
      for(unsigned iVector=0; iVector<inVector.size()-1; iVector++){
        outString += inVector[iVector] + delimiter;
      }
      outString += inVector[inVector.size()-1];
    }
    return 0;
  }
  
  std::string removeSpaces(std::string inString){
    inString.erase(std::remove(inString.begin(),inString.end(),' '),inString.end());
    return inString;
  }

  // based on pass_run2 with removal of non-existing branches
  const NamedFunc pass_2016("pass_2016", [](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_ =  b.pass_ra2_badmu() && (b.met()/b.met_calo()<5);
    if (b.type()<1000 && b.type()>0) { // Data
      pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
              b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu() && b.pass_eebadsc();
    } else {
      if (b.type()>=100e3) { // FastSim
        pass_ = pass_ && b.pass_fsjets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
      } else { //FullSim
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
      }
    }
    if (pass_) return 1.;
    else return 0.;
  });
  
  const NamedFunc weight_2016("weight_2016", [](const Baby &b) -> NamedFunc::ScalarType{
    // Data
    if (b.type()<1000 && b.type()>0) return 1.;
    else return b.weight()/b.w_btag()*b.w_bhig_deep();
  });

  TString nom2genmet(TString ibin){
    ibin.ReplaceAll("met", "met_tru");
    ibin.ReplaceAll("met_tru/met_tru_calo", "met/met_calo");
    return ibin;
  }
  
  TString nom2sys_bin(TString ibin, size_t shift_index){
    ibin.ReplaceAll("met", "sys_met["+to_string(shift_index)+"]");
    ibin.ReplaceAll("mt", "sys_mt["+to_string(shift_index)+"]");
    ibin.ReplaceAll("st", "sys_st["+to_string(shift_index)+"]");
    ibin.ReplaceAll("mj14", "sys_mj14["+to_string(shift_index)+"]");
    ibin.ReplaceAll("njets", "sys_njets["+to_string(shift_index)+"]");
    ibin.ReplaceAll("nbdm", "sys_nbdm["+to_string(shift_index)+"]");
    return ibin;
  }

  string getBaseFolder()
  {
    string baseFolder = "";
    string hostName = execute("echo $HOSTNAME");
    if(Contains(hostName, "cms") || Contains(hostName, "compute-")) baseFolder = "/net/cms2";
    return baseFolder;
  }

  void parseMassPoints(string mass_points_string, vector<pair<string, string> > & mass_points)
  {
    if (mass_points_string!="") 
    {
      size_t found = mass_points_string.find(",");
      size_t start = 0;
      string _tmp = "";
      while (found != string::npos) {
        _tmp = mass_points_string.substr(start, found-start);
        mass_points.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
        cout<<"Adding mass point: mgluino = "<<mass_points.back().first<<" mlsp = "<<mass_points.back().second<<endl;
        start = found+1;
        found = mass_points_string.find(",",start);
      } 
      _tmp = mass_points_string.substr(start, found-start);
      mass_points.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
      cout<<"Adding mass point: mgluino = "<<mass_points.back().first<<" mlsp = "<<mass_points.back().second<<endl;
    } 
  }
  
  void parseYears(string years_string, set<int> & years)
  {
    if (years_string!="")
    {
      size_t found = years_string.find(",");
      size_t start = 0;
      string _tmp = "";
      while (found != string::npos) {
        _tmp = years_string.substr(start, found-start);
        years.insert(atoi(_tmp.c_str()));
        start = found+1;
        found = years_string.find(",",start);
      } 
      _tmp = years_string.substr(start, found-start);
      years.insert(atoi(_tmp.c_str()));
    }
  }

  std::string setProcessName(std::string const & model, int const & mGluino, int const & mLSP)
  {
    return model+"_"+to_string(mGluino)+"_"+to_string(mLSP);
  }

  std::string setProcessNameLong(std::string const & model, int const & mGluino, int const & mLSP)
  {
    // Tune is added for utilites::parseMasses
    return model+"_mGluino-"+to_string(mGluino)+"_mLSP-"+to_string(mLSP)+"_Tune";
  }

  void getInfoFromProcessName(std::string const & processName, std::string & model, int & mGluino, int & mLSP)
  {
    vector<string> names;
    stringToVectorString(processName, names, "_");
    model = names[0];
    mGluino = atoi(names[1].c_str());
    mLSP = atoi(names[2].c_str());
  }

  void setABCDBins(map<string, string> xBins, map<string, string> yBins, map<string, vector<pair<string, string> > > dimensionBins, vector<pair<string, string> > & sampleBins)
  {
    // Combine dimensions to one list of dimensions.
    combineDimensionBins(dimensionBins);
  
    auto dimension = dimensionBins.begin();
    for (unsigned iBin = 0; iBin < dimension->second.size(); ++iBin)
    {
      for (auto yBin: yBins)
      {
        for (auto xBin : xBins)
        {
          string label = "x"+xBin.first+"_y"+yBin.first+"_"+(dimension->second)[iBin].first;
          string cut = xBin.second+"&&"+yBin.second+"&&"+(dimension->second)[iBin].second;
          //cout<<"bins: "<<label<<" "<<cut<<endl;
          sampleBins.push_back({label, cut});
        }
      }
    }
  }
  
  void combineDimensionBins(map<string, vector<pair<string, string> > > & dimensionBins)
  {
    unsigned nDimensions = dimensionBins.size();
    if (nDimensions == 1) return;
    auto dimension1 = dimensionBins.begin();
    auto dimension2 = dimension1++;
    // Combine dimensions
    dimensionBins[dimension1->first+"_"+dimension2->first];
    for (unsigned iDimension1 = 0; iDimension1 < dimension1->second.size(); ++iDimension1)
    {
      for (unsigned iDimension2 = 0; iDimension2 < dimension2->second.size(); ++iDimension2)
      {
        dimensionBins[dimension1->first+"_"+dimension2->first].push_back({dimension1->second[iDimension1].first+"_"+dimension2->second[iDimension2].first,dimension1->second[iDimension1].second+"&&"+dimension2->second[iDimension2].second});
        //cout<<dimension1->first+"_"+dimension2->first<<":"<<dimension1->second[iDimension1].first+"_"+dimension2->second[iDimension2].first<<" "<<dimension1->second[iDimension1].second+"&&"+dimension2->second[iDimension2].second<<endl;
      }
    }
    dimensionBins.erase(dimension1);
    dimensionBins.erase(dimension2);
    return combineDimensionBins(dimensionBins);
  }

  void setMcProcesses(set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // mcTag['ttx'] = {*.root,*.root}
    map<string, set<string> > mcTags;
  
    mcTags["t#bar{t]+X"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                       "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
    mcTags["V+jets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
    mcTags["QCD"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
    mcTags["Other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                       "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
    map<string, int> mcColors;
    Palette colors("txt/colors.txt", "default");
    mcColors["t#bar{t]+X"] = colors("tt_1l");
    mcColors["V+jets"] = kOrange+1;
    mcColors["QCD"] = colors("other");
    mcColors["Other"] = kGreen+1;
  
    // mcPaths['ttx'] = {folder/*.root,folder/*.root}
    map<string, set<string> > mcPaths;
    for (auto mcTag : mcTags)
    {
      for (auto file : mcTag.second)
      {
        for (auto year : years)
        {
          mcPaths[mcTag.first].insert(samplePaths["mc_"+to_string(year)]+file);
          cout<<"Adding "<<samplePaths["mc_"+to_string(year)]+file<<endl;
        }
      }
    }
  
    sampleProcesses["mc"] = vector<shared_ptr<Process> > ();
    for (auto mcPath : mcPaths)
    {
      sampleProcesses["mc"].push_back(Process::MakeShared<Baby_full>(mcPath.first, Process::Type::background, mcColors[mcPath.first], mcPath.second, filters));
    }
  }
  
  void setDataProcesses(set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // dataPaths['data'] = {folder/*.root,folder/*.root}
    map<string, set<string> > dataPaths;
    for (auto year : years)
    {
      dataPaths["Data"].insert(samplePaths["data_"+to_string(year)]+"*.root");
      cout<<"Adding "<<samplePaths["data_"+to_string(year)]+"*.root"<<endl;
    }
  
    sampleProcesses["data"] = vector<shared_ptr<Process> > ();
    for (auto dataPath : dataPaths)
    {
      sampleProcesses["data"].push_back(Process::MakeShared<Baby_full>(dataPath.first, Process::Type::data, 1, dataPath.second, filters));
    }
  }
  
  void setSignalProcesses(vector<pair<string, string> > massPoints, set<int> years, map<string, string> samplePaths, NamedFunc filters, map<string, vector<shared_ptr<Process> > > & sampleProcesses)
  {
    // mcPaths['ttx'] = {folder/*.root,folder/*.root}
    map<string, set<string> > signalPaths;
    for (auto massPoint : massPoints)
    {
     for (auto year : years)
     {
       signalPaths[HigUtilities::setProcessName("TChiHH",atoi(massPoint.first.c_str()),atoi(massPoint.second.c_str()))].insert(samplePaths["signal_"+to_string(year)]+"*mGluino-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root");
       cout<<"Adding "<<samplePaths["signal_"+to_string(year)]+"*mGluino-"+massPoint.first+"_mLSP-"+massPoint.second+"_*.root"<<endl;
     }
    }
  
    sampleProcesses["signal"] = vector<shared_ptr<Process> > ();
    for (auto signalPath : signalPaths)
    {
      sampleProcesses["signal"].push_back(Process::MakeShared<Baby_full>(signalPath.first, Process::Type::background, kBlack, signalPath.second, filters));
    }
  }

  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", baseline&&binCut.second,0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, TString (*replaceFunc)(TString), RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second),0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void addBinCuts(vector<pair<string, string> > sampleBins, string baseline, NamedFunc weight, string tag, TString (*replaceFunc)(TString, size_t) ,RowInformation & cutRows)
  {
    for (auto binCut : sampleBins)
    {
      cutRows.labels.push_back(tag +"_up_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second,1),0,0,weight));
      cutRows.labels.push_back(tag +"_down_" + binCut.first);
      cutRows.tableRows.push_back(TableRow("", replaceFunc(baseline+"&&"+binCut.second,2),0,0,weight));
      //cout<<cutRows.labels.back()<<" "<<cutRows.tableRows.back().cut_.Name()<<endl;
    }
  }
  
  void makePlots(map<string, RowInformation > & cutTable, map<string, vector<shared_ptr<Process> > > & sampleProcesses, float luminosity, PlotMaker & pm)
  {
    int tableIndex = 0;
    for (auto cutRow : cutTable)
    {
      // type = mc, data, signal
      string const & type = cutRow.first;
      pm.Push<Table>(type, cutTable[type].tableRows, sampleProcesses[type], true, false);
      cutTable[type].tableIndex= tableIndex;
      tableIndex++;
    }
    pm.multithreaded_ = true;
    pm.min_print_ = true; 
    pm.MakePlots(luminosity);  
  }

  void addToMapYields(string const & label, GammaParams & yield, TableRow & yieldMeta, map<string, pair<GammaParams, TableRow> > & mYields)
  {
      if (mYields.find(label) != mYields.end())
      {
        cout<<"[Warning] addYieldData(): "<<label<<" already exists in mYields. Will not add data."<<endl;
      }
      pair<GammaParams, TableRow> yieldData = {yield, yieldMeta};
      mYields.insert(pair<string, pair<GammaParams, TableRow> > (label, yieldData));
  }
  
  void fillDataYields(PlotMaker & pm, RowInformation & dataRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[dataRow.tableIndex].get());
    vector<GammaParams> yields = yieldTable->DataYield();
    for (unsigned ipar(0); ipar< yields.size(); ipar++) 
    {
      string const & label = dataRow.labels[ipar];
      addToMapYields(label, yields[ipar], dataRow.tableRows[ipar], mYields);
      if (verbose) cout<<label<<": "<<mYields.at(label).first.Yield()<<endl;
      //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
    }
  }
  
  void fillMcYields(PlotMaker & pm, float luminosity, RowInformation & mcRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[mcRow.tableIndex].get());
    vector<GammaParams> yields = yieldTable->BackgroundYield(luminosity);
    for (unsigned ipar(0); ipar< yields.size(); ipar++) 
    {
      string const & label = mcRow.labels[ipar];
      addToMapYields(label, yields[ipar], mcRow.tableRows[ipar], mYields);
      if (verbose) cout<<label<<": "<<mYields.at(label).first.Yield()<<endl;
      //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
    }
  }
  
  void fillSignalYieldsProcesses(PlotMaker & pm, float luminosity, vector<shared_ptr<Process> > & signalProcesses, RowInformation & signalRow, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    Table * yieldTable;
    yieldTable = static_cast<Table*>(pm.Figures()[signalRow.tableIndex].get());
    for (auto & process : signalProcesses)
    {
      vector<GammaParams> yields = yieldTable->Yield(process.get(), luminosity);
      for (unsigned ipar(0); ipar< yields.size(); ipar++) 
      {
        string const & label = process->name_ + "_" +signalRow.labels[ipar];
        addToMapYields(label, yields[ipar], signalRow.tableRows[ipar], mYields);
        if (verbose) cout<<label<<": "<<mYields.at(label).first.Yield()<<endl;
        //cout<<label<<": "<<mYields.at(label).first.Yield()<<" "<<mYields.at(label).second.cut_.Name()<<" "<<mYields.at(label).second.weight_.Name()<<endl;
      }
    }
  }
  
  void fillAverageGenMetYields(vector<shared_ptr<Process> > & processes, vector<pair<string, string> > sampleBins, string const & signalTag, string const & signalGenMetTag, string const & signalAverageGenMetTag, map<string, pair<GammaParams, TableRow> > & mYields, bool verbose)
  {
    for (auto & process : processes)
    {
      string const & processName = process->name_;
      for(auto & bin : sampleBins)
      {
        string yieldLabel = processName + "_" + signalTag + "_" + bin.first;
        string genMetYieldLabel = processName + "_" + signalGenMetTag + "_" + bin.first;
        //cout<<yieldLabel<<" "<<genMetYieldLabel<<endl;
        float averageValue = 0.5*(mYields.at(yieldLabel).first.Yield() + mYields.at(genMetYieldLabel).first.Yield());
        //float averageError = hypot(0.5*mYields.at(yieldLabel).first.Uncertainty(),0.5*mYields.at(genMetYieldLabel).first.Uncertainty());
        float averageError = max(mYields.at(yieldLabel).first.Uncertainty(),mYields.at(genMetYieldLabel).first.Uncertainty());
        GammaParams average;
        average.SetYieldAndUncertainty(averageValue,averageError);
        TableRow blank("0.5("+yieldLabel+"+"+genMetYieldLabel+")");
        addToMapYields(processName+"_"+signalAverageGenMetTag+"_"+bin.first, average, blank, mYields);
        if (verbose) cout<<processName+"_"+signalAverageGenMetTag+"_"+bin.first<<" "<<averageValue<<" "<<averageError<<endl;
      }
    }
  }
}
