
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <limits>

#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "core/utilities.hpp"
#include "core/cross_sections.hpp"

using namespace std;

double GetSignif(const std::string &filename);
std::string GetBaseName(const std::string &path);
void GetOptions(int argc, char *argv[]);

namespace{
  string pwd = getenv("PWD");
  string indir = "";
  string outdir = pwd;
  string datacard = "";
  bool do_signif = false;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  if(datacard == "") ERROR("Must supply an input file name");

  string model = "T1tttt";
  if(Contains(datacard, "T5tttt")) model = "T5tttt";
  if(Contains(datacard, "TChiHH")) model = "TChiHH";

  //// Parsing the gluino and LSP masses
  int mglu, mlsp;
  parseMasses(datacard, mglu, mlsp);
  double xsec, xsec_unc;
  if(model=="T1tttt" || model=="T5tttt") xsec::signalCrossSection(mglu, xsec, xsec_unc);
  else if(model == "TChiHH") xsec::higgsinoCrossSection(mglu, xsec, xsec_unc);
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
 
  ostringstream command;
  string done = "; ";
  command
    << "cd " << outdir << done;
  if (mglu < 801)
    command << "combine -M AsymptoticLimits --rMax 0.1 " << pwd << "/" << indir << "/" << datacard << done;
  else if (mglu < 901)
    command << "combine -M AsymptoticLimits --rMax 0.5 " << pwd << "/" << indir << "/" << datacard << done;
  else if (mglu < 1151)
    command << "combine -M AsymptoticLimits --rMax 2 " << pwd << "/" << indir << "/" << datacard << done;
  else 
    command << "combine -M AsymptoticLimits " << pwd << "/" << indir << "/" << datacard << done;
  if(do_signif){
    command
      << "combine -M ProfileLikelihood --significance --expectSignal=1 " 
      "--verbose=999999 --rMin=-10. --uncapped=1 " << pwd << "/" << indir << "/" << datacard
      << " < /dev/null &> signif_obs.log; "
      << "combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 "
      "--verbose=999999 --rMin=-10. --uncapped=1 --toysFreq " << pwd << "/" << indir << "/" << datacard
      << " < /dev/null &> signif_exp.log; ";
  }
  command << flush;
  execute(command.str());
  
  string limits_file_name = outdir+"/higgsCombineTest.AsymptoticLimits.mH120.root";
  TFile limits_file(limits_file_name.c_str(), "read");
  if(!limits_file.IsOpen()) ERROR("Could not open limits file "+limits_file_name);
  TTree *tree = static_cast<TTree*>(limits_file.Get("limit"));
  if(tree == nullptr) ERROR("Could not get limits tree");
  double limit;
  tree->SetBranchAddress("limit", &limit);
  int num_entries = tree->GetEntries();
  if(num_entries != 6) ERROR("Expected 6 tree entries. Saw "+to_string(num_entries));
  tree->GetEntry(0);
  double exp_2down = limit;
  tree->GetEntry(1);
  double exp_down = limit;
  tree->GetEntry(2);
  double exp = limit;
  tree->GetEntry(3);
  double exp_up = limit;
  tree->GetEntry(4);
  double exp_2up = limit;
  tree->GetEntry(5);
  double obs = limit;
  limits_file.Close();

  double sig_obs, sig_exp;
  if(do_signif){
    sig_obs = GetSignif(outdir+"/signif_obs.log");
    sig_exp = GetSignif(outdir+"/signif_exp.log");
  }

  cout << setprecision(numeric_limits<double>::max_digits10)
    << ' ' << mglu << ' ' << mlsp
    << ' ' << xsec << ' ' << xsec_unc
    << ' ' << obs << ' ' << obs << ' ' << obs 
    << ' ' << exp << ' ' << exp_up << ' ' << exp_down << ' ' << exp_2up << ' ' << exp_2down;
  if(do_signif)
    cout << ' ' << sig_obs << ' ' << sig_exp;
  cout << endl;
}

double GetSignif(const string &filename){
  double signif = 0.;
  ifstream file(filename);
  string line;
  while(getline(file, line)){
    auto pos = line.find("Significance: ");
    if(pos != 0) continue;
    string val = line.substr(14);
    signif = stod(val);
  }
  return signif;
}

string GetBaseName(const string &path){
  auto pos = path.rfind("/");
  if(pos == string::npos){
    return path;
  }else{
    return path.substr(pos+1);
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"datacard", required_argument, 0, 'd'},
      {"indir", required_argument, 0, 'i'},
      {"outdir", required_argument, 0, 'o'},
      {"signif", required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "d:i:o:s", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'd':
      datacard = optarg;
      break;
    case 'i':
      indir = optarg;
      break;
    case 'o':
      outdir = optarg;
      break;
    case 's':
      do_signif = false;
      break;
    default:
      cerr << "Bad option! getopt_long returned character code " << static_cast<int>(opt) << endl;
      break;
    }
  }
}
