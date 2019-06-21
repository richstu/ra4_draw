#include "ra4/scatter.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 1;

  string bfolder = "";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    bfolder = "/net/cms2";
  }
  Palette colors("txt/colors.txt", "default");

  TString baseline("st>500 && met>200 && njets>=6 && nbdm>=1 && nleps==1 && nveto==0");
  NamedFunc filters = Functions::hem_veto && Functions::pass_run2;
  NamedFunc nom_wgt = Functions::wgt_run2 * Functions::eff_trig_run2; 

  set<int> years = {2016, 2017, 2018};

  map<int, string> foldermc, folderdata, foldersig;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
  foldersig[2016] = bfolder+"/cms2r0/babymaker/babies/2019_05_16/T1tttt/skim_sys_abcd/";
  folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
  foldersig[2017] = bfolder+"/cms2r0/babymaker/babies/2019_05_17/T1tttt/skim_sys_abcd/";
  folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_stdnj5/";
  foldersig[2018] = bfolder+"/cms2r0/babymaker/babies/2019_05_18/T1tttt/skim_sys_abcd/";
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/merged_database_standard/";

  // Filling all other processes
  vector<string> vnames_other = {
    "_WJetsToLNu_HT", "_ST_", "_TTW","_TTZ", 
    "_DYJetsToLL_M-50_HT","_ZJet","_ttH",
     "_TTGJets","_TTTT","_WH_HToBB","_ZH_HToBB","_WWTo","_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"
    };

  set<string> tt1l_files, tt2l_files, other_files, data_files, sig_files;
  for (auto &yr: years) {
    tt1l_files.insert(foldermc[yr]+"*_TTJets*SingleLept*.root");
    tt2l_files.insert(foldermc[yr]+"*_TTJets*DiLept*.root");
    data_files.insert(folderdata[yr]+"*root");
    for(auto name : vnames_other)
      other_files.insert(foldermc[yr] + "*" + name + "*.root");

    sig_files.insert(foldersig[yr]+"*mGluino-2100_mLSP-100_*.root");
  }

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors.RGB(1,57,166),
    tt1l_files, baseline && filters && "stitch_met");
  tt1l->SetMarkerStyle(23);
  tt1l->SetMarkerSize(0.8);
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors.RGB(86,160,211),
    tt2l_files, baseline && filters && "stitch_met");
  tt2l->SetMarkerStyle(22);
  tt2l->SetMarkerSize(0.8);
  auto other = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    other_files, baseline && filters && "stitch_met");

  auto t1tttt = Process::MakeShared<Baby_full>("T1tttt(2100,100)", Process::Type::signal, colors("t1tttt"),
    sig_files, baseline && filters);
  t1tttt->SetMarkerStyle(21);
  t1tttt->SetMarkerSize(0.9);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    data_files, baseline && filters && Functions::trig_run2);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.);

  vector<shared_ptr<Process> > all_procs = {data, t1tttt, tt1l, tt2l, other};
  vector<shared_ptr<Process> > tt_sig = {tt1l, tt2l, t1tttt};

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm).Title(TitleType::data)};
  vector<PlotOpt> bkg_pts = {style().Stack(StackType::lumi_shapes).Title(TitleType::simulation)};

  vector<NamedFunc> met_bins = {"met>200 && njets>=7", "met>200&&met<=350 && njets>=7", "met>350&&met<=500 && njets>=7", "met>500 && njets>=6"};
  vector<set<double>> mj_lines = {{250, 400},{250, 400, 500},{250, 450, 650},{250, 500, 800}};
  vector<NamedFunc> nb_bins = {"nbdm>=2"};

  PlotMaker pm;
  for(unsigned imet(0); imet<met_bins.size(); imet++){
    for(const auto &nbdm_bin: nb_bins){
      NamedFunc cut = met_bins[imet] && nbdm_bin;
      if (imet==0) {
        pm.Push<Hist2D>(Axis(48, 0., 1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
          Axis(175, 0., 700., "mt", "m_{T} [GeV]", {140.}),
          cut, tt_sig, bkg_pts).Weight(nom_wgt);
      } else if (imet>1) {
        pm.Push<Hist2D>(Axis(48, 0., 1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
          Axis(25, 0., 700., "mt", "m_{T} [GeV]", {140.}),
          cut, all_procs, bkg_hist).Weight(nom_wgt);
      }
    }
  }

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_=true;
  pm.MakePlots(lumi);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
