///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 1.;
  int year = 2018;
  bool only_tt = false;
  bool debug = false;
  string tag = "nbd";
  pair<string, string> sig_nc = make_pair("2100","100");
  pair<string, string> sig_c = make_pair("1900","1250");
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years;
  if (year==0) years = {2016, 2017, 2018};
  else years = {year};

  map<int, string> foldermc, folderdata, foldersig;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
  foldersig[2016] = bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
  foldersig[2017] = "";//bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_01_18/mc/merged_mcbase_stdnj5/";
  foldersig[2018] = "";//bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/");
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_01_18/data/merged_database_standard/";
  
  set<string> vnames_other = {
    "_WJetsToLNu_HT","_ST_","_TTW","_TTZ", "_DYJetsToLL_M-50_HT","_ZJet","_ttH",
     "_TTGJets","_TTTT","_WH_HToBB","_ZH_HToBB","_WWTo","_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"
  };

  set<string> t1nc_files, t1c_files, tt1l_files, tt2l_files, tt_files, other_files, data_files;
  for (auto &yr: years) {
    t1nc_files.insert(foldersig[yr]+"*mGluino-"+sig_nc.first+"_mLSP-"+sig_nc.second+"_*.root");
    t1c_files.insert(foldersig[yr]+"*mGluino-"+sig_c.first+"_mLSP-"+sig_c.second+"_*.root");
    tt1l_files.insert(foldermc[yr]+"*_TTJets*SingleLept*.root");
    tt2l_files.insert(foldermc[yr]+"*_TTJets*DiLept*.root");
    tt_files.insert(foldermc[yr]+"*_TTJets*Lept*.root");
    data_files.insert(folderdata[yr]+"*root");
    for(auto name : vnames_other)
      other_files.insert(foldermc[yr] + "*" + name + "*.root");
      
  }

  NamedFunc baseline = "mj14>250 && st>500 && nleps==1 && nveto==0 && met>100 && njets>=5 && nbd>=1";
  baseline = baseline && Functions::hem_veto && "st<10000 && pass_ra2_badmu && met/met_calo<5";

  vector<shared_ptr<Process> > procs;
  // procs.push_back(Process::MakeShared<Baby_full>("t\\bar{t} (2\\ell)", Process::Type::background, 1, tt2l_files, 
  //                 baseline && "pass && stitch_met"));
  // procs.push_back(Process::MakeShared<Baby_full>("t\\bar{t} (1\\ell)", Process::Type::background, 1, tt1l_files, 
  //                 baseline && "pass && stitch_met"));
  // if (!only_tt)
  //   procs.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, 1, other_files, 
  //                 baseline && "pass && stitch_met"));
  // procs.push_back(Process::MakeShared<Baby_full>("("+sig_nc.first+","+sig_nc.second+")", Process::Type::signal, 1, t1nc_files, baseline));
  // procs.push_back(Process::MakeShared<Baby_full>("("+sig_c.first+","+sig_c.second+")", Process::Type::signal, 1, t1c_files, baseline));
  

  procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    data_files, baseline && Functions::trig_run2 && "pass"));
  vector<string> met, mjl, mjh;
  met.push_back("met>100 && met<=200"); mjl.push_back("350"); mjh.push_back("450");
  met.push_back("met>200 && met<=350"); mjl.push_back("400"); mjh.push_back("500");
  met.push_back("met>350 && met<=500"); mjl.push_back("450"); mjh.push_back("650");
  if (tag!="56j") {
    met.push_back("met>500");             mjl.push_back("500"); mjh.push_back("800");
  }

  vector<string> nbnj;
  if (tag=="56j"){
    nbnj.push_back("njets>=5 && njets<=6 && nbd==1");
    nbnj.push_back("njets>=5 && njets<=6 && nbd==2");
    nbnj.push_back("njets>=5 && njets<=6 && nbd>=3");
  } else {
    nbnj.push_back("nbd==1 && njets<=7");
    nbnj.push_back("nbd==1 && njets>=8");
    nbnj.push_back("nbd==2 && njets<=7");
    nbnj.push_back("nbd==2 && njets>=8");
    nbnj.push_back("nbd>=3 && njets<=7");
    nbnj.push_back("nbd>=3 && njets>=8");
  }

  vector<TString> abcd_lo  = {"mt<=140 && mj14<=MJ1X",
                              "mt<=140 && mj14> MJ1X && mj14<=MJ2X",
                              "mt>140  && mj14<=MJ1X",
                              "mt>140  && mj14> MJ1X && mj14<=MJ2X"};
  vector<TString> abcd_hi  = {"mt<=140 && mj14<=MJ1X",
                              "mt<=140 && mj14> MJ2X",
                              "mt>140  && mj14<=MJ1X",
                              "mt>140  && mj14> MJ2X"};
  vector<vector<TString>> abcds = {abcd_lo, abcd_hi};
  vector<string> abcd_names = {"lowmj","highmj"};

  NamedFunc wgt = Functions::wgt_run2 * Functions::eff_trig_run2;

  PlotMaker pm;
  vector<string> tabnames;
  for (size_t iabcd(0); iabcd<abcds.size(); iabcd++){
    tabnames.push_back("regions_"+tag+"_"+to_string(year)+"_"+abcd_names[iabcd]);
    vector<TableRow> table_rows;
    for (size_t imet(0); imet<met.size(); imet++){
      table_rows.push_back(TableRow("$"+CodeToLatex(met[imet])+"$"));
      for (size_t ir(0); ir<abcds[iabcd].size(); ir++){
        if (ir%2==1) {
          for (auto &ibj: nbnj) {
            TString _cut = met[imet] + " && " + abcds[iabcd][ir] + " && " + ibj;
            _cut.ReplaceAll("MJ1X", mjl[imet]).ReplaceAll("MJ2X", mjh[imet]);
            if (debug) cout<<"R"<<ir+1<<": "<<_cut<<endl;
            table_rows.push_back(TableRow("R"+to_string(ir+1), _cut,0,0, wgt));
          }
        } else {
          TString _cut = met[imet] + " && " + abcds[iabcd][ir];
          if (tag=="56j") _cut += "&& njets>=5 && njets<=6";
          _cut.ReplaceAll("MJ1X", mjl[imet]).ReplaceAll("MJ2X", mjh[imet]);
          if (debug) cout<<"R"<<ir+1<<": "<<_cut<<endl;
          table_rows.push_back(TableRow("R"+to_string(ir+1), _cut,1,1, wgt));
        }
      }
    }
    pm.Push<Table>(tabnames.back(), table_rows,procs,0);
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  for (auto &itab: tabnames) {
    TString _fname = itab+"_lumi_"+RoundNumber(lumi,1).ReplaceAll(".","p")+".tex";
    execute(("pdflatex -output-directory=tables tables/"+_fname+" > /dev/null").Data());
    _fname.ReplaceAll(".tex",".pdf");
    cout<<endl<<"open tables/"+_fname<<endl;
  }
  time(&endtime);
  cout<<endl<<"Making regions tables took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"debug",      no_argument, 0, 'd'}, 
      {"year",     required_argument, 0, 'y'},   
      {"tag",     required_argument, 0, 't'},   
      {"tt",      no_argument, 0, 0},   
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "dy:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'd':
      debug = true;
      break;
    case 'y':
      year = atoi(optarg);
      break;
    case 't':
      tag = optarg;
      break;
    case 0:
      if(optname == "tt"){
        only_tt = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

