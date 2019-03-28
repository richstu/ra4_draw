///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/plot_opt.hpp"

using namespace std;
namespace{
  string tag = "nom";
  int year = 2016;
}

NamedFunc offshellw("offshellw",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=24) continue;
      if (b.mc_mass()->at(i) > 140.) {
        return 1;
      }
    }
    return 0;
  });

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  set<int> years;
  if (year==0) years = {2016, 2017, 2018};
  else years = {year};

  map<int, string> foldermc;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_01_18/mc/merged_mcbase_stdnj5/"; 

  set<string> vnames_all = { "_TTJets_*Lept", "_WJetsToLNu_HT","_ST_","_TTW","_TTZ", 
  "_DYJetsToLL_M-50_HT","_ZJet","_ttH", "_TTGJets","_TTTT","_WH_HToBB","_ZH_HToBB","_WWTo",
  "_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"};
  set<string> vnames_other = {"_DYJetsToLL_M-50_HT","_ZJet","_ttH", "_WH_HToBB","_ZH_HToBB",
  "_WWTo","_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"};

  set<string> tt1lfiles, tt2lfiles, ttxfiles, ttttfiles, wfiles, stfiles, otherfiles, allfiles;
  for (auto &yr: years) {
    tt1lfiles.insert(foldermc[yr]+"*_TTJets_SingleLept*.root");
    tt2lfiles.insert(foldermc[yr]+"*_TTJets_DiLept*.root");
    ttxfiles.insert(foldermc[yr]+"*_TTW*.root");
    ttxfiles.insert(foldermc[yr]+"*_TTZ*.root");
    ttxfiles.insert(foldermc[yr]+"*_TTGJets*.root");
    ttttfiles.insert(foldermc[yr]+"*_TTTT*.root");
    wfiles.insert(foldermc[yr]+"*_WJetsToLNu_HT*.root");
    stfiles.insert(foldermc[yr]+"*_ST_*.root");
    for(auto &name : vnames_other) otherfiles.insert(foldermc[yr] + "*" + name + "*.root");
    for(auto &name : vnames_all) allfiles.insert(foldermc[yr] + "*" + name + "*.root");
  }

  string baseline = "pass && stitch_met && mj14>250 && st>500 && met>100 && njets>=5";
  NamedFunc baselinef = baseline && Functions::hem_veto && "st<10000 && pass_ra2_badmu && met/met_calo<5";

  NamedFunc multNeu = "(type==5000 || type==13000 || type==15000 || type==16000)";
  NamedFunc multNeu2l = "(ntruleps>=2 || (ntruleps<=1&&(type==5000 || type==13000 || type==15000 || type==16000)))";

  Process::Type bkg = Process::Type::background;
  map<string, vector<shared_ptr<Process> > > procs;
  procs["procs"] = vector<shared_ptr<Process> >();
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (l)", bkg, colors("tt_1l"),
    tt1lfiles,baselinef));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (ll)", bkg, colors("tt_2l"),
    tt2lfiles, baselinef && "ntrutaush==0"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (#tau_{h}l)", bkg, colors("tt_ltau"),
    tt2lfiles, baselinef && "ntrutaush>=1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("W+jets", bkg, colors("wjets"),
    wfiles, baselinef));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t", bkg, colors("single_t"),
    stfiles, baselinef));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+Z/W/#gamma", bkg, kOrange+7,
    ttxfiles, baselinef));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tttt", bkg, kBlue+2,
    ttttfiles, baselinef));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Other", bkg, kPink-2,
    otherfiles, baselinef));

  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>
    ("#geq2#nu^{prompt}", bkg, kCyan-3,
    allfiles, baselinef && multNeu2l));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
    ("#leq1#kern[.1]{#nu^{pr.}}, well-meas.", bkg, kAzure-7, 
    allfiles, baselinef && "mt<=140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
    ("#leq1#kern[.1]{#nu^{pr.}}, mismeas.", bkg, kRed-4, 
    allfiles, baselinef && "mt>140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
    ("#leq1#kern[.1]{#nu^{pr.}}, #geq1#kern[.1]{#nu}^{non-prompt}", bkg, kGreen-3, 
    allfiles, baselinef && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw==0.));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
    ("#leq1#kern[.1]{#nu^{pr.}}, off-shell W", bkg, kOrange, 
    allfiles, baselinef && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw>0.));

  PlotMaker pm;
  string smfile = "slide_pies_"+tag+".tex";
  SlideMaker sm(smfile,"1610");

  vector<string> metcuts;
  vector<vector<string>> mjcuts;

  if (Contains(tag, "1l")) {
    metcuts.push_back("met>100 && met<=150");
    mjcuts.push_back(vector<string>({"mj14>250 && mj14<=350", "mj14>350 && mj14<=450","mj14>450"}));

    metcuts.push_back("met>150 && met<=200");
    mjcuts.push_back(vector<string>({"mj14>250 && mj14<=350", "mj14>350 && mj14<=450","mj14>450"}));
  }
  metcuts.push_back("met>200 && met<=350");
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=400", "mj14>400 && mj14<=500","mj14>500"}));
  
  metcuts.push_back("met>350 && met<=500");
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=450", "mj14>450 && mj14<=650","mj14>650"}));
  
  metcuts.push_back("met>500");
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=500", "mj14>500 && mj14<=800","mj14>800"}));

  vector<string> nbcuts;
  nbcuts.push_back("nbd==1");
  nbcuts.push_back("nbd==2");
  nbcuts.push_back("nbd>=3");  

  vector<string> njcuts;
  if (Contains(tag, "1l")) {
    if (Contains(tag, "nj")) {
      njcuts.push_back("njets==5");
      njcuts.push_back("njets==6");
      njcuts.push_back("njets==7");
      njcuts.push_back("njets>=8");
    } else {
      njcuts.push_back("njets>=5");
    }
  } else if (Contains(tag, "2l")) {
      njcuts.push_back("njets==6");
      njcuts.push_back("njets>=7");
  } else if (Contains(tag, "veto")) {
      njcuts.push_back("njets==7");
      njcuts.push_back("njets>=8");
  }
  
  NamedFunc w = Functions::wgt_run2 * Functions::eff_trig_run2;
  
  string regcuts = "nleps==1 && nveto==0 && nbd>=1";
  if (Contains(tag, "2l")) regcuts = "nleps==2 && nbd<=1";
  if (Contains(tag, "veto")) regcuts = "nleps==1 && nveto==1 && nbd==1";

  //// nleps = 1      
  vector<TableRow> table_cuts_1l;
  vector<string> pnames_1l;
  for(unsigned inj(0); inj<njcuts.size(); inj++){ 
    for(unsigned imet(0); imet<metcuts.size(); imet++){ 
        string cuts = regcuts + " && mt>140 && "+njcuts[inj]+"&&"+metcuts[imet];
        table_cuts_1l.push_back(TableRow("", cuts,1,1,w));  
        pnames_1l.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(1,0)));  
    }
  }
  for(auto &ipr: procs){
    string label = "1l_abcd_"+ipr.first;
    pm.Push<Table>(label, table_cuts_1l, ipr.second, false, true, true, false);
    
    vector<string> col_labels = {""};
    for (auto &imet:metcuts) col_labels.push_back("$"+CodeToLatex(imet)+"$");
    vector<string> row_labels;
    for (auto &inj:njcuts) row_labels.push_back("$"+CodeToLatex(inj)+"$");
    string ttl = "1$\\ell$ sample";
    if (Contains(tag, "2l")) ttl = "2$\\ell$ sample";
    if (Contains(tag, "veto")) ttl = "1$\\ell$ + 1 veto track sample";
    sm.AddSlideWithReplace("XXX", label, pnames_1l, metcuts.size(), ttl, col_labels, row_labels);
  }
  


  // 2L CR pie charts
  // map<string, string> cuts_2lveto, njcuts_2lveto;
  
  // cuts_2lveto["2l"] = "nleps==2 && nbd<=1";
  // njcuts_2lveto["2l"] = {"njets==6", "njets>=7"};

  // cuts_2lveto["veto"] = "mt>140 && nleps==1 && nveto==1 && nbd==1";
  // njcuts_2lveto["veto"] = {"njets==7", "njets>=8"};
  
  // vector<TableRow> table_cuts_2l;
  // vector<string> pnames_2l;
  // if (Contains(tag, "2l")) {
  //   for(auto &icr: {"2l","veto"}){ 
  //     for (unsigned imet(0); imet<metcuts.size(); imet++) {
  //       for(unsigned inj(0); inj<njcuts_2lveto[icr].size(); inj++){ 
  //         string cuts = cuts_2lveto[icr]+"&&"+njcuts_2lveto[icr][inj];
  //         if (imet==(metcuts.size()-1)) cuts = cuts_2lveto[icr]+"&&"+njcuts_2lveto[icr][inj];
  //         table_cuts_2l.push_back(TableRow("", cuts,1,1,w));  
  //         pnames_2l.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(1,0)));  
  //       }
  //     }

  //     for(auto &ipr: procs){
  //       string label = "m2lveto_"+ipr.first;
  //       pm.Push<Table>(label, table_cuts_2l, ipr.second, false, true, true, true);

  //       vector<string> col_labels = {"CR"};
  //       for (auto &imet:metcuts) col_labels.push_back("$"+CodeToLatex(imet)+"$");
  //       vector<string> row_labels = {"Low $N_{jets}$", "High $N_{jets}$"};
  //       // for (auto &icr:cuts_2lveto) row_labels.push_back("$"+CodeToLatex(icr)+"$");
  //       sm.AddSlideWithReplace("XXX", label, pnames_2l, njcuts_2lveto[0].size(), 
  //                              "$2\\ell$ control region", 
  //                              col_labels, row_labels);
  //     }
  //   }
  // }
  pm.min_print_ = true;
  pm.MakePlots(1);
  sm.Close();

  time(&endtime);
  cout<<endl<<"Making piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"year", required_argument, 0, 'y'},
      {"tag",      no_argument, 0, 't'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "y:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'y':
      year = atoi(optarg);
      break;
    case 't':
      tag = optarg;
      break;      
    case 0:
      optname = long_options[option_index].name;
      if(false){

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
