///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
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
  float lumi = 135;
  string ntup_date = "2019_01_11";
  string filetag = "2016_new";
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

  string foldermc(bfolder+"/cms2r0/babymaker/babies/"+ntup_date+"/mc/merged_mcbase_stdnj5/");
  // string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");
  //if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  Palette colors("txt/colors.txt", "default");

  string ntupletag = "";
  set<string> allfiles = {foldermc+"*_TTJets_*Lept*",
         foldermc+"*_WJetsToLNu_HT*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
         foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
         foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
         foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",
         foldermc+"*DYJetsToLL_M-50_HT*"+ntupletag+"*.root",
         foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttH*"+ntupletag+"*.root",
         foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
         foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",foldermc+"*_ZZ_*"+ntupletag+"*.root"
       };

  // allfiles = set<string>({foldermc+"*_TTJets_Tune*"});

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch_met && mj14>250 && nleps>=1 && st>500 && met>100 && njets>=5 && st<10000 && pass_ra2_badmu && met/met_calo<5"; // Excluding one QCD event

  map<string, vector<shared_ptr<Process> > > procs;

  procs["procs"] = vector<shared_ptr<Process> >();
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root"},
    baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root"},
    baseline+" && ntrutaush==0"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (#tau_{h}l)", Process::Type::background, colors("tt_ltau"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntrutaush>=1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu_HT*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+Z/W/#gamma", Process::Type::background, kOrange+7,
    {foldermc+"*_TTZ*.root", foldermc+"*_TTW*.root", foldermc+"*_TTGJets*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tttt", Process::Type::background, kBlue+2,
    {foldermc+"*_TTTT*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kPink-2,
    {foldermc+"*DYJetsToLL_M-50_HT*.root",
    foldermc+"*_ZJet*.root",
    foldermc+"*_ttH*.root",
    foldermc+"*_WH_HToBB*.root",
    foldermc+"*_ZH_HToBB*.root",
    foldermc+"*_WWTo*.root",
    foldermc+"*_WZ*.root",
    foldermc+"*_ZZ_*.root",
    foldermc+"*QCD_HT*0_Tune*.root",
    foldermc+"*QCD_HT*Inf_Tune*.root"},
    baseline));


  NamedFunc multNeu = "(type==5000 || type==13000 || type==15000 || type==16000)";
  NamedFunc multNeu2l = "(ntruleps>=2 || (ntruleps<=1&&(type==5000 || type==13000 || type==15000 || type==16000)))";
  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#geq2#nu^{prompt}", Process::Type::background, kCyan-3,
           allfiles, baseline && multNeu2l));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#leq1#kern[.1]{#nu^{pr.}}, well-meas.", Process::Type::background, kAzure-7, allfiles, 
           baseline && "mt<=140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#leq1#kern[.1]{#nu^{pr.}}, mismeas.", Process::Type::background, kRed-4, allfiles, 
           baseline && "mt>140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, #geq1#kern[.1]{#nu}^{non-prompt}", Process::Type::background, kGreen-3, 
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw==0.));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, off-shell W", Process::Type::background, kOrange,
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw>0.));

  PlotMaker pm;
  string smfile = "slide_pies_"+filetag+".tex";
  SlideMaker sm(smfile,"1610");

  vector<string> metcuts;
  vector<vector<string>> mjcuts;

  metcuts.push_back("met>100 && met<=150");
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=350", "mj14>350 && mj14<=450","mj14>450"}));

  metcuts.push_back("met>150 && met<=200");
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=350", "mj14>350 && mj14<=450","mj14>450"}));
  
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
  njcuts.push_back("njets==5");
  njcuts.push_back("njets>=6 && njets<=7");
  njcuts.push_back("njets>=8");

  vector<string> mtcuts;
  mtcuts.push_back("mt>140");
  // mtcuts.push_back("mt<=140");
  
  //// nleps = 1      
  for(unsigned imj(0); imj<3; imj++){  
    vector<TableRow> table_cuts_1l;
    vector<string> pnames_1l;
    for(unsigned imt(0); imt<mtcuts.size(); imt++){  
      for(unsigned inj(0); inj<njcuts.size(); inj++){ 
        for(unsigned imet(0); imet<metcuts.size(); imet++){ 
            string cuts = "nleps==1 && nveto==0 && nbd>=1 && "
                          +njcuts[inj]+"&&"+metcuts[imet]+"&&"+mtcuts[imt]+"&&"+mjcuts[imet][imj];
            table_cuts_1l.push_back(TableRow("", cuts));  
            pnames_1l.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(lumi,0)));  
        }
      }
      for(auto &ipr: procs){
        string tag = "1l_abcd_"+ipr.first;
        pm.Push<Table>(tag, table_cuts_1l, ipr.second, true, true, true, true);
        
        vector<string> col_labels = {""};
        for (auto &imet:metcuts) col_labels.push_back("$"+CodeToLatex(imet)+"$");
        vector<string> row_labels;
        for (auto &inj:njcuts) row_labels.push_back("$"+CodeToLatex(inj)+"$");
        sm.AddSlideWithReplace("XXX", tag, pnames_1l, metcuts.size(), 
                               "5-jet control and signal regions, high $m_T$ composition, M$_J$ idx "+to_string(imj), 
                               col_labels, row_labels);
      }
    }
  }


  // 2L CR pie charts
  vector<string> cuts_2lveto; 
  vector<vector<string>> njcuts_2lveto;
  
  cuts_2lveto.push_back("nleps==2 && nbd<=2 && met<=500");
  njcuts_2lveto.push_back({"njets>=5 && njets<=6", "njets>=7"});

  cuts_2lveto.push_back("mt>140 && nleps==1 && nveto==1 && nbd>=1 && nbd<=2 && met<=500");
  njcuts_2lveto.push_back({"njets>=6 && njets<=7", "njets>=8"});
  
  vector<TableRow> table_cuts_2l;
  vector<string> pnames_2l;
  for(unsigned icr(0); icr<cuts_2lveto.size(); icr++){ 
    for(unsigned inj(0); inj<njcuts_2lveto[icr].size(); inj++){ 
      string cuts = cuts_2lveto[icr]+"&&"+njcuts_2lveto[icr][inj];
      table_cuts_2l.push_back(TableRow("", cuts));  
      pnames_2l.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(lumi,0)));  
    }
  }
  for(auto &ipr: procs){
    string tag = "m2lveto_"+ipr.first;
    pm.Push<Table>(tag, table_cuts_2l, ipr.second, true, true, true, true);

    vector<string> col_labels = {"CR","Low $N_{jets}$", "High $N_{jets}$"};
    vector<string> row_labels = {"2$\\ell$ only", "1$\\ell$ + 1 track"};
    // for (auto &icr:cuts_2lveto) row_labels.push_back("$"+CodeToLatex(icr)+"$");
    sm.AddSlideWithReplace("XXX", tag, pnames_2l, njcuts_2lveto[0].size(), 
                           "$2\\ell$ control region", 
                           col_labels, row_labels);
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);
  sm.Close();

  time(&endtime);
  cout<<endl<<"Making piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"ntup_date", required_argument, 0, 'n'},
      {"filetag",      no_argument, 0, 't'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "n:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      ntup_date = optarg;
      break;
    case 't':
      filetag = optarg;
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
