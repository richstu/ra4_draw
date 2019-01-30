///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

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

//// Lowers the njets requirements for dilepton bins
string lowerNjets(string &cut){
  string lowcut = cut;
  for(int nj=6; nj<=9; nj++){
    string nj_s = ""; nj_s += nj;
    string njlo_s = ""; njlo_s += nj-1;
    ReplaceAll(lowcut, nj_s, njlo_s);
  }
  return lowcut;
}

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/");
  // string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");
  //if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  Palette colors("txt/colors.txt", "default");

  string ntupletag = "";
  set<string> allfiles = {foldermc+"*_TTJets_*Lept*",
         foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
         foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
         foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
         foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",
         foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
         foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
         foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
         foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",foldermc+"*_ZZ_*"+ntupletag+"*.root"
       };

  // allfiles = set<string>({foldermc+"*_TTJets_Tune*"});

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch_met && mj14>250 && nleps>=1 && st>500 && met>100 && njets>=5"; // Excluding one QCD event

  map<string, vector<shared_ptr<Process> > > procs;

  procs["procs"] = vector<shared_ptr<Process> >();
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==2 && ntrutaush==0"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (#tau_{h}l)", Process::Type::background, colors("tt_ltau"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==2 && ntrutaush>=1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}Z", Process::Type::background, kOrange+7,
    {foldermc+"*_TTZ*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}W", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTW*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
  {foldermc+"*QCD_HT*0_Tune*.root",
    foldermc+"*QCD_HT*Inf_Tune*.root"},
    baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tt#gamma", Process::Type::background, kMagenta+3,
    {foldermc+"*_TTGJets*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tttt", Process::Type::background, kBlue+2,
    {foldermc+"*_TTTT*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kPink-2,
    {foldermc+"*DYJetsToLL*.root",
    foldermc+"*_ZJet*.root",
    foldermc+"*_ttHJetTobb*.root",
    foldermc+"*_WH_HToBB*.root",
    foldermc+"*_ZH_HToBB*.root",
    foldermc+"*_WWTo*.root",
    foldermc+"*_WZ*.root",
    foldermc+"*_ZZ_*.root"},
    baseline));


  NamedFunc multNeu = "(type==5000 || type==13000 || type==15000 || type==16000)";
  NamedFunc multNeu2l = "(ntruleps>=2 || (ntruleps<=1&&(type==5000 || type==13000 || type==15000 || type==16000)))";
  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#leq1#kern[.1]{#nu^{pr.}}, mismeas.", Process::Type::background, kRed-4, allfiles, 
           baseline && "mt>140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#leq1#kern[.1]{#nu^{pr.}}, well-meas.", Process::Type::background, kAzure-7, allfiles, 
           baseline && "mt<=140 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#geq2#nu^{prompt}", Process::Type::background, kCyan-3,
  			   allfiles, baseline && multNeu2l));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, #geq1#kern[.1]{#nu}^{non-prompt}", Process::Type::background, kGreen-3, 
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw==0.));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, off-shell W", Process::Type::background, kOrange,
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw>0.));

  PlotMaker pm;
  string smfile = "slide_pies.tex";
  SlideMaker sm(smfile,"1610");

  vector<string> metcuts;
  vector<vector<string>> mjcuts, mjcuts_old;
  // metcuts.push_back("met>100 && met<=150");
  // mjcuts.push_back(vector<string>({"mj14>250 && mj14<=400","mj14>400"}));

  metcuts.push_back("met>150 && met<=200");
  mjcuts_old.push_back(vector<string>({"mj14>250 && mj14<=400","mj14>400"}));
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=400", "mj14>400 && mj14<=500","mj14>500"})); //don't want to plot the same pie chart twice...
  
  metcuts.push_back("met>200 && met<=350");
  mjcuts_old.push_back(vector<string>({"mj14>250 && mj14<=400","mj14>400"}));
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=400", "mj14>400 && mj14<=500","mj14>500"}));
  
  metcuts.push_back("met>350 && met<=500");
  mjcuts_old.push_back(vector<string>({"mj14>250 && mj14<=400","mj14>400"}));
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=450", "mj14>450 && mj14<=650","mj14>650"}));
  
  metcuts.push_back("met>500");
  mjcuts_old.push_back(vector<string>({"mj14>250 && mj14<=400","mj14>400"}));
  mjcuts.push_back(vector<string>({"mj14>250 && mj14<=500", "mj14>500 && mj14<=800","mj14>800"}));

  vector<string> nbcuts;
  nbcuts.push_back("nbm==1");
  nbcuts.push_back("nbm==2");
  nbcuts.push_back("nbm>=3");  

  vector<string> njcuts, njcuts_old;
  njcuts_old.push_back("njets==5");
  njcuts_old.push_back("njets>=6 && njets<=8");
  njcuts_old.push_back("njets>=9");
  njcuts.push_back("njets==5");
  njcuts.push_back("njets>=6 && njets<=7");
  njcuts.push_back("njets>=8");

  vector<string> mtcuts;
  mtcuts.push_back("mt>140");
  // mtcuts.push_back("mt<=140");
  
  //// nleps = 1
  for(unsigned inj(0); inj<njcuts.size(); inj++){ 
    vector<TableRow> table_cuts;
    vector<string> pnames;
    for(unsigned imet(0); imet<metcuts.size(); imet++){ 
      for(unsigned imt(0); imt<mtcuts.size(); imt++){  
        for(unsigned imj(0); imj<mjcuts_old[imet].size(); imj++){  
          string cuts = "nleps==1 && nveto==0 && nbm>=1 && "+njcuts_old[inj]+"&&"+metcuts[imet]+"&&"+mtcuts[imt]+"&&"+mjcuts_old[imet][imj];
          table_cuts.push_back(TableRow("", cuts));  
          pnames.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(lumi,0)));  
        } 
        for(unsigned imj(0); imj<mjcuts[imet].size(); imj++){ 
          string cuts = "nleps==1 && nveto==0 && nbm>=1 && "+njcuts[inj]+"&&"+metcuts[imet]+"&&"+mtcuts[imt]+"&&"+mjcuts[imet][imj];
          table_cuts.push_back(TableRow("", cuts));  
          pnames.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(lumi,0)));  
        }
      }
    }
    for(auto &ipr: procs){
      string tag = "1l_abcd_"+ipr.first;
      pm.Push<Table>(tag, table_cuts, ipr.second, true, true, true, true);
      vector<string> col_labels = {"$p_{T}^{miss}$", "R3","Old R4","New low-$M_J$ R4","New high-$M_J$ R4"};
      vector<string> row_labels;
      for (auto &imet: metcuts) row_labels.push_back("$"+CodeToLatex(imet)+"$");
      sm.AddSlideWithReplace("XXX", tag, pnames, mjcuts[0].size()+mjcuts_old[0].size(), 
                             "$"+CodeToLatex(njcuts_old[inj])+"$", 
                             col_labels, row_labels);
    }
  }

  // alternate CR pie charts
  vector<string> crcuts;
  crcuts.push_back("nbm>=1 && njets==5");
  crcuts.push_back("nbm>=1 && njets==6");
  crcuts.push_back("nbm==1 && njets>=6");
  crcuts.push_back("nbm>=2 && njets>=7");
  
  for(unsigned imj(0); imj<mjcuts[0].size(); imj++){  
    vector<TableRow> table_cuts;
    vector<string> pnames;
    for(unsigned icr(0); icr<crcuts.size(); icr++){ 
      for(unsigned imet(0); imet<metcuts.size(); imet++){ 
        string cuts = "nleps==1 && nveto==0 && mt>140 && "+crcuts[icr]+"&&"+metcuts[imet]+"&&"+mjcuts[imet][imj];
        table_cuts.push_back(TableRow("", cuts));  
        pnames.push_back("pie_XXX_"+CodeToPlainText(cuts)+"_perc_lumi"+ToString(RoundNumber(lumi,0)));  
      } 
    }
    for(auto &ipr: procs){
      string tag = "altCRs_"+ipr.first;
      pm.Push<Table>(tag, table_cuts, ipr.second, true, true, true, true);

      vector<string> col_labels = {"CR"};
      for (auto &imet:metcuts) col_labels.push_back("$"+CodeToLatex(imet)+"$");
      vector<string> row_labels;
      for (auto &icr:crcuts) row_labels.push_back("$"+CodeToLatex(icr)+"$");
      sm.AddSlideWithReplace("XXX", tag, pnames, metcuts.size(), 
                             "$M_{J}$ region" + to_string(imj), 
                             col_labels, row_labels);
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);
  sm.Close();

  time(&endtime);
  cout<<endl<<"Making piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
