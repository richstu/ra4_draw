// Macro to make figure 22 from AN2016_187_v9

#include <cmath>
#include <stdio.h>
#include <chrono>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/table.hpp"
#include "core/slide_maker.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  double lumi =35.9;
  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string ntupletag=""; 
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/");

  Palette colors("txt/colors.txt", "default");


  NamedFunc baseline_1l = "nleps==1 && nveto==0 && njets>=7";
  NamedFunc baseline_2l = "nleps==2 && njets>=6";

  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root"};

  /// Study of ttbar 1l
  auto proc_1l_tt1l_lomt(Process::MakeShared<Baby_full>("t#bar{t} (1l), 1l m_{T}#leq140", 
    Process::Type::background, 1, ttfiles, 
    baseline_1l && "stitch_met && ntruleps<=1 && mt<=140"));
  auto proc_1l_ttltau(Process::MakeShared<Baby_full>("t#bar{t} (l#tau_{h}), 1l m_{T}>140", 
    Process::Type::background, kBlue-6, ttfiles, 
    baseline_1l && "stitch_met && ntruels+ntrumus+ntrutausl==1 && ntrutaush==1 && mt>140"));
  auto proc_1l_tt2l(Process::MakeShared<Baby_full>("t#bar{t} (2l), 1l", 
    Process::Type::background, colors("tt_2l"), ttfiles, 
    baseline_1l && "stitch_met && ntruels+ntrumus+ntrutausl==2"));
  auto proc_1l_tt1l_ghimt(Process::MakeShared<Baby_full>("t#bar{t} (1l, m_{T}^{tru}>140), 1l m_{T}>140", 
    Process::Type::background, kGreen-3, ttfiles, 
    baseline_1l && "stitch_met && ntruleps<=1 && mt>140&&mt_tru>140"));
  auto proc_1l_tt1l_bhimt(Process::MakeShared<Baby_full>("t#bar{t} (1l, m_{T}^{tru}#leq140), 1l m_{T}>140", 
    Process::Type::background, kRed-4, ttfiles, 
    baseline_1l && "stitch_met && ntruleps<=1 && mt>140&&mt_tru<140"));

  auto proc_2l_tt2l(Process::MakeShared<Baby_full>("t#bar{t} (2l), 1l", 
    Process::Type::background, colors("tt_2l"), ttfiles, 
    baseline_2l && "stitch_met && ntruels+ntrumus+ntrutausl==2"));

  // vector<shared_ptr<Process> > procs_1l = {proc_1l_tt1l_lomt, proc_1l_tt2l, proc_1l_ttltau, 
  //                                          proc_1l_tt1l_ghimt, proc_1l_tt1l_bhimt};
  
  vector<shared_ptr<Process> > procs_1l = {proc_1l_tt1l_lomt, proc_1l_tt2l};                                        
  vector<shared_ptr<Process> > procs_2l = {proc_1l_tt1l_lomt, proc_2l_tt2l};
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::shapes)
  .RatioMaximum(1.5)
  .RatioMinimum(0.5);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plots = {lin_shapes};//.PrintVals(true)};


  vector<string> njcuts;
  njcuts.push_back("njets==5");
  njcuts.push_back("njets==6");
  njcuts.push_back("njets==7");
  njcuts.push_back("njets>=8");

  vector<string> metcuts, mjcuts;
  metcuts.push_back("met>150 && met<=200"); mjcuts.push_back("mj14>450");
  metcuts.push_back("met>200 && met<=350"); mjcuts.push_back("mj14>500");
  metcuts.push_back("met>350 && met<=500"); mjcuts.push_back("mj14>650");
  metcuts.push_back("met>500");             mjcuts.push_back("mj14>800");

  vector<string> nbcuts;
  nbcuts.push_back("nbd==1");
  nbcuts.push_back("nbd==2");
  nbcuts.push_back("nbd>=3");

  PlotMaker pm;

  NamedFunc njets2("njets2",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nleps()==1) return b.njets();
    else return b.njets()+1;
  });

  // for (auto &imj:mjcuts) {
    // for (auto &inj:njcuts) {
    //   pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),inj && imj, procs_1l, plots).Tag("1l");
    //   pm.Push<Hist1D>(Axis(4,-0.5,3.5,Functions::ntrub,"True N_{b}"),inj && imj, procs_1l, plots).Tag("1l");
    // }
  //   for (auto &imet:metcuts) {
  //     pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),imet && imj, procs_1l, plots).Tag("1l");
  //     pm.Push<Hist1D>(Axis(4,-0.5,3.5,Functions::ntrub,"True N_{b}"),imet && imj, procs_1l, plots).Tag("1l");
  //   }
  // }
  for (size_t i(0); i<metcuts.size(); i++) {
    pm.Push<Hist1D>(Axis(5,6.5,11.5,"njets","N_{jets}"),metcuts[i] && mjcuts[i], procs_1l, plots).Tag("1l");
    pm.Push<Hist1D>(Axis(5,6.5,11.5,njets2,"N'_{jets}"),metcuts[i] && mjcuts[i], procs_2l, plots).Tag("2l");
    // pm.Push<Hist1D>(Axis(25,0,1200,"mj14","M_{J} [GeV]"),metcuts[i], procs_1l, plots).Tag("1l");
    // pm.Push<Hist1D>(Axis(25,0,1200,"mj14","M_{J} [GeV]"),metcuts[i], procs_2l, plots).Tag("2l");
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)
      <<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

