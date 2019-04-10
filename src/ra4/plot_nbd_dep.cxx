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


  NamedFunc baseline = "stitch_met && nleps==1 && nveto==0 && st>500 && met>100 && njets>=5 && weight<1";

  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root"};

  /// Study of ttbar 1l
  auto proc_tt1l_lomt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}#leq140", Process::Type::background, 
						       1, ttfiles, 
						       baseline && "ntruleps<=1 && mt<=140");
  auto proc_ttltau = Process::MakeShared<Baby_full>("t#bar{t} l#tau_{h}, m_{T}>140", Process::Type::background, 
						    kBlue-6, ttfiles, 
						    baseline && "ntruels+ntrumus+ntrutausl==1 && ntrutaush==1 && mt>140");
  auto proc_tt2l = Process::MakeShared<Baby_full>("t#bar{t} 2l, m_{T}>140", Process::Type::background, 
						  colors("tt_2l"), ttfiles, 
						  baseline && "ntruels+ntrumus+ntrutausl==2 && mt>140");
  auto proc_tt1l_ghimt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}>140, m_{T}^{tru}>140", Process::Type::background, 
						       kGreen-3, ttfiles, 
						       baseline && "ntruleps<=1 && mt>140&&mt_tru>140");
  auto proc_tt1l_bhimt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}>140, m_{T}^{tru}#leq140", Process::Type::background, 
						       kRed-4, ttfiles, 
						       baseline && "ntruleps<=1 && mt>140&&mt_tru<140");


  auto proc_1l = Process::MakeShared<Baby_full>("t#bar{t} 1l", Process::Type::background, 
                   colors("tt_1l"), {foldermc+"*_TTJets*SingleLept*"+ntupletag+"*.root"}, 
                   baseline && "mt<140");
  auto proc_2l = Process::MakeShared<Baby_full>("t#bar{t} 2l", Process::Type::data, 
              colors("tt_2l"), {foldermc+"*_TTJets*DiLept*"+ntupletag+"*.root"}, 
              baseline && "mt>140");

  vector<shared_ptr<Process> > tt_procs = {proc_tt1l_lomt, proc_tt2l, proc_ttltau, proc_tt1l_ghimt, proc_tt1l_bhimt};
  vector<shared_ptr<Process> > tt1l_procs = {proc_tt1l_lomt, proc_tt1l_ghimt, proc_tt1l_bhimt};
  vector<shared_ptr<Process> > ttnb_procs = {proc_1l, proc_2l};
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::data_norm)
  .RatioMaximum(1.5)
  .RatioMinimum(0.5);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {lin_shapes.PrintVals(true)};

  NamedFunc wgt = "weight";

  PlotMaker pm;
  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets==5 && met>200 && met<350 && mj14>250 && mj14<400", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");
  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets==5 && met>200 && met<350 && mj14>400 && mj14<500", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");
  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets==5 && met>200 && met<350 && mj14>500", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");

  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets>=7 && met>200 && met<350 && mj14>250 && mj14<400", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");
  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets>=7 && met>200 && met<350 && mj14>400 && mj14<500", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");
  pm.Push<Hist1D>(Axis(4,-0.5,3.5,"nbd","N_{b}"),"njets>=7 && met>200 && met<350 && mj14>500", ttnb_procs, plot_types).Weight(wgt).Tag("ttnb");


  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

