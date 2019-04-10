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
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/");

  Palette colors("txt/colors.txt", "default");


  NamedFunc baseline_1l = "nleps==1 && nveto==0 && njets>=7";
  NamedFunc baseline_2l = "nleps==2 && njets>=6";

  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root"};
  
  vector<string> metcuts;
  // metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=300");
  metcuts.push_back("met>300 && met<=400");
  metcuts.push_back("met>400 && met<=500");
  metcuts.push_back("met>500");

  vector<shared_ptr<Process> > procs_1l, procs_2l;
  for (size_t imet(0); imet<metcuts.size(); imet++) {
    procs_1l.push_back(Process::MakeShared<Baby_full>(metcuts[imet], 
      Process::Type::background, kCyan+imet, ttfiles, 
      "stitch_met &&"+metcuts[imet]));

    procs_2l.push_back(Process::MakeShared<Baby_full>(metcuts[imet], 
      Process::Type::background, kCyan+imet, ttfiles, 
      "stitch_met &&"+metcuts[imet]));
  }
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::shapes)
  .RatioMaximum(1.99)
  .RatioMinimum(0.);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plots = {lin_shapes};//.PrintVals(true)};

  PlotMaker pm;

  pm.Push<Hist1D>(Axis(30, 0., 1200.,"mj14","M_{J} [GeV]"), baseline_1l, procs_1l, plots);
  pm.Push<Hist1D>(Axis(30, 0., 1200.,"mj14","M_{J} [GeV]"), baseline_2l, procs_2l, plots);

  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)
      <<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

