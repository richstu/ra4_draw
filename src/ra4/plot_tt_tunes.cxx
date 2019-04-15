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


NamedFunc max_b_pt("max_b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=5) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });

NamedFunc max_t_pt("max_t_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=6) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });


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
  string folder_mc16old(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/");
  string folder_mc16new(bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/");
  string folder_mc17(bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/");

  Palette colors("txt/colors.txt", "default");
  NamedFunc baseline = "stitch_met && nleps==1 && nveto==0 && st>500 && met>100 && njets>=5 && nbdm>=1 && weight<1";

  auto proc_tt1l_mc16old = Process::MakeShared<Baby_full>("2016 MC bear", Process::Type::background, 
                   kBlue, {folder_mc16old+"*_TTJets*Lept*.root"}, baseline);
  auto proc_tt1l_mc16new = Process::MakeShared<Baby_full>("2016 MC quokka", Process::Type::background, 
                   kGreen+1, {folder_mc16new+"*_TTJets*Lept*.root"}, baseline);
  auto proc_tt1l_mc17 = Process::MakeShared<Baby_full>("2017 MC", Process::Type::background, 
						       kRed, {folder_mc17+"*_TTJets*Lept*.root"}, baseline);

  // vector<shared_ptr<Process> > tt1l_procs = {proc_tt1l_mc17, proc_tt1l_mc16old};
  vector<shared_ptr<Process> > tt1l_procs = {proc_tt1l_mc17, proc_tt1l_mc16new, proc_tt1l_mc16old};
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::shapes);
  // .RatioMaximum(2.);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {lin_shapes};

  PlotMaker pm;


    pm.Push<Hist1D>(Axis(16,0.,800.,"mj14","M_{J} [GeV]",{250., 400.}),"1", tt1l_procs, plot_types);
    pm.Push<Hist1D>(Axis(16,0.,800.,"mj14","M_{J} [GeV]",{250., 400.}),"ntruleps<=1 && mt>140&&mt_tru<140", tt1l_procs, plot_types);

    // pm.Push<Hist1D>(Axis(16,200.,1000.,"met","MET [GeV]",{200}),"1", tt1l_procs, plot_types).Weight(iweight);

    // pm.Push<Hist1D>(Axis(12,0.5,12.5,"njets","N_{jets} [GeV]",{5.5}),"1", tt1l_procs, plot_types).Weight(iweight);
    // pm.Push<Hist1D>(Axis(16,500.,1300.,"ht","H_{T} [GeV]"),"1", tt1l_procs, plot_types).Weight(iweight);
    // pm.Push<Hist1D>(Axis(16,0.,800.,"isr_tru_pt","t#bar{t} true p_{T} [GeV]"),"1", tt1l_procs, plot_types).Weight(iweight);
    // pm.Push<Hist1D>(Axis(16,0.,800.,max_t_pt,"Max t-quark p_{T} [GeV]"),"1", tt1l_procs, plot_types).Weight(iweight);

  
  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

