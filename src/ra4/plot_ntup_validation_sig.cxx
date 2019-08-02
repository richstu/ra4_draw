#include <iostream>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
  double lumi = 137;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string old_file = bfolder+"/cms2r0/babymaker/babies/2019_05_16/T1tttt/unskimmed/*-2100_mLSP-100_*.root"; 
  string new_file = bfolder+"/cms2r0/babymaker/babies/2019_07_16/T1tttt/unskimmed/*-2100_mLSP-100_*.root";

  NamedFunc filters = Functions::hem_veto;

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("New", Process::Type::signal, kRed+1, {new_file}, filters));
  procs.back()->SetLineStyle(2);
  procs.push_back(Process::MakeShared<Baby_full>("Old", Process::Type::signal, kBlack, {old_file}, filters));


  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info).YAxis(YAxisType::linear).Stack(StackType::lumi_shapes);
  // .Bottom(BottomType::ratio).YAxis(YAxisType::linear).Stack(StackType::lumi_shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> all_plot_types = {lin_lumi};

  PlotMaker pm;

  string baseline = "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbdm>=1";
  vector<NamedFunc> weights = {"w_lumi"};
  for(unsigned iw=0; iw<weights.size(); iw++){
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,500, 800}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {6.5, 8.5}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nleps", "N_{leps}", {5.5, 8.5}),
        "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "nbdm", "N_{b}", {1.5, 2.5,3.5}),
        "1", procs, all_plot_types).Weight(weights[iw]);

    pm.Push<Hist1D>(Axis(60, 0., 3., "w_btag_deep", "b-tag weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_lep", "lepton weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_fs_lep", "lepton weight (FastSim)", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_isr", "ISR weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(60, 0., 3., "w_pu", "Pileup weight", {}),
       "1", procs, all_plot_types).Weight(weights[iw]);
  }
   
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
