#include <iostream>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/palette.hpp"
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

  Palette colors("txt/colors.txt", "default");

  string tt_dir = bfolder+"/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_stdnj5/";  
  string t2tt_dir = bfolder+"/cms2r0/babymaker/babies/2019_07_18/T2tt/unskimmed/";
  string t54t_dir = bfolder+"/cms2r0/babymaker/babies/2019_07_18/T5tttt/unskimmed/";
  // string t14t_dir = bfolder+"/cms2r0/babymaker/babies/2019_07_18/T1tttt/unskimmed/";

  string baseline = "stitch_met && nleps==1 && st>500 && met>100 && njets>=5";

  vector<shared_ptr<Process> > procs;
  // procs.push_back(Process::MakeShared<Baby_full>("T1tttt(1650,1) 1L", Process::Type::signal, kGreen+3, 
  //                 {t14t_dir+"*-1650_mLSP-1_*.root"}, baseline + " && ntruleps<=1"));
  // procs.push_back(Process::MakeShared<Baby_full>("T1tttt(1650,1) 2L", Process::Type::signal, kGreen+1, 
  //                 {t14t_dir+"*-1650_mLSP-1_*.root"}, baseline + " && ntruleps>=2"));
  // procs.push_back(Process::MakeShared<Baby_full>("T5tttt(1650,1) 1L", Process::Type::signal, kGreen+3, 
  //                 {t54t_dir+"*-1650_mLSP-1_*.root"}, baseline + " && ntruleps<=1"));
  // procs.push_back(Process::MakeShared<Baby_full>("T5tttt(1650,1) 2L", Process::Type::signal, kGreen+1, 
  //                 {t54t_dir+"*-1650_mLSP-1_*.root"}, baseline + " && ntruleps>=2"));

  // procs.push_back(Process::MakeShared<Baby_full>("T2tt(175,1) 1L", Process::Type::signal, kRed+3, 
  //                 {t2tt_dir+"*mGluino-175_mLSP-1_*"}, baseline + " && ntruleps<=1"));
  // procs.back()->SetLineStyle(2);
  // procs.push_back(Process::MakeShared<Baby_full>("T2tt(175,1) 2L", Process::Type::signal, kRed+1, 
  //                 {t2tt_dir+"*mGluino-175_mLSP-1_*"}, baseline + " && ntruleps>=2"));
  // procs.back()->SetLineStyle(2);
  // procs.push_back(Process::MakeShared<Baby_full>("T2tt(200,25) 1L", Process::Type::signal, kGreen+3, 
  //                 {t2tt_dir+"*mGluino-200_mLSP-25_*"}, baseline + " && ntruleps<=1"));
  // procs.back()->SetLineStyle(2);
  // procs.push_back(Process::MakeShared<Baby_full>("T2tt(200,25) 2L", Process::Type::signal, kGreen+1, 
  //                 {t2tt_dir+"*mGluino-200_mLSP-25_*"}, baseline + " && ntruleps>=2"));
  // procs.back()->SetLineStyle(2);
  // procs.push_back(Process::MakeShared<Baby_full>("t#bar{t} 1L", Process::Type::signal, colors("tt_1l"), 
  //                 {tt1l_file}, baseline));
  // procs.push_back(Process::MakeShared<Baby_full>("t#bar{t} 2L", Process::Type::signal, colors("tt_2l"), 
  //                 {tt2l_file}, baseline));

  procs.push_back(Process::MakeShared<Baby_full>("T2tt(350,175)", Process::Type::signal, kRed, 
                  {t2tt_dir+"*mGluino-350_mLSP-175_*"}, baseline));
  procs.back()->SetLineStyle(2);
  procs.push_back(Process::MakeShared<Baby_full>("T2tt(375,200)", Process::Type::signal, kRed+2, 
                  {t2tt_dir+"*mGluino-375_mLSP-200_*"}, baseline));
  procs.back()->SetLineStyle(3);
  procs.push_back(Process::MakeShared<Baby_full>("T2tt(400,225)", Process::Type::signal, kRed+4, 
                  {t2tt_dir+"*mGluino-400_mLSP-225_*"}, baseline));
  procs.back()->SetLineStyle(3);
  // procs.push_back(Process::MakeShared<Baby_full>("T2tt(250,75)", Process::Type::signal, 1, 
  //                 {t2tt_dir+"*mGluino-250_mLSP-75_*"}, baseline));
  // procs.back()->SetLineStyle(3);
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::signal, colors("tt_1l"), 
                  {tt_dir+"*TTJets*SingleLept*.root", tt_dir+"*TTJets*DiLept*.root"}, baseline));
  // procs.push_back(Process::MakeShared<Baby_full>("T5tttt(1650,1)", Process::Type::signal, kGreen+1, 
  //                  {t54t_dir+"*-1650_mLSP-1_*.root"}, baseline));

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info).YAxis(YAxisType::linear);//.Stack(StackType::lumi_shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::log);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> all_plot_types = {lin_shapes};

  PlotMaker pm;

  vector<NamedFunc> weights = {"weight"};
  for(unsigned iw=0; iw<weights.size(); iw++){
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,500, 800}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(12, 200, 800., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(7, 4.5, 11.5, "njets", "N_{jets}", {6.5, 8.5}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "nbdm", "N_{b}", {1.5, 2.5,3.5}),
        baseline, procs, all_plot_types).Weight(weights[iw]);
  }
   
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
