#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
  double lumi = 35.9;

  Palette colors("txt/colors.txt", "default");

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string old_T1tttt_dir = bfolder+"/cms2r0/babymaker/babies/2019_01_11/T1tttt/unskimmed/"; 
  string new_T1tttt_dir = bfolder+"/cms2r0/babymaker/babies/2018_12_17/T1tttt/unskimmed/";
  string old_tag = "2016";
  string new_tag = "2017";

  NamedFunc filters = Functions::hem_veto;

  auto new_NC = Process::MakeShared<Baby_full>("(2100,100) "+new_tag, Process::Type::signal, colors("t1tttt"),
    {new_T1tttt_dir+"*T1tttt_mGluino-2100_mLSP-100_*.root"}, filters);

  auto old_NC = Process::MakeShared<Baby_full>("(2100,100) "+old_tag, Process::Type::signal, colors("t1tttt"),
    {old_T1tttt_dir+"*T1tttt_mGluino-2100_mLSP-100_*.root"}, filters);
  old_NC->SetLineStyle(2);

  auto new_C = Process::MakeShared<Baby_full>("(1900,1250) "+new_tag, Process::Type::signal, kAzure+2,
    {new_T1tttt_dir+"*T1tttt_mGluino-1900_mLSP-1250_*.root"}, filters);

  auto old_C = Process::MakeShared<Baby_full>("(1900,1250) "+old_tag, Process::Type::signal, kAzure+2,
    {old_T1tttt_dir+"*T1tttt_mGluino-1900_mLSP-1250_*.root"}, filters);
  old_C->SetLineStyle(2);


  vector<shared_ptr<Process> > procs = {new_NC,old_NC,new_C,old_C};
  
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    //    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::lumi_shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info};//, lin_lumi_info};

  PlotMaker pm;

  string baseline = "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1";
  vector<NamedFunc> weights = {"weight"};
  for(unsigned iw=0; iw<weights.size(); iw++){
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
        "nleps==1&&nveto==0&&met>200&&njets>=6&&nbm>=1",
       procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
        baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
        "nleps==1&&nveto==0&&st>500&&met>100&&njets>=6&&nbd>=1",
       procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
        baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
        "nleps==1&&nveto==0&&st>500&&met>200&&nbd>=1",
       procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    
    pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
        "st>500&&met>200&&njets>=6&&nbd>=1",
       procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
        "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbd", "N_{b}", {-999}),
    //     "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6",
    //    procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbd", "N_{b}", {-999}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    
    // pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
    //     baseline, procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");

    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
    //     "nleps==1&&st>500&&met>200&&njets>=6&&nbd>=1",
    //    procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
    //     "1", procs, all_plot_types).Weight(weights[iw]).Tag("sig");
    
  }
   
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
