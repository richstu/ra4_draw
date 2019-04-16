#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  double lumi = 1;
  bool paper = true;
}

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");
  
  string lsp = "{#lower[-0.1]{#tilde{#chi}}#lower[0.2]{#scale[0.95]{^{0}}}#kern[-1.3]{#scale[0.95]{_{1}}}}";
  string t1t_label = "#scale[0.95]{#tilde{g}#kern[0.2]{#tilde{g}}, #tilde{g}#rightarrowt#kern[0.18]{#bar{t}}#kern[0.18]"+lsp;
  string mt = "m#lower[-0.1]{_{T}}";
  string etmiss = "p#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.15]{#lower[0.15]{^{miss}}}}";
  string nj8 = "N#kern[-0.1]{#scale[1.15]{#lower[-0.15]{_{jets}}}}#kern[0.2]{#geq}#kern[0.2]{8}";
  string nj7 = "N#kern[-0.1]{#scale[1.15]{#lower[-0.15]{_{jets}}}}#kern[0.2]{#geq}#kern[0.2]{7}";
  string nj6 = "N#kern[-0.1]{#scale[1.15]{#lower[-0.15]{_{jets}}}}#kern[0.2]{#geq}#kern[0.2]{6}";
  string nb3 = "N#kern[-0.15]{#scale[1.15]{#lower[-0.15]{_{b}}}}#kern[0.2]{#geq}#kern[0.2]{3}";
  string nb2 = "N#kern[-0.15]{#scale[1.15]{#lower[-0.15]{_{b}}}}#kern[0.2]{#geq}#kern[0.2]{2}";

  TString baseline("st>500 && met>200 && njets>=6 && nbdm>=1 && nleps==1 && nveto==0");
  NamedFunc filters = Functions::hem_veto && Functions::pass_run2;
  NamedFunc nom_wgt = Functions::wgt_run2 * Functions::eff_trig_run2; 

  set<int> years = {2016, 2017, 2018};

  map<int, string> foldermc, folderdata, foldersig;
  foldersig[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/T1tttt/skim_sys_abcd/";
  folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

  foldersig[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/T1tttt/skim_sys_abcd/";
  folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

  foldersig[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/T1tttt/skim_sys_abcd/";
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/merged_database_standard/";

  set<string> data_files, sig_nc_files, sig_c_files;
  for (auto &yr: years) {
    data_files.insert(folderdata[yr]+"*root");
    sig_nc_files.insert(foldersig[yr]+"*mGluino-2100_mLSP-100_*.root");
    sig_c_files.insert(foldersig[yr]+"*mGluino-1900_mLSP-1250_*.root");
  }

  auto data_highmt = Process::MakeShared<Baby_full>(" Data, "+mt+" > 140 GeV", Process::Type::data, kBlack,
    data_files, baseline && filters && Functions::trig_run2 && "mt>140");
  auto data_lowmt = Process::MakeShared<Baby_full>(" Data, "+mt+" #leq 140 GeV", Process::Type::background, kAzure+10,
    data_files, baseline && filters && Functions::trig_run2 && "mt<=140");
  data_lowmt->SetFillColor(kWhite);
  data_lowmt->SetLineColor(kAzure+10);
  data_lowmt->SetLineWidth(2);

  auto t1tttt = Process::MakeShared<Baby_full>(" "+t1t_label+" (2100,100)}", Process::Type::signal, colors("t1tttt"),
    sig_nc_files, baseline && filters && "mt>140");
  t1tttt->SetLineWidth(2);
  auto t1ttttc = Process::MakeShared<Baby_full>(" "+t1t_label+" (1900,1250)}", Process::Type::signal, colors("t1tttt"),
    sig_c_files, baseline && filters && "mt>140");
  t1ttttc->SetLineWidth(2);
  t1ttttc->SetLineStyle(2);

  vector<shared_ptr<Process> > data1l_procs = {data_highmt,data_lowmt,t1tttt,t1ttttc};

  string style = "Preliminary";
  // if(paper) style = "PRLPaper";
  PlotOpt log_lumi("txt/plot_styles.txt", style);
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .RatioMaximum(1.86);
  if(!paper){
    log_lumi=log_lumi.Title(TitleType::preliminary);
  }
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  if(paper){
    log_shapes = log_lumi().Stack(StackType::shapes)
      .Bottom(BottomType::off)
      .ShowBackgroundError(false);
  }
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> log = {log_lumi_info,log_lumi};
  vector<PlotOpt> lin = {lin_lumi};


  PlotMaker pm;
  vector<string> metbins = { "met>200 && met<=350", "met>350 && met<=500", "met>500"};
  vector<set<double>> mj_lines = {{250, 400, 500}, {250, 450, 650}, {250, 500, 800}};
  for (unsigned imet(0); imet<metbins.size(); imet++){
    if (paper) {
      string metlabel = CodeToRootTex(metbins[imet]);
      ReplaceAll(metlabel,"E_{T}^{miss}",etmiss);
      if (metbins[imet]=="met>500") {
        pm.Push<Hist1D>(Axis(23, 50.,1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
          metbins[imet] + "&&nbdm>=2&&njets>=6", data1l_procs, lin).Tag("2b_6j").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
        .RightLabel({metlabel+" GeV, "+nj6+", "+nb2}).YAxisZoom(0.85).Weight(nom_wgt);        
      } else {
        // pm.Push<Hist1D>(Axis(23, 50.,1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
        //   metbins[imet] + "&&nbdm>=3&&njets==7", data1l_procs, lin).Tag("3b_7j").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
        // .RightLabel({metlabel+" GeV, "+nj7+", "+nb3}).YAxisZoom(0.85).Weight(nom_wgt);
         pm.Push<Hist1D>(Axis(23, 50.,1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
          metbins[imet] + "&&nbdm>=2&&njets>=7", data1l_procs, lin).Tag("2b_7j").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
        .RightLabel({metlabel+" GeV, "+nj7+", "+nb2}).YAxisZoom(0.85).Weight(nom_wgt);
        // pm.Push<Hist1D>(Axis(23, 50.,1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
        //   metbins[imet] + "&&nbdm>=3&&njets>=8", data1l_procs, lin).Tag("3b_8j").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
        // .RightLabel({metlabel+" GeV, "+nj8+", "+nb3}).YAxisZoom(0.85).Weight(nom_wgt);
      }
      pm.Push<Hist1D>(Axis(23, 50.,1200., "mj14", "M_{J} [GeV]", mj_lines[imet]),
        metbins[imet] + "&&nbdm>=2&&njets>=8", data1l_procs, lin).Tag("2b_8j").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
      .RightLabel({metlabel+" GeV, "+nj8+", "+nb2}).YAxisZoom(0.85).Weight(nom_wgt);

    } else {
      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]",{250.,400.}),
		      metbins[imet] + "&&nbdm==1", data1l_procs, lin).Tag("data1l1b").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
      .RightLabel({CodeToRootTex(metbins[imet])+" GeV",CodeToRootTex("njets>=6&&nbdm==1")}).YAxisZoom(0.93).Weight(nom_wgt);
      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]",{250.,400.}),
		      metbins[imet] + "&&nbdm>=2", data1l_procs, lin).Tag("data1l2b").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
      .RightLabel({CodeToRootTex(metbins[imet])+" GeV",CodeToRootTex("njets>=6&&nbdm>=2")}).YAxisZoom(0.93).Weight(nom_wgt);
    }
  } 


  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
