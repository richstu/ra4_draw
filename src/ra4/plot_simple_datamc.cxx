#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string.h>

#include "TMath.h"
#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

int main() {
  gErrorIgnoreLevel = 6000;

  Palette colors("txt/colors.txt", "default");
  Process::Type back = Process::Type::background;
  Process::Type data = Process::Type::data;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  map<int, string> foldermc, folderdata, foldersig;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
  foldersig[2016] = "";//bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
  foldersig[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/T1tttt/unskimmed/";
  folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_01_18/mc/merged_mcbase_stdnj5/";
  foldersig[2018] = "";//bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_01_18/data/merged_database_standard/";
  

  NamedFunc baseline = "mj14>400 && mj14<=500 && met>200 && met<=350  && mj14>250 && st>500";
  baseline = baseline && Functions::hem_veto && "st<10000 && pass && pass_ra2_badmu && met/met_calo<5";

  NamedFunc trigs = Functions::trig_run2;

  NamedFunc data_base = baseline && trigs;
  NamedFunc mc_base = baseline && "stitch_met";

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {folderdata[2016]+"*.root"},data_base);
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {folderdata[2017]+"*.root"},data_base);
	auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
	                 {folderdata[2018]+"*.root"},data_base);

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {foldermc[2016]+"*_TTJets*SingleLept*.root"}, mc_base);
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {foldermc[2016]+"*_TTJets*DiLept*.root"}, mc_base);
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {foldermc[2016]+"*_WJetsToLNu_HT*.root"}, mc_base);
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {foldermc[2016]+"*_ST_*.root"},mc_base);
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {foldermc[2016]+"*_TTWJets*.root", foldermc[2016]+"*_TTZ*.root", foldermc[2016]+"*_TTGJets*.root"}, mc_base);
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
	                     {foldermc[2016]+"*QCD_HT*0_Tune*.root", foldermc[2016]+"*QCD_HT*Inf_Tune*.root",
                        foldermc[2016]+"*_DYJetsToLL_M-50_HT*.root", 
                        foldermc[2016]+"*_ZJet*.root", foldermc[2016]+"*_ttHTobb_M125_*.root",
                        foldermc[2016]+"*_TTTT_*.root",
                        foldermc[2016]+"*_WH_HToBB*.root", foldermc[2016]+"*_ZH_HToBB*.root", 
                        foldermc[2016]+"*_WWTo*.root", foldermc[2016]+"*_WZ*.root",
                        foldermc[2016]+"_ZZ_*.root"}, mc_base);

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {foldermc[2017]+"*_TTJets*SingleLept*.root"}, mc_base);
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {foldermc[2017]+"*_TTJets*DiLept*.root"}, mc_base);
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {foldermc[2017]+"*_WJetsToLNu_HT*.root"}, mc_base);
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {foldermc[2017]+"*_ST_*.root"}, mc_base);
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {foldermc[2017]+"*_TTWJets*.root", foldermc[2017]+"*_TTZ*.root", foldermc[2017]+"*_TTGJets*.root"}, mc_base);
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
	                     {foldermc[2017]+"*QCD_HT*0_Tune*.root", foldermc[2017]+"*QCD_HT*Inf_Tune*.root",
                        foldermc[2017]+"*DYJetsToLL_M-50_HT*.root", 
                        foldermc[2017]+"*_ZJet*.root",              foldermc[2017]+"*_ttHTobb_M125_*.root",
                        foldermc[2017]+"*_TTTT_*.root",
                        foldermc[2017]+"*_WH_HToBB*.root",          foldermc[2017]+"*_ZH_HToBB*.root", 
                        foldermc[2017]+"*_WWTo*.root",           
                        foldermc[2017]+"*_WZ*.root",
                        foldermc[2017]+"_ZZ_*.root"}, mc_base);

	vector<shared_ptr<Process> > data16_mc16  = {data_2016, mc16_tt1l, mc16_tt2l, mc16_wjets, mc16_single_t, mc16_ttv, mc16_other};
	vector<shared_ptr<Process> > data17_mc17  = {data_2017, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<shared_ptr<Process> > data18_mc17  = {data_2018, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<vector<shared_ptr<Process> >> data_mc = {data16_mc16, data17_mc17, data18_mc17};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
	log_lumi.Title(TitleType::info)
					.YAxis(YAxisType::log)
					.Stack(StackType::data_norm)
	        .Bottom(BottomType::ratio)
					.FileExtensions({"pdf"});
  PlotOpt lin_stack_info = log_lumi().YAxis(YAxisType::linear); 
	vector<PlotOpt> lin_stack = {lin_stack_info};
	vector<PlotOpt> log_stack = {log_lumi};

	vector<string> tag = {"2016", "2017", "2018"};

  vector<string> cuts;
  cuts.push_back("nleps==1 && nveto==0 && mt>140  && njets==5");
  cuts.push_back("nleps==1 && nveto==0 && mt>140  && njets==6");
  cuts.push_back("nleps==1 && nveto==0 && mt<=140 && njets==5");
  cuts.push_back("nleps==1 && nveto==0 && mt<=140 && njets==6");
  cuts.push_back("nleps==1 && nveto==0 && mt<=140 && njets>=7");
  cuts.push_back("nleps==2 &&                        njets>=5 && njets<=6");
  cuts.push_back("nleps==2 &&                        njets>=7");

  PlotMaker pm;
  for (size_t ipr(0); ipr<data_mc.size(); ipr++) {
    for (size_t ic(0); ic<cuts.size(); ic++) {
      pm.Push<Hist1D>(Axis(5,-0.5,  4.5, "nbd",  "N_{b, Deep}",{}), cuts[ic],data_mc[ipr], lin_stack)
      .Weight(Functions::wgt_run2 * Functions::eff_trig_run2).Tag(tag[ipr]);
    }
  }
  pm.min_print_=true;
  pm.MakePlots(1);
}

