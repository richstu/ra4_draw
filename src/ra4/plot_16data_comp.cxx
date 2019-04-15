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
	Process::Type data = Process::Type::data;
	Process::Type back = Process::Type::background;

  // Data
  string data16_path_1el("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_ne1nj3/");
  string data16_path_1mu("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_nm1nj3/");
  string data16_path_unskimmed("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/unskimmed/");
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/");
  string data16_old_path("/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_standard/");
  string data16_old_path_nl1("/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_nl1nj3/");
  string data16_old_path_unskimmed("/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/unskimmed/");
	string q_cuts("(trig[19] || trig[20] || trig[21] || trig[24] || trig[40] || trig[41]) && pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {data16_path+"*.root"},q_cuts);
	auto data_2016_nl1nj3 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                       {data16_path_1el+"*.root",data16_path_1mu+"*.root"},q_cuts);
	auto data_2016_unskimmed = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                       {data16_path_unskimmed+"*.root"},q_cuts);
	auto data_2016_old = Process::MakeShared<Baby_full>("Old 2016 Data",back,kAzure,  
	                     {data16_old_path+"*.root"},q_cuts);
	auto data_2016_old_nl1nj3 = Process::MakeShared<Baby_full>("Old 2016 Data",back,kAzure,  
	                           {data16_old_path_nl1+"*.root"},q_cuts);
	auto data_2016_old_unskimmed = Process::MakeShared<Baby_full>("Old 2016 Data",back,kAzure,  
	                               {data16_old_path_unskimmed+"*.root"},q_cuts);


	vector<shared_ptr<Process> > data_16  = {data_2016_old, data_2016};
	vector<shared_ptr<Process> > data_16_nl1nj3  = {data_2016_old_nl1nj3, data_2016_nl1nj3};
	vector<shared_ptr<Process> > data_16_unskimmed  = {data_2016_old_unskimmed, data_2016_unskimmed};

	
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
	log_lumi.Title(TitleType::info)
					.YAxis(YAxisType::log)
					.FileExtensions({"pdf"});
  PlotOpt log_shapes_info = log_lumi().Stack(StackType::shapes)
	                                    .Bottom(BottomType::ratio)
// 																			.RatioMaximum(1.2)
// 																			.RatioMinimum(0.8)
																			.Title(TitleType::info);  
	vector<PlotOpt> log_shape = {log_shapes_info};
  PlotOpt lin_shapes_info = log_shapes_info().YAxis(YAxisType::linear);
	vector<PlotOpt> lin_shape = {lin_shapes_info};
	vector<PlotOpt> log_stack = {log_shapes_info.Stack(StackType::data_norm)};
  PlotMaker pm;
	TString baseline("nleps == 1 && leps_pt[0] > 50 && st > 500 && met > 200");
	string bef("16_ReReco_check");
	string cuts("");
  pm.Push<Hist1D>(Axis(30,0,600,    "met", "p_{T}^{miss} [GeV]",  {}),
									"nleps == 1 && leps_pt[0] > 50 && njets >= 3",  data_16_nl1nj3, log_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(16,200,1000,    "st", "S_{T} [GeV]",  {}),
									"nleps == 1 && leps_pt[0] > 50 && njets >= 3 && met > 100",  data_16_nl1nj3, log_shape).Tag(bef);
	pm.Push<Hist1D>(Axis(30,  0, 1200,  "mj14", "M_{J} [GeV]",      {}), baseline, data_16, log_shape).Tag(bef);
	pm.Push<Hist1D>(Axis(15,  0,  300,  "mt",   "m_{T} [GeV]",      {}), baseline, data_16, log_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(8,2.5,10.5,    "njets","N_{jets}",         {}), baseline, data_16, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(5,-0.5,  4.5,  "nbdm",  "N_{b, Deep}",      {}), baseline, data_16, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(16,0,800,  "jets_pt",  "Jet p_{T} [GeV]",  {}), baseline, data_16, log_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(20,0,1,  "jets_csvd",  "Jet Deep CSV",     {}), baseline, data_16, log_shape).Tag(bef);
  pm.min_print_=true;
  pm.MakePlots(1);
}

