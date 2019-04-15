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

NamedFunc BaselineCuts(string var = "", string extra = "1") {
	NamedFunc cuts[6] = {"nleps == 1", "st > 500", "met > 200", 
	                     "njets >= 6", "nbdm >= 1", "nveto == 0"};
  int num_cuts(6);
	if(extra == "2l") {
		cuts[0] = "nleps == 2";
		cuts[2] = "met > 200 && met < 500";
		cuts[3] = "njets >= 5";
		cuts[4] = "nbdm <= 2";
		num_cuts = 5;
	}
	else if(extra == "5j")
		cuts[3] = "njets == 5";
	int out(-1);
	if     (var == "nleps") cuts[0] = "nleps >= 1";
	else if(var == "met")   cuts[2] = "met > 100";
	else if(var == "njets") cuts[3] = "njets >= 4";
	else if(var == "nbm")   out = 4;
	NamedFunc baseline(cuts[0]);
	for(int i = 1; i < num_cuts; i++) 
		if(i != out) baseline = baseline && cuts[i];
	if(extra == "1" && var != "mt" && var != "mj") baseline = baseline && "(mj14 < 400 || mt < 140)";
	if(extra != "1" && extra != "2l" && extra != "5j") {
	  NamedFunc temp = extra;
	  baseline = baseline && temp;
	}
	return baseline;
}

int main() {
  gErrorIgnoreLevel = 6000;

  Palette colors("txt/colors.txt", "default");
	Process::Type back = Process::Type::background;
	Process::Type data = Process::Type::data;

  // Data
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/");
  string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/merged_database_standard/");
	string q_cuts("pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {data16_path+"*.root"},"trig_ra4&&"+q_cuts);
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {data17_path+"*.root"},"trig_ra4&&"+q_cuts);

	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_standard/");
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc16_path+"*_TTJets*SingleLept*.root"}, "stitch_met&&"+q_cuts);
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc16_path+"*_TTJets*DiLept*.root"}, "stitch_met&&"+q_cuts);
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {mc16_path+"*_WJetsToLNu*.root"}, "stitch_met&&"+q_cuts);
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {mc16_path+"*_ST_*.root"},"stitch_met&&"+q_cuts);
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root", mc16_path+"*_TTGJets*.root"}, "stitch_met&&"+q_cuts);
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
	                     {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
                        mc16_path+"*DYJetsToLL*.root", 
                        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
                        mc16_path+"*_TTTT_*.root",
                        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
                        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
                        mc16_path+"_ZZ_*.root"}, "stitch_met&&"+q_cuts);
												

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc17_path+"*_TTJets*SingleLept*.root"}, "stitch_met&&"+q_cuts);
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc17_path+"*_TTJets*DiLept*.root"}, "stitch_met&&"+q_cuts);
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc17_path+"*_WJetsToLNu_*.root"},
											 "pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000 && ((stitch_met && type!=2000) || (type==2000 && ht_isr_me<100))");
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc17_path+"*_ST_*.root"}, "stitch_met&&"+q_cuts);
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root"}, "stitch_met&&"+q_cuts);
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
	                     {mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
                        mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                        mc17_path+"*_TTTT_*.root",
                        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                        mc17_path+"*_WWTo*.root",           
                        mc17_path+"*_WZ*.root",
                        mc17_path+"_ZZ_*.root"}, "stitch_met&&"+q_cuts);

	vector<shared_ptr<Process> > data16_mc16  = {data_2016, mc16_tt1l, mc16_tt2l, mc16_wjets, mc16_single_t, mc16_ttv, mc16_other};
	vector<shared_ptr<Process> > data17_mc17  = {data_2017, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
	log_lumi.Title(TitleType::info)
					.YAxis(YAxisType::log)
					.Stack(StackType::data_norm)
	        .Bottom(BottomType::ratio)
					.FileExtensions({"pdf"});
	vector<PlotOpt> log_stack = {log_lumi};
	string tag("2016");
	string baseline("met > 200 && st > 500 && njets >= 5 && mt < 140");
  PlotMaker pm16;
	// Without w_prefire
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "jets_eta",  "Jet #eta",     {}), "nleps==1&&"+baseline, data16_mc16, log_stack).Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "leps_eta",  "Lepton #eta",  {}), "nleps==1&&"+baseline, data16_mc16, log_stack).Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "els_sceta", "electron #eta",{}), "nels==1&&" +baseline, data16_mc16, log_stack).Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "mus_eta",   "muon #eta",    {}), "nmus==1&&" +baseline, data16_mc16, log_stack).Tag(tag);
  // With w_prefire
	tag += "_prefire";
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "jets_eta",  "Jet #eta",     {}), "nleps==1&&"+baseline, data16_mc16, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "leps_eta",  "Lepton #eta",  {}), "nleps==1&&"+baseline, data16_mc16, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "els_sceta", "electron #eta",{}), "nels==1&&" +baseline, data16_mc16, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm16.Push<Hist1D>(Axis(25,-2.5,2.5, "mus_eta",   "muon #eta",    {}), "nmus==1&&" +baseline, data16_mc16, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm16.min_print_=true;
	pm16.MakePlots(35.9);

  PlotMaker pm17;
	// Without w_prefire
	tag = "2017";
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "jets_eta",  "Jet #eta",      {}), "nleps==1&&"+baseline, data17_mc17, log_stack).Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "leps_eta",  "Lepton #eta",   {}), "nleps==1&&"+baseline, data17_mc17, log_stack).Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "els_sceta", "electron #eta", {}), "nels==1&&" +baseline, data17_mc17, log_stack).Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "mus_eta",   "muon #eta",     {}), "nmus==1&&" +baseline, data17_mc17, log_stack).Tag(tag);
  // With w_prefire
	tag += "_prefire";
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "jets_eta",  "Jet #eta",      {}), "nleps==1&&"+baseline, data17_mc17, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "leps_eta",  "Lepton #eta",   {}), "nleps==1&&"+baseline, data17_mc17, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "els_sceta", "electron #eta", {}), "nels==1&&" +baseline, data17_mc17, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm17.Push<Hist1D>(Axis(25,-2.5,2.5, "mus_eta",   "muon #eta",     {}), "nmus==1&&" +baseline, data17_mc17, log_stack).Weight("weight*w_prefire").Tag(tag);
  pm17.min_print_=true;
	pm17.MakePlots(41.5);

}

