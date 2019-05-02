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
	                     "njets >= 7", "nbdm >= 1", "nveto == 0"};
  int num_cuts(6);
	if(extra == "2l") {
		cuts[0] = "nleps == 2";
		cuts[2] = "met > 200 && met < 500";
		cuts[3] = "njets >= 5";
		cuts[4] = "nbdm <= 2";
		num_cuts = 5;
	}
	else if(extra == "56j")
		cuts[3] = "njets >= 5 && njets <= 6";
	int out(-1);
	if     (var == "nleps") cuts[0] = "nleps >= 1";
	else if(var == "met")   cuts[2] = "met > 100";
	else if(var == "njets") out = 3;
// 	else if(var == "njets") cuts[3] = "njets >= 4";
	else if(var == "nbm")   out = 4;
	NamedFunc baseline(cuts[0]);
	for(int i = 1; i < num_cuts; i++) 
		if(i != out) baseline = baseline && cuts[i];
// 	if(extra == "1") { 
// 		if(var == "mt") baseline = baseline && "mj14 < 400";
// 		else if(var == "mt") baseline = baseline && "mt > 140";
// 		else baseline = baseline && "(mj14 < 400 || mt < 140)";
// 	}
	if(extra != "1" && extra != "2l" && extra != "56j") {
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

  const NamedFunc w_run2("w_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;
    double wgt = b.weight();
    if (b.SampleType()==2016)
      return wgt*b.w_prefire();
    else if (b.SampleType()==2017)
      return wgt*b.w_prefire()*Functions::wnpv2017(b);
    else 
      return wgt;
  });

  // Data
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/");
  string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/merged_database_standard/");
  string data18_path("/net/cms2/cms2r0/babymaker/babies/2019_03_30/data/merged_database_standard/");
// 	NamedFunc q_cuts("met/met_calo < 5 && pass_ra2_badmu && st < 10000");
	NamedFunc q_cuts(Functions::pass_run2);

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {data16_path+"*.root"},q_cuts && Functions::trig_run2);
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {data17_path+"*.root"},q_cuts && Functions::trig_run2);
	auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
	                 {data18_path+"*.root"},q_cuts && Functions::hem_veto && Functions::trig_run2);
	auto data_RunII = Process::MakeShared<Baby_full>("Run II Data",data,kBlack,
	                 {data16_path+"*.root",
	                  data17_path+"*.root",
	                  data18_path+"*.root"},q_cuts && Functions::hem_veto && Functions::trig_run2);

	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_standard/");
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");
  string mc18_path("/net/cms2/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_standard/");

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc16_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts);
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc16_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts);
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {mc16_path+"*_WJetsToLNu*.root"}, "stitch_met" && q_cuts);
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {mc16_path+"*_ST_*.root"},"stitch_met" && q_cuts);
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {mc16_path+"*_TTWJets*.root", 
											  mc16_path+"*_TTZ*.root", 
												mc16_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts);
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
	                     {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
                        mc16_path+"*DYJetsToLL*.root", 
                        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
                        mc16_path+"*_TTTT_*.root",
                        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
                        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
                        mc16_path+"_ZZ_*.root"}, "stitch_met" && q_cuts);
												

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc17_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts);
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc17_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts);
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc17_path+"*_WJetsToLNu_*.root"},
											 q_cuts && "((stitch_met && type!=2000) || (type==2000 && ht_isr_me<100))");
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc17_path+"*_ST_*.root"}, "stitch_met" && q_cuts);
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts);
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
	                     {mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
                        mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                        mc17_path+"*_TTTT_*.root",
                        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                        mc17_path+"*_WWTo*.root",           
                        mc17_path+"*_WZ*.root",
                        mc17_path+"_ZZ_*.root"}, "((stitch_met && type!=6000) || (type==6000 && ht_isr_me<100))" && q_cuts);


  auto mc18_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc18_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc18_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc18_path+"*_WJetsToLNu_*.root"},
				      					q_cuts && "(stitch_met && type!=2000) || (type==2000 && ht_isr_me<100)" && Functions::hem_veto);
  auto mc18_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc18_path+"*_ST_*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc18_path+"*_TTWJets*.root", mc18_path+"*_TTZ*.root", mc18_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
	                     {mc18_path+"*QCD_HT*0_Tune*.root", mc18_path+"*QCD_HT*Inf_Tune*.root",
                        mc18_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc18_path+"*_ZJet*.root",              mc18_path+"*_ttHTobb_M125_*.root",
                        mc18_path+"*_TTTT_*.root",
                        mc18_path+"*_WH_HToBB*.root",          mc18_path+"*_ZH_HToBB*.root", 
                        mc18_path+"*_WWTo*.root",           
                        mc18_path+"*_WZ*.root",
                        mc18_path+"_ZZ_*.root"}, "(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)" && q_cuts && Functions::hem_veto);


  auto mcRunII_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                        {mc16_path+"*_TTJets*SingleLept*.root",
	                         mc17_path+"*_TTJets*SingleLept*.root",
	                         mc18_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                        {mc16_path+"*_TTJets*DiLept*.root",
	                         mc17_path+"*_TTJets*DiLept*.root",
	                         mc18_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                        {mc16_path+"*_WJetsToLNu_*.root",
	                         mc17_path+"*_WJetsToLNu_*.root",
	                         mc18_path+"*_WJetsToLNu_*.root"},
				      					   q_cuts && "(stitch_met && type!=2000) || (type==2000 && ht_isr_me<100)" && Functions::hem_veto);
  auto mcRunII_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                        {mc16_path+"*_ST_*.root",
	                         mc17_path+"*_ST_*.root",
	                         mc18_path+"*_ST_*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                        {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root", mc16_path+"*_TTGJets*.root",
	                         mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root",
	                         mc18_path+"*_TTWJets*.root", mc18_path+"*_TTZ*.root", mc18_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
	                        {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
                           mc16_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc16_path+"*_ZJet*.root",              mc16_path+"*_ttHTobb_M125_*.root",
                           mc16_path+"*_TTTT_*.root",
                           mc16_path+"*_WH_HToBB*.root",          mc16_path+"*_ZH_HToBB*.root", 
                           mc16_path+"*_WWTo*.root",           
                           mc16_path+"*_WZ*.root",
                           mc16_path+"_ZZ_*.root",
	                         mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
                           mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                           mc17_path+"*_TTTT_*.root",
                           mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                           mc17_path+"*_WWTo*.root",           
                           mc17_path+"*_WZ*.root",
                           mc17_path+"_ZZ_*.root",
	                         mc18_path+"*QCD_HT*0_Tune*.root", mc18_path+"*QCD_HT*Inf_Tune*.root",
                           mc18_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc18_path+"*_ZJet*.root",              mc18_path+"*_ttHTobb_M125_*.root",
                           mc18_path+"*_TTTT_*.root",
                           mc18_path+"*_WH_HToBB*.root",          mc18_path+"*_ZH_HToBB*.root", 
                           mc18_path+"*_WWTo*.root",           
                           mc18_path+"*_WZ*.root",
                           mc18_path+"_ZZ_*.root"}, "(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)" && q_cuts && Functions::hem_veto);

	vector<shared_ptr<Process> > data16_mc16  = {data_2016, mc16_tt1l, mc16_tt2l, mc16_wjets, mc16_single_t, mc16_ttv, mc16_other};
	vector<shared_ptr<Process> > data17_mc17  = {data_2017, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<shared_ptr<Process> > data18_mc18  = {data_2018, mc18_tt1l, mc18_tt2l, mc18_wjets, mc18_single_t, mc18_ttv, mc18_other};
	vector<shared_ptr<Process> > data18_mc17  = {data_2018, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<shared_ptr<Process> > data_mc_RunII  = {data_RunII, mcRunII_tt1l, mcRunII_tt2l, mcRunII_wjets, mcRunII_single_t, mcRunII_ttv, mcRunII_other};
	vector<vector<shared_ptr<Process> >> data_mc = {data16_mc16, data17_mc17, data18_mc18, /*data18_mc17,*/ data_mc_RunII};

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

	vector<shared_ptr<Process> > temp;
	vector<string> sample_label = {"2016", "2017", "2018","RunII"};
	vector<double> sample_lumi = {35.9, 41.5, 60, 1};
	string tag;
  NamedFunc w_fr2(Functions::wgt_run2 * Functions::eff_trig_run2);
  w_fr2.Name("w_fr2");
  NamedFunc w_tot(w_run2 * Functions::eff_trig_run2);
	for(size_t i = 0; i < 4; i++) {
    PlotMaker pm;
	  temp = data_mc.at(i);
    if(i == 3) w_tot = w_fr2;
		tag = sample_label.at(i) + "_standard";
		// 1 lepton
		tag = sample_label.at(i) + "_1l";
    pm.Push<Hist1D>(Axis(10,200,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
		                BaselineCuts(),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,500,  2000, "st",  "S_{T} [GeV]",{}), 
		                BaselineCuts("st"),    temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                BaselineCuts("mt"),    temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0, 1500, "mj14", "M_{J} [GeV]",{}), 
		                BaselineCuts(),        temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(16,-0.5,15.5, "njets","N_{jets}",{}),    
		                BaselineCuts("njets"), temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(7,-0.5,  6.5, "nbdm",  "N_{b}",{}),
		                BaselineCuts("nbm"),   temp, lin_stack).Weight(w_tot).Tag(tag);
		// 2 lepton
		tag = sample_label.at(i) + "_2l";
    pm.Push<Hist1D>(Axis(10,200,  700,  "met", "p_{T}^{miss} [GeV]",{}), 
	                  BaselineCuts("","2l"),  temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  BaselineCuts("st","2l"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300,     "mt",  "m_{T} [GeV]",{}), 
		                BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0, 1500,  "mj14","M_{J} [GeV]",{}),
	                  BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(16,-0.5,15.5, "njets","N_{jets}",{}),
	                  BaselineCuts("njets","2l"),temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(7,-0.5,  6.5,  "nbdm","N_{b}",{}),
	                  BaselineCuts("nbm","2l"),  temp, lin_stack).Weight(w_tot).Tag(tag);
		// 5-6 jet
		tag = sample_label.at(i) + "_56jet";
    pm.Push<Hist1D>(Axis(10,200,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
	                  BaselineCuts("","56j"),  temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  BaselineCuts("st","56j"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                BaselineCuts("","56j"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0, 1500, "mj14", "M_{J} [GeV]",{}),
	                  BaselineCuts("","56j"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(10,-0.5,9.5, "njets","N_{jets}",{}),
	                  BaselineCuts("","56j"),  temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(7,-0.5,  6.5, "nbdm", "N_{b}",{}),
	                  BaselineCuts("nbm","56j"),  temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.min_print_=true;
    pm.MakePlots(sample_lumi.at(i));
  }
}

