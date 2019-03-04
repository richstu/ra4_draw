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
	                     "njets >= 6", "nbd >= 1", "nveto == 0"};
  int num_cuts(6);
	if(extra == "2l") {
		cuts[0] = "nleps == 2";
		cuts[2] = "met > 200 && met < 500";
		cuts[3] = "njets >= 5";
		cuts[4] = "nbd <= 2";
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

float dphi(double phi1, double phi2) {
  double abs_diff(abs(phi1-phi2));
	if(abs_diff > 6.283) return(abs_diff-6.283);
	else if(abs_diff > 3.142) return(6.283-abs_diff);
	else return(abs_diff);
}

int main() {
  gErrorIgnoreLevel = 6000;

  Palette colors("txt/colors.txt", "default");
	Process::Type back = Process::Type::background;
	Process::Type data = Process::Type::data;

	const NamedFunc hem_only("hem_only",[&](const Baby &b){
		if(b.SampleType() == -2018 && b.run() >= 319077) { 
			if(b.nels() > 0) {
				for(size_t i = 0; i < b.els_pt()->size(); i++) {
					if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
						return static_cast<float>(1);
				}
			}
			for(size_t i = 0; i < b.jets_pt()->size(); i++) {
				if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
					return static_cast<float>(1);
			}
		}
		else if(b.SampleType() == 2018) { 
			if(b.nels() > 0) {
				for(size_t i = 0; i < b.els_pt()->size(); i++) {
					if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
						return static_cast<float>(1);
				}
			}
			for(size_t i = 0; i < b.jets_pt()->size(); i++) {
				if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
					return static_cast<float>(1);
			}
		}
 		return static_cast<float>(0);
	});
	const NamedFunc wgt("wgt",[&](const Baby &b){
		return b.weight();
	});
  // Data
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2017_02_14/data/merged_database_standard/");
  string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/merged_database_standard/");
  string data18_path("/net/cms2/cms2r0/babymaker/babies/2019_01_18/data/merged_database_standard/");
	string q_cuts("pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {data16_path+"*.root"},"trig_ra4&&"+q_cuts);
	auto data_2016_back = Process::MakeShared<Baby_full>("2016 Data",back,kAzure,  
	                      {data16_path+"*.root"},"trig_ra4&&"+q_cuts);
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {data17_path+"*.root"},"trig_ra4&&"+q_cuts);
	auto data_2017_back = Process::MakeShared<Baby_full>("2017 Data",back,kAzure,
	                      {data17_path+"*.root"},"trig_ra4&&"+q_cuts);
	auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
	                 {data18_path+"*.root"},"trig_ra4&&"+q_cuts);

	auto data_2018_AC = Process::MakeShared<Baby_full>("2018 Data - Runs A-C",data,kBlack,
	                    {data18_path+"*Run2018A*.root",
	                     data18_path+"*Run2018B*.root",
	                     data18_path+"*Run2018C*.root" },"trig_ra4&&pass&&met/met_calo<5.0");
	auto data_2018_D  = Process::MakeShared<Baby_full>("2018 Data - Run D",data,kBlack,
	                    {data18_path+"*Run2018D*.root"},"trig_ra4&&pass&&met/met_calo<5.0");

	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/");
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
	vector<shared_ptr<Process> > data18_mc17  = {data_2018, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<shared_ptr<Process> > data16_data17 = {data_2017, data_2016_back};
	vector<shared_ptr<Process> > data17_data18 = {data_2018, data_2017_back};
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

	vector<shared_ptr<Process> > temp;
	vector<string> sample_label = {"2016", "2017", "2018"};
	vector<double> sample_lumi = {35.9, 41.5, 60};
	string tag;
	/*
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
                  BaselineCuts("mt"), data16_data17, log_stack).Tag("1617");
  pm.Push<Hist1D>(Axis(30,100, 700, "met",  "p_{T}^{miss} [GeV]",{}), 
                  BaselineCuts("met"), data16_data17, log_stack).Tag("1617");
  pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
                  BaselineCuts("mt"), data17_data18, log_stack).Tag("1718");
  pm.Push<Hist1D>(Axis(30,100, 700, "met",  "p_{T}^{miss} [GeV]",{}), 
                  BaselineCuts("met"), data17_data18, log_stack).Tag("1718");
  pm.min_print_=true;
	pm.MakePlots(1);
	*/
  NamedFunc w_tot("weight");;
	for(size_t i = 0; i < 1; i++) {
    PlotMaker pm;
	  temp = data_mc.at(i);
		tag = sample_label.at(i) + "_standard";
		if(i == 2) w_tot = wgt*Functions::hem_veto;
		// Standard
    pm.Push<Hist1D>(Axis(40,0, 80, "npv",  "N_{PV}",{}), 
                    "nleps>=1 && st>500 && met>100 && njets>=4", temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,100,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}),
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0,  300, "mt",   "m_{T} [GeV]",{}),        
	                  "nleps==1 && st>500 && met>100 && njets>=4", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(8,3.5,11.5,   "njets","N_{jets}",{}),
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(5,-0.5,  4.5, "nbd",  "N_{b, Deep}",{}),
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(25,0,  1, "jets_csvd",  "DeepCSV discriminator",{}),
	                  "nleps>=1 && st>500 && met>100 && njets>=4", temp, log_stack).Weight(w_tot).Tag(tag);
		// 1 lepton
		tag = sample_label.at(i) + "_1l";
    pm.Push<Hist1D>(Axis(30,100,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
		                BaselineCuts("met"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,500,  2000, "st",  "S_{T} [GeV]",{}), 
		                BaselineCuts("st"),    temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mj14 < 400", temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                BaselineCuts("mt"),        temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                BaselineCuts(),        temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}), 
		                BaselineCuts("mj"),        temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(8,3.5,11.5,   "njets","N_{jets}",{}),    
		                BaselineCuts("njets"), temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(5,-0.5,  4.5, "nbd",  "N_{b, Deep}",{}),
		                BaselineCuts("nbm"),   temp, lin_stack).Weight(w_tot).Tag(tag);
			// Electron variables
    pm.Push<Hist1D>(Axis(25,0, 500, "els_pt",  "electron p_{T} [GeV]",{}), 
		                BaselineCuts("","nels==1"), temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "els_sceta",  "electron #eta",{}), 
		                BaselineCuts("","nels==1"), temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(34,-3.4,3.4, "els_phi",  "electron #phi",{}), 
		                BaselineCuts("","nels==1"), temp, log_stack).Weight(w_tot).Tag(tag);
			// Muon variables
    pm.Push<Hist1D>(Axis(25,0, 500, "mus_pt",  "muon p_{T} [GeV]",{}), 
		                BaselineCuts("","nmus==1"), temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mus_eta",  "muon #eta",{}), 
		                BaselineCuts("","nmus==1"), temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(34,-3.4,3.4, "mus_phi",  "muon #phi",{}), 
		                BaselineCuts("","nmus==1"), temp, log_stack).Weight(w_tot).Tag(tag);
			// mT-MJ bins
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mt < 140", temp, log_stack).Weight(w_tot).Tag(tag+"_mt1_MJ");
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mt > 140", temp, log_stack).Weight(w_tot).Tag(tag+"_mt2_MJ");
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mj14>250 && mj14<400", temp, log_stack).Weight(w_tot).Tag(tag+"_MJ1_mt");
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mj14>400 && mj14<600", temp, log_stack).Weight(w_tot).Tag(tag+"_MJ2_mt");
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbd >= 1 && nveto == 0 && mj14>600", temp, log_stack).Weight(w_tot).Tag(tag+"_MJ3_mt");
		// 2 lepton
		tag = sample_label.at(i) + "_2l";
    pm.Push<Hist1D>(Axis(30,100,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
	                  BaselineCuts("met","2l"),  temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  BaselineCuts("st","2l"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}),
	                  BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(8,3.5,11.5,   "njets","N_{jets}",{}),
	                  BaselineCuts("njets","2l"),temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(5,-0.5,  4.5, "nbd", "N_{b, Deep}",{}),
	                  BaselineCuts("nbm","2l"),  temp, lin_stack).Weight(w_tot).Tag(tag);
		// 5 jet
		tag = sample_label.at(i) + "_5jet";
    pm.Push<Hist1D>(Axis(30,100,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
	                  BaselineCuts("met","5j"),  temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  BaselineCuts("st","5j"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                BaselineCuts("","5j"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(30,  0, 1200, "mj14", "M_{J} [GeV]",{}),
	                  BaselineCuts("","5j"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(5,-0.5,  4.5, "nbd", "N_{b, Deep}",{}),
	                  BaselineCuts("nbm","5j"),  temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.min_print_=true;
    pm.MakePlots(sample_lumi.at(i));
  }
}

