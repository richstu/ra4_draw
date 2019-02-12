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
	else if(var == "met")   cuts[2] = "met > 200";
	else if(var == "njets") out = 3;
// 	else if(var == "njets") cuts[3] = "njets >= 4";
	else if(var == "nbm")   out = 4;
	NamedFunc baseline(cuts[0]);
	for(int i = 1; i < num_cuts; i++) 
		if(i != out) baseline = baseline && cuts[i];
// 	if(extra == "1") baseline = baseline && "(mj14 < 400 || mt < 140)";
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

	const NamedFunc hem("hem_veto",[&](const Baby &b){
		if(b.SampleType() == 2018 && b.run() >= 319077) { 
			if(b.nels() > 0) {
				for(size_t i = 0; i < b.els_pt()->size(); i++) {
					if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
						return static_cast<float>(0);
				}
			}
			for(size_t i = 0; i < b.jets_pt()->size(); i++) {
				if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
					return static_cast<float>(0);
			}
		}
		else if(b.SampleType() == 2017 && (b.event()%1961) < 1296) { 
			if(b.nels() > 0) {
				for(size_t i = 0; i < b.els_pt()->size(); i++) {
					if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
						return static_cast<float>(0);
				}
			}
			for(size_t i = 0; i < b.jets_pt()->size(); i++) {
				if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
					return static_cast<float>(0);
			}
		}
 		return static_cast<float>(1);
	});
	const NamedFunc hem_only("hem_only",[&](const Baby &b){
		if(b.SampleType() == 2018 && b.run() >= 319077) { 
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
		else if(b.SampleType() == 2017) { 
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
  // Data
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/");
  string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/merged_database_standard/");
  string data18_path("/net/cms2/cms2r0/babymaker/babies/2019_01_18/data/merged_database_standard/");

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {data16_path+"*.root"},"trig_ra4&&pass&&met/met_calo<5.0");
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {data17_path+"*.root"},"trig_ra4&&pass&&met/met_calo<5.0");
	auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
	                 {data18_path+"*.root"},"trig_ra4&&pass&&met/met_calo<5.0");

	auto data_2018_AC = Process::MakeShared<Baby_full>("2018 Data - Runs A-C",data,kBlack,
	                    {data18_path+"*Run2018A*.root",
	                     data18_path+"*Run2018B*.root",
	                     data18_path+"*Run2018C*.root" },"trig_ra4&&pass&&met/met_calo<5.0");
	auto data_2018_D  = Process::MakeShared<Baby_full>("2018 Data - Run D",data,kBlack,
	                    {data18_path+"*Run2018D*.root"},"trig_ra4&&pass&&met/met_calo<5.0");

	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_standard/");
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc16_path+"*_TTJets*SingleLept*.root"}, "ntruleps<=1&&stitch_met");
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc16_path+"*_TTJets*DiLept*.root"}, "ntruleps>=2&&stitch_met");
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {mc16_path+"*_WJetsToLNu*.root"},"stitch_met");
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {mc16_path+"*_ST_*.root"});
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root"});
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
                       {mc16_path+"*DYJetsToLL*.root", mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
                        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
                        mc16_path+"*_TTGJets*.root", mc16_path+"*_TTTT_*.root",
                        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
                        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
                        mc16_path+"_ZZ_*.root"}, "stitch_met");

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc17_path+"*_TTJets*SingleLept*.root"}, "ntruleps<=1&&stitch_met");
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc17_path+"*_TTJets*DiLept*.root"}, "ntruleps>=2&&stitch_met");
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc17_path+"*_WJetsToLNu_*.root"},"(stitch_met && type!=2000) || (type==2000 && ht_isr_me<100)");
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc17_path+"*_ST_*.root"});
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root"});
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
                       {mc17_path+"*DYJetsToLL_M-50_HT*.root", mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
                        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                        mc17_path+"*_TTGJets*.root",           mc17_path+"*_TTTT_*.root",
                        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
//                         mc17_path+"*_WWTo*.root",           
                        mc17_path+"*_WZ*.root",
                        mc17_path+"_ZZ_*.root"}, "stitch_met");

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

	vector<shared_ptr<Process> > temp;
	vector<string> sample_label = {"2016", "2017", "2018"};
	vector<double> sample_lumi = {35.9, 41.5, 60};
	string tag;
  NamedFunc w_tot("weight");
	for(size_t i = 0; i < 3; i++) {
    PlotMaker pm;
	  temp = data_mc.at(i);
		tag = sample_label.at(i) + "_standard";
		if(i == 2) w_tot = w_tot*hem;
		// 1 lepton
		tag = sample_label.at(i) + "_1l";
    pm.Push<Hist1D>(Axis(10,200,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
		                BaselineCuts("met"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,500,  2000, "st",  "S_{T} [GeV]",{}), 
		                BaselineCuts("st"),    temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",  "m_{T} [GeV]",{}), 
		                BaselineCuts("mt"),    temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0, 1500, "mj14", "M_{J} [GeV]",{}), 
		                BaselineCuts(),        temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(16,-0.5,15.5,   "njets","N_{jets}",{}),    
		                BaselineCuts("njets"), temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(7,-0.5,  6.5, "nbd",  "N_{b, Deep}",{}),
		                BaselineCuts("nbm"),   temp, lin_stack).Weight(w_tot).Tag(tag);
		// 2 lepton
		tag = sample_label.at(i) + "_2l";
    pm.Push<Hist1D>(Axis(12,200,  700, "met",  "p_{T}^{miss} [GeV]",{}), 
	                  BaselineCuts("met","2l"),  temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,500,  2000, "st",  "S_{T} [GeV]",{}), 
	                  BaselineCuts("st","2l"),   temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,0, 300, "mt",      "m_{T} [GeV]",{}), 
		                BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(15,  0, 1500, "mj14", "M_{J} [GeV]",{}),
	                  BaselineCuts("","2l"),     temp, log_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(16,-0.5,15.5,   "njets","N_{jets}",{}),
	                  BaselineCuts("njets","2l"),temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.Push<Hist1D>(Axis(7,-0.5,  6.5, "nbd", "N_{b, Deep}",{}),
	                  BaselineCuts("nbm","2l"),  temp, lin_stack).Weight(w_tot).Tag(tag);
    pm.min_print_=true;
    pm.MakePlots(sample_lumi.at(i));
  }
}

