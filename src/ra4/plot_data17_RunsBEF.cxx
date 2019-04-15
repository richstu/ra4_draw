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

float dphi(double phi1, double phi2) {
  double abs_diff(abs(phi1-phi2));
	if(abs_diff > 6.283) return(abs_diff-6.283);
	else if(abs_diff > 3.142) return(6.283-abs_diff);
	else return(abs_diff);
}

// Run for RunB+RunC for mt plots
int main() {
  gErrorIgnoreLevel = 6000;
// If running on cms#, need to add /net/cms29/ to data_path names

  Palette colors("txt/colors.txt", "default");
	Process::Type data = Process::Type::data;
	Process::Type back = Process::Type::background;

	const NamedFunc dphi_lep("dphi_lep",[&](const Baby &b){
	  for(size_t i = 0; i < b.mus_pt()->size(); i++)
			if(b.mus_pt()->at(i) > 20 && abs(b.mus_eta()->at(i)) < 2.4 && b.mus_sigid()->at(i) && b.mus_miniso()->at(i) < 0.2)
			  return dphi(b.mus_phi()->at(0), b.met_phi());
	  for(size_t i = 0; i < b.els_pt()->size(); i++)
			if(b.els_pt()->at(i) > 20 && abs(b.els_sceta()->at(i)) < 2.5 && b.els_sigid()->at(i) && b.els_miniso()->at(i) < 0.1)
			  return dphi(b.els_phi()->at(0), b.met_phi());
		return static_cast<float>(3.4);
	});
	const NamedFunc dphi1("dphi1",[&](const Baby &b){
	    return dphi(b.met_phi(), b.jets_phi()->at(0));
	});
	const NamedFunc dphi2("dphi2",[&](const Baby &b){
	    return dphi(b.met_phi(), b.jets_phi()->at(1));
	});
  // Data
  string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/merged_database_standard/");
  string data17_prompt_path("/net/cms2/cms2r0/babymaker/babies/2018_01_30/data/merged_database_standard/");
// 	data17_path = data17_prompt_path;
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");
	string q_cuts("pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");
	bool doMC(true);

	auto data_2017_runBE = Process::MakeShared<Baby_full>("2017 Data - Runs B-E",data,kBlack,
	                      {data17_path+"*.root",
	                       data17_path+"*.root",
	                       data17_path+"*.root",
	                       data17_path+"*.root"},"run<=304797&&trig_ra4&&"+q_cuts);
	auto data_2017_runF = Process::MakeShared<Baby_full>("2017 Data - Run F",data,kBlack,
	                      {data17_path+"*.root"},"run>=305040&&trig_ra4&&"+q_cuts);
	auto data_2017_runF_back = Process::MakeShared<Baby_full>("2017 Data - Run F",back,kAzure,
	                           {data17_path+"*.root"},"run>=305040&&trig_ra4&&"+q_cuts);
	// MC											
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

	vector<shared_ptr<Process> > data17BE_data17F  = {data_2017_runBE, data_2017_runF_back};
	vector<shared_ptr<Process> > data17BE_mc17  = {data_2017_runBE, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<shared_ptr<Process> > data17F_mc17  = {data_2017_runF, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
	vector<vector<shared_ptr<Process> >> data17_BE_F = {data17BE_mc17, data17F_mc17};

	
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
	log_lumi.Title(TitleType::info)
					.YAxis(YAxisType::log)
					.FileExtensions({"pdf"});
  PlotOpt log_shapes_info = log_lumi().Stack(StackType::shapes)
	                                    .Bottom(BottomType::ratio)
																			.Title(TitleType::info);  
	vector<PlotOpt> log_shape = {log_shapes_info};
	vector<PlotOpt> log_stack = {log_shapes_info.Stack(StackType::data_norm)};
  PlotOpt lin_shapes_info = log_shapes_info().YAxis(YAxisType::linear);
	vector<PlotOpt> lin_shape = {lin_shapes_info};
	vector<PlotOpt> lin_stack = {lin_shapes_info.Stack(StackType::data_norm)};
  PlotMaker pm;
	TString mid_mT_1l("nleps == 1 && (mt > 100 && mt < 200) && st > 500 && met > 200 && njets >= 4 && nbdm >= 1 && nveto == 0");
	string bef("BEF");
  pm.Push<Hist1D>(Axis(15,0,300, "mt", "m_{T} [GeV]",{}), 
	                "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbdm >= 1 && nveto == 0", 
									data17BE_data17F, log_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(20,200,600, "met", "p_{T}^{miss} [GeV]",{}), 
	                mid_mT_1l, data17BE_data17F, log_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(20, 20, 420, "leps_pt[0]", "p_{T}^{lep} [GeV]",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(16,0,3.2, dphi_lep, "#Delta#phi(lep,MET)",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(16,0,3.2, dphi1, "#Delta#phi(Leading jet,MET)",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(16,0,3.2, dphi2, "#Delta#phi(Sub-Leading jet,MET)",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(20,80,880, "jets_pt[0]", "Leading jet p_{T} [GeV]",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.Push<Hist1D>(Axis(20,20,520, "jets_pt[1]", "Sub-Leading jet p_{T} [GeV]",{}), 
	                mid_mT_1l, data17BE_data17F, lin_shape).Tag(bef);
  pm.min_print_=true;
  pm.MakePlots(41.5);
	vector<double> lumis = {28.0, 13.5};
	vector<string> tags = {"BE_MC", "F_MC"};
	if(doMC) {
		for(int i = 0; i < 2; i++) {
	  	PlotMaker MCpm;
			vector<shared_ptr<Process> > temp = data17_BE_F.at(i);
			string tag = tags.at(i);
			MCpm.Push<Hist1D>(Axis(15,0,300, "mt", "m_{T} [GeV]",{}),
	                      "nleps == 1 && st > 500 && met > 200 && njets >= 4 && nbdm >= 1 && nveto == 0", 
												temp, log_stack).Tag(tag);
			MCpm.min_print_=true;
			MCpm.MakePlots(lumis.at(i));
		}
	}
}




