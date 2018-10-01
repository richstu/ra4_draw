#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string.h>

#include "TError.h"
#include "TColor.h"
#include "TVector3.h"

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

NamedFunc::ScalarType top_ijet(const Baby &b);
NamedFunc::ScalarType imax_NR(const Baby &b);
NamedFunc::ScalarType isub_NR(const Baby &b);
NamedFunc::ScalarType dRmin_lep(const Baby &b);

NamedFunc BaselineCuts(string var) {
	NamedFunc cuts[6] = {"nleps == 1", "st > 500", "met > 200", "nveto == 0", "njets >= 6", "nbm >= 1"};
	int out(-1);
	if     (var == "nleps") cuts[0] = "nleps >= 1";
	else if(var == "met")   out = 2;
	else if(var == "nveto") out = 3;
	else if(var == "njets") out = 4;
	else if(var == "nbm")   out = 5;
	NamedFunc baseline(cuts[0]);
	for(int i = 1; i < 6; i++) 
		if(i != out) baseline = baseline && cuts[i];
	return baseline;
}

string FuncName(NamedFunc var) {
	string name = var.Name();
	ReplaceAll(name,"_"," ");
	return name;
}

double dR(double deta, double dphi) {
	if      (dphi > 6.2830) dphi -= 6.2830;
	else if (dphi > 3.1415) dphi = 6.2830 - dphi;
	return(sqrt(pow(dphi,2)+pow(deta,2)));
}
// Run for RunB+RunC for mt plots
int main() {
  gErrorIgnoreLevel = 6000;
  bool quick(0);
// If running on cms#, need to add /net/cms29/ to data_path names
  string mc_standard_path("/net/cms29/cms29r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
  string ttbar1L(mc_standard_path+"*_TTJets*SingleLept*.root");
  string ttbar2L(mc_standard_path+"*_TTJets*DiLept*.root");
  string ttbarHT(mc_standard_path+"*_TTJets*HT*.root");
  string  signal1(mc_standard_path+"*mGluino-1200_mLSP-800*.root");
  string  signal2(mc_standard_path+"*mGluino-2000_mLSP-100*.root");
  if (quick) {
  	ttbar1L = mc_standard_path+"mergedbaby__TTJets_SingleLeptFromTbar_Tune_skim_standard_mcbase_nfiles_350.root";
  	ttbar2L = mc_standard_path+"mergedbaby__TTJets_DiLept_Tune_skim_standard_mcbase_nfiles_186.root";
	}
  Palette colors("txt/colors.txt", "default");
	Process::Type back   = Process::Type::background;
	Process::Type signal = Process::Type::signal;
	Process::Type data   = Process::Type::data;

  const NamedFunc top_pt("top_pt", [](const Baby &b) ->NamedFunc::ScalarType{
    int hadtop_id(0);
    TVector3 hadtop(0,0,0);
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find which top is hadronic (opposite sign as lepton)
      if(b.mc_id()->at(imc) == 11 || b.mc_id()->at(imc) == 13 || b.mc_id()->at(imc) == 15) hadtop_id = 6;
      else if(b.mc_id()->at(imc) == -11 || b.mc_id()->at(imc) == -13 || b.mc_id()->at(imc) == -15) hadtop_id = -6;
    }
    for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find hadronic top
      if(b.mc_id()->at(imc) == hadtop_id) hadtop.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
    }
    return hadtop.Pt();
    });
  const NamedFunc top_jet_mass_cut("top_jet_mass_cut", [](const Baby &b) ->NamedFunc::ScalarType{
		double mass(0);
		int itop = top_ijet(b);
		if(itop >= 0) mass = b.ak8jets_m()->at(itop);
    return(mass > 105 && mass < 210);
	});
  const NamedFunc top_jet_nom_score("top_jet_nom_score", [](const Baby &b) ->NamedFunc::ScalarType{
		double score(-1);
		int itop = top_ijet(b);
		if(itop >= 0) score = b.ak8jets_nom_bin_top()->at(itop);
    return score;
	});
  const NamedFunc top_jet_decor_score("top_jet_decor_score", [](const Baby &b) ->NamedFunc::ScalarType{
		double score(-1);
		int itop = top_ijet(b);
		if(itop >= 0) score = b.ak8jets_decor_bin_top()->at(itop);
    return score;
	});
  const NamedFunc is_300ak8("is_300ak8", [](const Baby &b) ->NamedFunc::ScalarType{
    int is300(0);
  	for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
	  	if(b.ak8jets_pt()->at(iak8) > 300) is300 = 1;
	  }
		return is300;
  });
	const NamedFunc deltaRmin_b("deltaRmin_b", [](const Baby &b) ->NamedFunc::ScalarType{
		double deltaR(10);
		double dphi(10), deta(10);
    for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++) { 
			if(b.ak8jets_pt()->at(ijet) > 300) {
				for(int ib(0); ib < b.njets(); ib++) {
					if(b.jets_csv()->at(ib) > 0.8484) {
						dphi = abs(b.ak8jets_phi()->at(ijet) - b.jets_phi()->at(ib));
						if      (dphi > 6.2830) dphi -= 6.2830;
						else if (dphi > 3.1415) dphi = 6.2830 - dphi;
						deta = b.ak8jets_eta()->at(ijet) - b.jets_eta()->at(ib);
						if(sqrt(pow(deta,2)+pow(dphi,2)) < deltaR || deltaR == -1) deltaR = sqrt(pow(deta,2)+pow(dphi,2));
						}
					}
				}
			}
		return deltaR;
	});
	const NamedFunc deltaRmin_lep("deltaRmin_lep", [](const Baby &b) ->NamedFunc::ScalarType{
		return dRmin_lep(b);
	});
	const NamedFunc max_NR("Highest_top_score",[&](const Baby &b){
		double score(-1);
		int imax = imax_NR(b);
		if(imax >= 0) score = b.ak8jets_nom_raw_top()->at(imax);
		return score;
	});
	const NamedFunc subleading_nom_raw("Subleading_NR_score",[&](const Baby &b){
		double score(-1);
		int isub = isub_NR(b);
		if(isub >= 0) score = b.ak8jets_nom_raw_top()->at(isub);
		return score;
	});
	const NamedFunc subleading_pt("Subleading_pt",[&](const Baby &b){
		double sub_pt(-1);
		if(b.nak8jets() > 1) sub_pt = b.ak8jets_pt()->at(1);
		return sub_pt;
	});
	const NamedFunc max_ak8_m("Highest_AK8_Mass",[&](const Baby &b){
		double max(0.), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			mass  = b.ak8jets_m()->at(ijet);
			if(mass > max && b.ak8jets_pt()->at(ijet) >= 300) max = mass;
		}
		return max;
	});
	const NamedFunc maxNR_ak8_m("MaxNR_AK8_Mass",[&](const Baby &b){
		double mass(-1);
		int imax = imax_NR(b);
		if(imax >= 0) mass = b.ak8jets_m()->at(imax);
		return mass;
	});
	const NamedFunc subleading_ak8_m("Subleading_AK8_mass",[&](const Baby &b){
		double max1(0.), max2(0.), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			mass = b.ak8jets_m()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 300);
			if(mass > max1) {
				max2 = max1;
				max1 = mass;
			}
			else if(mass > max2) max2 = mass;
		}
		return max2;
	});
	const NamedFunc ak8_deltaR("#DeltaR(Max NR, Sub NR)",[&](const Baby &b){
		int i1(imax_NR(b)), i2(isub_NR(b));
		if(i1 >= 0 && i2 >= 0) {
			double deta = abs(b.ak8jets_eta()->at(i1) - b.ak8jets_eta()->at(i2));
			double dphi = abs(b.ak8jets_phi()->at(i1) - b.ak8jets_phi()->at(i2));
			return dR(deta,dphi);
		}
		else return(9.);
	});
	const NamedFunc ak8_deltaPhi("#Delta#phi(Max NR, Sub NR)",[&](const Baby &b){
		int i1(imax_NR(b)), i2(isub_NR(b));
		if(i1 >= 0 && i2 >= 0) {
			double dphi = abs(b.ak8jets_phi()->at(i1) - b.ak8jets_phi()->at(i2));
			return dR(0,dphi);
		}
		else return(9.);
	});
	const NamedFunc ak8_met_deltaPhi("#Delta#phi(Max NR, MET)",[&](const Baby &b){
		int i1 = imax_NR(b);
		if(i1 >= 0) {
		double dphi = abs(b.ak8jets_phi()->at(i1) - b.met_phi());
		return dR(0,dphi);
		}
		else return(9.);
	});
// MC samples
  auto tt1l = Process::MakeShared<Baby_full>("1L t#bar{t}", back, kBlack,           {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met");
  auto tt2l = Process::MakeShared<Baby_full>("2L t#bar{t}", back, colors("tt_2l"),  {ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met");
  // MJ check
  auto tt1l_MJcheck = Process::MakeShared<Baby_full>("1L t#bar{t} - m_{T} #leq 140", back, kBlack,          {ttbar1L,ttbarHT}, "ntruleps==1&&stitch_met&&mt<=140");
  auto tt2l_MJcheck = Process::MakeShared<Baby_full>("2L t#bar{t} - m_{T} > 140",    back, colors("tt_2l"), {ttbar2L,ttbarHT}, "ntruleps==2&&stitch_met&&mt>140");
  // With atleast 1 ak8jet
  auto tt1l_wAK8 = Process::MakeShared<Baby_full>("1L t#bar{t} w/ AK8", back, colors("tt_1l"),  {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met" && "nak8jets>0" );
  auto tt2l_wAK8 = Process::MakeShared<Baby_full>("2L t#bar{t} w/ AK8", back, colors("tt_2l"),  {ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met" && "nak8jets>0" );
  // Nominal Scores (0.9828 bin, 0.7450 raw)
  vector<double>   nom_cuts = {0.1883 , .8511, 0.9377, 0.7450, 0.9897};
  vector<shared_ptr<Process> > ttbar1l_nom_samples = {};
  // Decorrelated Scores (0.6730 bin, 0.1850 raw)
  vector<double> decor_cuts = {0.04738, .4585, 0.6556, 0.1900, 0.8931};
  vector<shared_ptr<Process> > ttbar1l_decor_samples = {};
  for(unsigned i = 0; i < nom_cuts.size(); i++) {
		// Nominal Samples
    auto   nom_1l_temp = Process::MakeShared<Baby_full>("1L t#bar{t} - TMed & tagged AK8s", data, kBlue, {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met" &&  top_jet_nom_score > nom_cuts.at(i));
		// Decorellated Samples
    auto decor_1l_temp = Process::MakeShared<Baby_full>("1L t#bar{t} - TMed & tagged AK8s", data, kBlue, {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met" && (top_jet_decor_score > decor_cuts.at(i)) && top_jet_mass_cut == 1);
		ttbar1l_decor_samples.push_back(decor_1l_temp);
		ttbar1l_nom_samples.push_back(nom_1l_temp);
	}
  auto decor_2l_1Loose = Process::MakeShared<Baby_full>("2L t#bar{t} 1 D-Loose top", data, kBlue, {ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met && ntops_decor == 1");
  auto   nom_2l_1Loose = Process::MakeShared<Baby_full>("2L t#bar{t} 1 N-Loose top", data, kBlue, {ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met && ntops_nom == 1  ");
  auto   ttbar2l_NR    = Process::MakeShared<Baby_full>("2L t#bar{t} 1 NR top",      data, kBlue, {ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met" && max_NR > 0.40);
  vector<shared_ptr<Process> > ttbar2l_eff_decor_samples = {tt2l_wAK8, decor_2l_1Loose};
  vector<shared_ptr<Process> > ttbar2l_eff_nom_samples   = {tt2l_wAK8,   nom_2l_1Loose};
	vector<shared_ptr<Process> > ttbar2l_eff_NR = {tt2l_wAK8, ttbar2l_NR};

  auto tt1l_eff = Process::MakeShared<Baby_full>("1L t#bar{t} - ak8 jets truth-mathced to top",back , colors("tt_1l"), {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met" && top_jet_decor_score > -0.01);

  auto tt1l_wAK8_toppt0 = Process::MakeShared<Baby_full>("1L t#bar{t}, top p_{T} < 400",       back, colors("tt_1l"), {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met&&nak8jets>0" && top_pt < 400                );
  auto tt1l_wAK8_toppt1 = Process::MakeShared<Baby_full>("1L t#bar{t}, 400 < top p_{T} < 500", back, kMagenta-3,      {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met&&nak8jets>0" && top_pt < 500 && top_pt > 400);
  auto tt1l_wAK8_toppt2 = Process::MakeShared<Baby_full>("1L t#bar{t}, 500 < top p_{T} < 800", back, kPink-9,         {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met&&nak8jets>0" && top_pt < 800 && top_pt > 500);
  auto tt1l_wAK8_toppt3 = Process::MakeShared<Baby_full>("1L t#bar{t}, top p_{T} > 800",       back, kGreen-3,        {ttbar1L,ttbarHT}, "ntruleps<=1&&stitch_met&&nak8jets>0" &&                 top_pt > 800);

  vector<shared_ptr<Process> > ttbar_samples_pt_binned = {tt1l_wAK8_toppt0, tt1l_wAK8_toppt1, tt1l_wAK8_toppt2, tt1l_wAK8_toppt3};
	string sig_path("/net/cms29/cms29r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
  string MJcut_sig_path("/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/merged_mcbase_abcd/");
  auto T1tttt_2000_100  = Process::MakeShared<Baby_full>("S(2000,100)",  signal, colors("t1tttt"),{sig_path+"*2000_mLSP-100_*.root"}, "stitch_met" && is_300ak8);
  auto T1tttt_1200_800  = Process::MakeShared<Baby_full>("S(1200,800)",  signal, kMagenta        ,{sig_path+"*1200_mLSP-800_*.root"}, "stitch_met" && is_300ak8);
//   vector<shared_ptr<Process> > ttbar_samples =  {tt1l, tt2l, T1tttt_1800_1200, T1tttt_2000_100};
//   vector<shared_ptr<Process> > ttbar_samples_wAK8 =  {tt1l_wAK8, tt2l_wAK8, T1tttt_1800_1200, T1tttt_2000_100};
  vector<shared_ptr<Process> > ttbar_MJcheck_samples = {tt1l_MJcheck, tt2l_MJcheck};
	vector<shared_ptr<Process> > all_procs = {tt1l_wAK8, tt2l_wAK8, T1tttt_2000_100};
	vector<shared_ptr<Process> > ttbar_sig = {tt1l_wAK8, tt2l_wAK8, T1tttt_2000_100, T1tttt_1200_800};

	auto MJ_ttbar            = Process::MakeShared<Baby_full>("t#bar{t}",  back, kAzure, {"/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/*TTJets*.root"},"stitch_met");
  auto MJ_T1tttt_1800_1000 = Process::MakeShared<Baby_full>("S(1800,1000)", signal, kMagenta,      {MJcut_sig_path+"*1800_mLSP-1000_*.root"},"stitch_met");
  auto MJ_T1tttt_2100_100  = Process::MakeShared<Baby_full>("S(2100,100) low MET",  signal, kViolet-4,      {MJcut_sig_path+"*2100_mLSP-100_*.root"},"stitch_met");

	auto MJ_ttbar_lowMET  = Process::MakeShared<Baby_full>("t#bar{t} low MET",  back, kAzure, {"/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/*TTJets*.root"},"stitch_met && met>200 && met<=350");
	auto MJ_ttbar_highMET = Process::MakeShared<Baby_full>("t#bar{t} high MET", back, kAzure, {"/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/*TTJets*.root"},"stitch_met && met > 500");
  auto MJ_T1tttt_2100_100_lowMET  = Process::MakeShared<Baby_full>("S(2100,100) low MET",  signal, kViolet-4,      {MJcut_sig_path+"*2100_mLSP-100_*.root"},"stitch_met && met>200 && met<=350");
  auto MJ_T1tttt_2100_100_highMET = Process::MakeShared<Baby_full>("S(2100,100) high MET", signal, kViolet-4,      {MJcut_sig_path+"*2100_mLSP-100_*.root"},"stitch_met && met > 500");
	MJ_ttbar_highMET->SetLineStyle(2);
	MJ_T1tttt_2100_100_highMET->SetLineStyle(2);

	vector<shared_ptr<Process> > MJ_shape_MET = {MJ_ttbar_lowMET, MJ_ttbar_highMET, MJ_T1tttt_2100_100_lowMET, MJ_T1tttt_2100_100_highMET};
	vector<shared_ptr<Process> > MJ_shape     = {MJ_ttbar, MJ_T1tttt_2100_100, MJ_T1tttt_1800_1000};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
  log_ratios.YAxis(YAxisType::log)
    .Bottom(BottomType::ratio)
  	.BottomHeight(0.64)
		.RatioMinimum(0.0) // Set window ranges (default is 0.1-1.9)
		.RatioMaximum(0.5)
		.Title(TitleType::info);
	PlotOpt lin_ratios = log_ratios().YAxis(YAxisType::linear);
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .FileExtensions({"pdf"});
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::ratio);
	PlotOpt log_stacks_noRatio = log_lumi().Bottom(BottomType::off).Title(TitleType::info);
	PlotOpt lin_stacks_noRatio = log_stacks_noRatio().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info).PrintVals(false);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
	PlotOpt log_shapes_noRatio = log_shapes_info().Bottom(BottomType::off);
	PlotOpt lin_shapes_noRatio = log_shapes_noRatio().YAxis(YAxisType::linear);
  PlotOpt log_scores = log_lumi_info().CanvasWidth(800).RatioMaximum(2.4);
  PlotOpt log_scores_shape = log_shapes_info().CanvasWidth(800).RatioMaximum(2.4);
  PlotOpt lin_shapes_info = log_shapes_info().YAxis(YAxisType::linear);
  vector<PlotOpt> log_plot  = {log_lumi_info};
  vector<PlotOpt> log_shape = {log_shapes_info};
  vector<PlotOpt> lin_shape = {lin_shapes_info};
  vector<PlotOpt> log_shape_noRatio = {log_shapes_noRatio};
  vector<PlotOpt> lin_shape_noRatio = {lin_shapes_noRatio};
  vector<PlotOpt> log_stack_noRatio = {log_stacks_noRatio};
  vector<PlotOpt> lin_stack_noRatio = {lin_stacks_noRatio};
  vector<PlotOpt> log_score = {log_scores};
  vector<PlotOpt> log_ratio = {log_ratios};
  vector<PlotOpt> lin_ratio = {lin_ratios};
  vector<PlotOpt> log_score_shape = {log_scores_shape};
  string baseline(           "nleps == 1 && st > 500 && met > 200  && njets >= 6 && nbm >= 1 && nveto==0");
  string baseline_1l(        "nleps == 1 && st > 500 && met > 200  && njets >= 6 && nbm >= 1 && nveto==0 && mt < 140");
  string baseline_1l_nob(    "nleps == 1 && st > 500 && met > 200  && njets >= 6 &&             nveto==0 && mt < 140");
  string baseline_1l_MJcheck("nleps == 1 && st > 500 && met > 200     njets >= 6 && nbm >= 1 && nveto==0");
  string baseline_2l(        "nleps == 2 && st > 500 && mj14 > 200 && njets >= 5 && nbm <= 2             && mj14 < 500 ");
  string baseline_2l_nob(    "nleps == 2 && st > 500 && mj14 > 200 && njets >= 5                         && mj14 < 500 ");
  string baseline_mt(        "nleps == 1 && st > 500 && met > 200  && njets >= 6 && nbm >= 1 && nveto==0 && mt>140");
  string baseline_mt_mj(     "nleps == 1 && st > 500 && met > 200  && njets >= 6 && nbm >= 1 && nveto==0 && mt>140 && mj14>500");
//   vector<string> baseline     = {baseline_1l,     baseline_2l};
  vector<string> baseline_nob = {baseline_1l_nob, baseline_2l_nob};
  PlotMaker pm;
  TString var(""), tag(""), cuts("");
  vector<string> rtype = {"1Lep", "2Lep"};
  vector<string> mtype = {"nom", "decor"};
  vector<string> stype = {"raw", "bin"};
  vector<double> wp = {0.1883, 0.0473};

	pm.Push<Hist1D>(Axis(15, 0,1500, "mj14", "M_{J} [GeV]", {}), "st > 500 && nleps == 1 && njets >= 6 && met > 200 && met <= 350", MJ_shape,     lin_shape_noRatio).Tag("MJ_shape_lowMET");
	pm.Push<Hist1D>(Axis(15, 0,1500, "mj14", "M_{J} [GeV]", {}), "st > 500 && nleps == 1 && njets >= 6 && met > 350 && met <= 500", MJ_shape,     lin_shape_noRatio).Tag("MJ_shape_midMET");
	pm.Push<Hist1D>(Axis(15, 0,1500, "mj14", "M_{J} [GeV]", {}), "st > 500 && nleps == 1 && njets >= 6 && met > 500",               MJ_shape,     lin_shape_noRatio).Tag("MJ_shape_highMET");
	pm.Push<Hist1D>(Axis(15, 0,1500, "mj14", "M_{J} [GeV]", {}), "st > 500 && nleps == 1 && njets >= 6", 							              MJ_shape_MET, lin_shape_noRatio).Tag("MJ_shape");

//	// Baseline Distributions
//	pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "nleps", "N_{leps}",    {0.5,1.5}), BaselineCuts("nleps"), ttbar_sig, log_shape_noRatio).Tag("baseline");
//	pm.Push<Hist1D>(Axis(20, 500,2500, "st",    "S_{T} [GeV]",        {}), BaselineCuts("st"),    ttbar_sig, log_shape_noRatio).Tag("baseline");
//	pm.Push<Hist1D>(Axis(28, 100,1500, "met",   "p_{T}^{miss} [GeV]", {}), BaselineCuts("met"),   ttbar_sig, log_shape_noRatio).Tag("baseline");
//	pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "nveto", "N_{veto}",        {0.5}), BaselineCuts("nveto"), ttbar_sig, log_shape_noRatio).Tag("baseline");
//	pm.Push<Hist1D>(Axis(12,-0.5,11.5, "njets", "N_{jets}",        {5.5}), BaselineCuts("njets"), ttbar_sig, log_shape_noRatio).Tag("baseline");
//	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nbm",   "N_{b}",           {0.5}), BaselineCuts("nbm"),   ttbar_sig, log_shape_noRatio).Tag("baseline");
//	// Score distributions in R4
//	NamedFunc low_jet("njets >= 6 && njets <= 8"), high_jet("njets >= 9");
//	NamedFunc low_b("nbm == 1"), mid_b("nbm == 2"), high_b("nbm >= 3");
//	vector<NamedFunc> R24 = {low_jet && low_b, high_jet && low_b, low_jet && mid_b, high_jet && mid_b, low_jet && high_b, high_jet && high_b};
//	vector<NamedFunc> met_bins = {"met > 200 && met <= 350", "met > 350 && met <= 500", "met > 500"};
//	vector<NamedFunc> low_met_abcd  = {"mt<140 && mj14>250 && mj14<400","mt<140 && mj14>400","mt>140 && mj14>250 && mj14<400","mt>140 && mj14>400"}; 
//	vector<NamedFunc> mid_met_abcd  = {"mt<140 && mj14>250 && mj14<500","mt<140 && mj14>500","mt>140 && mj14>250 && mj14<500","mt>140 && mj14>500"}; 
//	vector<NamedFunc> high_met_abcd = {"mt<140 && mj14>250 && mj14<600","mt<140 && mj14>600","mt>140 && mj14>250 && mj14<600","mt>140 && mj14>600"}; 
//	vector<vector<NamedFunc> > abcds = {low_met_abcd, mid_met_abcd, high_met_abcd};
//	vector<string> mets = {"Low", "Mid", "High"};
//	for(size_t i = 0; i < met_bins.size(); i++) { // MET bins
//		pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR), {0.4}), met_bins.at(i) && BaselineCuts("met"), ttbar_sig, lin_stack_noRatio).Tag(mets.at(i)+"_met_cutflow_1");
//		for(size_t r = 0; r < 4; r++) { // ABCD regions
//			if(r == 1 || r == 3)
//				for(size_t j = 0; j < R24.size(); j++) // sub-ABCD regions for R2 and R4
//					pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR), {0.4}),met_bins.at(i) && BaselineCuts("met") && abcds.at(i).at(r) && R24.at(j),ttbar_sig,lin_stack_noRatio)
//						.Tag(mets.at(i)+"_met_cutflow_R"+to_string(r+1)+"_"+to_string(j+1));
//				pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR), {0.4}), met_bins.at(i) && BaselineCuts("met") && abcds.at(i).at(r),ttbar_sig,lin_stack_noRatio)
//					.Tag(mets.at(i)+"_met_cutflow_R"+to_string(r+1));
//		}
//	}
//
// 	pm.Push<Hist1D>(Axis(20,0,5,    ak8_deltaR,       ak8_deltaR.Name(),       {}), baseline && "nak8jets>1", all_procs, log_stack_noRatio);
// 	pm.Push<Hist1D>(Axis(32,0,3.2,  ak8_deltaPhi,     ak8_deltaPhi.Name(),     {}), baseline && "nak8jets>1", all_procs, log_stack_noRatio);
// 	pm.Push<Hist1D>(Axis(32,0,3.2,  ak8_met_deltaPhi, ak8_met_deltaPhi.Name(), {}), baseline && "nak8jets>0", all_procs, log_stack_noRatio);
// 	pm.Push<Hist1D>(Axis(5,0,5,"nak8jets", "N_{AK8}"        , {}), baseline, all_procs, lin_shape_noRatio);
	// Stack plots
   pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR),{}),    baseline && "mj14<=400",                    all_procs, log_stack_noRatio).Tag( "lowMJ_NR_stack");
   pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR),{}),    baseline && "mj14>400&&mj14<=800",          all_procs, log_stack_noRatio).Tag( "midMJ_NR_stack");
   pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR),{}),    baseline && "mj14>800",                     all_procs, log_stack_noRatio).Tag("highMJ_NR_stack");
   pm.Push<Hist1D>(Axis(20,0,1, max_NR, FuncName(max_NR),{0.4}), baseline,                                   all_procs, log_stack_noRatio).Tag( "allMJ_NR_stack");
   pm.Push<Hist1D>(Axis(24,200,1400, "mj14", "M_{J} [GeV]", {}), baseline && max_NR <= 0.4,                  all_procs, log_stack_noRatio).Tag( "lowNR_MJ_stack");
   pm.Push<Hist1D>(Axis(24,200,1400, "mj14", "M_{J} [GeV]", {}), baseline && max_NR >  0.4 && max_NR <= 0.8, all_procs, log_stack_noRatio).Tag( "midNR_MJ_stack");
   pm.Push<Hist1D>(Axis(24,200,1400, "mj14", "M_{J} [GeV]", {}), baseline && max_NR >  0.8,                  all_procs, log_stack_noRatio).Tag("highNR_MJ_stack");
   pm.Push<Hist1D>(Axis(24,200,1400, "mj14", "M_{J} [GeV]", {}), baseline,                                   all_procs, log_stack_noRatio).Tag( "allNR_MJ_stack");
// 	pm.Push<Hist1D>(Axis(10.,200.,1200., "mj14", "M_{J}",{400}),                " met > 350 && met < 500 && nleps == 1 && mt <= 140 && nbm >= 1 && njets >= 6 && nveto == 0", T1tttt_samples,log_shape_noRatio).Tag("sigMJshapes_midmet");
// 	pm.Push<Hist1D>(Axis(10.,200.,1200., "mj14", "M_{J}",{400}),                "              met > 500 && nleps == 1 && mt <= 140 && nbm >= 1 && njets >= 6 && nveto == 0", T1tttt_samples,log_shape_noRatio).Tag("sigMJshapes_highmet");
// 	pm.Push<Hist1D>(Axis(18., 100.,1000., "met" ,   "MET",{100,150,200,350,500}),"mj14 > 200              && nleps == 1 && mt <= 140 && nbm >= 1 && njets >= 6 && nveto == 0", T1tttt_samples,log_shape_noRatio).Tag("sigMETshapes");

  /*
  // akjets14_m
  for(unsigned r = 0; r < 2; r++) { // Single-lepton CR | Dilepton CR
	  for(unsigned k = 0; k < 2; k++) { // All ak8 masses | only leading jets
		var = "ak8jets_m";
		if(k == 1) var.ReplaceAll("_m","_m[0]"); 
	    pm.Push<Hist1D>(Axis(45, 0., 450., var, "Mass of AK8 jets", {}), baseline.at(r), ttbar_samples_wAK8, log_shape).Tag(rtype.at(r));
	  	for(unsigned i = 0; i < 2; i++) { // Nominal | Decorrelated
	  		for(unsigned j = 0; j < 2; j++) { // Raw | binarized
				cuts = "ak8jets_"+mtype.at(i)+"_"+stype.at(j)+"_top>"+to_string(wp.at(i));
				if(k == 1) cuts.ReplaceAll("_top","_top[0]"); 
	    		pm.Push<Hist1D>(Axis(45, 0., 450., var, "Mass of AK8 jets", {}), cuts && baseline.at(r), ttbar_samples_wAK8, log_shape).Tag(rtype.at(r));
		}}}}
  pm.Push<Hist1D>(Axis(45, 0., 450., "ak8jets_m", "Mass of AK8 jets", {}), baseline_1l, ttbar_samples_pt_binned, log_shape).Tag("top_pt_binned");
	*/
/*
  // MJ
  vector<TString> metcuts;
  //metcuts.push_back("met>100");
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");
  for(auto &imet: metcuts)
		pm.Push<Hist1D>(Axis(30.,0.,1500., "mj14", "M_{J}",{250,400}),"nleps == 1 && mt <= 140 && nbm >= 1 && njets >= 6 && nveto == 0 &&" + imet, ttbar_samples,log_shape).Tag("MJ_MET");
	pm.Push<Hist1D>(Axis(18.,100.,1000., "met" ,   "MET",{100,150,200,350,500}),"mj14 > 200 && nleps == 1 && mt <= 140 && nbm >= 1 && njets >= 6 && nveto == 0", ttbar_samples,log_shape).Tag("MJ_MET");
*/
/*
  vector<NamedFunc> nf_cuts = {"1", is_300ak8, "ntops_decor > 0", "ntops_nom > 0"};
  vector<string> mj_tags = {"","300ak8_","decor_","nom_"};
  for(unsigned r = 0; r < 2; r++) { // Single-lepton CR | Dilepton CR
  	for(unsigned j = 0; j < 4; j++) { // no cut | ak8 w/ pt > 300 | w/ decor top | w/ nom top
  		pm.Push<Hist1D>(Axis(30, 0., 1500., "mj14", "M_{J} [GeV]", {250.,400.}), nf_cuts.at(j) && baseline.at(r), ttbar_samples,log_shape).Tag(string(mj_tags.at(j))+rtype.at(r));
		}}
  // Top score
  double loose, medium, tight, tightest;
  string bcut("");
  for(unsigned r = 0; r < 2; r++) { // Single-lepton CR | Dilepton CR
	  for(unsigned i = 0; i < 2; i++) { // Nominal | Decorrelated scores
	  	for(unsigned j = 0; j < 2; j++) { // Raw | binarized scores
				// Only want working point markers on binarized plots 
				if(j == 1 && i == 0) { loose = 0.1883; medium = 0.8511; tight = 0.9377; tightest = 0.9897; }
				else if(j == 1)      { loose = 0.0474; medium = 0.4585; tight = 0.6556; tightest = 0.8931; }
				else                 { loose = -1;     medium = -1;     tight = -1;     tightest = -1; }
	  		for(unsigned k = 0; k < 5; k++) { // Different Nb bins (0,1,2,3+,1+/2-)
					if(k < 3)       bcut = "&& nbm == " + to_string(k);
	  			else if(k == 3) bcut = "&& nbm >= 3"; 
	  			else if(r == 0) bcut = "&& nbm >= 1"; 
					else            bcut = "&& nbm <= 2"; 
	  			for(unsigned t = 0; t < 2; t++) { // All ak8 jets | just leading
						if (t == 0) { // plot scores of all ak8 jets
							var = "ak8jets_"+mtype.at(i)+"_"+stype.at(j)+"_top";
							tag = mtype.at(i)+" "+stype.at(j)+" score"; 
						}
						else { // only plot scores of leading ak8 jet
							var = "ak8jets_"+mtype.at(i)+"_"+stype.at(j)+"_top[0]";
							tag = mtype.at(i)+" "+stype.at(j)+" score0"; 
						}
	  				pm.Push<Hist1D>(Axis(20, 0., 1., var, string(tag), {loose,medium,tight,tightest}), baseline_nob.at(r)+bcut, ttbar_samples_wAK8, log_score_shape).Tag(string(tag.ReplaceAll(" ","_")+"_nb"+to_string(k))+"_"+rtype.at(r));
					}	
				}
			}
		}
	}
  // Efficiency
  vector<string> tags = {"loose","medium","tight","eff50","tightest"};
  vector<shared_ptr<Process> > ttbar_eff_samples;
  for(unsigned r = 0; r < 2; r++) { // Single-lepton CR | Dilepton CR
		for(unsigned i = 0; i < tags.size(); i++) { // loose | medium | tight | 50% eff at pt 600 | tightest
			for(unsigned j = 0; j < 2; j++) { // Nominal | Decorrelated
				if(j == 0) ttbar_eff_samples = {ttbar1l_nom_samples.at(i),   tt1l_eff};
				else       ttbar_eff_samples = {ttbar1l_decor_samples.at(i), tt1l_eff};
	  			pm.Push<Hist1D>(Axis(11, 300., 1400, top_pt, "hadronic top p_{T} [GeV]", {}), baseline.at(r), ttbar_eff_samples, log_ratio).Tag(mtype.at(j)+"_"+tags.at(i)+"_top_pt_"+rtype.at(r));
			}
		}
	}
	*/
  string baseline_noNjets("nleps == 1 && st > 500 && met > 200  && nbm >= 1   && nveto==0");
  string baseline_noNbm(  "nleps == 1 && st > 500 && met > 200  && njets >= 6 && nveto==0");
//	pm.Push<Hist1D>(Axis(8, 4.5, 12.5,"njets", "N_{jets}",               {}), baseline_noNjets, ttbar2l_eff_NR, lin_ratio).Tag("NReff_2l_njets");
//	pm.Push<Hist1D>(Axis(5,-0.5, 4.5, "nbm",    "N_{b}",                 {}), baseline_noNbm,   ttbar2l_eff_NR, lin_ratio).Tag("NReff_2l_nbm");
//	pm.Push<Hist1D>(Axis(20,  0, 400, max_ak8_m,"Highest mass AK8 [GeV]",{}), baseline,         ttbar2l_eff_NR, log_ratio).Tag("NReff_2l_ak8m");
//	pm.Push<Hist1D>(Axis(20,  0, 400, maxNR_ak8_m,"MaxNR AK8 mass [GeV]",{}), baseline,         ttbar2l_eff_NR, log_ratio).Tag("NReff_2l_NRak8m");
//	pm.Push<Hist1D>(Axis(20,300,1100,"ak8jets_pt[0]","Max AK8 p_{T}[GeV]",{}), baseline,        ttbar2l_eff_NR, log_ratio).Tag("NReff_2l_ak8pt");
//	pm.Push<Hist1D>(Axis(30,  0, 1200,"mj14"   ,"M_{J} [GeV]",           {}), baseline,         ttbar2l_eff_NR, log_ratio).Tag("NReff_2l_mj");
//	pm.Push<Hist1D>(Axis(7,-0.5, 6.5, "nak8jets" ,"N_{AK8}",             {}), baseline,         ttbar2l_eff_NR, log_ratio).Tag("NReff_2l_nak8");
  pm.min_print_=true;
  pm.MakePlots(135);
}

NamedFunc::ScalarType hadtop_i(const Baby &b) {
  int hadtop_id(0), itop(-1);
  TVector3 hadtop(0,0,0);
  for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find which top is hadronic (opposite sign as lepton)
    if(b.mc_id()->at(imc) == 11 || b.mc_id()->at(imc) == 13 || b.mc_id()->at(imc) == 15) hadtop_id = 6;
    else if(b.mc_id()->at(imc) == -11 || b.mc_id()->at(imc) == -13 || b.mc_id()->at(imc) == -15) hadtop_id = -6;
  }
  for(size_t imc(0); imc < b.mc_pt()->size(); imc++) { // Find hadronic top
    if(b.mc_id()->at(imc) == hadtop_id) itop = imc;
  }
  return itop;
  }

NamedFunc::ScalarType top_ijet(const Baby &b) {
	int itop = hadtop_i(b);
	int itopjet(-1);
	TVector3 hadtop(0,0,0), jet(0,0,0);
	hadtop.SetPtEtaPhi(b.mc_pt()->at(itop),b.mc_eta()->at(itop),b.mc_phi()->at(itop));
  for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++) { // Truth-match jet to hadronic top
    jet.SetPtEtaPhi(b.ak8jets_pt()->at(ijet),b.ak8jets_eta()->at(ijet),b.ak8jets_phi()->at(ijet));
		if(jet.DeltaR(hadtop) < 0.6) itopjet = ijet;
  }
	return itopjet;
}

NamedFunc::ScalarType imax_NR(const Baby &b) {
	double max(0.), score(0);
	int imax(-1);
	for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
		score = b.ak8jets_nom_raw_top()->at(ijet);
// 		if(score > max && b.ak8jets_pt()->at(ijet) >= 300) {
		if(score > max && b.ak8jets_pt()->at(ijet) >= 30) {
			max = score;
			imax = ijet;
		}
	}
	return imax;
}

NamedFunc::ScalarType isub_NR(const Baby &b) {
	double max1(0.), max2(0.), score(0);
	int i1(-1), i2(-1);
	for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
// 		score = b.ak8jets_nom_raw_top()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 300);
		score = b.ak8jets_nom_raw_top()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 30);
		if(score > max1) {
			i2 = i1; i1 = ijet;
			max2 = max1; max1 = score;
		}
		else if(score > max2) {
			max2 = score;
			i2 = ijet;
		}
	}
	return i2;
}

NamedFunc::ScalarType dRmin_lep(const Baby &b) {
	double deltaR(10);
	double dphi(10), deta(10);
  for(size_t ijet(0); ijet < b.ak8jets_pt()->size(); ijet++) { 
		if(b.ak8jets_pt()->at(ijet) > 300) {
			for(int ilep(0); ilep < b.nleps(); ilep++) {
				if(b.leps_id()->at(ilep) == 1) {
					dphi = abs(b.ak8jets_phi()->at(ijet) - b.leps_phi()->at(ilep));
					if      (dphi > 6.2830) dphi -= 6.2830;
					else if (dphi > 3.1415) dphi = 6.2830 - dphi;
					deta = b.ak8jets_eta()->at(ijet) - b.leps_eta()->at(ilep);
					if(sqrt(pow(deta,2)+pow(dphi,2)) < deltaR || deltaR == -1) deltaR = sqrt(pow(deta,2)+pow(dphi,2));
					}
				}
			}
		}
	return deltaR;
}

