#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string.h>

#include "TError.h"
#include "TColor.h"
#include "TVector.h"
#include "TStyle.h"
#include "TString.h"

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
#include "core/styles.hpp"

using namespace std;
using namespace PlotOptTypes;

NamedFunc BaselineCuts(string var = "", string extra = "") {
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
	if(extra != "") baseline = baseline && extra;
	return baseline;
}

// NamedFunc::ScalarType ntops_nom(const Baby &b, double top_wp = 0.1883);
void MakeRoc(vector<Axis> vars, shared_ptr<Process> background, shared_ptr<Process> signal, string tag = "", NamedFunc baseline = "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0");
TString tex_name(TString name);
void setTitles(TH1 *h, TString xTitle, TString yTitle, TString Left, TString Right);
TGraph MarkWorkingPoints(Axis axis, TH1D signal, TH1D background);
string FuncName(NamedFunc var);
// Axis FuncAxis(NamedFunc var, std:set<double> WPs = {}, double xmin = 0, double xmax = 1);

int main() {
	gErrorIgnoreLevel = 6000;
	const NamedFunc max_nom_bin("Highest_Nominal_Binarized",[&](const Baby &b) {
		double max(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_nom_bin_top()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300) max = score;
		}
		return max;
	});
	const NamedFunc max_nom_raw("Highest_Nominal_Raw",[&](const Baby &b){
		double max(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_nom_raw_top()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300) max = score;
		}
		return max;
	});
	const NamedFunc max_decor_bin("Highest_Decorrelated_Binarized",[&](const Baby &b){
		double max(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_decor_bin_top()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300) max = score;
		}
		return max;
	});
	const NamedFunc max_decor_raw("Highest_Decorrelated_Raw",[&](const Baby &b){
		double max(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_decor_raw_top()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300) max = score;
		}
		return max;
	});
	const NamedFunc max_decor_bin_wM("Highest_Decorrelated_Binarized_w/_Mass_window",[&](const Baby &b){
		double max(0.), score(0), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_decor_bin_top()->at(ijet);
			mass  = b.ak8jets_m()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300 && mass > 105 && mass < 210) max = score;
		}
		return max;
	});
	const NamedFunc max_decor_raw_wM("Highest_Decorrelated_Raw_w/_Mass_window",[&](const Baby &b){
		double max(0.), score(0), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_decor_raw_top()->at(ijet);
			mass  = b.ak8jets_m()->at(ijet);
			if(score > max && b.ak8jets_pt()->at(ijet) >= 300 && mass > 105 && mass < 210) max = score;
		}
		return max;
	});
	const NamedFunc sum2_nom_raw("Sum_of_2_highest_NR_scores",[&](const Baby &b){
		double max1(0.), max2(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_nom_raw_top()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 300);
			if(score > max1) {
				max2 = max1;
				max1 = score;
			}
			else if(score > max2) max2 = score;
		}
		return max1+max2;
	});
	const NamedFunc subleading_nom_raw("Subleading_NR_score",[&](const Baby &b){
		double max1(0.), max2(0.), score(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			score = b.ak8jets_nom_raw_top()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 300);
			if(score > max1) {
				max2 = max1;
				max1 = score;
			}
			else if(score > max2) max2 = score;
		}
		return max2;
	});
	const NamedFunc max_ak8_m("Highest_AK8_Mass",[&](const Baby &b){
		double max(0.), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			mass  = b.ak8jets_m()->at(ijet);
			if(mass > max && b.ak8jets_pt()->at(ijet) >= 300) max = mass;
		}
		return max;
	});
	const NamedFunc sum2_ak8_m("Sum_of_2_highest_AK8_masses",[&](const Baby &b){
		double max1(0.), max2(0.), mass(0);
		for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
			mass = b.ak8jets_m()->at(ijet)*(b.ak8jets_pt()->at(ijet) >= 300);
			if(mass > max1) {
				max2 = max1;
				max1 = mass;
			}
			else if(mass > max2) max2 = mass;
		}
		return max1+max2;
	});
  const NamedFunc is_300ak8("is_300ak8", [](const Baby &b) ->NamedFunc::ScalarType{
    int is300(0);
  	for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
	  	if(b.ak8jets_pt()->at(iak8) > 300) is300 = 1;
	  }
		return is300;
  	});

// Samples
	Palette colors("txt/colors.txt", "default");
	string mc_standard_path("/net/cms2/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
	string ttbar1L(mc_standard_path+"*_TTJets*SingleLept*.root");
	string ttbar2L(mc_standard_path+"*_TTJets*DiLept*.root");
	string ttbarHT(mc_standard_path+"*_TTJets*HT*.root");
	string  signal1(mc_standard_path+"*mGluino-1200_mLSP-800*.root");
	string  signal2(mc_standard_path+"*mGluino-2000_mLSP-100*.root");
  
// Processes
	auto proc_tt1l             = Process::MakeShared<Baby_full>("tt_1l",    Process::Type::background, colors("tt_1l"), {ttbar1L,ttbarHT},         "ntruleps<=1&&stitch_met" && is_300ak8);
	auto proc_tt2l             = Process::MakeShared<Baby_full>("tt_2l",    Process::Type::background, colors("tt_2l"), {ttbar2L,ttbarHT},         "ntruleps>=2&&stitch_met" && is_300ak8);
	auto proc_ttbar            = Process::MakeShared<Baby_full>("tt",       Process::Type::background, colors("tt_1l"), {ttbar1L,ttbar2L,ttbarHT}, "ntruleps>=2&&stitch_met" && is_300ak8);
	auto proc_T1tttt_1200_800 =  Process::MakeShared<Baby_full>("sig_12_8", Process::Type::signal,      kRed,           {signal1},                              "stitch_met" && is_300ak8);
	auto proc_T1tttt_2000_100  = Process::MakeShared<Baby_full>("sig_20_1", Process::Type::signal,     kAzure,          {signal2},                              "stitch_met" && is_300ak8);

// Selections

// Options
	PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
	vector<PlotOpt> log_lumis = {log_lumi};

// Setup
	vector<NamedFunc> max_scores  = { max_nom_bin, max_nom_raw, max_decor_bin, max_decor_raw};
	vector<string> leading_scores = {"nom_bin",   "nom_raw",   "decor_bin",   "decor_raw"   };
	vector<double> wp_vals        = { 0.1883,      0.4,         0.0474,        -1           };
	Axis max_nr    = Axis(100,0,1, max_nom_raw, FuncName(max_nom_raw), {0.4});
	Axis max_nb    = Axis(100,0,1, max_nom_bin, FuncName(max_nom_bin), {0.1883});
	Axis max_db    = Axis(100,0,1, max_decor_bin, FuncName(max_decor_bin),{});
	Axis max_dr    = Axis(100,0,1, max_decor_raw, FuncName(max_decor_raw),{});
	Axis max_db_wM = Axis(100,0,1, max_decor_bin_wM, FuncName(max_decor_bin_wM),{});
	Axis max_dr_wM = Axis(100,0,1, max_decor_raw_wM, FuncName(max_decor_raw_wM),{});
	Axis MJ        = Axis(100,250,1250, "mj14", "M_{J}", {});
	Axis top_nr    = Axis(100,0,1,"ak8jets_nom_raw_top[0]",  "Leading p_{T} Nominal Raw",{});
	Axis top_nb    = Axis(100,0,1,"ak8jets_nom_bin_top[0]",  "Leading p_{T} Nominal Raw",{});
	Axis top_db    = Axis(100,0,1,"ak8jets_decor_raw_top[0]","Leading p_{T} Nominal Raw",{});
	Axis top_dr    = Axis(100,0,1,"ak8jets_decor_bin_top[0]","Leading p_{T} Nominal Raw",{});
	Axis sum2_nr   = Axis(100,0,2,  sum2_nom_raw, FuncName(sum2_nom_raw), {1});
	Axis subl_nr   = Axis(100,0,1,  subleading_nom_raw, FuncName(subleading_nom_raw), {});
	Axis max_m     = Axis(100,0,1000, max_ak8_m,  FuncName(max_ak8_m),  {150});
	Axis sum2_m    = Axis(100,0,1500, sum2_ak8_m, FuncName(sum2_ak8_m), {300});
	vector<Axis> score_axes = {max_nr, max_nb, max_dr, max_db};
	vector<Axis> leading_score_axes = {top_nr, top_nb, top_dr, top_db};
	vector<Axis> mj_comp = {MJ, max_nr, top_nr};
	vector<Axis> nr = {max_nr, sum2_nr};
	vector<Axis> mj_nr = {MJ, max_nr};
	vector<Axis> mj_nr_sum2nr = {MJ, max_nr, sum2_nr};
	vector<Axis> ak8m_nr = {max_nr, max_m, sum2_m};
	// Non-Compressed comparison
		// 1l score curves
	MakeRoc(score_axes, proc_tt1l, proc_T1tttt_2000_100, "1l_max_scores");
	MakeRoc(leading_score_axes, proc_tt1l, proc_T1tttt_2000_100, "1l_lead_scores");
	MakeRoc(mj_comp, proc_tt1l, proc_T1tttt_2000_100, "1l_mj_scores");
		// 2l score curves
	MakeRoc(score_axes, proc_tt2l, proc_T1tttt_2000_100, "2l_max_scores",BaselineCuts("","mj14 > 600"));
	MakeRoc(leading_score_axes, proc_tt2l, proc_T1tttt_2000_100, "2l_lead_scores");
	MakeRoc(mj_comp, proc_tt2l, proc_T1tttt_2000_100, "2l_mj_scores");
		// 2l score vs MJ
	MakeRoc(mj_nr, proc_tt2l, proc_T1tttt_2000_100, "2l_mj_nr");
		// 2l score at high MJ
	MakeRoc(nr, proc_tt2l, proc_T1tttt_2000_100, "SUSY_Talk_2l_nr",BaselineCuts());
	MakeRoc(mj_nr, proc_tt2l, proc_T1tttt_2000_100, "SUSY_Talk_2l_nr_mj",BaselineCuts("","mj14 > 600"));
		// 2l AK8 mass at high MJ
	MakeRoc(ak8m_nr, proc_tt2l, proc_T1tttt_2000_100, "2l_highMJ_ak8m",BaselineCuts());
	return(0);
}

void MakeRoc(vector<Axis> vars, shared_ptr<Process> background, shared_ptr<Process> signal, string tag, NamedFunc baseline) {
	PlotOpt style("txt/plot_styles.txt", "CMSPaper"); setPlotStyle(style);
	TCanvas can;
	vector<Int_t> colz = {kRed,kBlue,kGreen-2,kMagenta, kViolet-8};
	vector<PlotOpt> styles = {style};
	TH1D base_histo("base","",1,0.03,1.0);
	base_histo.SetMinimum(0.0); base_histo.SetMaximum(1.0);
	base_histo.SetDirectory(0); base_histo.Draw();
  setTitles(&base_histo, tex_name(signal->name_)+" efficiency", tex_name(background->name_)+" efficiency","#scale[0.4]{"+CodeToRootTex(baseline.Name())+"}","");
	base_histo.GetYaxis()->SetTitleOffset(1.6);
	double legX = style.LeftMargin()+0.03, legY = 1-style.TopMargin()-0.02, legSingle = 0.05;
	double legW = 0.2, legH = legSingle*vars.size();
	TLegend leg(legX, legY-legH, legX+legW, legY);
	leg.SetTextSize(0.036); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextFont(42);
	PlotMaker pm;
	vector<shared_ptr<Process> > back = {background}, sig = {signal};
			for(size_t i = 0; i < vars.size(); i++) {
				pm.Push<Hist1D>(vars.at(i), baseline, back, styles);
				pm.Push<Hist1D>(vars.at(i), baseline, sig,  styles);
			}
	pm.min_print_=true;
	pm.MakePlots(1);
	TString leg_title;
	vector<TGraph> rocks;
	vector<TGraph> rocks_wps;
	for(size_t a = 0; a < vars.size(); a++) {
		// Background hist
		Hist1D * back_h1d = static_cast<Hist1D*>(pm.Figures()[2*a].get());
		Hist1D::SingleHist1D* back_sh = static_cast<Hist1D::SingleHist1D*>(back_h1d->GetComponent(back[0].get()));
		TH1D back_h = back_sh->scaled_hist_;
		// Signal hist
		Hist1D * sig_h1d = static_cast<Hist1D*>(pm.Figures()[2*a+1].get());
		Hist1D::SingleHist1D* sig_sh = static_cast<Hist1D::SingleHist1D*>(sig_h1d->GetComponent(sig[0].get()));
		TH1D sig_h = sig_sh->scaled_hist_;
		TGraph roc(0), wp(0);
		for(size_t p = 0; p < vars.at(a).Nbins(); p++) 
			roc.SetPoint(roc.GetN(),sig_h.Integral(p,-1)/sig_h.Integral(),back_h.Integral(p,-1)/back_h.Integral());
		wp = MarkWorkingPoints(vars[a], sig_h, back_h);
		roc.SetPoint(roc.GetN(),0,0);
		roc.SetLineColor(colz.at(a));
		if(vars.size() > 3) { 
			roc.SetLineWidth(3);
			wp.SetMarkerSize(2.4);
		}
		else {
			roc.SetLineWidth(4);
			wp.SetMarkerSize(3.4);
		}
		rocks.push_back(roc);
		wp.SetMarkerStyle(33);
		wp.SetMarkerColor(colz.at(a));
		rocks_wps.push_back(wp);
	}
	TString fname("plots/roc_"+tag+"_"+background->name_+"_"+signal->name_+"_"+CodeToPlainText(baseline.Name())+".pdf");
	for(size_t a = 0; a < vars.size(); a++) {
		rocks[a].Draw("Csame");
		leg_title = tex_name(vars.at(a).Title());
		leg.AddEntry(&(rocks[a]),leg_title,"l");
	}
	for(size_t a = 0; a < vars.size(); a++) 
		rocks_wps[a].Draw("Psame");
	leg.Draw("same");
	can.Print(fname);
	return;
}	

TString tex_name(TString name) {
	if     (name == "tt_1l") return "1-Lepton t#bar{t}";
	else if(name == "tt_2l") return "2-Lepton t#bar{t}";
	else if(name == "tt")    return "t#bar{t}";
	else if(name.Contains("sig")) return name.ReplaceAll("sig_","T1tttt(").ReplaceAll("_","00,") + "00)"; 
	if(name.Contains("nom"))
		name.ReplaceAll("nom_","Nominal ").ReplaceAll("bin","Binarized").ReplaceAll("raw","Raw");
	else if(name.Contains("decor"))
		name.ReplaceAll("decor_","Decorrelated ").ReplaceAll("bin","Binarized").ReplaceAll("raw","Raw");
	return name;
}

void setTitles(TH1 *h, TString xTitle, TString yTitle, TString Left, TString Right) {
	PlotOpt cms("txt/plot_styles.txt", "CMSPaper");
  if (0==h) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(xTitle); h->SetYTitle(yTitle);
    TLatex label; label.SetNDC(kTRUE);
    label.SetTextSize(0.06);
    label.SetTextAlign(11);
    label.DrawLatex(cms.LeftMargin(),1-cms.TopMargin()+0.02,Left);
    label.SetTextAlign(31);
    label.DrawLatex(1-cms.RightMargin(),1-cms.TopMargin()+0.02,Right);
  }
}

TGraph MarkWorkingPoints(Axis axis, TH1D signal, TH1D background) {
	const vector<double> bins = axis.Bins();
	int nbins = bins.size()-1;
	double xmin(bins[0]), xmax(bins[nbins]);
	int wpbin;
	TGraph wp(0);
	for(double cut: axis.cut_vals_) {
		if(cut > xmin && cut < xmax) {
			wpbin = round((cut - xmin)*nbins/(xmax-xmin));
			wp.SetPoint(wp.GetN(),signal.Integral(wpbin,-1)/signal.Integral(),background.Integral(wpbin,-1)/background.Integral());
		}
	}
	return wp;
}

string FuncName(NamedFunc var) {
	string name = var.Name();
	ReplaceAll(name,"_"," ");
	return name;
}
