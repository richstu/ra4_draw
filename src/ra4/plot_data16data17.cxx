#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string.h>

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

string cuts(std::string selection, std::string var = "") {
  string cut("");
  if(strstr(selection.c_str(),"standard") != NULL) {
	if(var == "miniso" || var == "reliso" || var == "mus_pt" || var == "els_pt") cut = "st>500&&met>200&&njets>=4";
	else if(var == "njets") cut = "nleps>=1&&st>500&&met>200";
	else cut = "nleps>=1&&st>500&&met>200&&njets>=4";
	}
  else if(strstr(selection.c_str(),"1Lepton") != NULL) {
  	if(var == "njets")                               cut = "nleps==1&&st>500&&met>200&&(mt<140||mj14<400) && nbmReReco >= 1 && nveto==0";
		else if(var == "nleps")                          cut = "st>500&&met>200&&(mt<140||mj14<400)&&njets>=6 && nbmReReco >= 1 && nveto==0";
  	else if(var == "nbm") cut = "nleps==1&&st>500&&met>200&&(mt<140||mj14<400)&&njets>=6&&nveto==0";
  	else if(var == "mt")   											     cut = "nleps==1&&st>500&&mj14<400&&met>200&&njets>=6 && nbmReReco >= 1 && nveto==0";
		else if(var == "mj14") 											     cut = "nleps==1&&st>500&&met>200&&(mt<140)&&njets>=6 && nbmReReco >= 1 && nveto==0";
		else if(var == "miniso" || var == "reliso")      cut = "st>500&&met>200&&(mt<140||mj14<400)&&njets>=6 && nbmReReco >= 1";
		else if(var == "nmus" || var == "nels") cut = "nleps>0&&st>500&&met>200&&(mt<140||mj14<400)&&njets>=6 && nbmReReco >= 1 && nveto==0";
		else if(var == "jets_csv" || var == "jets_csvd") cut = "nleps==1&&st>500&&met>200&&(mt<140||mj14<400)&&njets>=6&&nveto==0";
		else 																	 cut = "nleps==1&&st>500&&met>200&&(mt<140||mj14<400)&&njets>=6 && nbmReReco >= 1 && nveto==0";
  	}
  else if(strstr(selection.c_str(),"2Lepton") != NULL) {
		if(var == "njets") cut = "st>500&&nleps==2&&nbm<=2&&met>200&&met<500";
		else if(var == "nbm") cut = "st>500&&nleps==2&&njets>=4&&met>200&&met<500";
		else if(var == "met") cut = "st>500&&nleps==2&&njets>=4&&nbm<=2&&met>200";
		else cut = "st>500&&nleps==2&&njets>=4&&nbm<=2&&met>200&&met<500";
	}
  return cut;
}

//Run for RunB+RunC for mt plots
int main(){
  gErrorIgnoreLevel = 6000;
  string data17_cuts(""), data17new_cuts(""), data17old_cuts("");
  string data17_Prompt_path(       "/net/cms2/cms2r0/babymaker/babies/2017_10_20/data/merged_database_standard/*.root");
  string data17_Prompt_met100_path("/net/cms2/cms2r0/babymaker/babies/2017_10_20/data/merged_database_met100/*.root");
  string data17_ReReco_path(       "/net/cms2/cms2r0/babymaker/babies/2018_01_30/data/merged_database_standard/*.root");
  string data17_ReReco_met100_path("/net/cms2/cms2r0/babymaker/babies/2018_01_30/data/merged_database_met100/*.root");
  string data17_ReReco_nl1nj3_path("/net/cms2/cms2r0/babymaker/babies/2018_01_30/data/merged_database_n*nj3/*.root");
//   double lumi17(0), lumi16(0);
//If running on cms#, need to add /net/cms2/ to data_path names
//   lumi16 = 18.83;
  string data16_met100_path = "/cms2r0/babymaker/babies/2017_02_14/data/merged_database_met100/*.root";
  string data16_standard_path = "/net/cms2/cms2r0/babymaker/babies/2017_02_14/data/merged_database_standard/*.root";
  string data16_cuts = "(met/met_calo<5.0)&&pass&&trig_ra4&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell&&pass_cschalo)";
//   lumi17 = 32.5;
  data17_cuts = "(met/met_calo<5.0)&&(trig[3]||trig[7]||trig[9]||trig[15]||trig[20]||trig[21]||trig[23]||trig[24])&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell&&pass_cschalo)&&pass";
//  mc16_path = "/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/*.root";

	string data18_path("/net/cms2/cms2r0/babymaker/babies/2018_08_04/data/merged_database_standard/*.root");

  
  Palette colors("txt/colors.txt", "default");
  //
  // 2016 data
  //
	//Standard Selection
  auto data_16 = Process::MakeShared<Baby_full>("2016 Data - 36.2 ifb", Process::Type::background, kBlack,
    {data16_standard_path},data16_cuts);
  data_16->SetFillColor(kWhite);
  data_16->SetLineColor(kAzure-2);//kBlue-7);
  data_16->SetLineWidth(3);
	//met > 100 Selection
  auto data_16_met100 = Process::MakeShared<Baby_full>("2016 Data - 36.2 ifb", Process::Type::background, kBlack,
    {data16_met100_path},data16_cuts);
  data_16_met100->SetFillColor(kWhite);
  data_16_met100->SetLineColor(kAzure-2);//kBlue-7);
  data_16_met100->SetLineWidth(3);
  //
  // 2017 data
  //
  //Prompt Reco
	//Standard Selection
  auto data_17_Prompt = Process::MakeShared<Baby_full>("2017 Prompt Reco Data - 18.8 ifb", Process::Type::background, kBlack,
    {data17_Prompt_path},data17_cuts);
  data_17_Prompt->SetFillColor(kWhite);
  data_17_Prompt->SetLineColor(kAzure-2);//kBlue-7);
  data_17_Prompt->SetLineWidth(3);
	//Standard Selection
  auto data_17_Prompt_runs_BD = Process::MakeShared<Baby_full>("2017 Prompt Reco Data - 18.8 ifb", Process::Type::background, kBlack,
    {data17_Prompt_path},data17_cuts+"&&run<302663");
  data_17_Prompt_runs_BD->SetFillColor(kWhite);
  data_17_Prompt_runs_BD->SetLineColor(kAzure-2);//kBlue-7);
  data_17_Prompt_runs_BD->SetLineWidth(3);
	//met100 Selection
  auto data_17_Prompt_met100 = Process::MakeShared<Baby_full>("2017 Prompt Reco Data - 18.8 ifb", Process::Type::background, kBlack,
    {data17_Prompt_met100_path},data17_cuts);
  data_17_Prompt_met100->SetFillColor(kWhite);
  data_17_Prompt_met100->SetLineColor(kAzure-2);//kBlue-7);
  data_17_Prompt_met100->SetLineWidth(3);
  
  // 2017 ReReco data
	//Standard Selection
  auto data_17_ReReco = Process::MakeShared<Baby_full>("2017 ReReco Data - 42.0 ifb", Process::Type::data, kBlack,
    {data17_ReReco_path},data17_cuts);
  // 2017 ReReco data - Runs B-D
  auto data_17_ReReco_runs_BD_data = Process::MakeShared<Baby_full>("2017 ReReco Runs B-D - 18.8 ifb", Process::Type::data, kBlack,
    {data17_ReReco_path},data17_cuts+"&&run<303434");
  auto data_17_ReReco_runs_BD = Process::MakeShared<Baby_full>("2017 ReReco Runs B-D - 18.8 ifb", Process::Type::background, kBlack,
    {data17_ReReco_path},data17_cuts+"&&run<303434");
  data_17_ReReco_runs_BD->SetFillColor(kWhite);
  data_17_ReReco_runs_BD->SetLineColor(kAzure-2);//kBlue-7);
  data_17_ReReco_runs_BD->SetLineWidth(3);
  // 2017 ReReco data - Run F
  auto data_17_ReReco_run_F = Process::MakeShared<Baby_full>("2017 ReReco Run F - 13.7 ifb", Process::Type::data, kBlack,
    {data17_ReReco_path},data17_cuts+"&&run>305040");
  // 2017 ReReco - met100 sample
  auto data_17_ReReco_met100 = Process::MakeShared<Baby_full>("2017 ReReco Data - 42.0 ifb", Process::Type::data, kBlack,
    {data17_ReReco_met100_path},data17_cuts);
  // 2017 ReReco data - Runs B-D met100 sample
  auto data_17_ReReco_met100_runs_BD = Process::MakeShared<Baby_full>("2017 ReReco Runs B-D - 18.8 ifb", Process::Type::background, kBlack,
    {data17_ReReco_met100_path},data17_cuts+"&&run<303434");
  data_17_ReReco_met100_runs_BD->SetFillColor(kWhite);
  data_17_ReReco_met100_runs_BD->SetLineColor(kAzure-2);//kBlue-7);
  data_17_ReReco_met100_runs_BD->SetLineWidth(3);
  // 2017 ReReco data - Run F met100
  auto data_17_ReReco_met100_run_F = Process::MakeShared<Baby_full>("2017 ReReco Run F - 13.7 ifb", Process::Type::data, kBlack,
    {data17_ReReco_met100_path},data17_cuts+"&&run>305040");
	
	// 2018 data
	auto data_18 = Process::MakeShared<Baby_full>("2018 Data - 20.0 ifb", Process::Type::data, kBlack, {data18_path},data17_cuts);


//  auto mc_16 = Process::MakeShared<Baby_full>("2016 MC", Process::Type::background, kBlack,
  //  {mc16_path},mc16_cuts);

  vector<shared_ptr<Process> > full_data_16_17_compare = {data_17_ReReco, data_16};
  vector<shared_ptr<Process> > full_data_ReReco_BD_F_compare = {data_17_ReReco_runs_BD, data_17_ReReco_run_F};
  vector<shared_ptr<Process> > full_data_Prompt_ReReco_compare = {data_17_ReReco_runs_BD_data, data_17_Prompt_runs_BD};
  vector<shared_ptr<Process> > full_data_16_17_met100_compare = {data_17_ReReco_met100, data_16_met100};
  vector<shared_ptr<Process> > full_data_ReReco_BD_F_met100_compare = {data_17_ReReco_met100_run_F,data_17_ReReco_met100_runs_BD};
  vector<shared_ptr<Process> > full_data_Prompt_ReReco_met100_compare = {data_17_ReReco_met100, data_17_Prompt_met100};

	vector<shared_ptr<Process> > data_16_18_compare = {data_16, data_18};
  

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info).PrintVals(false);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info).PrintVals(false);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info};
  vector<PlotOpt> lin_plot = {lin_lumi_info};
  vector<PlotOpt> log_plot = {log_lumi_info};
  string bjets_16_17("&&((run<297047&&jets_csv>.8484)||(run>297047&&jets_csv>.8838))");
  string bjets("&&jets_csv>.8838");
  string standard("nleps>=1&&st>500&&met>100");
  string single_lepton("nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm>=1&&nveto==0");
  string dilepton("st>500&&nleps==2&&njets>=4&&nbm<=2&&met>200&&met<500");
  string cut("2Lepton");
  string sample("2016_2017");
  vector<shared_ptr<Process> > data_samples;
  if(strstr(sample.c_str(),"2016_2017") != NULL) {
	if(strstr(cut.c_str(),"met100") != NULL) data_samples = full_data_16_17_met100_compare;
	else data_samples = full_data_16_17_compare;
	bjets = bjets_16_17;
	}
  else if(strstr(sample.c_str(),"2017_Prompt_ReReco") != NULL) {
	if(strstr(cut.c_str(),"met100") != NULL) data_samples = full_data_Prompt_ReReco_met100_compare;
	else data_samples = full_data_Prompt_ReReco_compare;
	}
  else if(strstr(sample.c_str(),"2017-ReReco_runsBD_runF") != NULL) {
	if(strstr(cut.c_str(),"met100") != NULL) data_samples = full_data_ReReco_BD_F_met100_compare;
	else data_samples = full_data_ReReco_BD_F_compare;
	}
  PlotMaker pm;
//  pm.Push<Hist1D>(Axis(100, 0., 100., "npv", "Number of Vertices", {}),"met>100&&nleps>=1",full_data_16_17_met100_compare, all_plot_types).Weight(Functions::npv_weights_16_to_17);
//  pm.Push<Hist1D>(Axis(50, 0., 100., "npv", "Number of Vertices", {}), "met>100",full_data_16_17_met100_compare, lin_plot);
//  pm.Push<Hist1D>(Axis(50, 0., 100., "npv", "Number of Vertices", {}), "met>100",full_data_16_17_met100_compare, lin_plot);
//  pm.Push<Hist1D>(Axis(50, 0., 100., "npv", "Number of Vertices", {}), "met>100&&nleps>=1",full_data_ReReco_BD_F_met100_compare, lin_plot);
//  pm.Push<Hist1D>(Axis(50, 0., 100., "npv", "Number of Vertices", {}), "met>100",full_data_ReReco_BD_F_met100_compare, lin_plot);
/*
//Pileup
  pm.Push<Hist1D>(Axis(50, 0., 100., "npv", "Number of  Vertices", {}), "met>100&&nleps>=1",data_samples, lin_plot);

//Electrons
  pm.Push<Hist1D>(Axis(32,-3.2,3.2, "els_sceta", "#eta_{e}", {}), cuts(cut), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(34,-3.4,3.4, "els_phi", "#phi_{e}", {}), cuts(cut), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(28,0,140, "els_pt", "Electron p_{T} [GeV]", {10}), cuts(cut), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nels", "N_{e}", {0.5}), cuts(cut,"nels"), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(20,0,0.1, "els_miniso*(els_sigid==1)*(els_pt>20)*(els_sceta<2.5&&els_sceta>-2.5)", "Electron Miniso", {0.1}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(28,10,150,"els_pt*(els_sigid==1)*(els_miniso<0.1)*(els_sceta<2.5&&els_sceta>-2.5)", "Electron p_{T} [GeV]", {20}), cuts(cut,"miniso"), data_samples, log_plot);
//Muons
  pm.Push<Hist1D>(Axis(26,-2.6,2.6, "mus_eta", "#eta_{#mu}", {}), cuts(cut), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(34,-3.4,3.4, "mus_phi", "#phi_{#mu}", {}), cuts(cut), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(28,0,140, "mus_pt", "Muon p_{T} [GeV]", {10}), cuts(cut), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nmus", "N_{#mu}", {0.5}), cuts(cut,"nmus"), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(20,0,0.2, "mus_miniso*(mus_sigid==1)*(mus_pt>20)*(mus_eta<2.4&&mus_eta>-2.4)", "Muon Miniso", {0.2}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(28,10,150,"mus_pt*(mus_sigid==1)*(mus_miniso<0.2)*(mus_eta<2.4&&mus_eta>-2.4)", "Muon p_{T} [GeV]", {20}), cuts(cut,"miniso"), data_samples, log_plot);
//Jets
  pm.Push<Hist1D>(Axis(30,-3,3, "jets_eta", "#eta_{jet}", {}), cuts(cut), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(34,-3.4,3.4, "jets_phi", "#phi_{jet}", {}), cuts(cut), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(20,0,1000, "jets_pt", "jet p_{T} [GeV]", {20}), cuts(cut), data_samples, log_plot);
//B-jets
  pm.Push<Hist1D>(Axis(30,-3,3, "jets_eta", "#eta_{b-jet}", {}), cuts(cut)+bjets, data_samples, log_plot);
  pm.Push<Hist1D>(Axis(34,-3.4,3.4, "jets_phi", "#phi_{b-jet}", {}), cuts(cut)+bjets, data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(20,0,1000, "jets_pt", "b-jet p_{T} [GeV]", {20}), cuts(cut)+bjets, data_samples, log_plot);
  pm.Push<Hist1D>(Axis(6,-0.5,5.5, "nbm", "N_{b}", {0.5}), cuts(cut,"nbm"), data_samples, lin_plot);
  pm.Push<Hist1D>(Axis(25,0,1, "jets_csv", "Jet CSV", {.8484,.8838}), cuts(cut,"jets_csv"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(25,0,1, "jets_csvd", "Jet Deep CSV", {.4942,.6324}), cuts(cut,"jets_csvd"), data_samples, log_plot);
	*/
//SUSY vars
	cut = "1Lepton";
  pm.Push<Hist1D>(Axis(10,200,700, "met", "E_{T}^{miss} [GeV]", {200,500}), "nleps == 1 && st > 500 && met > 200 && (mt < 140 || mj14 < 400) && njets >=6 && nbmReReco >= 1 && nveto ==0", data_16_18_compare, log_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(12,-0.5,11.5, "njets", "N_{jets}", {5.5}), cuts(cut,"njets"), data_16_18_compare, lin_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(15,500,2000, "st", "S_{T} [GeV]", {}),     cuts(cut,"st"),    data_16_18_compare, log_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(17,0,340, "mt", "m_{T} [GeV]", {140}),     cuts(cut,"mt"),    data_16_18_compare, log_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(15,0,1500, "mj14", "M_{J} [GeV]", {400}),  cuts(cut,"mj14"),  data_16_18_compare, log_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(7,-0.5,6.5, "nbmReReco", "N_{b}", {0.5}),        cuts(cut,"nbm"),   data_16_18_compare, lin_plot).Tag("2018");
  pm.Push<Hist1D>(Axis(6,-0.5,5.5, "nleps", "N_{leps}", {0.5}),   cuts(cut,"nleps"), data_16_18_compare, lin_plot).Tag("2018");

  pm.Push<Hist1D>(Axis(10,200,700, "met", "E_{T}^{miss} [GeV]", {200,500}), "nleps == 1 && st > 500 && met > 200 && (mt < 140 || mj14 < 400) && njets >=6 && nbm >=1 && nveto ==0", 
																																										 full_data_16_17_compare, log_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(12,-0.5,11.5, "njets", "N_{jets}", {5.5}), cuts(cut,"njets"), full_data_16_17_compare, lin_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(15,500,2000, "st", "S_{T} [GeV]", {}),     cuts(cut,"st"),    full_data_16_17_compare, log_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(17,0,340, "mt", "m_{T} [GeV]", {140}),     cuts(cut,"mt"),    full_data_16_17_compare, log_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(15,0,1500, "mj14", "M_{J} [GeV]", {400}),  cuts(cut,"mj14"),  full_data_16_17_compare, log_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(7,-0.5,6.5, "nbmReReco", "N_{b}", {0.5}),        cuts(cut,"nbm"),   full_data_16_17_compare, lin_plot).Tag("2017");
  pm.Push<Hist1D>(Axis(6,-0.5,5.5, "nleps", "N_{leps}", {0.5}),   cuts(cut,"nleps"), full_data_16_17_compare, lin_plot).Tag("2017");
	/*
//Plots that need met100 sample
  //Electrons
  pm.Push<Hist1D>(Axis(28,10,150,"els_pt*(els_sigid==1)*(els_miniso<0.1)*(els_sceta<2.5&&els_sceta>-2.5)", "Electron p_{T}", {20}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(50,0,1, "els_miniso*(els_sigid==1)*(els_pt>20)*(els_sceta<2.5&&els_sceta>-2.5)", "Electron Miniso", {0.1}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(50,0,1, "el_miniso", "Electron Miniso", {0.1}), cuts(cut,"miniso"), data_samples, log_plot);  
//  pm.Push<Hist1D>(Axis(50,0,1, "els_reliso*(els_sigid==1)*(els_pt>20)", "Electron Reliso", {}), cuts(cut,"miniso"), data_samples, log_plot);
  //Muons
  pm.Push<Hist1D>(Axis(28,10,150,"mus_pt*(mus_sigid==1)*(mus_miniso<0.2)*(mus_eta<2.4&&mus_eta>-2.4)", "Muon p_{T}", {20}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(50,0,1, "mus_miniso*(mus_sigid==1)*(mus_pt>20)*(mus_eta<2.4&&mus_eta>-2.4)", "Muon Miniso", {0.2}), cuts(cut,"miniso"), data_samples, log_plot);
  pm.Push<Hist1D>(Axis(50,0,1, "mu_miniso", "Muon Miniso", {0.2}), cuts(cut,"miniso"), data_samples, log_plot);
//  pm.Push<Hist1D>(Axis(50,0,1, "mus_reliso*(mus_sigid==1)*(mus_pt>20)", "Muon Reliso", {}), cuts(cut,"miniso"), data_samples, log_plot);
*/
  pm.min_print_=true;
  pm.MakePlots(42.0);

}

