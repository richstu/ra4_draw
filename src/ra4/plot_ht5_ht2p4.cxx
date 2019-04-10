#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string.h>

#include "TMath.h"
#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
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
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/skim_met150/";
  

  NamedFunc baseline = "met>300";
  baseline = baseline && Functions::hem_veto && "st<10000 && pass_ra2_badmu";// && met/met_calo<5";

  NamedFunc trigs = Functions::trig_run2;

  NamedFunc data_base = baseline && trigs;
  NamedFunc mc_base = baseline && "stitch_met";

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
	                 {folderdata[2016]+"*.root"},data_base);
	auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
	                 {folderdata[2017]+"*.root"},data_base);
  auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
                   {folderdata[2018]+"*.root"},data_base);
	auto data_2018A = Process::MakeShared<Baby_full>("2018 Data",back,kBlack,
	                 {folderdata[2018]+"*Run2018A_0_Run2018*.root"},data_base);
  auto data_2018D = Process::MakeShared<Baby_full>("2018 Data",back,kBlack,
                   {folderdata[2018]+"*Run2018D_0_Run2018*.root"},data_base);

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
  vector<shared_ptr<Process> > data18A  = {data_2018A};
	vector<shared_ptr<Process> > data18D  = {data_2018D};
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

  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style2D().Stack(StackType::data_norm).Title(TitleType::preliminary)};
  vector<PlotOpt> bkg_pts = {style2D().Stack(StackType::lumi_shapes).Title(TitleType::info)};

	vector<string> tag = {"2016", "2017", "2018"};



  vector<string> cuts;

  cuts.push_back("nleps>=1 && met>100 && st>500 && njets>=5 && nbd>=1");

  const NamedFunc dphi_leadje24_met("dphi_leadje24_met", [](const Baby &b) ->NamedFunc::ScalarType{
    double dphi = -1.;
    if (b.jets_pt()->size()>0) {
      dphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(0)));
    }
    return dphi;
  });

  const NamedFunc dphi_leadje5_met("dphi_leadje5_met", [](const Baby &b) ->NamedFunc::ScalarType{
    double dphi = -1.;
    if (b.jets_pt()->size()>0) {
      if (b.ejets_pt()->size()>0) {
        if (b.jets_pt()->at(0)>b.ejets_pt()->at(0)) 
          dphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(0)));
        else 
          dphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.ejets_phi()->at(0)));
      } else {
        dphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(0)));
      }
    } else if (b.ejets_pt()->size()>0) {
      dphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.ejets_phi()->at(0)));
    }
    return dphi;
  });

  PlotMaker pm;
  pm.Push<Hist2D>(Axis(50, 0, 3.3, dphi_leadje24_met,  "#Delta#phi(j_{1}^{2.4}, MET)",{}), 
                  Axis(50, 1., 3., "ht_ejets/ht",  "H_{T}^{5}/H_{T}^{2.4}",{}), 
                  "met>300",
                  data18A, bkg_hist).Weight("1").Tag("Run2018A");
  pm.Push<Hist2D>(Axis(50, 0, 3.3, dphi_leadje5_met,  "#Delta#phi(j_{1}^{5}, MET)",{}), 
                  Axis(50, 1., 3., "ht_ejets/ht",  "H_{T}^{5}/H_{T}^{2.4}",{}), 
                  "met>300",
                  data18A, bkg_hist).Weight("1").Tag("Run2018A");
  pm.Push<Hist2D>(Axis(50, 0, 3.3, dphi_leadje24_met,  "#Delta#phi(j_{1}^{2.4}, MET)",{}), 
                  Axis(50, 1., 3., "ht_ejets/ht",  "H_{T}^{5}/H_{T}^{2.4}",{}), 
                  "met>300",
                  data18D, bkg_hist).Weight("1").Tag("Run2018D");
  pm.Push<Hist2D>(Axis(50, 0, 3.3, dphi_leadje5_met,  "#Delta#phi(j_{1}^{5}, MET)",{}), 
                  Axis(50, 1., 3., "ht_ejets/ht",  "H_{T}^{5}/H_{T}^{2.4}",{}), 
                  "met>300",
                  data18D, bkg_hist).Weight("1").Tag("Run2018D");

  pm.min_print_=true;
  pm.MakePlots(1);
}

