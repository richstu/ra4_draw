#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 135.0;

  string mc_folder_old = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/"; // 2016 BG MC
  string sig_folder_old = "/net/cms2/cms2r0/babymaker/babies/2017_02_07/T1tttt/unskimmed/"; // 2016 Signal MC
  string mc_folder_new = "/net/cms2/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/";
  string sig_folder_new = "/net/cms2/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/"; // same as new BG folder
 
  Palette colors("txt/colors.txt", "default");

  NamedFunc nhighpt_jets("nhighpt_jets", [](const Baby &b) ->NamedFunc::ScalarType{
      int _njets = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
        if(b.ak8jets_pt()->at(ijet)>300) _njets++;
      } 
      return _njets;
    });

  NamedFunc ntopmass_jets("ntopmass_jets", [](const Baby &b) ->NamedFunc::ScalarType{
      int _njets = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
        if(b.ak8jets_m()->at(ijet)>105. && b.ak8jets_m()->at(ijet)<210. && b.ak8jets_pt()->at(ijet)>300.) _njets++;
      } 
      return _njets;
    });

  NamedFunc ntop_loose_nom_raw("ntop_loose_nom_raw", [](const Baby &b) ->NamedFunc::ScalarType{
      int _ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
        if(b.ak8jets_pt()->at(ijet)>300. && b.ak8jets_nom_raw_top()->at(ijet)>0.4) _ntop++;
      } 
      return _ntop;
    });


  string c_ps = "pass && stitch_met";
  // Background samples
  auto tt_1l = Process::MakeShared<Baby_full>("t#bar{t} #, (1l)", Process::Type::background, colors("tt_1l"),
    {mc_folder_new+"*TTJets*SingleLept*"},
    c_ps+"&&ntruleps==1");

  auto tt_2l = Process::MakeShared<Baby_full>("t#bar{t} #, (2l)", Process::Type::background, colors("tt_1l"),
    {mc_folder_new+"*TTJets*DiLept*"},
    c_ps+"&&ntruleps==2");

  //Signal Samples
  // (1800,100) and (1400,1000) for comparison to paper
  auto T1tttt1800_100 = Process::MakeShared<Baby_full>("T1tttt (1800,100)", Process::Type::signal, colors("tt_1l"),
    {sig_folder_old+"*T1tttt*mGluino-1800_mLSP-100_*"});

  auto T1tttt1400_1000 = Process::MakeShared<Baby_full>("T1tttt (1400,1000)", Process::Type::signal, colors("tt_1l"),
    {sig_folder_old+"*T1tttt*mGluino-1400_mLSP-1000_*"});

  // (1200,800), (1500,100), (2000,100) working points for new studies
  auto T1tttt1200_800 = Process::MakeShared<Baby_full>("T1tttt (1200,800)", Process::Type::signal, colors("tt_1l"),
    {sig_folder_new+"*T1tttt*mGluino-1200_mLSP-800_*"});

  auto T1tttt1500_100 = Process::MakeShared<Baby_full>("T1tttt (1500,100)", Process::Type::signal, colors("tt_1l"),
    {sig_folder_new+"*T1tttt*mGluino-1500_mLSP-100_*"});

  auto T1tttt2000_100 = Process::MakeShared<Baby_full>("T1tttt (2000,100)", Process::Type::signal, colors("tt_1l"),
    {sig_folder_new+"*T1tttt*mGluino-2000_mLSP-100_*"});
  
  vector<shared_ptr<Process> > samples = {tt_1l, tt_2l, T1tttt1200_800, T1tttt2000_100};

  string baseline_s = "mt>140&&mj14>400&&nleps==1&&st>500&&met>200&&nveto==0&&njets>=6&&nbm>=1&&pass_ra2_badmu&&met/met_calo<5&&nak8jets>=1";

  PlotMaker pm;
  pm.Push<Table>("rpv_regions", vector<TableRow>{
      TableRow("Baseline", baseline_s,0,0,"weight"),
	TableRow("Baseline + 1 Leading ak8 jet$ p_{T} > 300$ GeV, mass$ > 105$ GeV", baseline_s && "ak8jets_pt[0]>300&&ak8jets_m[0]>105" ,0,0,"weight"),
	TableRow("$200$ GeV$ < E_{T}^{miss} \\leq 350$ GeV"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=1$", baseline_s && "met<=350&&njets<=8&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=1$", baseline_s && "met<=350&&njets>=9&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=2$", baseline_s && "met<=350&&njets<=8&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=2$", baseline_s && "met<=350&&njets>=9&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b} \\geq 3$", baseline_s && "met<=350&&njets<=8&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b} \\geq 3$", baseline_s && "met<=350&&njets>=9&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$350$ GeV$ < E_{T}^{miss} \\leq 500$ GeV"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=1$", baseline_s && "met>350&&met<=500&&njets<=8&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=1$", baseline_s && "met>350&&met<=500&&njets>=9&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=2$", baseline_s && "met>350&&met<=500&&njets<=8&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=2$", baseline_s && "met>350&&met<=500&&njets>=9&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b} \\geq 3$", baseline_s && "met>350&&met<=500&&njets<=8&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b} \\geq 3$", baseline_s && "met>350&&met<=500&&njets>=9&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$E_{T}^{miss} > 500$ GeV"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=1$", baseline_s && "met>500&&njets<=8&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=1$", baseline_s && "met>500&&njets>=9&&nbm==1&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b}=2$", baseline_s && "met>500&&njets<=8&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b}=2$", baseline_s && "met>500&&njets>=9&&nbm==2&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$6 \\leq N_{jets} \\leq 8, N_{b} \\geq 3$", baseline_s && "met>500&&njets<=8&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	TableRow("$N_{jets} \\geq 9, N_{b} \\geq 3$", baseline_s && "met>500&&njets>=9&&nbm>=3&&ak8jets_pt[0]>300&&ak8jets_m[0]>105",0,0,"weight"),
	},samples,true);
  
  pm.min_print_ = true;
  pm.MakePlots(lumi);
}
