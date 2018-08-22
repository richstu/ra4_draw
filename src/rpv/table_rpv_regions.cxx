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

  string c_ps = "pass && stitch_met";
  // Background samples
  auto tt_1l = Process::MakeShared<Baby_full>("t#bar{t} #, (1l)", Process::Type::background, colors("tt_1l"),
    {mc_folder_new+"*TTJets*SingleLept*"},
    c_ps+"&&ntruleps>=1");

  auto tt_2l = Process::MakeShared<Baby_full>("t#bar{t} #, (2l)", Process::Type::background, colors("tt_1l"),
    {mc_folder_new+"*TTJets*DiLept*"},
    c_ps+"&&ntruleps>=1");

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
  
  vector<shared_ptr<Process> > samples = {tt_1l, tt_2l, T1tttt1200_800, T1tttt1500_100, T1tttt2000_100};

  PlotMaker pm;
  pm.Push<Table>("rpv_regions", vector<TableRow>{
      TableRow("Baseline, $m_{T}>140$ GeV, $M_{J}>400$ GeV", "mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&nveto==0&&njets>=6&&nbm>=1",0,0,"weight"),
	TableRow("Baseline, $m_{T}>140$ GeV, $M_{J}>400$ GeV, 1 T top", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&nveto==0&&njets>=6&&nbm>=1",0,0,"weight"),
	TableRow("$200 < E_{T}^{miss}\\leq 350 GeV$"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=6&&njets<=8&&nbm==1",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=9&&nbm==1",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=6&&njets<=8&&nbm==2",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=9&&nbm==2",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=6&&njets<=8&&nbm>=3",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>200&&met<=350&&nveto==0&&njets>=9&&nbm>=3",0,0,"weight"),
	TableRow("$350 < E_{T}^{miss}\\leq 500 GeV$"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=6&&njets<=8&&nbm==1",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=9&&nbm==1",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=6&&njets<=8&&nbm==2",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=9&&nbm==2",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=6&&njets<=8&&nbm>=3",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>350&&met<=500&&nveto==0&&njets>=9&&nbm>=3",0,0,"weight"),
	TableRow("$E_{T}^{miss} > 500 GeV$"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=6&&njets<=8&&nbm==1",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=1$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=9&&nbm==1",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=6&&njets<=8&&nbm==2",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}=2$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=9&&nbm==2",0,0,"weight"),
	TableRow("$6\\leq N_{jets}\\leq 8, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=6&&njets<=8&&nbm>=3",0,0,"weight"),
	TableRow("$N_{jets}\\geq 9, N_{b}\\geq 3$", "ntop_tight_decor>=1&&mt>140&&mj14>400&&st>500&&nleps==1&&met>500&&nveto==0&&njets>=9&&nbm>=3",0,0,"weight"),
	},samples,true);
  
  pm.min_print_ = true;
  pm.MakePlots(lumi);
}
