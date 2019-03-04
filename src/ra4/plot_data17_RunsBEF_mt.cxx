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

  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  // Data
  string data17_path(bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/");
  string mc17_path(bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/");
  string q_cuts("pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");

  auto data_2017_runBE = Process::MakeShared<Baby_full>("2017 Data - Runs B-E",data,kBlack,
                        {data17_path+"*Run2017B*.root",
                         data17_path+"*Run2017C*.root",
                         data17_path+"*Run2017D*.root",
                         data17_path+"*Run2017E*.root"}, Functions::trig_run2 && "run<=304797" && q_cuts);
  auto data_2017_runF = Process::MakeShared<Baby_full>("2017 Data - Run F",data,kBlack,
                        {data17_path+"*Run2017F*.root"}, Functions::trig_run2 && "run>=305040" && q_cuts);
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

  vector<shared_ptr<Process> > data17BE_mc17  = {data_2017_runBE, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
  vector<shared_ptr<Process> > data17F_mc17  = {data_2017_runF, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
  .YAxis(YAxisType::log)
  .Stack(StackType::data_norm)
  .Bottom(BottomType::ratio);
  vector<PlotOpt> log_stack = {log_lumi};

  NamedFunc baseline = "nleps == 1 && st > 500 && met > 200 && njets >= 5 && nbd >= 1 && nveto == 0 && mj14<400";
  const NamedFunc wgtBE("wgtBE", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;
    return b.weight()*b.w_prefire();
  });
  const NamedFunc wgtF("wgtF", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;
    return b.weight()*b.w_prefire()*Functions::wnpv2017(b);
  });

  PlotMaker pm;
  pm.Push<Hist1D>(Axis(15,0,300, "mt", "m_{T} [GeV]",{}), 
    baseline, data17BE_mc17, log_stack).Weight(wgtBE).Tag("RunBE");
  pm.Push<Hist1D>(Axis(15,0,300, "mt", "m_{T} [GeV]",{}), 
    baseline, data17F_mc17, log_stack).Weight(wgtBE).Tag("RunF");
  pm.Push<Hist1D>(Axis(15,0,300, "mt", "m_{T} [GeV]",{}), 
    baseline, data17F_mc17, log_stack).Weight(wgtF).Tag("RunFnpv");


  pm.Push<Hist1D>(Axis(50,0,100, "npv", "NPV",{}), 
    baseline, data17BE_mc17, log_stack).Weight(wgtBE).Tag("RunBE");
  pm.Push<Hist1D>(Axis(50,0,100, "npv", "NPV",{}), 
    baseline, data17F_mc17, log_stack).Weight(wgtBE).Tag("RunF");
  pm.Push<Hist1D>(Axis(50,0,100, "npv", "NPV",{}), 
    baseline, data17F_mc17, log_stack).Weight(wgtF).Tag("RunFnpv");

  pm.min_print_=true;
  pm.MakePlots(41.5);

}




