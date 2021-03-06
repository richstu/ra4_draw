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
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
	Process::Type back = Process::Type::background;
	Process::Type sig = Process::Type::signal;

  Palette colors("txt/colors.txt", "default");
	bool old16(false);
	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_standard/");
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");
  string mc18_path("/net/cms2/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_standard/");
	NamedFunc q_cuts(Functions::pass_run2);
  const NamedFunc w_run2("w_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;
    double wgt = b.weight();
    if (b.SampleType()==2016)
      return wgt*b.w_prefire();
    else if (b.SampleType()==2017)
      return wgt*b.w_prefire()*Functions::wnpv2017(b);
    else 
      return wgt;
  });
	if(old16) mc16_path = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/";

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc16_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts);
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc16_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts);
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {mc16_path+"*_WJetsToLNu*.root"}, "stitch_met" && q_cuts);
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {mc16_path+"*_ST_*.root"},"stitch_met" && q_cuts);
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root", mc16_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts);
  auto mc16_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                     {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root"},"stitch_met" && q_cuts);
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
                       {mc16_path+"*DYJetsToLL*.root", 
                        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
                        mc16_path+"*_TTTT_*.root",
                        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
                        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
                        mc16_path+"_ZZ_*.root"}, "stitch_met" && q_cuts);
												

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc17_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts);
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc17_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts);
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc17_path+"*_WJetsToLNu_*.root"},
											 "pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000 && ((stitch_met && type!=2000) || (type==2000 && ht_isr_me<100))");
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc17_path+"*_ST_*.root"}, "stitch_met" && q_cuts);
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts);
  auto mc17_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                     {mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root"},"stitch_met" && q_cuts);
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
                       {mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                        mc17_path+"*_TTTT_*.root",
                        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                        mc17_path+"*_WWTo*.root",           
                        mc17_path+"*_WZ*.root",
                        mc17_path+"_ZZ_*.root"}, "stitch_met" && q_cuts);

  auto mc18_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc18_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc18_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc18_path+"*_WJetsToLNu_*.root"},
				      					q_cuts && "(stitch_met && type!=2000) || (type==2000 && ht_isr_me<100)" && Functions::hem_veto);
  auto mc18_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc18_path+"*_ST_*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc18_path+"*_TTWJets*.root", mc18_path+"*_TTZ*.root", mc18_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mc18_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                     {mc18_path+"*QCD_HT*0_Tune*.root", mc18_path+"*QCD_HT*Inf_Tune*.root"},"stitch_met" && q_cuts);
  auto mc18_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
                       {mc18_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc18_path+"*_ZJet*.root",              mc18_path+"*_ttHTobb_M125_*.root",
                        mc18_path+"*_TTTT_*.root",
                        mc18_path+"*_WH_HToBB*.root",          mc18_path+"*_ZH_HToBB*.root", 
                        mc18_path+"*_WWTo*.root",           
                        mc18_path+"*_WZ*.root",
                        mc18_path+"_ZZ_*.root"}, "(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)" && q_cuts && Functions::hem_veto);

  auto mcRunII_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                        {mc16_path+"*_TTJets*SingleLept*.root",
	                         mc17_path+"*_TTJets*SingleLept*.root",
	                         mc18_path+"*_TTJets*SingleLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                        {mc16_path+"*_TTJets*DiLept*.root",
	                         mc17_path+"*_TTJets*DiLept*.root",
	                         mc18_path+"*_TTJets*DiLept*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                        {mc16_path+"*_WJetsToLNu_*.root",
	                         mc17_path+"*_WJetsToLNu_*.root",
	                         mc18_path+"*_WJetsToLNu_*.root"},
				      					   q_cuts && "(stitch_met && type!=2000) || (type==2000 && ht_isr_me<100)" && Functions::hem_veto);
  auto mcRunII_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                        {mc16_path+"*_ST_*.root",
	                         mc17_path+"*_ST_*.root",
	                         mc18_path+"*_ST_*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                        {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root", mc16_path+"*_TTGJets*.root",
	                         mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root",
	                         mc18_path+"*_TTWJets*.root", mc18_path+"*_TTZ*.root", mc18_path+"*_TTGJets*.root"}, "stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                        {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
	                         mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
	                         mc18_path+"*QCD_HT*0_Tune*.root", mc18_path+"*QCD_HT*Inf_Tune*.root"},"stitch_met" && q_cuts && Functions::hem_veto);
  auto mcRunII_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
                          {mc16_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc16_path+"*_ZJet*.root",              mc16_path+"*_ttHTobb_M125_*.root",
                           mc16_path+"*_TTTT_*.root",
                           mc16_path+"*_WH_HToBB*.root",          mc16_path+"*_ZH_HToBB*.root", 
                           mc16_path+"*_WWTo*.root",           
                           mc16_path+"*_WZ*.root",
                           mc16_path+"_ZZ_*.root",
                           mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                           mc17_path+"*_TTTT_*.root",
                           mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                           mc17_path+"*_WWTo*.root",           
                           mc17_path+"*_WZ*.root",
                           mc17_path+"_ZZ_*.root",
                           mc18_path+"*DYJetsToLL_M-50_HT*.root", 
                           mc18_path+"*_ZJet*.root",              mc18_path+"*_ttHTobb_M125_*.root",
                           mc18_path+"*_TTTT_*.root",
                           mc18_path+"*_WH_HToBB*.root",          mc18_path+"*_ZH_HToBB*.root", 
                           mc18_path+"*_WWTo*.root",           
                           mc18_path+"*_WZ*.root",
                           mc18_path+"_ZZ_*.root"}, "(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)" && q_cuts && Functions::hem_veto);
  string sig16_path("/net/cms27/cms27r0/babymaker/babies/2019_01_11/T1tttt/unskimmed/");
  string sig17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/T1tttt/unskimmed/");
  string sig18_path("/net/cms2/cms2r0/babymaker/babies/2019_03_30/T1tttt/unskimmed/");

  auto sig16_NC = Process::MakeShared<Baby_full>("NC",sig, colors("t1tttt"), {sig16_path+"*mGluino-2100_mLSP-100_*.root"},q_cuts);
  auto sig16_C  = Process::MakeShared<Baby_full>( "C",sig, colors("t1tttt"), {sig16_path+"*mGluino-1900_mLSP-1250*.root"},q_cuts);
  auto sig17_NC = Process::MakeShared<Baby_full>("NC",sig, colors("t1tttt"), {sig17_path+"*mGluino-2100_mLSP-100_*.root"},q_cuts);
  auto sig17_C  = Process::MakeShared<Baby_full>( "C",sig, colors("t1tttt"), {sig17_path+"*mGluino-1900_mLSP-1250*.root"},q_cuts);
  auto sig18_NC = Process::MakeShared<Baby_full>("NC",sig, colors("t1tttt"), {sig18_path+"*mGluino-2100_mLSP-100_*.root"},q_cuts && Functions::hem_veto);
  auto sig18_C  = Process::MakeShared<Baby_full>( "C",sig, colors("t1tttt"), {sig18_path+"*mGluino-1900_mLSP-1250*.root"},q_cuts && Functions::hem_veto);
  auto sigRunII_NC = Process::MakeShared<Baby_full>("NC",sig, colors("t1tttt"), 
                     {sig16_path+"*mGluino-2100_mLSP-100_*.root",
                      sig17_path+"*mGluino-2100_mLSP-100_*.root",
                      sig18_path+"*mGluino-2100_mLSP-100_*.root"},q_cuts && Functions::hem_veto);
  auto sigRunII_C  = Process::MakeShared<Baby_full>( "C",sig, colors("t1tttt"), 
                     {sig16_path+"*mGluino-1900_mLSP-1250*.root",
                      sig17_path+"*mGluino-1900_mLSP-1250*.root",
                      sig18_path+"*mGluino-1900_mLSP-1250*.root"},q_cuts && Functions::hem_veto);

	vector<shared_ptr<Process> > samples_16  = {mc16_other, mc16_qcd, mc16_ttv, mc16_single_t, mc16_wjets, mc16_tt1l, mc16_tt2l, sig16_NC, sig16_C};
	vector<shared_ptr<Process> > samples_17  = {mc17_other, mc17_qcd, mc17_ttv, mc17_single_t, mc17_wjets, mc17_tt1l, mc17_tt2l, sig17_NC, sig17_C};
	vector<shared_ptr<Process> > samples_18  = {mc18_other, mc18_qcd, mc18_ttv, mc18_single_t, mc18_wjets, mc18_tt1l, mc18_tt2l, sig18_NC, sig18_C};
	vector<shared_ptr<Process> > samples_RunII  = {mcRunII_other, mcRunII_qcd, mcRunII_ttv, mcRunII_single_t, mcRunII_wjets, mcRunII_tt1l, mcRunII_tt2l, sigRunII_NC, sigRunII_C};

	vector<string> cuts= {"nleps==1 && st>500 && met>200",
                        "nveto==0",
                        "(njets>=7 && met<500) || (njets>=6 && met>500)",
                        "nbdm>=1",
	                      "mj14>250",
                        "mt>140",
                        "mj14>400",
                        "nbdm>=2",
                        "met>350 && mj14>450",
                        "met>500 && mj14>500",
                        "njets>=8"};
	vector<string> cuts_2l = {"nleps>=1 && st>500 && met>200", "nleps==2", "njets>=6", "met<500", "nbdm<=2"};
	vector<NamedFunc> cutflow;
	vector<NamedFunc> cutflow_2l;
	NamedFunc cut("1");
	for(size_t i = 0; i < cuts.size(); i++) {
	  cut = cut && cuts.at(i);
		cutflow.push_back(cut);
	}
	cut = "1";
	for(size_t i = 0; i < cuts_2l.size(); i++) {
	  cut = cut && cuts_2l.at(i);
		cutflow_2l.push_back(cut);
	}
  vector<string> label = {"RunII", "2018", "2017", "2016"};
  vector<double> lumi  = {137,       59.6,   41.5,   35.9};
  vector<vector<shared_ptr<Process> > > sample = {samples_RunII, samples_18, samples_17, samples_16};
  NamedFunc w_fr2(Functions::wgt_run2 * Functions::eff_trig_run2);
  for(int i = 0; i < 4; i++) {
    PlotMaker pm;
    pm.Push<Table>(label.at(i) + "_cutflow", vector<TableRow>{
        TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV",               cutflow.at(0) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("Track veto",                                               cutflow.at(1) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm jets} \\geq$ 7(6 if \\ptmiss $>$500) GeV",        cutflow.at(2) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm b} \\geq 1$",                                     cutflow.at(3) ,0,1,w_fr2/lumi.at(i)),
  	    TableRow("$M_{J} > 250$ GeV",                                        cutflow.at(4) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$m_{T} > 140$ GeV",                                        cutflow.at(5) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$M_{J} > 400$ GeV",                                        cutflow.at(6) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm b} \\geq 2$",                                     cutflow.at(7) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("MET $> 350$ GeV and $M_{J} > 450$ GeV",                    cutflow.at(8) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("MET $> 500$ GeV and $M_{J} > 500$ GeV",                    cutflow.at(9) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm jets} \\geq 8$",                                  cutflow.at(10),0,0,w_fr2/lumi.at(i))
  	    },sample.at(i),false,true,false,false);
    pm.Push<Table>(label.at(i)+"_cutflow_eff", vector<TableRow>{
        TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV",               cutflow.at(0) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("Track veto",                                               cutflow.at(1) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm jets} \\geq$ 7(6 if \\ptmiss $>$500) GeV",        cutflow.at(2) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm b} \\geq 1$",                                     cutflow.at(3) ,0,1,w_fr2/lumi.at(i)),
  	    TableRow("$M_{J} > 250$ GeV",                                        cutflow.at(4) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$m_{T} > 140$ GeV",                                        cutflow.at(5) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$M_{J} > 400$ GeV",                                        cutflow.at(6) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm b} \\geq 2$",                                     cutflow.at(7) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("MET $> 350$ GeV and $M_{J} > 450$ GeV",                    cutflow.at(8) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("MET $> 500$ GeV and $M_{J} > 500$ GeV",                    cutflow.at(9) ,0,0,w_fr2/lumi.at(i)),
  	    TableRow("$N_{\\rm jets} \\geq 8$",                                  cutflow.at(10),0,0,w_fr2/lumi.at(i))
  	    },sample.at(i),false,true,false,false,true,true);
    pm.min_print_ = true;
    pm.MakePlots(lumi.at(i));
  }
}
