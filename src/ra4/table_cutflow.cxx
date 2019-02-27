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

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
	Process::Type back = Process::Type::background;

  Palette colors("txt/colors.txt", "default");
	bool old16(false);
	// MC samples
  string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_standard/");
  string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_standard/");
	string q_cuts("stitch_met && pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");
	if(old16) mc16_path = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/";

  auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc16_path+"*_TTJets*SingleLept*.root"}, q_cuts);
  auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc16_path+"*_TTJets*DiLept*.root"}, q_cuts);
  auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets",       back, colors("wjets"), 
	                     {mc16_path+"*_WJetsToLNu*.root"}, q_cuts);
  auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
	                     {mc16_path+"*_ST_*.root"},q_cuts);
  auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
	                     {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root", mc16_path+"*_TTGJets*.root"}, q_cuts);
  auto mc16_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                     {mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root"},q_cuts);
  auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
                       {mc16_path+"*DYJetsToLL*.root", 
                        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
                        mc16_path+"*_TTTT_*.root",
                        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
                        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
                        mc16_path+"_ZZ_*.root"}, q_cuts);
												

  auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
	                     {mc17_path+"*_TTJets*SingleLept*.root"}, q_cuts);
  auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
	                     {mc17_path+"*_TTJets*DiLept*.root"}, q_cuts);
  auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
	                     {mc17_path+"*_WJetsToLNu_*.root"},
											 "pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000 && ((stitch_met && type!=2000) || (type==2000 && ht_isr_me<100))");
  auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
	                     {mc17_path+"*_ST_*.root"}, q_cuts);
  auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
	                     {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root", mc17_path+"*_TTGJets*.root"}, q_cuts);
  auto mc17_qcd      = Process::MakeShared<Baby_full>("QCD", back, colors("qcd"), 
	                     {mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root"},q_cuts);
  auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
                       {mc17_path+"*DYJetsToLL_M-50_HT*.root", 
                        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
                        mc17_path+"*_TTTT_*.root",
                        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
                        mc17_path+"*_WWTo*.root",           
                        mc17_path+"*_WZ*.root",
                        mc17_path+"_ZZ_*.root"}, q_cuts);

	vector<shared_ptr<Process> > samples_16  = {mc16_other, mc16_qcd, mc16_ttv, mc16_single_t, mc16_wjets, mc16_tt1l, mc16_tt2l};
	vector<shared_ptr<Process> > samples_17  = {mc17_other, mc17_qcd, mc17_ttv, mc17_single_t, mc17_wjets, mc17_tt1l, mc17_tt2l};

	vector<string> cuts= {"nleps==1 && st>500 && met>200", "nveto==0", "njets>=6", "nbd>=1",
	                          "mj14>250","mt>140","mj14>400","nbd>=2","met>350","met>500","njets>=9"};
	vector<string> cuts_2l = {"nleps>=1 && st>500 && met>200", "nleps==2", "njets>=6", "met<500", "nbd<=2"};
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
  PlotMaker pm16;
  pm16.Push<Table>("2016_cutflow", vector<TableRow>{
      TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow.at(0) ,0,0),
	    TableRow("Track veto",                                 cutflow.at(1) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 6$",                    cutflow.at(2) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 1$",                       cutflow.at(3) ,0,1),
	    TableRow("$M_{J} > 250$ GeV",                          cutflow.at(4) ,0,0),
	    TableRow("$m_{T} > 140$ GeV",                          cutflow.at(5) ,0,0),
	    TableRow("$M_{J} > 400$ GeV",                          cutflow.at(6) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 2$",                       cutflow.at(7) ,0,0),
	    TableRow("MET $> 350$ GeV",                            cutflow.at(8) ,0,0),
	    TableRow("MET $> 500$ GeV",                            cutflow.at(9) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 9$",                    cutflow.at(10),0,0),
	    },samples_16,false);
  pm16.Push<Table>("2016_cutflow_eff", vector<TableRow>{
      TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow.at(0) ,0,0),
	    TableRow("Track veto",                                 cutflow.at(1) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 6$",                    cutflow.at(2) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 1$",                       cutflow.at(3) ,0,1),
	    TableRow("$M_{J} > 250$ GeV",                          cutflow.at(4) ,0,0),
	    TableRow("$m_{T} > 140$ GeV",                          cutflow.at(5) ,0,0),
	    TableRow("$M_{J} > 400$ GeV",                          cutflow.at(6) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 2$",                       cutflow.at(7) ,0,0),
	    TableRow("MET $> 350$ GeV",                            cutflow.at(8) ,0,0),
	    TableRow("MET $> 500$ GeV",                            cutflow.at(9) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 9$",                    cutflow.at(10),0,0),
	    },samples_16,false,true,true);
  pm16.Push<Table>("2016_2l_cutflow", vector<TableRow>{
      TableRow("$\\geq 1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow_2l.at(0),0,0),
	    TableRow("$2\\ell$",                                         cutflow_2l.at(1),0,0),
	    TableRow("$N_{\\rm jets} \\geq 5$",                          cutflow_2l.at(2),0,0),
	    TableRow("$200 <$ MET $< 500$ GeV",                          cutflow_2l.at(3),0,0),
	    TableRow("$N_{\\rm b} \\leq 2$",                             cutflow_2l.at(4),0,0),
	    },samples_16,false);
  pm16.min_print_ = true;
  pm16.MakePlots(35.9);
  PlotMaker pm17;
  pm17.Push<Table>("2017_cutflow", vector<TableRow>{
      TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow.at(0) ,0,0),
	    TableRow("Track veto",                                 cutflow.at(1) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 6$",                    cutflow.at(2) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 1$",                       cutflow.at(3) ,0,1),
	    TableRow("$M_{J} > 250$ GeV",                          cutflow.at(4) ,0,0),
	    TableRow("$m_{T} > 140$ GeV",                          cutflow.at(5) ,0,0),
	    TableRow("$M_{J} > 400$ GeV",                          cutflow.at(6) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 2$",                       cutflow.at(7) ,0,0),
	    TableRow("MET $> 350$ GeV",                            cutflow.at(8) ,0,0),
	    TableRow("MET $> 500$ GeV",                            cutflow.at(9) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 9$",                    cutflow.at(10),0,0),
	    },samples_17,false);
  pm17.Push<Table>("2017_cutflow_eff", vector<TableRow>{
      TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow.at(0) ,0,0),
	    TableRow("Track veto",                                 cutflow.at(1) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 6$",                    cutflow.at(2) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 1$",                       cutflow.at(3) ,0,1),
	    TableRow("$M_{J} > 250$ GeV",                          cutflow.at(4) ,0,0),
	    TableRow("$m_{T} > 140$ GeV",                          cutflow.at(5) ,0,0),
	    TableRow("$M_{J} > 400$ GeV",                          cutflow.at(6) ,0,0),
	    TableRow("$N_{\\rm b} \\geq 2$",                       cutflow.at(7) ,0,0),
	    TableRow("MET $> 350$ GeV",                            cutflow.at(8) ,0,0),
	    TableRow("MET $> 500$ GeV",                            cutflow.at(9) ,0,0),
	    TableRow("$N_{\\rm jets} \\geq 9$",                    cutflow.at(10),0,0),
	    },samples_17,false,true,true);
  pm17.Push<Table>("2017_2l_cutflow", vector<TableRow>{
      TableRow("$\\geq 1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow_2l.at(0),0,0),
	    TableRow("$2\\ell$",                                         cutflow_2l.at(1),0,0),
	    TableRow("$N_{\\rm jets} \\geq 5$",                          cutflow_2l.at(2),0,0),
	    TableRow("$200 <$ MET $< 500$ GeV",                          cutflow_2l.at(3),0,0),
	    TableRow("$N_{\\rm b} \\leq 2$",                             cutflow_2l.at(4),0,0),
	    },samples_17,false);
  pm17.min_print_ = true;
  pm17.MakePlots(41.5);
}
