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
// 	Process::Type data = Process::Type::data;

  Palette colors("txt/colors.txt", "default");
  // Data
  string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/");
  string data16_old_path("/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_standard/");
	string q_cuts("trig[30] && pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");

	auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",back,kBlack,  
	                 {data16_path+"*.root"},q_cuts);
	auto data_2016_old = Process::MakeShared<Baby_full>("Old 2016 Data",back,kAzure,  
	                     {data16_old_path+"*.root"},q_cuts);

	// MC samples
	vector<shared_ptr<Process> > samples = {data_2016_old, data_2016};

	vector<string> cuts= {"nleps==1 && st>500 && met>200", "nveto==0", "njets>=6", "nbdm>=1",
	                          "mj14>250","mt>140","mj14>400","nbdm>=2","met>350","met>500","njets>=9"};
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
  PlotMaker pm;
  pm.Push<Table>("data_cutflow_trig30", vector<TableRow>{
      TableRow("$1\\ell, S_{T} > 500$ GeV, MET $> 200$ GeV", cutflow.at(0) ,0,0,"1"),
	    TableRow("Track veto",                                 cutflow.at(1) ,0,0,"1"),
	    TableRow("$N_{\\rm jets} \\geq 6$",                    cutflow.at(2) ,0,0,"1"),
	    TableRow("$N_{\\rm b} \\geq 1$",                       cutflow.at(3) ,0,1,"1"),
	    TableRow("$M_{J} > 250$ GeV",                          cutflow.at(4) ,0,0,"1"),
	    TableRow("$m_{T} > 140$ GeV",                          cutflow.at(5) ,0,0,"1"),
	    TableRow("$M_{J} > 400$ GeV",                          cutflow.at(6) ,0,0,"1"),
	    TableRow("$N_{\\rm b} \\geq 2$",                       cutflow.at(7) ,0,0,"1"),
	    TableRow("MET $> 350$ GeV",                            cutflow.at(8) ,0,0,"1"),
	    TableRow("MET $> 500$ GeV",                            cutflow.at(9) ,0,0,"1"),
	    TableRow("$N_{\\rm jets} \\geq 9$",                    cutflow.at(10),0,0,"1"),
	    },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(1);
}
