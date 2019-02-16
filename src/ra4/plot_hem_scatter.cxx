#include "ra4/scatter.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/functions.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt", "default");
	Process::Type back = Process::Type::background;
  // Data
  string data18_path("/net/cms2/cms2r0/babymaker/babies/2019_01_18/data/merged_database_standard/");
	string q_cuts("pass && met/met_calo < 5 && pass_ra2_badmu && st < 10000");
	auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",back,kBlack, {data18_path+"*.root"},"trig_ra4&&"+q_cuts);
  vector<shared_ptr<Process> > data18 = {data_2018};
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist10 = {style().Stack(StackType::data_norm)
	                                     .LogMinimum(10)
	                                     .CanvasWidth(600)
	                                     .Title(TitleType::info)};
  vector<PlotOpt> bkg_hist1k = {style().Stack(StackType::data_norm)
	                                     .LogMinimum(2000)
	                                     .CanvasWidth(600)
	                                     .Title(TitleType::info)};
  NamedFunc baseline = "nleps==1&&st>500&&met>150&&njets>=6&&nbd>=1&&nveto==0";
  PlotMaker pm18;
	pm18.Push<Hist2D>(Axis(60, -3, 3,   "leps_eta", "Lepton #eta", {}), 
	                  Axis(64,-3.2, 3.2,"leps_phi", "Lepton #phi", {}),     baseline, data18, bkg_hist10).Tag("lep_2018");
	pm18.Push<Hist2D>(Axis(60, -3, 3,   "els_sceta", "Electron #eta", {}), 
	                  Axis(64,-3.2, 3.2,"els_phi",   "Electron #phi", {}),  baseline, data18, bkg_hist10).Tag("el_2018");
	pm18.Push<Hist2D>(Axis(60, -3, 3,   "mus_eta", "Muon #eta", {}), 
	                  Axis(64,-3.2, 3.2,"mus_phi", "Muon #phi", {}),        baseline, data18, bkg_hist10).Tag("mu_2018");
	pm18.Push<Hist2D>(Axis(60, -3, 3,   "jets_eta", "jet #eta", {}), 
	                  Axis(64,-3.2, 3.2,"jets_phi", "jet #phi", {}),        baseline, data18, bkg_hist1k).Tag("jet_2018");
	pm18.min_print_=true;
	pm18.MakePlots(60.0);
}
