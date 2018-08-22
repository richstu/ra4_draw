///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 150;
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  else if (Contains(hostname, "cms"))
    bfolder = "cms2";

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2018_07_25/mc/merged_mcbase_stdnj5/");

  Palette colors("txt/colors.txt", "default");
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm).LogMinimum(100);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes};
  vector<PlotOpt> shape_types = {lin_shapes};
  vector<PlotOpt> lumi_types = {log_lumi, lin_lumi};

  NamedFunc nhighpt_jets("nhighpt_jets", [](const Baby &b) ->NamedFunc::ScalarType{
      int _njets = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
        if(b.ak8jets_pt()->at(ijet)>300) _njets++;
      } 
      return _njets;
    });

  // Define processes to be used in the plots
  string c_ps = "pass && stitch_met";
  auto tt_1l = Process::MakeShared<Baby_full>("tt 1l", Process::Type::background, colors("tt_1l"),
   {foldermc+"*TTJets_SingleLeptFromT_Tune*", foldermc+"*TTJets_SingleLeptFromT_genMET*", 
   foldermc+"*TTJets_SingleLeptFromTbar_Tune*", foldermc+"*TTJets_SingleLeptFromTbar_genMET*"},c_ps);
  auto tt_1l_top = Process::MakeShared<Baby_full>("tt 1l top", Process::Type::background, kBlack,
   {foldermc+"*TTJets_SingleLeptFromT_Tune*", foldermc+"*TTJets_SingleLeptFromT_genMET*", 
   foldermc+"*TTJets_SingleLeptFromTbar_Tune*", foldermc+"*TTJets_SingleLeptFromTbar_genMET*"},"ntop_loose_nom>=2 &&" + c_ps);

  auto tt_2l = Process::MakeShared<Baby_full>("tt 2l", Process::Type::background, colors("tt_2l"),
   {foldermc+"*TTJets_DiLept_Tune*", foldermc+"*TTJets_DiLept_genMET*.root"},c_ps);
  
  auto sig2000 = Process::MakeShared<Baby_full>("T1tttt(2000,100)", 
    Process::Type::signal, kRed, {foldermc+"*T1tttt*_mGluino-2000*.root"}, "1");
  auto sig1200 = Process::MakeShared<Baby_full>("T1tttt(1200,800)", 
    Process::Type::signal, kGreen, {foldermc+"*T1tttt*_mGluino-1200*.root"}, "1");
  
  // define combinations of processes for various types of plots
  vector<shared_ptr<Process> > procs_bkg = {tt_1l, tt_2l};
  vector<shared_ptr<Process> > procs_tt_1l = {tt_1l_top, tt_1l};
  vector<shared_ptr<Process> > procs_all = {tt_1l, tt_2l, sig2000, sig1200};

  // predefine some cuts
  string baseline_s = "mj14>200 && nleps>=1 && met>100 && njets>=5 && st<10000 && pass_ra2_badmu && met/met_calo<5";

    
  // add plot definitions to the plot maker instance
  PlotMaker pm;
  // one dimensional stack plot
  pm.Push<Hist1D>(Axis(12, 0, 300, "mt", "m_{T} [GeV]", {140.}),
                  baseline_s, procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(12, 0, 300, "mt", "m_{T} [GeV]", {140.}),
                  baseline_s+"&& nak8jets>1 && ak8jets_pt[0]>300 && ak8jets_pt[1]>300", procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(20, 100, 600, "met", "E_{T}^{miss}", {200., 350., 500.}),
          baseline_s, procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(20, 100, 600, "met", "E_{T}^{miss}", {200., 350., 500.}),
          baseline_s+"&& nak8jets>1 && ak8jets_pt[0]>300 && ak8jets_pt[1]>300", procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(16, 200, 1000, "mj14", "M_{J}", {400.}),
          baseline_s, procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(16, 200, 1000, "mj14", "M_{J}", {400.}),
          baseline_s+"&& nak8jets>1 && ak8jets_pt[0]>300 && ak8jets_pt[1]>300", procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(5, 5, 10, "njets", "n_{jets}", {6, 9}),
          baseline_s, procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(5, 5, 10, "njets", "n_{jets}", {6, 9}),
          baseline_s+"&& nak8jets>1 && ak8jets_pt[0]>300 && ak8jets_pt[1]>300", procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(5, 0, 5, "nbm", "n_{b}", {1, 2, 3}),
          baseline_s, procs_tt_1l, shape_types);

  pm.Push<Hist1D>(Axis(5, 0, 5, "nbm", "n_{b}", {1, 2, 3}),
          baseline_s+"&& nak8jets>1 && ak8jets_pt[0]>300 && ak8jets_pt[1]>300", procs_tt_1l, shape_types);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
