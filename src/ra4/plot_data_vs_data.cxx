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

int main() {
    gErrorIgnoreLevel = 6000;

    Palette colors("txt/colors.txt", "default");
    Process::Type data = Process::Type::data;
    Process::Type back = Process::Type::background;

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

    string bfolder("");
    string hostname = execute("echo $HOSTNAME");
    if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
      bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

    string fold(bfolder+"/cms2r0/babymaker/babies/2019_01_18/data/skim_met100/");
    string fnew(bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/skim_met100/");

    TString baseline("nleps>=1");
    NamedFunc filters = Functions::hem_veto && "st<10000 && pass_ra2_badmu && met/met_calo<5 && pass";

    auto data_old = Process::MakeShared<Baby_full>("Old 2018",back,kAzure,  
        {fold+"*.root"},Functions::trig_run2 && baseline && filters);
    data_old->SetFillColor(kWhite);
    data_old->SetLineColor(kBlue-7);
    data_old->SetLineWidth(2);
    auto data_new = Process::MakeShared<Baby_full>("New 2018",data,kBlack,  
        {fnew+"*.root"},Functions::trig_run2 && baseline && filters);

    vector<shared_ptr<Process> > procs  = {data_old, data_new};

    PlotMaker pm;

    pm.Push<Hist1D>(Axis(30,0,600,    "met", "p_{T}^{miss} [GeV]",  {}),"njets>=5 && mj14>250 && nbd>=1",  procs, log_stack);
    pm.Push<Hist1D>(Axis(20,0,1500,    "st", "S_{T} [GeV]",  {}),    "met>200",  procs, log_stack);
    pm.Push<Hist1D>(Axis(30,  0, 1200,  "mj14", "M_{J} [GeV]",      {}), "1", procs, log_stack);
    pm.Push<Hist1D>(Axis(15,  0,  300,  "mt",   "m_{T} [GeV]",      {}), "met>100 && njets>=5 && mj14>250 && nbd>=1", procs, log_stack);
    pm.Push<Hist1D>(Axis(11,-0.5,10.5,    "njets","N_{jets}",         {}), "1", procs, log_stack);
    pm.Push<Hist1D>(Axis(5,-0.5,  4.5,  "nbd",  "N_{b, Deep}",      {}), "1", procs, log_stack);
    pm.Push<Hist1D>(Axis(16,0,800,  "jets_pt[0]",  "Jet p_{T} [GeV]",  {}), "met>200 && njets>=5 && mj14>250 && nbd>=1", procs, log_stack);
    pm.Push<Hist1D>(Axis(16,0,800,  "jets_pt[1]",  "Jet p_{T} [GeV]",  {}), "met>200 && njets>=5 && mj14>250 && nbd>=1", procs, log_stack);
    pm.Push<Hist1D>(Axis(16,0,800,  "jets_pt[2]",  "Jet p_{T} [GeV]",  {}), "met>200 && njets>=5 && mj14>250 && nbd>=1", procs, log_stack);
    pm.Push<Hist1D>(Axis(16,0,800,  "jets_pt[3]",  "Jet p_{T} [GeV]",  {}), "met>200 && njets>=5 && mj14>250 && nbd>=1", procs, log_stack);
    pm.min_print_=true;
    pm.MakePlots(1);
}

