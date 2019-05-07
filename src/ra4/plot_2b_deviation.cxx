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

    string bfolder("");
    string hostname = execute("echo $HOSTNAME");
    if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
       bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

    Palette colors("txt/colors.txt", "default");
    Process::Type back = Process::Type::background;
    Process::Type data = Process::Type::data;

    TString baseline("st>500 && nleps==1 && nveto==0");
    NamedFunc filters = Functions::hem_veto && Functions::pass_run2;
    NamedFunc nom_wgt = Functions::wgt_run2 * Functions::eff_trig_run2; 

    map<int, string> foldermc, folderdata, foldersig;
    foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
    foldersig[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/T1tttt/skim_sys_abcd/";
    folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

    foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
    foldersig[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/T1tttt/skim_sys_abcd/";
    folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

    foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_stdnj5/";
    foldersig[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/T1tttt/skim_sys_abcd/";
    folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/merged_database_standard/";

    auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
        {folderdata[2016]+"*.root"},baseline && filters && Functions::trig_run2);
    auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
        {folderdata[2017]+"*.root"},baseline && filters && Functions::trig_run2);
    auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
        {folderdata[2018]+"*.root"},baseline && filters && Functions::trig_run2);

    auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
        {foldermc[2016]+"*_TTJets*SingleLept*.root"}, baseline && filters && "stitch_met");
    auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
        {foldermc[2016]+"*_TTJets*DiLept*.root"},baseline && filters && "stitch_met");
    auto mc16_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
        {foldermc[2016]+"*_WJetsToLNu_*.root"}, baseline && filters && "stitch_met");
    auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
        {foldermc[2016]+"*_ST_*.root"}, baseline && filters && "stitch_met");
    auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
        {foldermc[2016]+"*_TTWJets*.root", foldermc[2016]+"*_TTZ*.root"}, baseline && filters && "stitch_met");
    auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
        {foldermc[2016]+"*QCD_HT*0_Tune*.root", foldermc[2016]+"*QCD_HT*Inf_Tune*.root",
        foldermc[2016]+"*_ZJet*.root", foldermc[2016]+"*_ttHTobb_M125_*.root", foldermc[2016]+"*DYJetsToLL_M-50_*.root",
        foldermc[2016]+"*_TTGJets*.root", foldermc[2016]+"*_TTTT_*.root",
        foldermc[2016]+"*_WH_HToBB*.root", foldermc[2016]+"*_ZH_HToBB*.root", 
        foldermc[2016]+"*_WWTo*.root", foldermc[2016]+"*_WZ*.root",
        foldermc[2016]+"_ZZ_*.root"}, baseline && filters && "stitch_met");

    auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
        {foldermc[2017]+"*_TTJets*SingleLept*.root"}, baseline && filters && "stitch_met");
    auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
        {foldermc[2017]+"*_TTJets*DiLept*.root"}, baseline && filters && "stitch_met");
    auto mc17_wjets    = Process::MakeShared<Baby_full>("W+jets", back, colors("wjets"), 
        {foldermc[2017]+"*_WJetsToLNu_*.root"}, baseline && filters && "stitch_met");
    auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
        {foldermc[2017]+"*_ST_*.root"}, baseline && filters && "stitch_met");
    auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
        {foldermc[2017]+"*_TTWJets*.root", foldermc[2017]+"*_TTZ*.root"}, baseline && filters && "stitch_met");
    auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
        {foldermc[2017]+"*QCD_HT*0_Tune*.root", foldermc[2017]+"*QCD_HT*Inf_Tune*.root",
        foldermc[2017]+"*_ZJet*.root",              foldermc[2017]+"*_ttHTobb_M125_*.root", foldermc[2017]+"*DYJetsToLL_M-50_*.root",
        foldermc[2017]+"*_TTGJets*.root",           foldermc[2017]+"*_TTTT_*.root",
        foldermc[2017]+"*_WH_HToBB*.root",          foldermc[2017]+"*_ZH_HToBB*.root", 
        foldermc[2017]+"*_WWTo*.root",           
        foldermc[2017]+"*_WZ*.root",
        foldermc[2017]+"_ZZ_*.root"}, baseline && filters && "stitch_met");

    vector<shared_ptr<Process> > data16_mc16  = {data_2016, mc16_tt1l, mc16_tt2l, mc16_wjets, mc16_single_t, mc16_ttv, mc16_other};
    vector<shared_ptr<Process> > data17_mc17  = {data_2017, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};
    vector<shared_ptr<Process> > data18_mc17  = {data_2018, mc17_tt1l, mc17_tt2l, mc17_wjets, mc17_single_t, mc17_ttv, mc17_other};

    PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
    PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
    log_lumi.Title(TitleType::info)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .Bottom(BottomType::ratio);
    PlotOpt lin_stack_info = log_lumi().YAxis(YAxisType::linear); 
    vector<PlotOpt> lin_stack = {lin_stack_info};

    // PlotMaker pm16;
    // pm16.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
    //     baseline,   data16_mc16, log_stack).Weight(Functions::wgt_run2).Tag("2016");
    // pm16.min_print_=true;
    // pm16.MakePlots(35.9);

    // const NamedFunc wgt("wgt", [](const Baby &b) -> NamedFunc::ScalarType{
    //   if (b.type()<1000) return 1.;
    //   return b.weight()*b.w_prefire();//*wnpv(b);
    // });

    PlotMaker pm17;
    pm17.Push<Hist1D>(Axis(6,-0.5,  5.5, "nbdm",  "N_{b}",{}),
                    "met>200 && met<=350 && njets>=5 && njets<=6 && mj14>400 && mj14<=500 && mt<=140",   
                    data17_mc17, lin_stack).Weight(nom_wgt);

    pm17.Push<Hist1D>(Axis(6,-0.5,  5.5, "nbdm",  "N_{b}",{}),
                    "met>200 && met<=350 && njets>=5 && njets<=6 && mj14>400 && mj14<=500 && mt>140",   
                    data17_mc17, lin_stack).Weight(nom_wgt);
    pm17.min_print_=true;
    pm17.MakePlots(1.);


    // PlotMaker pm18;
    // pm18.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
    //     baseline && hem,   data18_mc17, log_stack).Weight(Functions::wgt_run2).Tag("2018");
    // pm18.min_print_=true;
    // pm18.MakePlots(60.0);
    
}
