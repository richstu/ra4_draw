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

bool isRelIsoEl(const Baby &b, size_t iel);
NamedFunc::ScalarType nRelIsoEls(const Baby &b);
NamedFunc::ScalarType maxRelIsoElsPt(const Baby &b);

bool isRelIsoMu(const Baby &b, size_t imu);
NamedFunc::ScalarType nRelIsoMus(const Baby &b);
NamedFunc::ScalarType maxRelIsoMusPt(const Baby &b);

int main() {
    gErrorIgnoreLevel = 6000;

    const NamedFunc wnpv_runF("wnpv_RunF",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.SampleType()<0) return 1;
      int _npv = b.npv();
      if(_npv <  2) return  0.621; // +- 0.058
      else if(_npv <  4) return  0.708; // +- 0.033
      else if(_npv <  6) return  0.802; // +- 0.024
      else if(_npv <  8) return  0.900; // +- 0.017
      else if(_npv < 10) return  0.848; // +- 0.010
      else if(_npv < 12) return  0.686; // +- 0.006
      else if(_npv < 14) return  0.568; // +- 0.004
      else if(_npv < 16) return  0.517; // +- 0.003
      else if(_npv < 18) return  0.493; // +- 0.003
      else if(_npv < 20) return  0.498; // +- 0.002
      else if(_npv < 22) return  0.525; // +- 0.002
      else if(_npv < 24) return  0.561; // +- 0.002
      else if(_npv < 26) return  0.624; // +- 0.002
      else if(_npv < 28) return  0.710; // +- 0.003
      else if(_npv < 30) return  0.825; // +- 0.003
      else if(_npv < 32) return  0.966; // +- 0.003
      else if(_npv < 34) return  1.166; // +- 0.004
      else if(_npv < 36) return  1.372; // +- 0.005
      else if(_npv < 38) return  1.586; // +- 0.005
      else if(_npv < 40) return  1.765; // +- 0.006
      else if(_npv < 42) return  1.961; // +- 0.007
      else if(_npv < 44) return  2.196; // +- 0.009
      else if(_npv < 46) return  2.354; // +- 0.010
      else if(_npv < 48) return  2.551; // +- 0.012
      else if(_npv < 50) return  2.773; // +- 0.015
      else if(_npv < 52) return  3.093; // +- 0.019
      else if(_npv < 54) return  3.429; // +- 0.024
      else if(_npv < 56) return  3.750; // +- 0.030
      else if(_npv < 58) return  4.065; // +- 0.038
      else if(_npv < 60) return  4.376; // +- 0.048
      else if(_npv < 62) return  5.309; // +- 0.067
      else if(_npv < 64) return  5.848; // +- 0.086
      else if(_npv < 66) return  6.571; // +- 0.112
      else if(_npv < 68) return  6.839; // +- 0.139
      else if(_npv < 70) return  7.045; // +- 0.165
      else if(_npv < 72) return  8.197; // +- 0.226
      else if(_npv < 74) return  8.906; // +- 0.290
      else if(_npv < 76) return  9.676; // +- 0.355
      else if(_npv < 78) return 10.002; // +- 0.437
      else if(_npv < 80) return 11.745; // +- 0.594
      else if(_npv < 82) return 12.839; // +- 0.751
      else if(_npv < 84) return 11.304; // +- 0.788
      else if(_npv < 86) return 11.312; // +- 0.963
      else if(_npv < 88) return 12.399; // +- 1.228
      else if(_npv < 90) return 15.106; // +- 1.575
      else if(_npv < 92) return 15.361; // +- 1.863
      else if(_npv < 94) return 11.858; // +- 1.950
      else if(_npv < 96) return 12.607; // +- 2.073
      else if(_npv < 98) return  7.499; // +- 1.875
      else return 11.104; // +- 1.434
    });

    const NamedFunc wnpv_runBE("wnpv_RunBE",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.SampleType()<0) return 1;
      int _npv = b.npv();
      if(_npv <  2) return  0.226;// +- 0.023
      else if(_npv <  4) return  0.280;// +- 0.013
      else if(_npv <  6) return  0.395;// +- 0.011
      else if(_npv <  8) return  0.554;// +- 0.009
      else if(_npv < 10) return  0.740;// +- 0.006
      else if(_npv < 12) return  0.839;// +- 0.004
      else if(_npv < 14) return  0.907;// +- 0.003
      else if(_npv < 16) return  0.957;// +- 0.003
      else if(_npv < 18) return  0.998;// +- 0.002
      else if(_npv < 20) return  1.020;// +- 0.002
      else if(_npv < 22) return  1.058;// +- 0.002
      else if(_npv < 24) return  1.080;// +- 0.002
      else if(_npv < 26) return  1.096;// +- 0.002
      else if(_npv < 28) return  1.111;// +- 0.002
      else if(_npv < 30) return  1.101;// +- 0.002
      else if(_npv < 32) return  1.068;// +- 0.002
      else if(_npv < 34) return  1.036;// +- 0.002
      else if(_npv < 36) return  0.979;// +- 0.003
      else if(_npv < 38) return  0.901;// +- 0.003
      else if(_npv < 40) return  0.826;// +- 0.003
      else if(_npv < 42) return  0.784;// +- 0.003
      else if(_npv < 44) return  0.761;// +- 0.003
      else if(_npv < 46) return  0.726;// +- 0.004
      else if(_npv < 48) return  0.743;// +- 0.004
      else if(_npv < 50) return  0.757;// +- 0.005
      else if(_npv < 52) return  0.803;// +- 0.006
      else if(_npv < 54) return  0.866;// +- 0.008
      else if(_npv < 56) return  0.933;// +- 0.010
      else if(_npv < 58) return  1.006;// +- 0.012
      else if(_npv < 60) return  1.075;// +- 0.016
      else if(_npv < 62) return  1.266;// +- 0.021
      else if(_npv < 64) return  1.448;// +- 0.028
      else if(_npv < 66) return  1.560;// +- 0.036
      else if(_npv < 68) return  1.601;// +- 0.044
      else if(_npv < 70) return  1.678;// +- 0.053
      else if(_npv < 72) return  2.004;// +- 0.073
      else if(_npv < 74) return  2.323;// +- 0.097
      else if(_npv < 76) return  2.266;// +- 0.113
      else if(_npv < 78) return  2.383;// +- 0.140
      else if(_npv < 80) return  2.852;// +- 0.192
      else if(_npv < 82) return  2.871;// +- 0.233
      else if(_npv < 84) return  3.112;// +- 0.271
      else if(_npv < 86) return  3.486;// +- 0.350
      else if(_npv < 88) return  2.663;// +- 0.373
      else if(_npv < 90) return  3.245;// +- 0.478
      else if(_npv < 92) return  3.494;// +- 0.582
      else if(_npv < 94) return  3.855;// +- 0.729
      else if(_npv < 96) return  2.488;// +- 0.604
      else if(_npv < 98) return  1.611;// +- 0.569
      else return 2.783;// +- 0.470
    });

    const NamedFunc hem("hem",[](const Baby &b) -> NamedFunc::ScalarType{
        if (b.type()<1000) {
            if (b.run() >= 319077) { 
                if(b.nels() > 0) {
                    for(size_t i = 0; i < b.els_pt()->size(); i++) {
                        if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
                            return static_cast<float>(0);
                    }
                }
                for(size_t i = 0; i < b.jets_pt()->size(); i++) {
                    if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
                        return static_cast<float>(0);
                }
            }
        } else {
            if ((b.event()%1961) < 1296) { 
                if(b.nels() > 0) {
                    for(size_t i = 0; i < b.els_pt()->size(); i++) {
                        if(b.els_pt()->at(i) > 20 && b.els_sceta()->at(i) < -1.5 && (b.els_phi()->at(i) > -1.6 && b.els_phi()->at(i) < -0.8) && b.els_sigid()->at(i)) 
                            return static_cast<float>(0);
                    }
                }
                for(size_t i = 0; i < b.jets_pt()->size(); i++) {
                    if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
                        return static_cast<float>(0);
                }
            }
        }
        return static_cast<float>(1);
    });

    Palette colors("txt/colors.txt", "default");
    Process::Type back = Process::Type::background;
    Process::Type data = Process::Type::data;

    string data16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/data/skim_zcandnb0/");
    string data17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/data/skim_zcandnb0/");
    string data17old_path("/net/cms27/cms27r0/babymaker/babies/2018_01_30/data/skim_zcandnb0/");
    string data18_path("/net/cms2/cms2r0/babymaker/babies/2019_01_18/data/skim_zcandnb0/");


  // 2016 trigs: 19 - IsoMu24, 20 - IsoMu27, 21 - Mu50, 40 - Ele27, 41 - Ele115
  string trig16 = "(trig[19] || trig[20] || trig[21] || trig[24] || trig[40] || trig[41])";
  // 2017 trigs: 19 - IsoMu24, 20 - IsoMu27, 21 - Mu50, 23 - Ele35, 24 - Ele115
  string trig17 = "(trig[19] || trig[20] || trig[21] || trig[23] || trig[24])";

    auto data_2016 = Process::MakeShared<Baby_full>("2016 Data",data,kBlack,  
        {data16_path+"*.root"},trig16);
    auto data_2017 = Process::MakeShared<Baby_full>("2017 Data",data,kBlack,
        {data17_path+"*.root"},trig17);
    auto data_2017old = Process::MakeShared<Baby_full>("2017 Data (Old)",data,kBlack,
        {data17old_path+"*.root"},trig17);
    auto data_2018 = Process::MakeShared<Baby_full>("2018 Data",data,kBlack,
        {data18_path+"*.root"},trig17);

    string mc16_path("/net/cms2/cms2r0/babymaker/babies/2019_01_11/mc/skim_zcandnb0/");
    string mc17_path("/net/cms2/cms2r0/babymaker/babies/2018_12_17/mc/skim_zcandnb0/");

    NamedFunc nreliso_els("nreliso_els",nRelIsoEls);
    NamedFunc nreliso_mus("nreliso_mus",nRelIsoMus);
    NamedFunc baseline("pass && nleps==2&&leps_pt[0]>60&&((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && nbd==0");
    baseline = baseline && (nreliso_els+nreliso_mus>=1);

    auto mc16_dyjets    = Process::MakeShared<Baby_full>("DY+jets", back, colors("dy"), 
        {mc16_path+"*DYJetsToLL_M-50_*.root"},"(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)");
    auto mc16_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
        {mc16_path+"*_TTJets*DiLept*.root"}, "ntruleps>=2 && stitch_met");
    auto mc16_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
        {mc16_path+"*_TTJets*SingleLept*.root"}, "ntruleps<=1 && stitch_met");
    auto mc16_single_t = Process::MakeShared<Baby_full>("Single t",  back, colors("single_t"), 
        {mc16_path+"*_ST_*.root"});
    auto mc16_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V",      back, colors("ttv"), 
        {mc16_path+"*_TTWJets*.root", mc16_path+"*_TTZ*.root"});
    auto mc16_other    = Process::MakeShared<Baby_full>("Other",        back, colors("other"),
        {mc16_path+"*_WJetsToLNu_Tune*.root", mc16_path+"*QCD_HT*0_Tune*.root", mc16_path+"*QCD_HT*Inf_Tune*.root",
        mc16_path+"*_ZJet*.root", mc16_path+"*_ttHTobb_M125_*.root",
        mc16_path+"*_TTGJets*.root", mc16_path+"*_TTTT_*.root",
        mc16_path+"*_WH_HToBB*.root", mc16_path+"*_ZH_HToBB*.root", 
        mc16_path+"*_WWTo*.root", mc16_path+"*_WZ*.root",
        mc16_path+"_ZZ_*.root"}, "1");

    auto mc17_dyjets    = Process::MakeShared<Baby_full>("DY+jets", back, colors("dy"), 
        {mc17_path+"*DYJetsToLL_M-50_*.root"},"(stitch_met && type!=6000) || (type==6000 && ht_isr_me<100)");
    auto mc17_tt2l     = Process::MakeShared<Baby_full>("t#bar{t} (2l)", back, colors("tt_2l"), 
        {mc17_path+"*_TTJets*DiLept*.root"}, "ntruleps>=2 && stitch_met");
    auto mc17_tt1l     = Process::MakeShared<Baby_full>("t#bar{t} (1l)", back, colors("tt_1l"), 
        {mc17_path+"*_TTJets*SingleLept*.root"}, "ntruleps<=1 && stitch_met");
    auto mc17_single_t = Process::MakeShared<Baby_full>("Single t", back, colors("single_t"), 
        {mc17_path+"*_ST_*.root"});
    auto mc17_ttv      = Process::MakeShared<Baby_full>("t#bar{t}V", back, colors("ttv"), 
        {mc17_path+"*_TTWJets*.root", mc17_path+"*_TTZ*.root"});
    auto mc17_other    = Process::MakeShared<Baby_full>("Other", back, colors("other"),
        {mc17_path+"*_WJetsToLNu_Tune*.root", mc17_path+"*QCD_HT*0_Tune*.root", mc17_path+"*QCD_HT*Inf_Tune*.root",
        mc17_path+"*_ZJet*.root",              mc17_path+"*_ttHTobb_M125_*.root",
        mc17_path+"*_TTGJets*.root",           mc17_path+"*_TTTT_*.root",
        mc17_path+"*_WH_HToBB*.root",          mc17_path+"*_ZH_HToBB*.root", 
        mc16_path+"*_WWTo*.root",           
        mc17_path+"*_WZ*.root",
        mc17_path+"_ZZ_*.root"}, "1");

    vector<shared_ptr<Process> > data16_mc16  = {data_2016, mc16_dyjets, mc16_tt2l, mc16_tt1l, mc16_single_t, mc16_ttv, mc16_other};
    vector<shared_ptr<Process> > data17_mc17  = {data_2017, mc17_dyjets, mc17_tt2l, mc17_tt1l, mc17_single_t, mc17_ttv, mc17_other};
    vector<shared_ptr<Process> > data17old_mc17  = {data_2017old, mc17_dyjets, mc17_tt2l, mc17_tt1l, mc17_single_t, mc17_ttv, mc17_other};
    vector<shared_ptr<Process> > data18_mc17  = {data_2018, mc17_dyjets, mc17_tt2l, mc17_tt1l, mc17_single_t, mc17_ttv, mc17_other};

    PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
    PlotOpt log_ratios("txt/plot_styles.txt", "CMSPaper");
    log_lumi.Title(TitleType::info)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .Bottom(BottomType::ratio);
    PlotOpt lin_stack_info = log_lumi().YAxis(YAxisType::linear); 
    vector<PlotOpt> lin_stack = {lin_stack_info};
    vector<PlotOpt> log_stack = {log_lumi.PrintVals(true)};
    NamedFunc Run2017F("(type>=1000 || run>=304911)");
    NamedFunc Run2017BCDE("(type>=1000 || run<304911)");

    PlotMaker pm16;
    pm16.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
        baseline,   data16_mc16, log_stack).Weight(Functions::wgt_run2).Tag("2016");
    pm16.min_print_=true;
    pm16.MakePlots(35.9);

    const NamedFunc wgt("wgt", [](const Baby &b) -> NamedFunc::ScalarType{
      if (b.type()<1000) return 1.;
      return b.weight()*b.w_prefire();//*wnpv(b);
    });

    PlotMaker pm17;
    pm17.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
        baseline,   data17_mc17, log_stack).Weight(Functions::wgt_run2).Tag("2017");
    pm17.Push<Hist1D>(Axis(50,0, 100, "npv",  "NPV",{}), 
        baseline,   data17_mc17, log_stack).Weight(Functions::wgt_run2).Tag("2017");
    pm17.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
        baseline && Run2017BCDE,   data17_mc17, log_stack).Weight(wnpv_runBE*"weight").Tag("2017");
    pm17.Push<Hist1D>(Axis(50,0, 100, "npv",  "NPV",{}), 
        baseline && Run2017F,   data17_mc17, log_stack).Weight(wnpv_runF*"weight").Tag("2017");
    pm17.min_print_=true;
    pm17.MakePlots(1.);


    PlotMaker pm18;
    pm18.Push<Hist1D>(Axis(25,0, 250, "met",  "p_{T}^{miss} [GeV]",{}), 
        baseline && hem,   data18_mc17, log_stack).Weight(Functions::wgt_run2).Tag("2018");
    pm18.min_print_=true;
    pm18.MakePlots(60.0);
    
}


bool isRelIsoEl(const Baby &b, size_t iel){
  return iel<b.els_pt()->size()
      && b.els_pt()->at(iel)>30.
      && fabs(b.els_sceta()->at(iel))<2.
      && b.els_tight()->at(iel)
      && b.els_reliso()->at(iel) < 0.1;
}

NamedFunc::ScalarType nRelIsoEls(const Baby &b){
  int nels = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (isRelIsoEl(b,iel)) nels++;
  }
  return nels;
}

NamedFunc::ScalarType maxRelIsoElsPt(const Baby &b){
  double max_pt = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (isRelIsoEl(b,iel) && b.els_pt()->at(iel)>max_pt) max_pt = b.els_pt()->at(iel);
  }
  return max_pt;
}

bool isRelIsoMu(const Baby &b, size_t imu){
  return imu<b.mus_pt()->size()
      && b.mus_pt()->at(imu)>30.
      && fabs(b.mus_eta()->at(imu))<2.
      && b.mus_tight()->at(imu)
      && b.mus_reliso()->at(imu) < 0.1;
}

NamedFunc::ScalarType nRelIsoMus(const Baby &b){
  int nmus = 0;
  for (size_t imu(0); imu<b.mus_pt()->size(); imu++){
    if (isRelIsoMu(b,imu)) nmus++;
  }
  return nmus;
}

NamedFunc::ScalarType maxRelIsoMusPt(const Baby &b){
  double max_pt = 0;
  for (size_t imu(0); imu<b.mus_pt()->size(); imu++){
    if (isRelIsoMu(b,imu) && b.mus_pt()->at(imu)>max_pt) max_pt = b.mus_pt()->at(imu);
  }
  return max_pt;
}