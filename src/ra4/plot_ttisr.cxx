// Macro for deriving ISR weights for strong production

#include <cmath>
#include <algorithm>

#include "TError.h"
#include "TVector2.h"
#include "TString.h"

#include <getopt.h>

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  const string isrtype = "ttisr";
  bool do_tt1l = false;
  int year = 2016;

  // double CSVMedium = 0.800;
  double DeepCSVMedium_2016 = 0.6324;
  double DeepCSVMedium_2017 = 0.4941;

  bool single_thread = false;
  bool quick = false; 
}

void addSlices(PlotMaker &pm, const vector<double> slices, NamedFunc svar,
               const vector<double> xbins, NamedFunc xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits=0);

NamedFunc::ScalarType zCandidate(const Baby &b);
bool isRelIsoEl(const Baby &b, size_t iel);
NamedFunc::ScalarType nRelIsoEls(const Baby &b);
NamedFunc::ScalarType maxRelIsoElsPt(const Baby &b);

bool isRelIsoMu(const Baby &b, size_t imu);
NamedFunc::ScalarType nRelIsoMus(const Baby &b);
NamedFunc::ScalarType maxRelIsoMusPt(const Baby &b);

bool isGoodJet(const Baby &b, size_t ijet);
NamedFunc::VectorType isrJetsPt(const Baby &b, float ptThresh=30.);
NamedFunc::ScalarType isrSystemPt(const Baby &b);

NamedFunc::ScalarType nJetsWeights_ttisr(const Baby &b, bool use_baby_nisr, string wgtopt);
NamedFunc::ScalarType nJetsWeights_visr(const Baby &b);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  float lumi = 41.5;
  if (year==2016) lumi = 35.9;
  cout<<"Running on "<<year<<" data."<<endl;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  //// Processes for ISR skims
  string dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_isrmc_"+isrtype+"/";
  string dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_isrdata_"+isrtype+"/";
  if (year==2016) {
    dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_isrmc_"+isrtype+"/";
    dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2017_01_27/data/merged_isrdata_"+isrtype+"/";
  }
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {dir_mc_isr+"*_TTJets_SingleLeptFromT_Tune*.root", 
    dir_mc_isr+"*_TTJets_SingleLeptFromTbar_Tune*.root", 
    dir_mc_isr+"*_TTJets_HT*.root"}, "ntruleps<=1 && stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {dir_mc_isr+"*_TTJets_DiLept_Tune*.root", dir_mc_isr+"*_TTJets_HT*.root"}, "ntruleps>=2 && stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_isr+"*_ST_*.root"});
  auto dyjets = Process::MakeShared<Baby_full>("DY+jets", Process::Type::background, kOrange+1,
    {dir_mc_isr+"*DYJetsToLL_M-50_Tune*.root"}); // Inclusive DY only since stitch is wrong
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, kTeal-8,
    {dir_mc_isr+"*_TTWJets*.root", dir_mc_isr+"*_TTZTo*.root", dir_mc_isr+"*_TTGJets*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*_ZJet*.root", 
    dir_mc_isr+"*QCD_HT*0_Tune*.root", dir_mc_isr+"*QCD_HT*Inf_Tune*.root",
    dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*_ttHJetTobb*.root",
    dir_mc_isr+"*_TTTT*.root", dir_mc_isr+"*_WWTo*.root",
    dir_mc_isr+"*_WH_HToBB*.root",dir_mc_isr+"*_ZH_HToBB*.root",
    dir_mc_isr+"*_WZ*.root",dir_mc_isr+"*_ZZ_*.root"},"stitch");
  
  //// Only used for W+jets
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {dir_mc_isr+"*_WJetsToLNu*.root"},"stitch");
  auto other_w = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*DYJetsToLL_M-50_Tu*.root",
        dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WWTo*.root",
        dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*ttHJetTobb*.root",
        dir_mc_isr+"*_TTTT_*.root",
        dir_mc_isr+"*_WH_HToBB*.root", dir_mc_isr+"*_WZTo*.root",
        dir_mc_isr+"*_ZH_HToBB*.root", dir_mc_isr+"_ZZ_*.root"});

  // 2016 trigs: 19 - IsoMu24, 20 - IsoMu27, 21 - Mu50, 40 - Ele27, 41 - Ele115
  string trig = "(trig[19] || trig[20] || trig[21] || trig[24] || trig[40] || trig[41])";
  // 2017 trigs: 19 - IsoMu24, 20 - IsoMu27, 21 - Mu50, 23 - Ele35, 24 - Ele115
  if (year!=2016) trig = "(trig[19] || trig[20] || trig[21] || trig[23] || trig[24])";
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Process::MakeShared<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {dir_data_isr+"*.root"},
    "pass &&" + trig);

  vector<shared_ptr<Process> > procs;
  if (isrtype=="zisr") procs = {data, dyjets, tt2l, tt1l, single_t, ttv, other};
  else if (isrtype=="ttisr") {
    procs = {data, tt2l, tt1l, dyjets, single_t, ttv, other};
    if (quick) procs = {data,tt1l};
  } else if (isrtype=="wisr") procs = {data, wjets, tt1l, tt2l, single_t, ttv, other_w};
  else {cout<<isrtype<<" not supported, exiting"<<endl<<endl; return 0;}

  //// Processes for 1l ttbar closure
  string dir_mc_std(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/");
  string dir_data_std(bfolder+"/cms2r0/babymaker/babies/2017_01_27/data/merged_database_standard/");

  auto std_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {dir_data_std+"/*.root"},
    "pass && trig_ra4");

  auto std_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {dir_mc_std+"*_TTJets*SingleLept*.root", dir_mc_std+"*_TTJets_HT*.root"},
    "ntruleps==1 && stitch");
  auto std_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {dir_mc_std+"*_TTJets*DiLept*.root", dir_mc_std+"*_TTJets_HT*.root"},
    "ntruleps==2 && stitch");
  auto std_wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {dir_mc_std+"*_WJetsToLNu*.root"}, "stitch");
  auto std_singlet = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_std+"*_ST_*.root"});
  auto std_ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {dir_mc_std+"*_TTWJets*.root", dir_mc_std+"*_TTZTo*.root"});
  auto std_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_std+"*DYJetsToLL*.root", 
    dir_mc_std+"*_ZJet*.root", 
    dir_mc_std+"*QCD_HT*0_Tune*.root", dir_mc_std+"*QCD_HT*Inf_Tune*.root",
    dir_mc_std+"*ggZH_HToBB*.root", dir_mc_std+"*_ttHJetTobb*.root",
    dir_mc_std+"*_TTTT*.root", dir_mc_std+"*_WWTo*.root",
    dir_mc_std+"*_WH_HToBB*.root",dir_mc_std+"*_ZH_HToBB*.root",
    dir_mc_std+"*_WZ*.root",dir_mc_std+"*_ZZ_*.root"},"stitch");

  vector<shared_ptr<Process> > procs_1l = {std_data, std_tt1l, std_tt2l, std_wjets, std_singlet, std_ttv, std_other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {log_lumi};//, lin_lumi};
  vector<PlotOpt> plot_vals = {lin_lumi().PrintVals(true)};
  PlotMaker pm;

  // tt_isr skim def:
  // "nvleps==2 && nleps>=1 && nbm==2 &&
  // max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30"
  NamedFunc nreliso_els("nreliso_els",nRelIsoEls);
  NamedFunc nreliso_mus("nreliso_mus",nRelIsoMus);
  NamedFunc zcand("zcand", zCandidate);
  NamedFunc baseline = "nleps==2" && zcand<1; // && nreliso_els+nreliso_mus>=1
  NamedFunc baseline_w = nreliso_els+nreliso_mus==1 && "ht>200&&met>100&&nbl==0";
  NamedFunc baseline_1l = "nleps==1 && ht>500 && met>200 && nbm>=1 && mt<140 && pass";
  if(isrtype=="wisr") baseline = baseline_w;

  NamedFunc max_reliso_elspt("max_reliso_elspt",maxRelIsoElsPt);
  NamedFunc max_reliso_muspt("max_reliso_muspt",maxRelIsoMusPt);
  NamedFunc isr_jetspt("isr_jetspt",[&](const Baby &b){
      return isrJetsPt(b, 30.);
    });
  NamedFunc nisrjets("nisrjets", [&](const Baby &b){
      return isrJetsPt(b, 30.).size();
    });
  NamedFunc isr_ht("isr_ht", [&](const Baby &b){
      vector<double> jets_pt = isrJetsPt(b, 30.);
      double ht = 0;
      for (auto &jpt: jets_pt) ht += jpt;
      return ht;
    });
  NamedFunc nisrjets50("nisrjets50", [&](const Baby &b){
      return isrJetsPt(b, 50).size();
    });
  NamedFunc nisrjets75("nisrjets75", [&](const Baby &b){
      return isrJetsPt(b, 75).size();
    });
  NamedFunc isr_syspt("isr_syspt", isrSystemPt);

  // definitions for njets in slices of ISR pT
  const vector<double> isr_syspt_slices = {0, 50, 100,150,200, 300,400};
  const vector<double> nisrjet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5};

  // definitions for ISR pT in slices of njets
  vector<double> nisrjet_slices = {0,1,2,3,4,5};
  const vector<double> isr_syspt_bins = {0, 50, 100, 150, 200, 300, 400, 600, 800};

  const vector<double> isr_ht_slices = {0,100, 200, 300, 400, 500, 2000};
  const vector<double> isr_ht_bins = {0,100, 200, 300};

  vector<double> ptbins = {30,40,50,75,100,150,200,300,400,600};
  vector<double> lepptbins = {20,30,40,50,75,100,150,200,300,400,600};
  vector<double> ptbins_zoom = {20,25,30,35,40,50,75,100,150,200};

  vector<NamedFunc> weight_opts;
  // weight_opts.push_back(NamedFunc("ichep", 
  //   [](const Baby &b) -> NamedFunc::ScalarType{return b.weight();}));
  if (isrtype=="ttisr"){
    //reweight TTJets only according to nisr = b.nisr()
    weight_opts.push_back(NamedFunc("none", 
      [](const Baby &b) -> NamedFunc::ScalarType{return nJetsWeights_ttisr(b, true,"none");}));
    // weight_opts.push_back(NamedFunc("is2016", 
    //   [](const Baby &b) -> NamedFunc::ScalarType{return nJetsWeights_ttisr(b, true,"is2016");}));
    //reweight TTJets only according to nisr = b.njets()-2 
    // weight_opts.push_back(NamedFunc("check", 
    //   [](const Baby &b) -> NamedFunc::ScalarType{return nJetsWeights_ttisr(b, false,"moriond");}));
  } else if (isrtype=="zisr"){
    weight_opts.push_back(NamedFunc("w_visr", nJetsWeights_visr));
  }
  for (const auto &iweight: weight_opts){
    if(do_tt1l) {//// 1l ttbar closure
      pm.Push<Hist1D>(Axis(13, -0.5, 12.5, "njets", "Number of jets"), 
        baseline_1l, procs_1l, plot_types).Weight(iweight).Tag("tt1l");
      pm.Push<Hist1D>(Axis(20,200.,700., "met", "MET [GeV]"), 
        baseline_1l, procs_1l, plot_types).Weight(iweight).Tag("tt1l");
      pm.Push<Hist1D>(Axis(15,500.,2000., "ht", "H_{T} [GeV]"), 
        baseline_1l, procs_1l, plot_types).Weight(iweight).Tag("tt1l");
      pm.Push<Hist1D>(Axis(15,0.,1500., "mj14", "M_{J} [GeV]"), 
        baseline_1l, procs_1l, plot_types).Weight(iweight).Tag("tt1l");       
      pm.Push<Hist1D>(Axis(lepptbins, "leps_pt[0]", "Lepton p_{T} [GeV]"), 
        baseline_1l, procs_1l, plot_types).Weight(iweight).Tag("tt1l");
    } else {
      addSlices(pm, isr_syspt_slices, isr_syspt, nisrjet_bins, nisrjets, "ISR jet multiplicity", 
        baseline, iweight, procs, plot_types);
      addSlices(pm, nisrjet_slices, nisrjets, isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", 
        baseline, iweight, procs, {log_lumi});
      addSlices(pm, nisrjet_slices, nisrjets, isr_syspt_bins, isr_ht, "ISR H_{T} [GeV]", 
        baseline, iweight, procs, {log_lumi});

      //print ratio
      pm.Push<Hist1D>(Axis(nisrjet_bins, nisrjets, "ISR jet multiplicity"), 
        baseline, procs, plot_vals).Weight(iweight).Tag(isrtype+"_vals");

      pm.Push<Hist1D>(Axis(isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]"), 
        baseline, procs, vector<PlotOpt>{log_lumi}).Weight(iweight).Tag(isrtype);
      pm.Push<Hist1D>(Axis(isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]"), 
        baseline && "ht>300", procs, vector<PlotOpt>{log_lumi}).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(ptbins, isr_jetspt[0.], "Leading ISR jet p_{T} [GeV]"), 
      //   baseline && nisrjets>0., procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(ptbins, isr_jetspt[1], "2^{nd} ISR jet p_{T} [GeV]"), 
      //   baseline && nisrjets>1, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(ptbins, isr_jetspt[2], "3^{rd} ISR jet p_{T} [GeV]"), 
      //   baseline && nisrjets>2, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(ptbins, isr_jetspt[3], "4^{th} ISR jet p_{T} [GeV]"), 
      //   baseline && nisrjets>3, procs, plot_types).Weight(iweight).Tag(isrtype);

      // pm.Push<Hist1D>(Axis(lepptbins, "leps_pt[0]", "Lepton p_{T} [GeV]"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);

      // pm.Push<Hist1D>(Axis(20,0.,500., "met", "MET [GeV]"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(15,0.,1500., "ht", "H_{T} [GeV]"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(15,0.,1500., "mj14", "M_{J} [GeV]"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(nisrjet_bins, nisrjets50, "Number of 50 GeV ISR jets"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);
      // pm.Push<Hist1D>(Axis(nisrjet_bins, nisrjets75, "Number of 75 GeV ISR jets"), 
      //   baseline, procs, plot_types).Weight(iweight).Tag(isrtype);
      
      if(isrtype=="zisr"){
        pm.Push<Hist1D>(Axis(nisrjet_bins, nisrjets, "ISR jet multiplicity"), 
          baseline && "ht>200", procs, plot_types).Weight(iweight).Tag(isrtype);
        pm.Push<Hist1D>(Axis(isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]"), 
          baseline && "ht>200", procs, plot_types).Weight(iweight).Tag(isrtype);
      }
    } // if not wjets_tt1l
  } // Loop over weights

  if (single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);
}

NamedFunc::ScalarType zCandidate(const Baby &b){
  if (b.nels()==2) {
    if (b.elel_m()>80 && b.elel_m()<100) return 1;
  } else if (b.nmus()==2) {
    if (b.mumu_m()>80 && b.mumu_m()<100) return 1;
  } 
  return -1;
}

void addSlices(PlotMaker &pm, const vector<double> slices, NamedFunc svar,
               const vector<double> xbins, NamedFunc xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits){

  //add the inclusive version first
  pm.Push<Hist1D>(Axis(xbins, xvar, xlabel), baseline, proc, plot_types).Weight(weight).Tag(isrtype+"_incl");
  for(unsigned i(0); i<slices.size(); i++){
    NamedFunc cut = baseline && svar>=slices[i];
    if (i<(slices.size()-1)) cut = cut && svar<slices[i+1];

    string tag = CodeToPlainText(isrtype+"_"+svar.Name()+RoundNumber(slices[i],tag_digits).Data());
    pm.Push<Hist1D>(Axis(xbins, xvar, xlabel), cut, proc, plot_types).Weight(weight).Tag(tag);
  }
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

bool isGoodJet(const Baby &b, size_t ijet){
  return ijet<b.jets_pt()->size()
      && fabs(b.jets_eta()->at(ijet))<2.4
      && !b.jets_islep()->at(ijet);
}

NamedFunc::VectorType isrJetsPt(const Baby &b, float ptThresh){
  vector<double> isr_jetspt;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if (!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<ptThresh) continue;
    float dcsv = DeepCSVMedium_2017;
    if (year==2016) dcsv = DeepCSVMedium_2016;
    if (isrtype=="ttisr" && b.jets_csvd()->at(ijet)>dcsv) continue;
    isr_jetspt.push_back(b.jets_pt()->at(ijet));
  }
  std::sort(isr_jetspt.begin(), isr_jetspt.end(), std::greater<double>());
  return isr_jetspt;
}

NamedFunc::ScalarType isrSystemPt(const Baby &b){
    if (isrtype=="ttisr") return b.jetsys_nobd_pt();
    else return b.jetsys_pt();
}

NamedFunc::ScalarType nJetsWeights_ttisr(const Baby &b, bool use_baby_nisr, string wgtopt){
  if (b.w_isr()<0) return 1.; // Do not reweight Data

  if (wgtopt=="none") return b.weight()/b.w_isr();
  else if (b.type()<1000 || b.type()>=2000) return b.weight();

  int nisrjets = b.njets()-2;
  if (use_baby_nisr) nisrjets = b.nisr();
  
  float isr_wgt = 1.;
  if (wgtopt=="is2016") {
         if (nisrjets==1) isr_wgt = 0.920; //  +- 0.014
    else if (nisrjets==2) isr_wgt = 0.821; //  +- 0.020
    else if (nisrjets==3) isr_wgt = 0.715; //  +- 0.031
    else if (nisrjets==4) isr_wgt = 0.662; //  +- 0.051
    else if (nisrjets==5) isr_wgt = 0.561; //  +- 0.088
    else if (nisrjets>=6) isr_wgt = 0.511; //  +- 0.133
    isr_wgt *= 1.090; //normalization factor
  } else if (wgtopt=="is2017"){
         if (nisrjets==1) isr_wgt = 0.882; //  +- 0.014
    else if (nisrjets==2) isr_wgt = 0.792; //  +- 0.020
    else if (nisrjets==3) isr_wgt = 0.702; //  +- 0.031
    else if (nisrjets==4) isr_wgt = 0.648; //  +- 0.051
    else if (nisrjets==5) isr_wgt = 0.601; //  +- 0.088
    else if (nisrjets>=6) isr_wgt = 0.515; //  +- 0.133
    isr_wgt *= 1.090; //normalization factor
  }

  // if (use_baby_nisr) cout<<"Weight "<<b.type()<<" baby = "<<RoundNumber(b.weight(),6)<<" alg = "<<RoundNumber(wgt,6)<<endl;
  return b.weight()/b.w_isr()*isr_wgt;
}

NamedFunc::ScalarType nJetsWeights_visr(const Baby &b){
  if (b.ntrupv()<0) return 1.; // Do not reweight Data

  float wgt = b.weight()/b.w_isr();
  if((b.type()>=2000 && b.type()<3000) ||          //wjets
    (b.type()>=6000 && b.type()<7000)) return wgt; //dyjets

  int nisrjets(b.njets());
  // weights derived in DY+jets
  if      (nisrjets==0) return 0.981*wgt; //  +- 0.001
  else if (nisrjets==1) return 1.071*wgt; //  +- 0.001
  else if (nisrjets==2) return 1.169*wgt; //  +- 0.003
  else if (nisrjets==3) return 1.157*wgt; //  +- 0.007
  else if (nisrjets==4) return 1.014*wgt; //  +- 0.013
  else if (nisrjets==5) return 0.920*wgt; //  +- 0.025
  else if (nisrjets==6) return 0.867*wgt; //  +- 0.048
  else if (nisrjets>=7) return 0.935*wgt; //  +- 0.088
  else return wgt;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"year", required_argument, 0, 'y'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "y:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'y':
      year = atoi(optarg);
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}