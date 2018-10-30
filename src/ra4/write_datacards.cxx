#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>
#include <dirent.h>

#include "TSystem.h"
#include "TString.h"
#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"

using namespace std;
namespace {
  bool fake_PU = false;
  bool unblind = false;
  bool debug = false;
  double lumi = 135;
  int ibatch = -1;
  int nfiles = 10;
  string model = "T1tttt";
  string tag = "nominal";
  string outfolder = getenv("PWD");
  string mass_pts_str = "";
  TString nom_wgt = "weight*eff_trig"; // nominal weight to use
  enum SysType {kConst, kWeight, kSmear, kCorr, kMetSwap, kPU};
  TString syst = "all";

  const vector<double> v_data_npv{6.540e-06, 2.294e-05, 6.322e-05, 8.558e-05, 1.226e-04, 1.642e-04, 1.917e-04, 3.531e-04, 9.657e-04, 2.155e-03, 4.846e-03, 9.862e-03, 1.651e-02, 2.401e-02, 3.217e-02, 4.078e-02, 4.818e-02, 5.324e-02, 5.612e-02, 5.756e-02, 5.841e-02, 5.886e-02, 5.831e-02, 5.649e-02, 5.376e-02, 5.044e-02, 4.667e-02, 4.257e-02, 3.833e-02, 3.406e-02, 2.982e-02, 2.567e-02, 2.169e-02, 1.799e-02, 1.464e-02, 1.170e-02, 9.178e-03, 7.058e-03, 5.306e-03, 3.884e-03, 2.757e-03, 1.890e-03, 1.247e-03, 7.901e-04, 4.795e-04, 2.783e-04, 1.544e-04, 8.181e-05, 4.141e-05, 2.004e-05, 9.307e-06, 4.178e-06, 1.846e-06, 8.350e-07, 4.150e-07, 2.458e-07, 1.779e-07, 1.488e-07, 1.339e-07, 1.238e-07, 1.153e-07, 1.071e-07, 9.899e-08, 9.095e-08, 8.301e-08, 7.527e-08, 6.778e-08, 6.063e-08, 5.387e-08, 4.753e-08, 4.166e-08, 3.627e-08, 3.136e-08, 2.693e-08, 2.297e-08};
  TH1D h_data_npv("h_data_npv", "Data;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  TH1D h_mc_npv("h_mc_npv", "MC;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  double pu_low = 0.;
  double pu_high = 0.;
  bool do_syst = true;
}

class bindef {
public:
  bindef(string itag, string icut): tag(itag), cut(icut){};
  string tag, cut;
};

class sysdef {
public:
  sysdef(TString ilabel, TString itag, SysType isystype): label(ilabel), tag(itag), sys_type(isystype) {
    v_wgts = vector<TString>();
  }
  // as will appear in the AN latex table
  TString label;
  // as will appear in the file handed to ra4 stats
  TString tag;
  // Is it a const, a weight or does it actually change the analysis variables, e.g. like lumi, b_tag or JEC
  SysType sys_type;
  // if sys_type = kSmear, what is the index in e.g. sys_met, check in babymaker:
  // https://github.com/manuelfs/babymaker/blob/2c0d9b2bde517b0bb129b8b3afffa77a581123e1/bmaker/interface/utilities.hh#L17 
  // if sys_type = kCorr, what is the index in e.g. sys_met, where the shifted *Up* value is stored, assuming Down is Up+1
  size_t shift_index; 
  // if sys_type = kWeight, add all weights to be used
  vector<TString> v_wgts;
  // here, we will store where this systematic begins in the big yields & entires vectors that we get from getYields()
  // from there, indices are order like: nVariations*iBin + iVariation 
  size_t ind;
};

TString nom2sys_bin(TString ibin, size_t shift_index);
TString nom2genmet(TString ibin);
void GetOptions(int argc, char *argv[]);
void fillTtbarSys(ofstream &fsys);

vector<double> getYields(Baby_full &baby, const NamedFunc &baseline, const vector<NamedFunc> &bincuts,
                         vector<double> &yield, vector<double> &w2, double lumi, const TString &flag = "");

const NamedFunc ntop("ntop", [](const Baby &b) -> NamedFunc::ScalarType{
    int _ntop = 0;
    for (unsigned i(0); i<b.ak8jets_pt()->size(); i++) {
      if (b.ak8jets_pt()->at(i)<=300) continue;
      if (b.ak8jets_nom_raw_top()->at(i)>0.4) _ntop++;
    }
    return _ntop;
  });

const NamedFunc masstop("masstop", [](const Baby &b) -> NamedFunc::ScalarType{
    int _masstop = 0;
    for (unsigned i(0); i<b.ak8jets_pt()->size(); i++) {
      if (b.ak8jets_pt()->at(i)<=300) continue;
      if (b.ak8jets_m()->at(i)>150) _masstop++;
    }
    return _masstop;
  });

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  // Bear
  string foldersig = bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  string foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_abcd/";
  string folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/";
  // Shrew
  // string foldersig = bfolder+"/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_stdnj5/";
  // string foldermc = bfolder+"/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_stdnj5/";
  // string folderdata = "";

  GetOptions(argc, argv);

  if (do_syst && (Contains(foldersig, "merged_") || Contains(foldersig, "skim_"))) {
    cout<<"Systematics cannot be derived from skim!!"<<endl;
    exit(1);
  }
  gSystem->mkdir(outfolder.c_str(), kTRUE);


  vector<pair<string, string>> mass_pts;
  if (mass_pts_str!="") {
    size_t found = mass_pts_str.find(",");
    size_t start = 0;
    string _tmp = "";
    while (found != string::npos) {
      _tmp = mass_pts_str.substr(start, found-start);
      mass_pts.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
      start = found+1;
      found = mass_pts_str.find(",",start);
    } 
    _tmp = mass_pts_str.substr(start, found-start);
    mass_pts.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
  } else {
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (foldersig.c_str())) != NULL) {
      int ifile = 0;
      while ((ent = readdir (dir)) != NULL) {
        if (Contains(ent->d_name,".root")){
          if (ibatch==-1 || (ifile>=ibatch*nfiles && ifile<(ibatch+1)*nfiles)) {
            string _tmp = ent->d_name;
            string mglu = _tmp.substr(_tmp.find("ino-")+4,_tmp.find("_mLSP")-_tmp.find("ino-")-4);
            string mlsp = _tmp.substr(_tmp.find("LSP-")+4,_tmp.find("_Tune")-_tmp.find("LSP-")-4);
            mass_pts.push_back(make_pair(mglu, mlsp));
            cout<<"Including mass point: "<<mglu<<" "<<mlsp<<endl;
          }
          ifile++;
        }
      }
      closedir (dir);
    } else {
      cout<<"Could not find directory with signal samples"<<endl;
      exit(1);
    }
  }

  TString baseline("st>500 && met>200 && mj14>250 && njets>=6 && nbm>=1 && nleps==1 && nveto==0");
  NamedFunc filters = "st<10000 && pass_ra2_badmu && met/met_calo<5";
  if (Contains(tag, "masstop")) filters = filters && masstop>=1;
  else if (Contains(tag, "top")) filters = filters && ntop>=1;

  // --------------------------------------
  //            Processes
  //---------------------------------------
  set<string> other_files = {
    foldermc+"*_WJetsToLNu*",
     foldermc+"*_ST_*", 
    foldermc+"*_TTW*", foldermc+"*_TTZ*", foldermc+"*_TTGJets*", foldermc+"*_TTTT*",
    foldermc+"*DYJetsToLL*", foldermc+"*_ZJet*", foldermc+"*_ttHTobb_M125_*", 
    foldermc+"*_WH_HToBB*", foldermc+"*_ZH_HToBB*", foldermc+"*_WWTo*", 
    foldermc+"*_WZ*", foldermc+"*_ZZ_*",
    foldermc+"*QCD_HT*0_Tune*",foldermc+"*QCD_HT*Inf_Tune*"
  };

  auto proc_tt1l = Process::MakeShared<Baby_full>("tt_1l",    Process::Type::background, 1,  
    {foldermc+ "*_TTJets*SingleLept*.root"}, filters && "stitch_met && pass");  
  auto proc_tt2l = Process::MakeShared<Baby_full>("tt_2l",    Process::Type::background, 1,  
    {foldermc+ "*_TTJets*DiLept*.root"}, filters && "stitch_met && pass");  
  auto proc_other = Process::MakeShared<Baby_full>("other",    Process::Type::background, 1,  
    other_files, filters && "stitch_met && pass");

  vector<shared_ptr<Process> > bkg_procs = {proc_tt1l, proc_tt2l, proc_other};

  vector<shared_ptr<Process> > data_procs;
  data_procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    set<string>{folderdata+"*root"}, filters && "trig_ra4 && pass"));

  vector<shared_ptr<Process> > sig_procs;
  for (auto &imass: mass_pts)
    sig_procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::signal, kBlack,
      set<string>{foldersig+"/*mGluino-"+imass.first+"_mLSP-"+imass.second+"_*.root"}, filters));

  // --------------------------------------
  //            Binning
  //---------------------------------------

  string c_mt_r1 = "mt<=140";

  // N.B.: currently, changing naming convention would break writing the closure systematics derived from CRs
  vector<string> vl_met, vc_met, vc_mj_r1;
  vl_met.push_back("met200to350"); vc_met.push_back("met>200 && met<=350"); vc_mj_r1.push_back("mj14<=400");
  vl_met.push_back("met350to500"); vc_met.push_back("met>350 && met<=500"); vc_mj_r1.push_back(Contains(tag,"vary_mj")? "mj14<=450" : "mj14<=400");
  vl_met.push_back("met500");      vc_met.push_back("met>500");             vc_mj_r1.push_back(Contains(tag,"vary_mj")? "mj14<=500" : "mj14<=400");
  unsigned nbins_met(vl_met.size());

  vector<string> vl_nb, vc_nb;
  vl_nb.push_back("1b");     vc_nb.push_back("nbm==1");
  vl_nb.push_back("2b");     vc_nb.push_back("nbm==2");
  vl_nb.push_back("ge3b");   vc_nb.push_back("nbm>=3");
  unsigned nbins_nb(vl_nb.size());

  vector<string> vl_nj, vc_nj;
  vl_nj.push_back("6to8j");  vc_nj.push_back("njets>=6 && njets<=7");
  vl_nj.push_back("ge9j");   vc_nj.push_back("njets>=8");
  unsigned nbins_nj(vl_nj.size());

  vector<string> vl_mj, vc_mj;
  vl_mj.push_back("mj400toXXX");  vc_mj.push_back("mj14<=XXX");
  vl_mj.push_back("mjXXX");  vc_mj.push_back("mj14>XXX");
  unsigned nbins_mj(vl_mj.size());
  vector<string> v_metdep_midmj = {"500","650","800"};
  // vector<string> v_metdep_midmj = {"400","650"};
  if (v_metdep_midmj.size()!=vc_met.size()) {
    cout<<"ERROR: Intermediate MJ thresholds not specified for each MET bin"<<endl;
    exit(1);
  }

  vector<bindef> vbins; 
  for (unsigned imet(0); imet<nbins_met; imet++){
    string _label(""), _cut("");
    if (!Contains(tag,"noabcd")) {
      _label = "r1_"+vl_met[imet];
      _cut = vc_met[imet]+"&&"+c_mt_r1+"&&"+vc_mj_r1[imet];
      vbins.push_back(bindef(_label, _cut));
      for (unsigned inb(0); inb<nbins_nb; inb++){
        for (unsigned inj(0); inj<nbins_nj; inj++){
          for (unsigned imj(0); imj<nbins_mj; imj++){
            _label = "r2_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
            _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
            _cut = vc_met[imet]+"&&"
                   +c_mt_r1+"&&"+CopyReplaceAll(vc_mj_r1[imet],"<=",">")+"&&"
                   +vc_nb[inb]+"&&"+vc_nj[inj];
            _cut += "&&"+CopyReplaceAll(vc_mj[imj],"XXX",v_metdep_midmj[imet]);
            vbins.push_back(bindef(_label, _cut));
          }
        }
      }
      _label = "r3_"+vl_met[imet];
      _cut = vc_met[imet]+"&&"+CopyReplaceAll(c_mt_r1,"<=",">")+"&&"+vc_mj_r1[imet];
      vbins.push_back(bindef(_label, _cut));
    }
    for (unsigned inb(0); inb<nbins_nb; inb++){
      for (unsigned inj(0); inj<nbins_nj; inj++){
        for (unsigned imj(0); imj<nbins_mj; imj++){
          _label = "r4_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
          _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
          _cut = vc_met[imet]+"&&"
                 +CopyReplaceAll(c_mt_r1,"<=",">")+"&&"+CopyReplaceAll(vc_mj_r1[imet],"<=",">")+"&&"
                 +vc_nb[inb]+"&&"+vc_nj[inj];
          _cut += "&&"+CopyReplaceAll(vc_mj[imj],"XXX",v_metdep_midmj[imet]);
          vbins.push_back(bindef(_label, _cut));
        }
      }
    }
  }
  unsigned nbins(vbins.size());

  // ---------------------------------------------------
  //     Vector of signal systematic variations
  //----------------------------------------------------
  vector<sysdef> v_sys;
  // order as they will appear in latex table
  // *Nominal must stay in the first spot!!* (will be skipped in table)
  v_sys.push_back(sysdef("Nominal", "nominal", kWeight)); 
  v_sys.back().v_wgts.push_back("1.");

  v_sys.push_back(sysdef("Gen vs reco MET FS", "fs_genmet",kMetSwap));

  if (do_syst) {
    v_sys.push_back(sysdef("Lepton efficiency", "lepeff", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_lep["+to_string(i)+"]/w_lep");
    v_sys.push_back(sysdef("Lepton efficiency FS", "fs_lepeff", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_lep["+to_string(i)+"]/w_fs_lep");
    v_sys.push_back(sysdef("Trigger efficiency", "trig", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_trig["+to_string(i)+"]/eff_trig"); 
    v_sys.push_back(sysdef("B-tag efficiency", "bctag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_bctag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("B-tag efficiency FS", "fs_bctag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_bctag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Mistag efficiency", "udsgtag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_udsgtag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Mistag efficiency FS", "fs_udsgtag",kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_udsgtag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Jet energy corrections", "jec", kCorr));
    v_sys.back().shift_index = 1; // JEC Up index in sys_met, etc.
    // v_sys.push_back(sysdef("Jet energy resolution", "jer", kSmear));
    // v_sys.back().shift_index = 0; // JER index in sys_met, etc.
    // v_sys.push_back(sysdef("PDFs", "pdf", kWeight));
    // for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_pdf["+to_string(i)+"]");
    // v_sys.push_back(sysdef("RMS PDFs", "rms_pdf", kWeight));
    // for (size_t i = 0; i<100; ++i) v_sys.back().v_wgts.push_back("w_pdf["+to_string(i)+"]");
    v_sys.push_back(sysdef("QCD scales", "murf",kWeight));
    for (size_t i = 0; i<2; ++i) {
      v_sys.back().v_wgts.push_back("sys_mur["+to_string(i)+"]");
      v_sys.back().v_wgts.push_back("sys_muf["+to_string(i)+"]");
      v_sys.back().v_wgts.push_back("sys_murf["+to_string(i)+"]");
    }
    v_sys.push_back(sysdef("ISR", "isr", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_isr["+to_string(i)+"]/w_isr");
    v_sys.push_back(sysdef("Jet ID FS", "jetid", kConst));
    v_sys.back().v_wgts.push_back("0.01");
    // v_sys.push_back(sysdef("Pile up", "pu", kPU));
    v_sys.push_back(sysdef("Luminosity", "lumi", kConst));
    v_sys.back().v_wgts.push_back("0.026");
  } 

  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  //  No further change beyond this point needed when just adding systematics
  

  //------------------------------------------
  //       Vector of all cuts
  //------------------------------------------
  sysdef nom = v_sys[0];
  if (nom.tag != "nominal"){
    cerr<<" The first entry in the v_sys vector must be the nominal"<<endl;
    exit(1);
  }
  vector<TableRow> cuts_nosys;
  for (auto &bin: vbins) 
    cuts_nosys.emplace_back(TableRow("", baseline+"&&"+bin.cut,0,0,nom_wgt));

  vector<TableRow> cuts;
  for (auto &sys: v_sys) {
    sys.ind = cuts.size(); 
    if (sys.sys_type == kConst){
      continue;
    } else if (sys.sys_type == kWeight) {
      for (auto &bin: vbins) {
        for (auto &wgt: sys.v_wgts) {
          cuts.emplace_back(TableRow("", baseline+"&&"+bin.cut,0,0,nom_wgt+"*"+wgt));
        }
      }
    } else if (sys.sys_type == kCorr || sys.sys_type == kSmear) {
      for (auto &bin: vbins) {
        cuts.emplace_back(TableRow("", nom2sys_bin(baseline+"&&"+bin.cut, sys.shift_index),0,0,nom_wgt));
        if (sys.sys_type == kCorr) { //if it is a correction, need to push the 'down' variation as well
          cuts.emplace_back(TableRow("", nom2sys_bin(baseline+"&&"+bin.cut, sys.shift_index+1),0,0,nom_wgt));
        }
      }
    } else if (sys.sys_type == kMetSwap){
      for (auto &bin: vbins) {
        cuts.emplace_back(TableRow("", nom2genmet(baseline+"&&"+bin.cut),0,0,nom_wgt));
      }
    } else if (sys.sys_type == kPU) {
      for(const auto &bin: vbins){
        cuts.emplace_back(TableRow("", baseline+"&&"+bin.cut+"&&"+"ntrupv<=20",0,0,nom_wgt));
        cuts.emplace_back(TableRow("", "ntrupv<=20",0,0,nom_wgt));
        cuts.emplace_back(TableRow("", baseline+"&&"+bin.cut+"&&"+"ntrupv>=21",0,0,nom_wgt));
        cuts.emplace_back(TableRow("", "ntrupv>=21",0,0,nom_wgt));
      }
    }
  }

  //------------------------------------------
  //       Get BKG yields
  //------------------------------------------
  PlotMaker pm;
  pm.Push<Table>("tmc",  cuts_nosys, bkg_procs, true, false);
  pm.Push<Table>("tsig",  cuts, sig_procs, true, false);
  if (unblind) pm.Push<Table>("tdata",  cuts_nosys, data_procs, true, false);
  pm.multithreaded_ = true;
  pm.min_print_ = true; 
  pm.MakePlots(lumi);  

  vector<GammaParams> bkg_params;
  Table * yield_table;
  yield_table = static_cast<Table*>(pm.Figures()[0].get());
  bkg_params = yield_table->BackgroundYield(lumi);
  vector<float> bkg_yields;
  for (unsigned ipar(0); ipar< bkg_params.size(); ipar++)
    bkg_yields.push_back(bkg_params[ipar].Yield());  

  // calculate kappas
  vector<float> vkappas, vkappas_unc;
  if (!Contains(tag,"noabcd")) {
    vector<float> powers = {1, -1, -1, 1};
    // order of for loops must match the one used to define cut vector, i.e. met, nb, nj, mj!
    for (unsigned imet(0); imet<nbins_met; imet++){
      for (unsigned inb(0); inb<nbins_nb; inb++){
        for (unsigned inj(0); inj<nbins_nj; inj++){
          for (unsigned imj(0); imj<nbins_mj; imj++){
            unsigned r1_idx = imet*(2*nbins_nb*nbins_nj*nbins_mj+2);  
            unsigned r2_idx = r1_idx + inb*nbins_nj*nbins_mj + inj*nbins_mj + imj + 1;
            unsigned r3_idx = r1_idx + nbins_nb*nbins_nj*nbins_mj + 1; 
            unsigned r4_idx = r3_idx + inb*nbins_nj*nbins_mj + inj*nbins_mj + imj + 1;  
            // repackage yields as required for input to calcKappa
            vector<vector<float>> _entries, _weights;
            for (auto &_idx: vector<unsigned>{r1_idx, r2_idx, r3_idx, r4_idx}){
              _entries.push_back(vector<float>{static_cast<float>(bkg_params[_idx].NEffective())});
              _weights.push_back(vector<float>{static_cast<float>(bkg_params[_idx].Weight())});
            }
            float kappa_up(0), kappa_dn(0);
            double kappa = calcKappa(_entries, _weights, powers, kappa_dn, kappa_up);
            vkappas.push_back(kappa);
            vkappas_unc.push_back(kappa_up>kappa_dn ? kappa_up:kappa_dn);
          } //mj
        } //nj
      } // nb
    } // met bin loop
  } // if abcd

  //------------------------------------------
  //       Get SIG yields
  //------------------------------------------
  vector<vector<GammaParams>> sig_params;
  yield_table = static_cast<Table*>(pm.Figures()[1].get());
  for (auto &isig: sig_procs) {
    sig_params.push_back(yield_table->Yield(isig.get(), lumi));
  }

  //------------------------------------------
  //       Get DATA yields
  //------------------------------------------
  vector<float> data_yields;
  if(unblind){
    yield_table = static_cast<Table*>(pm.Figures()[2].get());
    vector<GammaParams> data_params = yield_table->DataYield();
    for (auto &ipar: data_params)
      data_yields.push_back(ipar.Yield());    
  } else {
      data_yields = bkg_yields;   
  }

  for (unsigned isig(0); isig<sig_params.size(); isig++) {
    //    calculate average of yields with GEN and RECO MET for signal
    // -----------------------------------------------------------------------
    vector<GammaParams> nom_met_avg;
    for (auto &sys: v_sys) {
      if (sys.sys_type == kMetSwap) {
        for (size_t ibin = 0; ibin<nbins; ++ibin) {
          GammaParams tmp_gps;
          tmp_gps.SetYieldAndUncertainty(0.5*(sig_params[isig][sys.ind + ibin].Yield()+sig_params[isig][ibin].Yield()),
                        max(sig_params[isig][sys.ind + ibin].Uncertainty(), sig_params[isig][ibin].Uncertainty()));
          nom_met_avg.push_back(tmp_gps);
          // nom_met_avg.push_back(0.5*(sig_params[isig][sys.ind + ibin]+sig_params[isig][ibin]));
        }
      }
    }
    //   Writing datacard header
    //---------------------------------------
    TString outpath = outfolder+"/datacard_SMS-"+model;
    outpath += "_mGluino-"+mass_pts[isig].first+"_mLSP-"+mass_pts[isig].second;
    outpath += "_"+RoundNumber(lumi,1).ReplaceAll(".","p")+"ifb_"+tag+".txt";
    if (!do_syst)  outpath.ReplaceAll(".txt*","_nosys.txt");
    cout<<"open "<<outpath<<endl;
    unsigned wname(28), wdist(7), wbin(31);
    for (size_t ibin(0); ibin<nbins; ibin++) 
      if(vbins[ibin].tag.length() > wbin) wbin = vbins[ibin].tag.length();
    wbin+=1;
    unsigned digit = 2;
    if (unblind) digit = 0;

    // --------- write header
    ofstream fcard(outpath);
    fcard<<"imax "<<nbins<<"  number of channels\n";
    fcard<<"jmax 1  number of backgrounds\n";
    fcard<<"kmax *  number of nuisance parameters\n";
    fcard<<"shapes * * FAKE\n";
    fcard<<endl<<setw(wname)<<"bin"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<" "<<setw(wbin)<<vbins[ibin].tag;
    fcard<<endl<<setw(wname)<<"Observation"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<" "<<setw(wbin)<<RoundNumber(data_yields[ibin],digit);

    fcard<<endl<<endl<<setw(wname)<<"bin"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<vbins[ibin].tag<<setw(wbin)<<vbins[ibin].tag;
    fcard<<endl<<setw(wname)<<"process"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<"sig"<<setw(wbin)<<"bkg";
    fcard<<endl<<setw(wname)<<"process"<<setw(wdist)<<" ";
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<"0"<<setw(wbin)<<"1";
    fcard<<endl<<setw(wname)<<"rate"<<setw(wdist)<<" ";
    if (!Contains(tag,"noabcd")) {
      for (size_t ibin(0); ibin<nbins; ibin++) 
        fcard<<setw(wbin)<<Form("%.2f",nom_met_avg[ibin].Yield())<<setw(wbin)<<"1";
    } else {
      for (size_t ibin(0); ibin<nbins; ibin++) 
        fcard<<setw(wbin)<<Form("%.2f",nom_met_avg[ibin].Yield())<<setw(wbin)<<bkg_yields[ibin];     
    }
    fcard<<endl<<endl;
    cout<<"Wrote headers"<<endl;

    //--------- Signal statistical uncertainties ----------------------------
    for (size_t ibin(0); ibin<nbins; ibin++) {
      fcard<<setw(wname)<<"statsig_"+vbins[ibin].tag<<setw(wdist)<<"lnN";
      TString sig_stat = Form("%.2f",1+nom_met_avg[ibin].Uncertainty()/nom_met_avg[ibin].Yield());
      for (size_t jbin(0); jbin<nbins; jbin++) {
        if (ibin==jbin) fcard<<setw(wbin)<<sig_stat<<setw(wbin)<<"-";
        else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      fcard<<endl;
    }
    cout<<"Wrote signal stat. uncertainties"<<endl;

    // ------------ Closure uncertainties
    // ordered as: 2l lownj, 2l high nj, 5j lowmet, 5j mid met
    vector<float> cr_unc = {1.06, 1.16, 1.16, 1.4}; // nominal 35.9 ifb, based on data CR 
    if (lumi>130) cr_unc = {1.03, 1.07, 1.08, 1.17}; // nominal 135 ifb, based on data CR 
    if (Contains(tag,"vary_mj")) cr_unc = {1.03, 1.07, 1.08, 1.21}; // nominal 35.9, based on data CR 

    if (do_syst){
      fcard<<endl<<setw(wname)<<"dilep_lownj"<<setw(wdist)<<"lnN";
      for (size_t ibin(0); ibin<nbins; ibin++) {
        if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,"6to8j")) 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc[0],2);
        else 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      fcard<<endl<<setw(wname)<<"dilep_highnj"<<setw(wdist)<<"lnN";
      for (size_t ibin(0); ibin<nbins; ibin++) {
        if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,"ge9j")) 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc[1],2);
        else 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      fcard<<endl<<setw(wname)<<"fivejet_lowmet"<<setw(wdist)<<"lnN";
      for (size_t ibin(0); ibin<nbins; ibin++) {
        if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,"met200to350")) 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc[2],2);
        else 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      fcard<<endl<<setw(wname)<<"fivejet_highmet"<<setw(wdist)<<"lnN";
      for (size_t ibin(0); ibin<nbins; ibin++) {
        if(Contains(vbins[ibin].tag,"r2_") && 
          (Contains(vbins[ibin].tag,"met350to500") || Contains(vbins[ibin].tag,"met500"))) 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc[3],2);
        else 
          fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      cout<<"Wrote CR-based closure uncertainties"<<endl;
    }

    //calculate uncertainties and write results to three files
    ofstream fsys(outpath.ReplaceAll("datacard_","sys_"));
    cout<<"Writing to "<<outpath<<endl;
    fillTtbarSys(fsys);
    if(do_syst){
      for (auto &sys: v_sys) {
        if (sys.tag != "nominal") {
          if (sys.tag != "rms_pdf") fsys<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
          fcard<<setw(wname)<<sys.tag<<setw(wdist)<<"lnN";
        }
        for (size_t ibin = 0; ibin<nbins; ++ibin) {
          const double nom_yield(sig_params[isig][ibin].Yield());
          // size_t nwgts = sys.v_wgts.size();
          double up(0.), dn(0.); 
          if (sys.sys_type == kConst) {
            up = stod(sys.v_wgts[0].Data());
            dn = -up;
          } else if (sys.sys_type == kWeight) {
            if (sys.tag == "nominal") {
              continue;
            } else if (sys.tag == "rms_pdf") { 
              // double sumw2(0), mean(0);
              // for (size_t iwgt = 0; iwgt<nwgts; ++iwgt) {
              //   sumw2 += pow(sig_params[isig][sys.ind + nwgts*ibin + iwgt],2);
              //   mean += sig_params[isig][sys.ind + nwgts*ibin + iwgt];
              // }
              // mean = mean/nwgts;
              // up = sqrt((sumw2-nwgts*pow(mean,2))/(nwgts-1))/nom_yield;  // RMS
              // dn = -up;
            } else if (sys.tag == "murf") {
              //max/min of all weights mur_up, muf_up and murf_up
              // up = *max_element(sig_params[isig].begin() + sys.ind + nwgts*ibin, 
              //                   sig_params[isig].begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; 
              // dn = *min_element(sig_params[isig].begin() + sys.ind + nwgts*ibin, 
              //                   sig_params[isig].begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; 
            } else {
              up = sig_params[isig][sys.ind + 2*ibin].Yield()/nom_yield - 1;
              dn = sig_params[isig][sys.ind + 2*ibin + 1].Yield()/nom_yield - 1;
            }
          } else if (sys.sys_type == kSmear) {
            up = sig_params[isig][sys.ind + ibin].Yield()/nom_yield - 1;
            dn = -up;
          } else if (sys.sys_type == kMetSwap) {
            //Use average of met yield and met_tru yield as central value
            dn = sig_params[isig][sys.ind + ibin].Yield();
            dn = dn/(0.5*(sig_params[isig][sys.ind + ibin].Yield()+nom_yield)) - 1;
            up = -dn;
          } else if (sys.sys_type == kCorr) {
            up = sig_params[isig][sys.ind + 2*ibin].Yield()/nom_yield - 1;
            dn = sig_params[isig][sys.ind + 2*ibin + 1].Yield()/nom_yield - 1;
          } else if (sys.sys_type == kPU ) {
            double eff_low  = sig_params[isig][sys.ind+4*ibin+0].Yield()/sig_params[isig][sys.ind+4*ibin+1].Yield();
            double eff_high = sig_params[isig][sys.ind+4*ibin+2].Yield()/sig_params[isig][sys.ind+4*ibin+3].Yield();
            double m = (eff_high-eff_low)/(pu_high-pu_low);
            double b = (eff_low*pu_high-eff_high*pu_low)/(pu_high-pu_low);
            double eff_data = 0., eff_mc = 0.;
            for(size_t i = 0; i < v_data_npv.size(); ++i){
              double fx = m*static_cast<double>(i)+b;
              eff_data += fx*h_data_npv.GetBinContent(i+1);
              eff_mc += fx*h_mc_npv.GetBinContent(i+1);
            }
            up = (eff_data-eff_mc)/eff_mc;
            dn = -up;

            //Temporary stand-in until better method available
            if(fake_PU){
              if((Contains(vbins.at(ibin).tag,"lowmet") || Contains(vbins.at(ibin).tag,"medmet"))
                 && (Contains(vbins.at(ibin).tag,"lownj") || Contains(vbins.at(ibin).tag,"r1_") || Contains(vbins.at(ibin).tag,"r3_"))){
                up = 0.1;
                dn = 0.1;
              }else{
                up = 0.15;
                dn = 0.15;
              }
            }
          }
          // convert to ra4_stats input and write to file
          double ln = (up>0 ? 1:-1)*max(up>0 ? up : (1/(1+up)-1), dn>0 ? dn : (1/(1+dn)-1));
          if (sys.sys_type == kConst) ln = up;
          if (sys.tag !="rms_pdf") {
            if(sys.tag.Contains("trig")){
              fsys<<"    " <<left<<setw(25)<<vbins[ibin].tag <<" "<<right<<setw(10)<<Form("%.3f",ln) <<endl;
            } else {
              fsys<<"    " <<left<<setw(25)<<vbins[ibin].tag <<" "<<right<<setw(10)<<Form("%.2f",ln) <<endl;
            }
          }

          // write systematics to datacard
          ln = max(up>0 ? 1+up : 1/(up+1), dn>0 ? 1+dn : 1/(dn+1));
          if (std::isnan(ln) || std::isinf(ln)) {
            cout <<" Found bad unc. set to 0 -> "<<left<<setw(10)<<sys.tag <<left<<setw(10)<<vbins[ibin].tag 
                 <<" "<<right<<setprecision(0)<<setw(25)<<sig_params[isig][ibin].NEffective() 
                 <<" "<<setprecision(5)<<setw(15)<<sig_params[isig][ibin].Yield() 
                 <<" "<<setprecision(10)<<setw(15)<<sig_params[isig][ibin].Weight()<<endl;  
            ln = 0;
          } 
          if (sys.sys_type == kConst) ln = 1+up;
          fcard<<setw(wbin)<<Form("%.2f",ln)<<setw(wbin)<<"-";
        } // loop over bins
        fcard<<endl;
      } // loop over systematics
    }
    fsys.close();

    //   Writing datacard param section
    //---------------------------------------
    if (!Contains(tag,"noabcd")) {
      unsigned iyield(0);
      wbin +=3; // to allow the label to fit even with "rp_" in front
      // order of for loops must match the one used to define cut vector, i.e. met, nb, nj!
      if(debug) cout <<endl<<left<<setw(25)<<"bin"
                     <<right<<setw(15)<<"Sig."
                     <<right<<setw(15)<<"Sig. (MET avg)"
                     <<right<<setw(15)<<"Bkg. MC"
                     <<right<<setw(10)<<"Data"
                     <<endl;
      for (unsigned imet(0); imet<nbins_met; imet++){
        ostringstream kappas; kappas.clear();
        string _label = "r1_"+vl_met[imet];
        fcard<<"rp_"<<left<<setw(wbin)<<_label<<setw(10)<<"rateParam"<<left<<setw(wbin)<<_label<<"bkg "
             <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
             <<(data_yields[iyield]<50 ? " [0,100]":"")<<endl;
        if(debug) cout <<left<<setw(wbin)<<_label 
                       <<right<<setw(15)<<RoundNumber(sig_params[0][iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(nom_met_avg[iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(1+nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(bkg_yields[iyield],2)
                       <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
                       <<endl;
        iyield++;
        for (unsigned inb(0); inb<nbins_nb; inb++){
          for (unsigned inj(0); inj<nbins_nj; inj++){
            for (unsigned imj(0); imj<nbins_mj; imj++){
              _label = "r2_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
              _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
              fcard<<"rp_"<<left<<setw(wbin)<<_label<<setw(10)<<"rateParam"<<left<<setw(wbin)<<_label<<"bkg "
                   <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
                   <<(data_yields[iyield]<50 ? " [0,100]":"")<<endl;
              if(debug) cout <<left<<setw(wbin)<<_label
                             <<right<<setw(15)<<RoundNumber(sig_params[0][iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(nom_met_avg[iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(1+nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(bkg_yields[iyield],2)
                             <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
                             <<endl;
              iyield++;
            }
          }
        }
        _label = "r3_"+vl_met[imet];
        fcard<<"rp_"<<left<<setw(wbin)<<_label<<setw(10)<<"rateParam"<<left<<setw(wbin)<<_label<<"bkg "
             <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
             <<(data_yields[iyield]<50 ? " [0,100]":"")<<endl;
        if(debug) cout <<left<<setw(wbin)<<_label
                       <<right<<setw(15)<<RoundNumber(sig_params[0][iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(nom_met_avg[iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(1+nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),2)
                       <<right<<setw(15)<<RoundNumber(bkg_yields[iyield],2)
                       <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
                       <<endl;
        iyield++;
        for (unsigned inb(0); inb<nbins_nb; inb++){
          for (unsigned inj(0); inj<nbins_nj; inj++){
            for (unsigned imj(0); imj<nbins_mj; imj++){
              _label = "r4_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
              _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
              fcard<<"rp_"<<left<<setw(wbin)<<_label<<setw(10)<<"rateParam"<<left<<setw(wbin)<<_label<<"bkg "
              <<"(@0*@1/@2)*@3 "<<"rp_r3_"+vl_met[imet]<<",rp_"+CopyReplaceAll(_label,"r4_","r2_")
              <<",rp_r1_"+vl_met[imet]<<","+CopyReplaceAll(_label,"r4_","kappa_")
              <<endl;
              if(debug) cout <<left<<setw(wbin)<<_label
                             <<right<<setw(15)<<RoundNumber(sig_params[0][iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(nom_met_avg[iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(1+nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),2)
                             <<right<<setw(15)<<RoundNumber(bkg_yields[iyield],2)
                             <<right<<setw(10)<<RoundNumber(data_yields[iyield],digit)
                             <<endl;
              unsigned ikappa = imet*nbins_nb*nbins_nj*nbins_mj +inb*nbins_nj*nbins_mj +inj*nbins_mj+imj;
              kappas<<left<<setw(wbin)<<CopyReplaceAll(_label,"r4_","kappa_")<<setw(10)<<"param "
                   <<RoundNumber(vkappas[ikappa], 2)<<"  "<<RoundNumber(vkappas_unc[ikappa], 2)<<endl;
              iyield++;
            }
          }
        } 
        fcard<<kappas.str();
        if (debug) cout<<kappas.str();
      } // met bin loop
    } // if doing an ABCD

    fcard.close();
  } // loop over signal points

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

TString nom2sys_bin(TString ibin, size_t shift_index){
  ibin.ReplaceAll("met", "sys_met["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mt", "sys_mt["+to_string(shift_index)+"]");
  ibin.ReplaceAll("st", "sys_st["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mj14", "sys_mj14["+to_string(shift_index)+"]");
  ibin.ReplaceAll("njets", "sys_njets["+to_string(shift_index)+"]");
  ibin.ReplaceAll("nbm", "sys_nbm["+to_string(shift_index)+"]");
  return ibin;
}

TString nom2genmet(TString ibin){
  ibin.ReplaceAll("met", "met_tru");
  ibin.ReplaceAll("met_tru/met_tru_calo", "met/met_calo");
  return ibin;

}

void fillTtbarSys(ofstream &fsys){
    fsys << "SYSTEMATIC dilep_lownj" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_lownj_1b    0.06" << endl;
    fsys << "  r2_lowmet_lownj_2b    0.06" << endl;
    fsys << "  r2_lowmet_lownj_3b    0.06" << endl;
    fsys << "  r2_medmet_lownj_1b    0.06" << endl;
    fsys << "  r2_medmet_lownj_2b    0.06" << endl;
    fsys << "  r2_medmet_lownj_3b    0.06" << endl;
    fsys << "  r2_highmet_lownj_1b   0.06" << endl;
    fsys << "  r2_highmet_lownj_2b   0.06" << endl;
    fsys << "  r2_highmet_lownj_3b   0.06" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC dilep_highnj" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_highnj_1b   0.16" << endl;
    fsys << "  r2_lowmet_highnj_2b   0.16" << endl;
    fsys << "  r2_lowmet_highnj_3b   0.16" << endl;
    fsys << "  r2_medmet_highnj_1b   0.16" << endl;
    fsys << "  r2_medmet_highnj_2b   0.16" << endl;
    fsys << "  r2_medmet_highnj_3b   0.16" << endl;
    fsys << "  r2_highmet_highnj_1b  0.16" << endl;
    fsys << "  r2_highmet_highnj_2b  0.16" << endl;
    fsys << "  r2_highmet_highnj_3b  0.16" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC fivejet_lowmet" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_lownj_1b    0.15" << endl;
    fsys << "  r2_lowmet_lownj_2b    0.15" << endl;
    fsys << "  r2_lowmet_lownj_3b    0.15" << endl;
    fsys << "  r2_lowmet_highnj_1b   0.15" << endl;
    fsys << "  r2_lowmet_highnj_2b   0.15" << endl;
    fsys << "  r2_lowmet_highnj_3b   0.15" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC fivejet_highmet" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_medmet_lownj_1b    0.37" << endl;
    fsys << "  r2_medmet_lownj_2b    0.37" << endl;
    fsys << "  r2_medmet_lownj_3b    0.37" << endl;
    fsys << "  r2_medmet_highnj_1b   0.37" << endl;
    fsys << "  r2_medmet_highnj_2b   0.37" << endl;
    fsys << "  r2_medmet_highnj_3b   0.37" << endl;
    fsys << "  r2_highmet_lownj_1b   0.37" << endl;
    fsys << "  r2_highmet_lownj_2b   0.37" << endl;
    fsys << "  r2_highmet_lownj_3b   0.37" << endl;
    fsys << "  r2_highmet_highnj_1b  0.37" << endl;
    fsys << "  r2_highmet_highnj_2b  0.37" << endl;
    fsys << "  r2_highmet_highnj_3b  0.37" << endl;
    fsys << endl;
}

vector<double> getYields(Baby_full &baby, const NamedFunc &/*baseline*/, const vector<NamedFunc> &bincuts,
                         vector<double> &yield, vector<double> &w2, const TString &flag){
  for(size_t i = 0; i <v_data_npv.size(); ++i){
    h_data_npv.SetBinContent(i+1, v_data_npv.at(i));
    h_data_npv.SetBinError(i+1, 0.);
    h_mc_npv.SetBinContent(i+1, 0.);
    h_mc_npv.SetBinError(i+1, 0.);
  }
  
  vector<double> entries = vector<double>(bincuts.size(), 0);
  yield = vector<double>(bincuts.size(), 0);
  w2 = yield;
  long nentries = baby.GetEntries();

  bool isData(false), isBkg(false);
  baby.GetEntry(0);
  if (baby.type()<1000) isData = true;
  else if (baby.type()<100e3) isBkg = true;
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    
    if (isBkg || isData) {
      if(!baby.pass()) continue;
    } else {
      h_mc_npv.Fill(baby.ntrupv(), baby.weight()*baby.eff_trig());
    }

    if (isBkg && !baby.stitch_met()) continue;

    if (isData && !baby.trig_ra4()) continue;
    
    for(size_t ind = 0; ind<bincuts.size(); ++ind){
      float wgt = bincuts.at(ind).GetScalar(baby);
      if(wgt != 0.){
        ++entries.at(ind);
        
        if(flag=="jer_tail"){
          double jet_res_min = *min_element(baby.jets_pt_res()->begin(), baby.jets_pt_res()->end());
          double jet_res_max = *max_element(baby.jets_pt_res()->begin(), baby.jets_pt_res()->end());
          if((jet_res_min>0&&jet_res_min<0.675) || jet_res_max>1.391)
            wgt *= 1.5;
        }

        yield.at(ind) += wgt;
        w2.at(ind) += wgt*wgt;
      }
    }
  } // Loop over entries
  if (!isData) {
    for(size_t ind = 0; ind<bincuts.size(); ++ind){ 
      yield.at(ind) *= lumi;
      w2.at(ind) *= pow(lumi, 2);
    }
  }
  if (!isData && !isBkg) {
    h_data_npv.Scale(1./h_data_npv.Integral());
    h_mc_npv.Scale(1./h_mc_npv.Integral());
    pu_low = 0.;
    pu_high = 0.;
    double norm = 0.;
    for(size_t npv = 0; npv <= 20 && npv < v_data_npv.size(); ++npv){
      pu_low += npv*h_mc_npv.GetBinContent(npv+1);
      norm += h_mc_npv.GetBinContent(npv+1);
    }
    pu_low /= norm;
    norm = 0.;
    for(size_t npv = 21; npv < v_data_npv.size(); ++npv){
      pu_high += npv*h_mc_npv.GetBinContent(npv+1);
      norm += h_mc_npv.GetBinContent(npv+1);
    }
    pu_high /= norm;
  }
  return entries;
}

void GetOptions(int argc, char *argv[]){
  string blah;
  while(true){
    static struct option long_options[] = {
      {"syst", required_argument, 0, 's'},
      {"model", required_argument, 0, 'm'},
      {"tag", required_argument, 0, 't'},
      {"mass_pts", required_argument, 0, 'p'},
      {"lumi", required_argument, 0, 'l'},
      {"ibatch", required_argument, 0, 'b'},
      {"nfiles", required_argument, 0, 'n'},
      {"fake_PU", no_argument, 0, 0},
      {"no_syst", no_argument, 0, 0},
      {"unblind", no_argument, 0, 'u'},
      {"debug", no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:m:t:o:p:l:b:n:ud", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's': syst = optarg; break;
    case 'm': model = optarg; break;
    case 't': tag = optarg; break;
    case 'o': outfolder = optarg; break;
    case 'p': mass_pts_str = optarg; break;
    case 'l': lumi = atof(optarg); break;
    case 'b': ibatch = atoi(optarg); break;
    case 'n': nfiles = atoi(optarg); break;
    case 'u': unblind = true; break;
    case 'd': debug = true; break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "fake_PU"){
        fake_PU = true;
      }else if(optname == "no_syst"){
        do_syst = false;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default: printf("Bad option! getopt_long returned character code 0%o\n", opt); break;
    }
  }
}