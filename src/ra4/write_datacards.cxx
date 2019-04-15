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
  bool unblind = false;
  bool debug = false;
  int year = 0;
  string model = "T1tttt";
  string xoption = "nom";
  string outfolder = getenv("PWD");
  string mass_pts_str = "";
  enum SysType {kConst, kWeight, kSmear, kCorr, kMetSwap, kPU};
  bool do_syst = true;
  float lumi = 1.; // should not be changed, actual lumi assigned by the weight
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

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  if (do_syst) {
    cout<<"INFO:: Systematics are ON. Make sure to run on appropriate babies, i.e. unskimmed or skim_sys_abcd!!"<<endl;
  }
  gSystem->mkdir(outfolder.c_str(), kTRUE);


  // -----------------------------------------------------------
  //        Determine signal mass points
  //------------------------------------------------------------
  vector<pair<string, string>> mass_pts;
  if (mass_pts_str!="") {
    size_t found = mass_pts_str.find(",");
    size_t start = 0;
    string _tmp = "";
    while (found != string::npos) {
      _tmp = mass_pts_str.substr(start, found-start);
      mass_pts.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
      cout<<"Adding mass point: mgluino = "<<mass_pts.back().first<<" mlsp = "<<mass_pts.back().second<<endl;
      start = found+1;
      found = mass_pts_str.find(",",start);
    } 
    _tmp = mass_pts_str.substr(start, found-start);
    mass_pts.push_back(make_pair(_tmp.substr(0,_tmp.find("_")), _tmp.substr(_tmp.find("_")+1)));
    cout<<"Adding mass point: mgluino = "<<mass_pts.back().first<<" mlsp = "<<mass_pts.back().second<<endl;
  } 

  // --------------------------------------
  //            Processes
  //---------------------------------------
  TString baseline("st>500 && met>200 && mj14>250 && njets>=6 && nbdm>=1 && nleps==1 && nveto==0");
  NamedFunc filters = Functions::hem_veto && Functions::pass_run2;
  NamedFunc nom_wgt = Functions::wgt_run2 * Functions::eff_trig_run2; 

  set<int> years;
  if (year==0) years = {2016, 2017, 2018};
  else years = {year};

  map<int, string> foldermc, folderdata, foldersig;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_abcd/";
  foldersig[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/T1tttt/skim_sys_abcd/";
  folderdata[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/data/merged_database_standard/";

  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_abcd/";
  foldersig[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/T1tttt/skim_sys_abcd/";
  folderdata[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/data/merged_database_stdnj5/";

  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_abcd/";
  foldersig[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/T1tttt/skim_sys_abcd/";
  folderdata[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/data/merged_database_standard/";

  // Filling all other processes
  vector<string> vnames_other = {
    "_WJetsToLNu_HT", "_ST_", "_TTW","_TTZ", 
    "_DYJetsToLL_M-50_HT","_ZJet","_ttH",
     "_TTGJets","_TTTT","_WH_HToBB","_ZH_HToBB","_WWTo","_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"
    };

  set<string> tt1l_files, tt2l_files, other_files, data_files;
  for (auto &yr: years) {
    tt1l_files.insert(foldermc[yr]+"*_TTJets*SingleLept*.root");
    tt2l_files.insert(foldermc[yr]+"*_TTJets*DiLept*.root");
    data_files.insert(folderdata[yr]+"*root");
    for(auto name : vnames_other)
      other_files.insert(foldermc[yr] + "*" + name + "*.root");
      
  }

  auto proc_tt1l = Process::MakeShared<Baby_full>("tt_1l",    Process::Type::background, 1,  
    tt1l_files, filters && "stitch_met");  
  auto proc_tt2l = Process::MakeShared<Baby_full>("tt_2l",    Process::Type::background, 1,  
    tt2l_files, filters && "stitch_met");  
  auto proc_other = Process::MakeShared<Baby_full>("other",    Process::Type::background, 1,  
    other_files, filters && "stitch_met");

  vector<shared_ptr<Process> > bkg_procs = {proc_tt1l, proc_tt2l, proc_other};

  vector<shared_ptr<Process> > data_procs;
  data_procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    data_files, filters && Functions::trig_run2));

  vector<shared_ptr<Process> > sig_procs;
  for (auto &imass: mass_pts) {
    set<string> sig_files;
    for (auto &yr: years) 
      sig_files.insert(foldersig[yr]+"*mGluino-"+imass.first+"_mLSP-"+imass.second+"_*.root");
    
    sig_procs.push_back(Process::MakeShared<Baby_full>("T1tttt", Process::Type::signal, kBlack,
      sig_files, filters));
  }

  // --------------------------------------
  //            Binning
  //---------------------------------------

  string c_mt_r1 = "mt<=140";

  // N.B.: currently, changing naming convention would break writing the closure systematics derived from CRs
  vector<string> vl_met, vc_met, vc_mj_r1;
  vl_met.push_back("lmet"); vc_met.push_back("met>200 && met<=350"); vc_mj_r1.push_back("mj14<=400");
  vl_met.push_back("mmet"); vc_met.push_back("met>350 && met<=500"); vc_mj_r1.push_back("mj14<=450");
  vl_met.push_back("hmet"); vc_met.push_back("met>500");             vc_mj_r1.push_back("mj14<=500");
  unsigned nbins_met(vl_met.size());

  vector<string> vl_nb, vc_nb;
  vl_nb.push_back("lnb");     vc_nb.push_back("nbdm==1");
  vl_nb.push_back("mnb");     vc_nb.push_back("nbdm==2");
  vl_nb.push_back("hnb");   vc_nb.push_back("nbdm>=3");
  unsigned nbins_nb(vl_nb.size());

  vector<string> vl_nj; vector<vector<string>> vc_nj;
  vl_nj.push_back("lnj");  vc_nj.push_back({"njets==7", "njets==7","njets>=6 && njets<=7"});
  vl_nj.push_back("hnj");   vc_nj.push_back({"njets>=8", "njets>=8", "njets>=8"});
  unsigned nbins_nj(vl_nj.size());

  vector<string> vl_mj, vc_mj;
  vl_mj.push_back("lmj");  vc_mj.push_back("mj14<=XXX");
  vl_mj.push_back("hmj");  vc_mj.push_back("mj14>XXX");
  unsigned nbins_mj(vl_mj.size());
  vector<string> v_metdep_midmj = {"500","650","800"};
  if (v_metdep_midmj.size()<vc_met.size()) {
    cout<<"ERROR: Intermediate MJ thresholds not specified for each MET bin"<<endl;
    exit(1);
  }

  vector<bindef> vbins; 
  for (unsigned imet(0); imet<nbins_met; imet++){
    string _label(""), _cut("");
    if (!Contains(xoption,"noabcd")) {
      _label = "r1_"+vl_met[imet];
      _cut = vc_met[imet]+"&&"+c_mt_r1+"&&"+vc_mj_r1[imet]+"&&(";
      for (unsigned inb(0); inb<nbins_nb; inb++){
        if (inb!=0) _cut +="||";
        for (unsigned inj(0); inj<nbins_nj; inj++){
          if (inj!=0) _cut +="||";
          _cut +="("+vc_nb[inb]+"&&"+vc_nj[inj][imet]+")";
        }
      }
      _cut+=")";
      vbins.push_back(bindef(_label, _cut));
      if (debug) cout<<_cut<<endl;

      for (unsigned inb(0); inb<nbins_nb; inb++){
        for (unsigned inj(0); inj<nbins_nj; inj++){
          for (unsigned imj(0); imj<nbins_mj; imj++){
            _label = "r2_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
            _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
            _cut = vc_met[imet]+"&&"
                   +c_mt_r1+"&&"+CopyReplaceAll(vc_mj_r1[imet],"<=",">")+"&&"
                   +vc_nb[inb]+"&&"+vc_nj[inj][imet];
            _cut += "&&"+CopyReplaceAll(vc_mj[imj],"XXX",v_metdep_midmj[imet]);
            vbins.push_back(bindef(_label, _cut));
            if (debug) cout<<_cut<<endl;
          }
        }
      }

      _label = "r3_"+vl_met[imet];
      _cut = vc_met[imet]+"&&"+CopyReplaceAll(c_mt_r1,"<=",">")+"&&"+vc_mj_r1[imet]+"&&(";
      for (unsigned inb(0); inb<nbins_nb; inb++){
        if (inb!=0) _cut +="||";
        for (unsigned inj(0); inj<nbins_nj; inj++){
          if (inj!=0) _cut +="||";
          _cut +="("+vc_nb[inb]+"&&"+vc_nj[inj][imet]+")";
        }
      }
      _cut+=")";
      vbins.push_back(bindef(_label, _cut));
      if (debug) cout<<_cut<<endl;
    }

    for (unsigned inb(0); inb<nbins_nb; inb++){
      for (unsigned inj(0); inj<nbins_nj; inj++){
        for (unsigned imj(0); imj<nbins_mj; imj++){
          _label = "r4_"+vl_met[imet]+'_'+vl_nb[inb]+'_'+vl_nj[inj];
          _label += '_'+CopyReplaceAll(vl_mj[imj],"XXX",v_metdep_midmj[imet]);
          _cut = vc_met[imet]+"&&"
                 +CopyReplaceAll(c_mt_r1,"<=",">")+"&&"+CopyReplaceAll(vc_mj_r1[imet],"<=",">")+"&&"
                 +vc_nb[inb]+"&&"+vc_nj[inj][imet];
          _cut += "&&"+CopyReplaceAll(vc_mj[imj],"XXX",v_metdep_midmj[imet]);
          vbins.push_back(bindef(_label, _cut));
          if (debug) cout<<_cut<<endl;
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
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_bctag_deep["+to_string(i)+"]/w_btag_deep");
    v_sys.push_back(sysdef("B-tag efficiency FS", "fs_bctag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_bctag_deep["+to_string(i)+"]/w_btag_deep");
    v_sys.push_back(sysdef("Mistag efficiency", "udsgtag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_udsgtag_deep["+to_string(i)+"]/w_btag_deep");
    v_sys.push_back(sysdef("Mistag efficiency FS", "fs_udsgtag",kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_udsgtag_deep["+to_string(i)+"]/w_btag_deep");
    v_sys.push_back(sysdef("Jet energy corrections", "jec", kCorr));
    v_sys.back().shift_index = 1; // JEC Up index in sys_met, etc.
    // v_sys.push_back(sysdef("QCD scales", "murf",kWeight));
    // for (size_t i = 0; i<2; ++i) {
    //   v_sys.back().v_wgts.push_back("sys_mur["+to_string(i)+"]");
    //   v_sys.back().v_wgts.push_back("sys_muf["+to_string(i)+"]");
    //   v_sys.back().v_wgts.push_back("sys_murf["+to_string(i)+"]");
    // }
    v_sys.push_back(sysdef("ISR", "isr", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_isr["+to_string(i)+"]/w_isr");
    v_sys.push_back(sysdef("Jet ID FS", "jetid", kConst));
    v_sys.back().v_wgts.push_back("0.01");
    v_sys.push_back(sysdef("Pileup", "pu", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_pu["+to_string(i)+"]/w_pu");
    v_sys.push_back(sysdef("Luminosity", "lumi", kConst));
    v_sys.back().v_wgts.push_back("0.025");
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
          cuts.emplace_back(TableRow("", baseline+"&&"+bin.cut,0,0,nom_wgt * wgt));
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
  if (debug) cout<<"calculating kappas..."<<endl;
  vector<float> vkappas, vkappas_unc;
  if (!Contains(xoption,"noabcd")) {
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
  if (debug) cout<<"Getting signal yields..."<<endl;
  vector<vector<GammaParams>> sig_params;
  yield_table = static_cast<Table*>(pm.Figures()[1].get());
  for (auto &isig: sig_procs) {
    sig_params.push_back(yield_table->Yield(isig.get(), lumi));
  }

  //------------------------------------------
  //       Get DATA yields
  //------------------------------------------
  if (debug) cout<<"Getting data yields..."<<endl;
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
    if (debug) cout<<"Going over sig_params..."<<endl;
    vector<GammaParams> nom_met_avg;
    for (auto &sys: v_sys) {
      if (sys.sys_type == kMetSwap) {
        for (size_t ibin = 0; ibin<nbins; ++ibin) {
          GammaParams tmp_gps;
          tmp_gps.SetYieldAndUncertainty(0.5*(sig_params[isig][sys.ind + ibin].Yield()+sig_params[isig][ibin].Yield()),
                        max(sig_params[isig][sys.ind + ibin].Uncertainty(), sig_params[isig][ibin].Uncertainty()));
          nom_met_avg.push_back(tmp_gps);
        }
      }
    }
    //   Writing datacard header
    //---------------------------------------
    if (debug) cout<<"Writing data cards..."<<endl;
    TString outpath = outfolder+"/datacard_SMS-"+model;
    outpath += "_mGluino-"+mass_pts[isig].first+"_mLSP-"+mass_pts[isig].second+"_";
    outpath += year;
    outpath += "_"+xoption;
    if (!unblind) outpath +="_blind";
    outpath +=".txt";
    if (!do_syst)  outpath.ReplaceAll(".txt*","_nosys.txt");
    cout<<"open "<<outpath<<endl;
    unsigned wname(23), wdist(7), wbin(23);
    for (size_t ibin(0); ibin<nbins; ibin++) 
      if(vbins[ibin].tag.length() > wbin) wbin = vbins[ibin].tag.length();
    wbin+=1;
    unsigned digit = 2;
    if (unblind) digit = 0;

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
    for (size_t ibin(0); ibin<nbins; ibin++) 
      fcard<<setw(wbin)<<Form("%.2f",nom_met_avg[ibin].Yield())<<setw(wbin)<<"1";
    fcard<<endl<<endl;
    cout<<"Wrote headers"<<endl;

    //--------- Signal statistical uncertainties ----------------------------
    ofstream fsys(outpath.ReplaceAll("datacard_","sys_")); //start writing also in a sys file
    fsys<<"\nSYSTEMATIC MC_stat\n  PROCESSES signal\n";
    for (size_t ibin(0); ibin<nbins; ibin++) {
      TString sig_stat = Form("%.2f",1+nom_met_avg[ibin].Uncertainty()/nom_met_avg[ibin].Yield());
      fsys<<"    " <<left<<setw(25)<<vbins[ibin].tag <<" "<<right<<setw(10)<<nom_met_avg[ibin].Yield()
          <<" "<<right<<setw(10)<<sig_stat<<endl;
      fcard<<setw(wname)<<"stat_"+vbins[ibin].tag<<setw(wdist)<<"lnN";
      for (size_t jbin(0); jbin<nbins; jbin++) {
        if (ibin==jbin) fcard<<setw(wbin)<<sig_stat<<setw(wbin)<<"-";
        else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
      }
      fcard<<endl;
    }
    cout<<"Wrote signal stat. uncertainties"<<endl;
    

    // ------------ Closure uncertainties
    if (do_syst){
      vector<vector<float>> cr_unc_nj = {{1.09, 1.09}, {1.09, 1.09}}; 
      vector<vector<float>> cr_unc_nb = {{1.10, 1.20, 1.25}, {1.10, 1.20, 1.25}}; 
      vector<vector<float>> cr_unc_met = {{1.15, 1.21}, {1.18, 1.30}}; //with 9% unc on Njets already folded in 
      for (unsigned imj(0); imj<vl_mj.size(); imj++){
        fcard<<endl<<setw(wname)<<"cr_unc_nj_"+vl_mj[imj]<<setw(wdist)<<"lnN";
        for (size_t ibin(0); ibin<nbins; ibin++) {
          if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,vl_mj[imj])) {
            if(Contains(vbins[ibin].tag,"lnj")) 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_nj[imj][0],2);
            else 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_nj[imj][1],2);
          } else {
            fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          }
        }
      }
      for (unsigned imj(0); imj<vl_mj.size(); imj++){
        fcard<<endl<<setw(wname)<<"cr_unc_nb_"+vl_mj[imj]<<setw(wdist)<<"lnN";
        for (size_t ibin(0); ibin<nbins; ibin++) {
          if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,vl_mj[imj])) {
            if (Contains(vbins[ibin].tag,"lnb")) 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_nb[imj][0],2);
            else if (Contains(vbins[ibin].tag,"mnb")) 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_nb[imj][1],2);
            else
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_nb[imj][2],2);
          } else {
            fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          }
        }
      }
      for (unsigned imj(0); imj<vl_mj.size(); imj++){
        fcard<<endl<<setw(wname)<<"cr_unc_met_"+vl_mj[imj]<<setw(wdist)<<"lnN";
        for (size_t ibin(0); ibin<nbins; ibin++) {
          if(Contains(vbins[ibin].tag,"r2_") && Contains(vbins[ibin].tag,vl_mj[imj])) {
            if (Contains(vbins[ibin].tag,"lmet")) 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
            else if (Contains(vbins[ibin].tag,"mmet")) 
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_met[imj][0],2);
            else
              fcard<<setw(wbin)<<"-"<<setw(wbin)<<RoundNumber(cr_unc_met[imj][1],2);
          } else {
            fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          }
        }
      }
      fcard<<endl;
      cout<<"Wrote CR-based closure uncertainties"<<endl;
    }

    //calculate uncertainties and write results to three files
    cout<<"Writing to "<<outpath<<endl;
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
          } 
          // convert to ra4_stats input and write to file

          // write systematics to datacard
          double unc = fabs(up)>fabs(dn) ? fabs(up) : fabs(dn);
          unc = (up>0 ? 1:-1)*unc + 1.;

          if (std::isnan(unc) || std::isinf(unc)) {
            cout <<" Found bad unc. set to 0 -> "<<left<<setw(10)<<sys.tag <<left<<setw(10)<<vbins[ibin].tag 
                 <<" "<<right<<setprecision(0)<<setw(25)<<sig_params[isig][ibin].NEffective() 
                 <<" "<<setprecision(5)<<setw(15)<<sig_params[isig][ibin].Yield() 
                 <<" "<<setprecision(10)<<setw(15)<<sig_params[isig][ibin].Weight()<<endl;  
            unc = 2; // put a 100% uncertainty, this is probably very poor stats bin for signal, i.e. 0
          } 
          fsys<<"    " <<left<<setw(25)<<vbins[ibin].tag <<" "<<right<<setw(10)<<nom_met_avg[ibin].Yield()
              <<" "<<right<<setw(10)<<Form("%.2f",unc-1) <<endl;
          fcard<<setw(wbin)<<Form("%.2f",unc)<<setw(wbin)<<"-";
        } // loop over bins
        fcard<<endl;
      } // loop over systematics
    }
    fsys.close();

    //   Writing datacard param section
    //---------------------------------------
    if (!Contains(xoption,"noabcd")) {
      unsigned iyield(0);
      wbin +=3; // to allow the label to fit even with "rp_" in front
      // order of for loops must match the one used to define cut vector, i.e. met, nb, nj!
      if(debug) cout <<endl<<left<<setw(wbin)<<"bin"
                     <<right<<setw(15)<<"Sig."
                     <<right<<setw(15)<<"MET avg"
                     <<right<<setw(15)<<"MET unc"
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
                       <<right<<setw(15)<<RoundNumber(100*nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),1)
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
                             <<right<<setw(15)<<RoundNumber(100*nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),1)
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
                       <<right<<setw(15)<<RoundNumber(100*nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),1)
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
                             <<right<<setw(15)<<RoundNumber(100*nom_met_avg[iyield].Uncertainty()/nom_met_avg[iyield].Yield(),1)
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
  ibin.ReplaceAll("nbdm", "sys_nbdm["+to_string(shift_index)+"]");
  return ibin;
}

TString nom2genmet(TString ibin){
  ibin.ReplaceAll("met", "met_tru");
  ibin.ReplaceAll("met_tru/met_tru_calo", "met/met_calo");
  return ibin;

}

void GetOptions(int argc, char *argv[]){
  string blah;
  while(true){
    static struct option long_options[] = {
      {"model", required_argument, 0, 'm'},
      {"xoption", required_argument, 0, 'x'},
      {"mass_pts", required_argument, 0, 'p'},
      {"no_syst", no_argument, 0, 0},
      {"unblind", no_argument, 0, 'u'},
      {"year", no_argument, 0, 'y'},
      {"debug", no_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:x:o:p:uy:d", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'm': model = optarg; break;
    case 'x': xoption = optarg; break;
    case 'o': outfolder = optarg; break;
    case 'p': mass_pts_str = optarg; break;
    case 'u': unblind = true; break;
    case 'y': year = atoi(optarg); break;
    case 'd': debug = true; break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "no_syst"){
        do_syst = false;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default: printf("Bad option! getopt_long returned character code 0%o\n", opt); break;
    }
  }
}