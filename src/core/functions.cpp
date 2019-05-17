#include "core/functions.hpp"

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Functions{

  float wnpv2017(const Baby &b){
    if (b.SampleType()<0) return 1;
    if (b.event()%1000<675) return 1;
    int _npv = b.npv();
    if(_npv <  2) return  0.621;// +- 0.058  
    else if(_npv <  4) return  0.706;// +- 0.032  
    else if(_npv <  6) return  0.803;// +- 0.024  
    else if(_npv <  8) return  0.896;// +- 0.017  
    else if(_npv < 10) return  0.847;// +- 0.010  
    else if(_npv < 12) return  0.686;// +- 0.006  
    else if(_npv < 14) return  0.569;// +- 0.004  
    else if(_npv < 16) return  0.518;// +- 0.003  
    else if(_npv < 18) return  0.494;// +- 0.003  
    else if(_npv < 20) return  0.499;// +- 0.002  
    else if(_npv < 22) return  0.525;// +- 0.002  
    else if(_npv < 24) return  0.561;// +- 0.002  
    else if(_npv < 26) return  0.623;// +- 0.002  
    else if(_npv < 28) return  0.711;// +- 0.003  
    else if(_npv < 30) return  0.826;// +- 0.003  
    else if(_npv < 32) return  0.966;// +- 0.003  
    else if(_npv < 34) return  1.166;// +- 0.004  
    else if(_npv < 36) return  1.373;// +- 0.005  
    else if(_npv < 38) return  1.587;// +- 0.005  
    else if(_npv < 40) return  1.766;// +- 0.006  
    else if(_npv < 42) return  1.962;// +- 0.007  
    else if(_npv < 44) return  2.192;// +- 0.009  
    else if(_npv < 46) return  2.352;// +- 0.010  
    else if(_npv < 48) return  2.543;// +- 0.012  
    else if(_npv < 50) return  2.758;// +- 0.015  
    else if(_npv < 52) return  3.084;// +- 0.019  
    else if(_npv < 54) return  3.423;// +- 0.024  
    else if(_npv < 56) return  3.698;// +- 0.030  
    else if(_npv < 58) return  4.008;// +- 0.037  
    else if(_npv < 60) return  4.309;// +- 0.047  
    else if(_npv < 62) return  5.114;// +- 0.064  
    else if(_npv < 64) return  5.274;// +- 0.077  
    else if(_npv < 66) return  6.153;// +- 0.105  
    else if(_npv < 68) return  6.396;// +- 0.130  
    else if(_npv < 70) return  7.101;// +- 0.167  
    else if(_npv < 72) return  7.078;// +- 0.195  
    else if(_npv < 74) return  8.941;// +- 0.291  
    else if(_npv < 76) return  9.725;// +- 0.357  
    else if(_npv < 78) return 10.044;// +- 0.439  
    else if(_npv < 80) return 11.836;// +- 0.599  
    else if(_npv < 82) return 12.909;// +- 0.755  
    else if(_npv < 84) return 11.381;// +- 0.793  
    else if(_npv < 86) return 11.387;// +- 0.969  
    else if(_npv < 88) return 12.642;// +- 1.252  
    else if(_npv < 90) return 15.268;// +- 1.592  
    else if(_npv < 92) return 15.476;// +- 1.877  
    else if(_npv < 94) return 11.923;// +- 1.960  
    else if(_npv < 96) return 12.651;// +- 2.080  
    else if(_npv < 98) return  7.552;// +- 1.888  
    else return 11.175;// +- 1.443  
  }

  // based on r136 of https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  const NamedFunc pass_run2("pass_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_ = (b.st()<10000) && b.pass_ra2_badmu() && (b.met()/b.met_calo()<5);
    if (b.SampleType()<0) { // Data
      pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_cschalo_tight() && b.pass_hbhe() && 
              b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu() && b.pass_eebadsc();
      if (b.SampleType()!=-2016) pass_ = pass_ && b.pass_badcalib();
    } else {
      if (b.type()>=100e3) { // FastSim
        pass_ = pass_ && b.pass_fsjets() && b.pass_goodv() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        if (b.SampleType()!=2016) pass_ = pass_ && b.pass_badcalib();
      } else { //FullSim
        pass_ = pass_ && b.pass_jets() && b.pass_goodv() && b.pass_cschalo_tight() && b.pass_hbhe() && 
                b.pass_hbheiso() && b.pass_ecaldeadcell() && b.pass_badpfmu();
        if (b.SampleType()!=2016) pass_ = pass_ && b.pass_badcalib();
      }
    }
    if (pass_) return 1.;
    else return 0.;
  });

  const NamedFunc wgt_run2("wgt_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;

    double wgt = b.weight();
    if (b.SampleType()==2016){
      return wgt*b.w_prefire()*35.9;
    } else if (b.SampleType()==2017){
      return wgt*b.w_prefire()*wnpv2017(b)*41.5;
    } else {
      return wgt*59.6;
    }
  });

  const NamedFunc mht_ratio("mht_ratio", [](const Baby &b) -> NamedFunc::ScalarType{
    double mht_px(0.), mht_py(0.);
    for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
      if (b.jets_pt()->at(ijet) <= 30. && fabs(b.jets_eta()->at(ijet))>2.4) continue;
      mht_px -= b.jets_pt()->at(ijet)*cos(b.jets_phi()->at(ijet));
      mht_py -= b.jets_pt()->at(ijet)*sin(b.jets_phi()->at(ijet));
    }
    return b.mht()/hypot(mht_px, mht_py);
  });

  const NamedFunc fake_met("fake_met", [](const Baby &b) -> NamedFunc::ScalarType{
    return 10*wnpv2017(b);
  });

  const NamedFunc trig_run2("trig_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()>0) return 1.;
    
    int _pass(0);
    if (b.SampleType()==-2016){
      _pass = b.trig()->at(3) || b.trig()->at(4) || b.trig()->at(34) || b.trig()->at(7) || 
              b.trig()->at(8) || b.trig()->at(36) || b.trig()->at(14) || b.trig()->at(15) || 
              b.trig()->at(30) || b.trig()->at(31) || b.trig()->at(19) || b.trig()->at(55) || 
              b.trig()->at(20) || b.trig()->at(21) || b.trig()->at(40) || b.trig()->at(41);

    } else { // 2017 or 2018
      _pass = b.trig()->at(2) || b.trig()->at(3) || b.trig()->at(22) || b.trig()->at(6) || 
              b.trig()->at(7) || b.trig()->at(30) || b.trig()->at(9) || b.trig()->at(10) || 
              b.trig()->at(15) || b.trig()->at(13) || b.trig()->at(19) || b.trig()->at(20) || 
              b.trig()->at(21) || b.trig()->at(26) || b.trig()->at(23) || b.trig()->at(24);
    }
    return _pass;
  });

  const NamedFunc eff_trig_run2("eff_trig_run2", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.;
    double _met = b.met();
    if (b.SampleType()==2016) {
      if (b.nmus()>=1) {
        if (_met<20) return 0.948;
        else if (_met<40) return 0.949;
        else if (_met<60) return 0.951;
        else if (_met<80) return 0.955;
        else if (_met<100) return 0.960;
        else if (_met<120) return 0.970;
        else if (_met<140) return 0.975;
        else if (_met<160) return 0.984;
        else if (_met<180) return 0.993;
        else return 1.;
      } else if (b.nels()>=1) {
        if (_met<20) return 0.84109;
        else if (_met<40) return 0.833;
        else if (_met<60) return 0.850;
        else if (_met<80) return 0.862;
        else if (_met<100) return 0.890;
        else if (_met<120) return 0.901;
        else if (_met<140) return 0.937;
        else if (_met<160) return 0.947;
        else if (_met<180) return 0.967;
        else if (_met<200) return 0.976;
        else if (_met<220) return 0.993;
        else return 1.;
      } else {
        return 0.;
      }
    } else if (b.SampleType()==2017) {
      if (b.nmus()>=1) {
        if (_met<20) return 0.933;
        else if (_met<40) return 0.923;
        else if (_met<60) return 0.924;
        else if (_met<80) return 0.934;
        else if (_met<100) return 0.944;
        else if (_met<120) return 0.957;
        else if (_met<140) return 0.969;
        else if (_met<160) return 0.978;
        else if (_met<180) return 0.986;
        else if (_met<200) return 0.993;
        else return 1.;
      } else if (b.nels()>=1) {
        if (_met<20) return 0.791;
        else if (_met<40) return 0.779;
        else if (_met<60) return 0.791;
        else if (_met<80) return 0.794;
        else if (_met<100) return 0.828;
        else if (_met<120) return 0.854;
        else if (_met<140) return 0.886;
        else if (_met<160) return 0.917;
        else if (_met<180) return 0.945;
        else if (_met<200) return 0.968;
        else if (_met<220) return 0.978;
        else if (_met<240) return 0.990;
        else if (_met<260) return 0.991;
        else return 1.;
      } else {
        return 0.;
      }
    } else { //2018
      if (b.nmus()>=1) {
        if (_met<20) return 0.956;
        else if (_met<40) return 0.964;
        else if (_met<60) return 0.964;
        else if (_met<80) return 0.969;
        else if (_met<100) return 0.969;
        else if (_met<120) return 0.976;
        else if (_met<140) return 0.979;
        else if (_met<160) return 0.995;
        else if (_met<180) return 0.993;
        else return 1.;
      } else if (b.nels()>=1) {
        if (_met<20) return 0.827;
        else if (_met<40) return 0.814;
        else if (_met<60) return 0.820;
        else if (_met<80) return 0.817;
        else if (_met<100) return 0.834;
        else if (_met<120) return 0.851;
        else if (_met<140) return 0.872;
        else if (_met<160) return 0.909;
        else if (_met<180) return 0.935;
        else if (_met<200) return 0.965;
        else if (_met<220) return 0.976;
        else if (_met<240) return 0.991;
        else if (_met<260) return 0.993;
        else return 1.;
      } else {
        return 0.;
      }
    }
  });

  const NamedFunc hem_veto("hem_veto",[](const Baby &b) -> NamedFunc::ScalarType{
    if(abs(b.SampleType()) == 2018) {
      if (b.SampleType()<0) {
        if (b.run() >= 319077) { 
          if(b.nels() > 0) {
            for(size_t i = 0; i < b.leps_pt()->size(); i++) {
              if (abs(b.leps_id()->at(i))==13) continue;
              if(b.leps_eta()->at(i) < -1.5 && (b.leps_phi()->at(i) > -1.6 && b.leps_phi()->at(i) < -0.8)) 
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
            for(size_t i = 0; i < b.leps_pt()->size(); i++) {
              if (abs(b.leps_id()->at(i))==13) continue;
              if(b.leps_eta()->at(i) < -1.5 && (b.leps_phi()->at(i) > -1.6 && b.leps_phi()->at(i) < -0.8)) 
                return static_cast<float>(0);
            }
          }
          for(size_t i = 0; i < b.jets_pt()->size(); i++) {
            if(Functions::IsGoodJet(b,i) && b.jets_eta()->at(i) < -1.5 && (b.jets_phi()->at(i) > -1.6 && b.jets_phi()->at(i) < -0.8)) 
              return static_cast<float>(0);
          }
        }
      }
    }
    return static_cast<float>(1);
  });

  std::pair<double, double> calc_adj_met(const Baby &b){
    if (b.SampleType() > 0) {
      int rj;
      string run;
      vector<double> corrs, fracs;
      if(b.SampleType() == 2016) {
        rj = b.event()%3590;
        if(rj <= 505)       { corrs = {0.90, 0.95, 1.15}; fracs = {0.672, 0.994}; } // 2016B
        else if(rj <=  745) { corrs = {0.96, 1.00, 1.17}; fracs = {0.621, 0.993}; } // 2016C
        else if(rj <= 1171) { corrs = {0.96, 1.00, 1.17}; fracs = {0.655, 0.994}; } // 2016D
        else if(rj <= 1576) { corrs = {1.01, 1.06, 1.19}; fracs = {0.646, 0.994}; } // 2016E
        else if(rj <= 1887) { corrs = {1.04, 1.07, 1.15}; fracs = {0.624, 0.994}; } // 2016F
        else if(rj <= 2641) { corrs = {1.00, 1.02, 1.20}; fracs = {0.672, 0.994}; } // 2016G
        else                { corrs = {1.03, 1.07, 1.33}; fracs = {0.642, 0.990}; } // 2016H
      }
      else if(b.SampleType() == 2017) {
        rj = b.event()%4150;
        if(rj <= 482)       { corrs = {1.01, 0.98, 1.09}; fracs = {0.426, 0.983}; } // 2017B
        else if(rj <= 1448) { corrs = {1.00, 0.95, 1.15}; fracs = {0.485, 0.988}; } // 2017C
        else if(rj <= 1873) { corrs = {1.00, 0.96, 1.14}; fracs = {0.454, 0.987}; } // 2017D
        else if(rj <= 2801) { corrs = {1.09, 1.07, 1.16}; fracs = {0.318, 0.978}; } // 2017E
        else                { corrs = {1.15, 1.15, 1.19}; fracs = {0.257, 0.970}; } // 2017F
      }
      else {
        rj = b.event()%6000;
        if(rj <= 1400)      { corrs = {1.12, 1.12, 1.26}; fracs = {0.331, 0.983}; } // 2018A
        else if(rj <= 2110) { corrs = {1.10, 1.13, 1.18}; fracs = {0.354, 0.988}; } // 2018B
        else if(rj <= 2804) { corrs = {1.15, 1.15, 1.17}; fracs = {0.298, 0.985}; } // 2018C
        else                { corrs = {1.10, 1.13, 1.20}; fracs = {0.361, 0.987}; } // 2018D
      }
      TRandom3 rng(b.event());
      double ri(rng.Rndm());
      int i = (ri < fracs.at(0)) + (ri < fracs.at(1)) + (ri >= fracs.at(1));
      double metx(b.met()*cos(b.met_phi())), mety(b.met()*sin(b.met_phi()));
      double met_trux(b.met_tru()*cos(b.met_tru_phi())), met_truy(b.met_tru()*sin(b.met_tru_phi()));
      double n_metx = (metx-met_trux)*corrs.at(i) + met_trux;
      double n_mety = (mety-met_truy)*corrs.at(i) + met_truy;
      double n_metphi = atan2(n_mety,n_metx);
      return std::make_pair(hypot(n_metx,n_mety), n_metphi);
    }
    return std::make_pair(b.met(), b.met_phi());
  };

  const NamedFunc adj_met("adj_met", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType() > 0) {
      std::pair<double, double> met_val_phi = calc_adj_met(b);
      return met_val_phi.first;
    }
    return b.met();
  });

  const NamedFunc adj_mt("adj_mt", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType() > 0) {
      std::pair<double, double> met_val_phi = calc_adj_met(b);
      return sqrt(2*b.leps_pt()->at(0)*met_val_phi.first*(1-cos(met_val_phi.second-b.leps_phi()->at(0))));
    }
    return b.mt();
  });

  const NamedFunc ntrub("ntrub", [](const Baby &b){
      int ntrub_(0);
      for(size_t i = 0; i < b.jets_pt()->size(); i++) {
        if(IsGoodJet(b,i) && b.jets_hflavor()->at(i)==5)
          ntrub_++;
      }
      return ntrub_;
    });

  const NamedFunc n_mus_bad("n_mus_bad", [](const Baby &b) -> NamedFunc::ScalarType{
      int n=0;
      for(unsigned int i=0; i< b.mus_pt()->size(); i++){
	if(b.mus_bad()->at(i)) n++;
      }
      return n;
    });

  const NamedFunc n_mus_bad_dupl("n_mus_bad_dupl", [](const Baby &b) -> NamedFunc::ScalarType{
      int n=0;
      for(unsigned int i=0; i< b.mus_pt()->size(); i++){
	if(b.mus_bad_dupl()->at(i)) n++;
      }
      return n;
    });

  const NamedFunc n_mus_bad_trkmu("n_mus_bad_trkmu", [](const Baby &b) -> NamedFunc::ScalarType{
      int n=0;
      for(unsigned int i=0; i< b.mus_pt()->size(); i++){
	if(b.mus_bad_trkmu()->at(i)) n++;
      }
      return n;
    });

  const NamedFunc n_isr_match("n_isr_match", NISRMatch);

  const NamedFunc njets_weights_ttisr("njets_weights_ttisr", [](const Baby &b){
      return NJetsWeights_ttISR(b, false);
    });

  const NamedFunc njets_weights_visr("njets_weights_visr", NJetsWeights_vISR);

  const NamedFunc min_dphi_lep_met("min_dphi_lep_met", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
      double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
      if(phi1 != -999 && phi2 != -999){
        return std::min(dphi1,dphi2);
      }else if(phi1 != -999 && phi2 == -999){
        return dphi1;
      }else if(phi1 == -999 && phi2 != -999){
        return dphi2;
      }
      return -1;
    });

  const NamedFunc max_dphi_lep_met("max_dphi_lep_met", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
      double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
      if(phi1 != -999 && phi2 != -999){
        return std::max(dphi1,dphi2);
      }else if(phi1 != -999 && phi2 == -999){
        return dphi1;
      }else if(phi1 == -999 && phi2 != -999){
        return dphi2;
      }
      return -1;
    });

  const NamedFunc min_dphi_lep_jet("min_dphi_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double minphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
        double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
        double thisdphi = -1;
        if(phi1 != -999 && phi2 != -999){
          thisdphi = std::min(dphi1, dphi2);
        }else if(phi1 != -999 && phi2 == -999){
          thisdphi = dphi1;
        }else if(phi1 == -999 && phi2 != -999){
          thisdphi = dphi2;
        }
        if(minphi < 0. || thisdphi < minphi){
          minphi = thisdphi;
        }
       }
      return minphi;
    });

  const NamedFunc nbm_moriond("nbm_moriond", [](const Baby &b) ->NamedFunc::ScalarType{
      int nbm = 0;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
	if(b.jets_csv()->at(ijet) > 0.8484) nbm++;
	//if(b.jets_csv()->at(ijet) > 0.800) nbm++;
      } // Loop over jets
      return nbm;
    });

  const NamedFunc ntop_loose_nom("ntop_loose_nom", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
      	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.1883) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc ntop_med_nom("ntop_med_nom", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.8511) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc ntop_tight_nom("ntop_tight_nom", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.9377) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc ntop_loose_decor("ntop_loose_decor", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
      	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.04738 && b.ak8jets_m()->at(ijet)>105 && b.ak8jets_m()->at(ijet)<210) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc ntop_med_decor("ntop_med_decor", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.4585 && b.ak8jets_m()->at(ijet)>105 && b.ak8jets_m()->at(ijet)<210) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc ntop_tight_decor("ntop_tight_decor", [](const Baby &b) ->NamedFunc::ScalarType{
      int ntop = 0;
      for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ++ijet){
	if(!IsGoodak8Jet(b,ijet)) continue;
	if(b.ak8jets_decor_bin_top()->at(ijet) > 0.6556 && b.ak8jets_m()->at(ijet)>105 && b.ak8jets_m()->at(ijet)<210) ntop++;
      } // Loop over all ak8 jets
      return ntop;
    });

  const NamedFunc max_dphi_lep_jet("max_dphi_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double maxphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
        double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
        double thisdphi = -1;
        if(phi1 != -999 && phi2 != -999){
          thisdphi = std::max(dphi1, dphi2);
        }else if(phi1 != -999 && phi2 == -999){
          thisdphi = dphi1;
        }else if(phi1 == -999 && phi2 != -999){
          thisdphi = dphi2;
        }
        if(maxphi < 0. || thisdphi > maxphi){
          maxphi = thisdphi;
        }
      }
      return maxphi;
    });

  const NamedFunc min_dphi_met_jet("min_dphi_met_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double minphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
        if(minphi < 0. || thisdphi < minphi){
          minphi = thisdphi;
        }
      }
      return minphi;
    });

  const NamedFunc max_dphi_met_jet("max_dphi_met_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double maxphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
        if(maxphi < 0. || thisdphi > maxphi){
          maxphi = thisdphi;
        }
      }
      return maxphi;
    });

  const NamedFunc min_dr_lep_jet("min_dr_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double minr = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
        double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
        double thisdr = -1;
        if(phi1 != -999 && phi2 != -999){
          thisdr = std::min(dr1, dr2);
        }else if(phi1 != -999 && phi2 == -999){
          thisdr = dr1;
        }else if(phi1 == -999 && phi2 != -999){
          thisdr = dr2;
        }
        if(minr < 0. || thisdr < minr){
          minr = thisdr;
        }
      }
      return minr;
    });

  const NamedFunc max_dr_lep_jet("max_dr_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
      double maxr = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
        double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
        double thisdr = -1;
        if(phi1 != -999 && phi2 != -999){
          thisdr = std::max(dr1, dr2);
        }else if(phi1 != -999 && phi2 == -999){
          thisdr = dr1;
        }else if(phi1 == -999 && phi2 != -999){
          thisdr = dr2;
        }
        if(maxr < 0. || thisdr > maxr){
          maxr = thisdr;
        }
      }
      return maxr;
    });

  const NamedFunc offshellw("offshellw",[](const Baby &b) -> NamedFunc::ScalarType{
      for (unsigned i(0); i<b.mc_pt()->size(); i++){
	if (abs(b.mc_id()->at(i))!=24) continue;
	if (b.mc_mass()->at(i) > 140.) {
	  return 1;
	}
      }
      return 0;
    });

  bool IsGoodJet(const Baby &b, size_t ijet){
    return ijet<b.jets_pt()->size()
      && b.jets_pt()->at(ijet) > 30.
      && fabs(b.jets_eta()->at(ijet))<2.4
      && !b.jets_islep()->at(ijet);
  }

  bool IsGoodak8Jet(const Baby &b, size_t ijet){
    return ijet<b.ak8jets_pt()->size()
      && b.ak8jets_pt()->at(ijet) > 300.
      && fabs(b.ak8jets_eta()->at(ijet))<2.4;
  }

  bool IsGoodElectron(const Baby &b, size_t iel){
    return iel<b.els_pt()->size()
      && b.els_pt()->at(iel)>20.
      && fabs(b.els_sceta()->at(iel))<2.5
      && b.els_sigid()->at(iel)
      && b.els_miniso()->at(iel) >= 0.
      && b.els_miniso()->at(iel) < 0.1;
  }

  bool IsGoodMuon(const Baby &b, size_t imu){
    return imu<b.mus_pt()->size()
      && b.mus_pt()->at(imu)>20.
      && fabs(b.mus_eta()->at(imu))<2.4
      && b.mus_sigid()->at(imu)
      && b.mus_miniso()->at(imu) >= 0.
      && b.mus_miniso()->at(imu) < 0.2;
  }

  bool IsGoodTrack(const Baby &b, size_t itk){
    if(itk >= b.tks_pt()->size()) return false;

    if(abs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 && b.tks_mt2()->at(itk)<60 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 ){
      return true;
    }else if(abs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
      return true;
    }else if(abs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
      return true;
    }else{
      return false;
    }
  }

  NamedFunc::ScalarType NJetsWeights_ttISR(const Baby &b, bool use_baby_nisr){
    if (b.ntrupv()<0) return 1.; // Do not reweight Data

    float wgt = b.weight()/b.eff_trig()/b.w_toppt();

    int nisrjets = use_baby_nisr ? b.nisr() : NISRMatch(b);

    // weights derived in TTJets and applied using the nisr calculation algorithm
    if      (nisrjets==0) return 1.099*wgt; // ;// +- 0.012
    else if (nisrjets==1) return 0.969*wgt; // ;// +- 0.014
    else if (nisrjets==2) return 0.870*wgt; // ;// +- 0.020
    else if (nisrjets==3) return 0.772*wgt; // ;// +- 0.031
    else if (nisrjets==4) return 0.712*wgt; // ;// +- 0.051
    else if (nisrjets==5) return 0.661*wgt; // ;// +- 0.088
    else if (nisrjets>=6) return 0.566*wgt; // ;// +- 0.133
    else return wgt;
  }

  NamedFunc::ScalarType NJetsWeights_vISR(const Baby &b){
    if (b.ntrupv()<0) return 1.; // Do not reweight Data
    
    float wgt = b.weight()/b.eff_trig()/b.w_toppt();
    if(b.SampleType()<30 && b.SampleType()>=60) return wgt;
    
    int nisrjets(b.njets());
    // weights derived in DY+jets
    if      (nisrjets==0) return 0.981*wgt; // ;// +- 0.001
    else if (nisrjets==1) return 1.071*wgt; // ;// +- 0.001
    else if (nisrjets==2) return 1.169*wgt; // ;// +- 0.003
    else if (nisrjets==3) return 1.157*wgt; // ;// +- 0.007
    else if (nisrjets==4) return 1.014*wgt; // ;// +- 0.013
    else if (nisrjets==5) return 0.920*wgt; // ;// +- 0.025
    else if (nisrjets==6) return 0.867*wgt; // ;// +- 0.048
    else if (nisrjets>=7) return 0.935*wgt; // ;// +- 0.088
    else return wgt;
  }

  int NISRMatch(const Baby &b){
    int Nisr=0;
    for (size_t ijet(0); ijet<b.jets_pt()->size(); ++ijet){
      if(!IsGoodJet(b, ijet)) continue;
      bool matched=false;
      for (size_t imc(0); imc<b.mc_pt()->size(); ++imc){
	if(b.mc_status()->at(imc)!=23 || abs(b.mc_id()->at(imc))>5) continue;
	if(!(abs(b.mc_mom()->at(imc))==6 || abs(b.mc_mom()->at(imc))==23 ||
	     abs(b.mc_mom()->at(imc))==24 || abs(b.mc_mom()->at(imc))==15)) continue; // In our ntuples where all taus come from W
	float dR = deltaR(b.jets_eta()->at(ijet), b.jets_phi()->at(ijet), b.mc_eta()->at(imc), b.mc_phi()->at(imc));
	if(dR<0.4){
	  matched = true;
	  break;
	}
      } // Loop over MC particles
      if(!matched) ++Nisr;
    } // Loop over jets

    return Nisr;
  }

  

  void DileptonAngles(const Baby &b,
                      NamedFunc::ScalarType &eta1, NamedFunc::ScalarType &phi1,
                      NamedFunc::ScalarType &eta2, NamedFunc::ScalarType &phi2){
    phi1 = -999.; eta1 = -999.;
    phi2 = -999.; eta2 = -999.;
    bool h1=false, h2=false;
    if(b.nels()==2 && b.nmus()==0){
      for(size_t iel = 0; iel < b.els_pt()->size() && !(h1&&h2); ++iel){
        if(IsGoodElectron(b, iel)){
          if(!h1){
            phi1 = b.els_phi()->at(iel);
            eta1 = b.els_sceta()->at(iel);
            h1 = true;
          }else if(!h2){
            phi2 = b.els_phi()->at(iel);
            eta2 = b.els_sceta()->at(iel);
            h2 = true;
          }
        }
      }
    }else if(b.nels()==1 && b.nmus()==1){
      for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
        if(IsGoodElectron(b, iel)){
          if(!h1){
            phi1 = b.els_phi()->at(iel);
            eta1 = b.els_sceta()->at(iel);
            h1 = true;
          }else if(!h2){
            phi2 = b.els_phi()->at(iel);
            eta2 = b.els_sceta()->at(iel);
            h2 = true;
          }
        }
      }
      for(size_t imu = 0; imu < b.mus_pt()->size() && h1 && !h2; ++imu){
        if(IsGoodMuon(b, imu)){
          if(!h1){
            phi1 = b.mus_phi()->at(imu);
            eta1 = b.mus_eta()->at(imu);
            h1 = true;
          }else if(!h2){
            phi2 = b.mus_phi()->at(imu);
            eta2 = b.mus_eta()->at(imu);
            h2 = true;
          }
        }
      }
    }else if(b.nels()==0 && b.nmus()==2){
      for(size_t imu = 0; imu < b.mus_pt()->size() && !(h1&&h2); ++imu){
        if(IsGoodMuon(b, imu)){
          if(!h1){
            phi1 = b.mus_phi()->at(imu);
            eta1 = b.mus_eta()->at(imu);
            h1 = true;
          }else if(!h2){
            phi2 = b.mus_phi()->at(imu);
            eta2 = b.mus_eta()->at(imu);
            h2 = true;
          }
        }
      }
    }else if(b.nels()==1 && b.nmus()==0 && b.nveto()==1){
      for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
        if(IsGoodElectron(b, iel)){
          if(!h1){
            phi1 = b.els_phi()->at(iel);
            eta1 = b.els_sceta()->at(iel);
            h1 = true;
          }else if(!h2){
            phi2 = b.els_phi()->at(iel);
            eta2 = b.els_sceta()->at(iel);
            h2 = true;
          }
        }
      }
      for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
        if(IsGoodTrack(b, itk)){
          if(!h1){
            phi1 = b.tks_phi()->at(itk);
            eta1 = b.tks_eta()->at(itk);
            h1 = true;
          }else if(!h2){
            phi2 = b.tks_phi()->at(itk);
            eta2 = b.tks_eta()->at(itk);
            h2 = true;
          }
        }
      }
    }else if(b.nels()==0 && b.nmus()==1 && b.nveto()==1){
      for(size_t imu = 0; imu < b.mus_pt()->size() && !h1; ++imu){
        if(IsGoodMuon(b, imu)){
          if(!h1){
            phi1 = b.mus_phi()->at(imu);
            eta1 = b.mus_eta()->at(imu);
            h1 = true;
          }else if(!h2){
            phi2 = b.mus_phi()->at(imu);
            eta2 = b.mus_eta()->at(imu);
            h2 = true;
          }
        }
      }
      for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
        if(IsGoodTrack(b, itk)){
          if(!h1){
            phi1 = b.tks_phi()->at(itk);
            eta1 = b.tks_eta()->at(itk);
            h1 = true;
          }else if(!h2){
            phi2 = b.tks_phi()->at(itk);
            eta2 = b.tks_eta()->at(itk);
            h2 = true;
          }
        }
      }
    }
  }

  NamedFunc MismeasurementCorrection(const string &file_path,
                                     const string &mismeas_scenario,
                                     Variation variation){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);

    string name = "w_mismeas";
    NamedFunc reweight_cut = cp.GetOpt<std::string>("reweight_cut");
    double wgt = 1.;
    switch(variation){
    case Variation::central:
      wgt = cp.GetOpt<double>("w_central");
      break;
    case Variation::up:
      wgt = cp.GetOpt<double>("w_up");
      break;
    case Variation::down:
      wgt = cp.GetOpt<double>("w_down");
      break;
    default:
      wgt = 1.;
      ERROR("Invalid variation: "+to_string(static_cast<int>(variation)));
      break;
    }

    if(reweight_cut.IsScalar()){
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          return reweight_cut.GetScalar(b) ? wgt : 1.;
        });
    }else{
      return NamedFunc(name, [reweight_cut, wgt](const Baby &b){
          NamedFunc::VectorType cuts = reweight_cut.GetVector(b);
          NamedFunc::VectorType wgts(cuts.size());
          for(size_t i = 0; i < cuts.size(); ++i){
            wgts.at(i) = cuts.at(i) ? wgt : 1.;
          }
          return wgts;
        });
    }
  }

  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario){
    ConfigParser cp;
    cp.Load(file_path, mismeas_scenario);
    NamedFunc cut = cp.GetOpt<std::string>("mismeas_cut");
    NamedFunc wgt = cp.GetOpt<std::string>("mismeas_wgt");
    string name = "w_sys_mm";
    if(cut.IsScalar()){
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            return cut.GetScalar(b) ? wgt.GetScalar(b) : 1.;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            return cut.GetScalar(b) ? wgts : NamedFunc::VectorType(wgts.size(), 1.);
          });
      }
    }else{
      if(wgt.IsScalar()){
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts(cuts.size());
            NamedFunc::ScalarType scalar_wgt = wgt.GetScalar(b);
            for(size_t i = 0; i < wgts.size(); ++i){
              wgts.at(i) = cuts.at(i) ? scalar_wgt : 1.;
            }
            return wgts;
          });
      }else{
        return NamedFunc(name, [cut, wgt](const Baby &b){
            NamedFunc::VectorType cuts = cut.GetVector(b);
            NamedFunc::VectorType wgts = wgt.GetVector(b);
            size_t min_size = min(cuts.size(), wgts.size());
            NamedFunc::VectorType out(min_size);
            for(size_t i = 0; i < min_size; ++i){
              out.at(i) = cuts.at(i) ? wgts.at(i) : 1.;
            }
            return out;
          });
      }
    }
  }
}
