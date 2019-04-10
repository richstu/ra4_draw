#ifndef H_FUNCTIONS
#define H_FUNCTIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
 
namespace Functions{
  extern const NamedFunc n_isr_match;
  extern const NamedFunc ntrub;
  extern const NamedFunc n_mus_bad;
  extern const NamedFunc n_mus_bad_dupl;
  extern const NamedFunc n_mus_bad_trkmu;
  extern const NamedFunc njets_weights_ttisr;
  extern const NamedFunc njets_weights_visr;
  extern const NamedFunc min_dphi_lep_met;
  extern const NamedFunc max_dphi_lep_met;
  extern const NamedFunc min_dphi_lep_jet;
  extern const NamedFunc max_dphi_lep_jet;
  extern const NamedFunc min_dphi_met_jet;
  extern const NamedFunc max_dphi_met_jet;
  extern const NamedFunc min_dr_lep_jet;
  extern const NamedFunc max_dr_lep_jet;
  extern const NamedFunc nbm_moriond;
  extern const NamedFunc ntop_loose_decor;
  extern const NamedFunc ntop_med_decor;
  extern const NamedFunc ntop_tight_decor;
  extern const NamedFunc ntop_loose_nom;
  extern const NamedFunc ntop_med_nom;
  extern const NamedFunc ntop_tight_nom;
  extern const NamedFunc offshellw;
  extern const NamedFunc nbd;
  extern const NamedFunc hem_veto;
  extern const NamedFunc wgt_run2;
  extern const NamedFunc wgt_run2_partial;
  extern const NamedFunc wgt_run2_partial_2017BCDE;
  extern const NamedFunc wgt_run2_partial_2017F;
  extern const NamedFunc wgt_run2_partial_2017;
  extern const NamedFunc wgt_run2_partial_2018AB;
  extern const NamedFunc wgt_run2_partial_2018D;
  extern const NamedFunc wgt_run2_partial_2018;
  extern const NamedFunc mht_ratio;
  extern const NamedFunc fake_met;
  extern const NamedFunc trig_run2;
  extern const NamedFunc eff_trig_run2;

  float wnpv2017(const Baby &b);
  bool IsGoodJet(const Baby &b, std::size_t ijet);
  bool IsGoodak8Jet(const Baby &b, std::size_t ijet);
  bool IsGoodElectron(const Baby &b, std::size_t iel);
  bool IsGoodMuon(const Baby &b, std::size_t imu);
  bool IsGoodTrack(const Baby &b, std::size_t itk);

  NamedFunc::ScalarType NJetsWeights_ttISR(const Baby &b, bool use_baby_nisr);
  NamedFunc::ScalarType NJetsWeights_vISR(const Baby &b);
  int NISRMatch(const Baby &b);

  void DileptonAngles(const Baby &b,
		      NamedFunc::ScalarType &eta1, NamedFunc::ScalarType &phi1,
		      NamedFunc::ScalarType &eta2, NamedFunc::ScalarType &phi2);

  enum class Variation{central, up, down};
  NamedFunc MismeasurementCorrection(const std::string &file_path,
                                     const std::string &mismeas_scenario,
                                     Variation variation = Variation::central);
  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario);
}

#endif
