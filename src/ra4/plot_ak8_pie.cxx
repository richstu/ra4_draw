///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/plot_opt.hpp"

namespace{
  float lumi = 135.;
  bool do_allplots = false;
}

using namespace std;

NamedFunc::ScalarType imax_NR(const Baby &b);

NamedFunc offshellw("offshellw",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=24) continue;
      if (b.mc_mass()->at(i) > 140.) {
        return 1;
      }
    }
    return 0;
  });

//// Lowers the njets requirements for dilepton bins
TString lowerNjets(TString &cut){
  TString lowcut = cut;
  for(int nj=6; nj<=9; nj++){
    TString nj_s = ""; nj_s += nj;
    TString njlo_s = ""; njlo_s += nj-1;
    lowcut.ReplaceAll(nj_s, njlo_s);
  }
  return lowcut;
}

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

	const NamedFunc max_NR("Highest_top_score",[&](const Baby &b){
		double score(-1);
		int imax = imax_NR(b);
		if(imax >= 0) score = b.ak8jets_nom_raw_top()->at(imax);
		return score;
	});
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
//   if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  Palette colors("txt/colors.txt", "default");
  string mc_standard_path("/net/cms2/cms2r0/babymaker/babies/2018_08_03/mc/merged_mcbase_standard/");
  string ttbar1L(mc_standard_path+"*_TTJets*SingleLept*.root");
  string ttbar2L(mc_standard_path+"*_TTJets*DiLept*.root");
  string ttbarHT(mc_standard_path+"*_TTJets*HT*.root");
  string  signal1(mc_standard_path+"*mGluino-1200_mLSP-800*.root");
  string  signal2(mc_standard_path+"*mGluino-2000_mLSP-100*.root");

  string ntupletag = "";
  set<string> allfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
         foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
         foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
         foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
         foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",
         foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
         foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
         foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
         foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",foldermc+"*_ZZ_*"+ntupletag+"*.root"
       };

  const NamedFunc ak8_Ntops_nom("ak8_Ntops_nom", [](const Baby &b) ->NamedFunc::ScalarType{
    int tops(0);
	if(b.ak8jets_pt()->size() > 0) {
	  	for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
			if(b.ak8jets_nom_bin_top()->at(iak8) > 0.1883) tops++;
	}	}
	return tops;
  	});
  const NamedFunc ak8_Ntops_decor("ak8_Ntops_decor", [](const Baby &b) ->NamedFunc::ScalarType{
    int tops(0);
	if(b.ak8jets_pt()->size() > 0) {
	  	for(size_t iak8(0); iak8 < b.ak8jets_pt()->size(); iak8++) {
			if(b.ak8jets_decor_bin_top()->at(iak8) > 0.04738) tops++;
	} 	}
	return tops;
  	});

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch_met && mj14>250 && nleps==1 && st>500 && met>100 && njets>=5 && nveto==0";
//   string baseline_s = "mj14>250 && nleps>=1 && met>100 && njets>=5 && st<10000 && pass_ra2_badmu && met/met_calo<5";

  map<string, vector<shared_ptr<Process> > > procs;
  vector<string> cats_lab = {"#geq2#nu^{prompt}", "#leq1#kern[.1]{#nu^{pr.}}, mismeas.", "#leq1#kern[.1]{#nu^{pr.}}, #geq1#kern[.1]{#nu}^{non-prompt}","#leq1#kern[.1]{#nu^{pr.}}, off-shell W"};
//   vector<string> samples = {"base","ak8","ak8_1LooseNom","ak8_1LooseDecor"};
//   vector<NamedFunc> sample_cuts = {"1","nak8jets > 0",ak8_Ntops_nom >= 1,ak8_Ntops_decor >= 1};
  NamedFunc multNeu = "(type==5000 || type==13000 || type==15000 || type==16000)";
  NamedFunc multNeu2l = "(ntruleps>=2 || (ntruleps<=1&&(type==5000 || type==13000 || type==15000 || type==16000)))";

  procs["base"] = vector<shared_ptr<Process> >();
  procs["base"].push_back(Process::MakeShared<Baby_full>("1L t#bar{t}", Process::Type::background, colors("tt_1l"),  
								{ttbar1L,ttbarHT}, baseline && "ntruleps<=1&&stitch_met"));
  procs["base"].push_back(Process::MakeShared<Baby_full>("2L t#bar{t}", Process::Type::background, colors("tt_2l"),  
								{ttbar2L,ttbarHT}, baseline && "ntruleps>=2&&stitch_met"));
  // With atleast 1 ak8
  procs["ak8"] = vector<shared_ptr<Process> >();
  procs["ak8"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(0), Process::Type::background, kCyan-3, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "nak8jets > 0" && multNeu2l));
  procs["ak8"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(1), Process::Type::background, kRed-4,  
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "nak8jets > 0 && ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["ak8"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(2), Process::Type::background, kGreen-3,
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "nak8jets > 0 && ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw==0.));
  procs["ak8"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(3), Process::Type::background, kOrange, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "nak8jets > 0 && ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw>0.));
  // With atleast 1 Loose Nominal top
  procs["ak8_1LNom"] = vector<shared_ptr<Process> >();
  procs["ak8_1LNom"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(0), Process::Type::background, kCyan-3, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_nom >= 1 && multNeu2l));
  procs["ak8_1LNom"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(1), Process::Type::background, kRed-4,  
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_nom >= 1 && "ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["ak8_1LNom"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(2), Process::Type::background, kGreen-3,
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_nom >= 1 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw==0.));
  procs["ak8_1LNom"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(3), Process::Type::background, kOrange, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_nom >= 1 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw>0.));
  // With atleast 1 Loose Decorrelated top
  procs["ak8_1LDecor"] = vector<shared_ptr<Process> >();
  procs["ak8_1LDecor"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(0), Process::Type::background, kCyan-3, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_decor >= 1 && multNeu2l));
  procs["ak8_1LDecor"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(1), Process::Type::background, kRed-4,  
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_decor >= 1 && "ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["ak8_1LDecor"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(2), Process::Type::background, kGreen-3,
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_decor >= 1 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw==0.));
  procs["ak8_1LDecor"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(3), Process::Type::background, kOrange, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && ak8_Ntops_decor >= 1 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw>0.));
  // With atleast 1 NR top
  procs["ak8_1NR"] = vector<shared_ptr<Process> >();
  procs["ak8_1NR"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(0), Process::Type::background, kCyan-3, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && max_NR > 0.4 && multNeu2l));
  procs["ak8_1NR"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(1), Process::Type::background, kRed-4,  
								{ttbar1L,ttbarHT,ttbar2L}, baseline && max_NR > 0.4 && "ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["ak8_1NR"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(2), Process::Type::background, kGreen-3,
								{ttbar1L,ttbarHT,ttbar2L}, baseline && max_NR > 0.4 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw==0.));
  procs["ak8_1NR"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(3), Process::Type::background, kOrange, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && max_NR > 0.4 && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw>0.));

  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(0), Process::Type::background, kCyan-3, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && multNeu2l));
  procs["cats"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(1), Process::Type::background, kRed-4,  
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(2), Process::Type::background, kGreen-3,
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw==0.));
  procs["cats"].push_back(Process::MakeShared<Baby_full>(cats_lab.at(3), Process::Type::background, kOrange, 
								{ttbar1L,ttbarHT,ttbar2L}, baseline && "ntruleps<=1 && mt_tru>140"  && !multNeu && offshellw>0.));

  PlotMaker pm;

  vector<TString> metcuts;
  //metcuts.push_back("met>100");
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");

  vector<TString> nbcuts;
  //nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm>=1");
  if(do_allplots){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
    nbcuts.push_back("nbm>=3");
  }

  vector<TString> njcuts;
   njcuts.push_back("njets==5");
//   njcuts.push_back("njets>=5");
//   njcuts.push_back("njets>=6");
  njcuts.push_back("njets>=6 && njets<=8");
//   njcuts.push_back("njets>=9");
  
  vector<TString> mtcuts({"mt>140"});
  
  // Adding nleps==1 cuts
  vector<TString> cuts;
  vector<TableRow> table_cuts;

  //// nleps = 1
  for(auto &imet: metcuts) 
    for(auto &inb: nbcuts) 
      for(auto &inj: njcuts) 
	 			for(auto &imt: mtcuts) 
		  		cuts.push_back("nleps==1 && nveto==0 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);

  //// nleps = 2 and ll+lv
//  for(auto &imet: metcuts) 
//    for(size_t ind=0; ind<njcuts.size(); ind++) {
//      if(njcuts[ind]=="njets==5") continue;
//      cuts.push_back("nleps==2 && nbm<=2 && "+imet+"&&"+lowerNjets(njcuts[ind]));
//      if(do_allplots) 
//		cuts.push_back("((nleps==2 && nbm<=2 && "+lowerNjets(njcuts[ind])+") || (nleps==1 && nveto==1 && nbm>=1 && mt>140 && "+njcuts[ind]+")) && "+imet);
//    }

  for(size_t icut=0; icut<cuts.size(); icut++)
		table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
  for(auto &ipr: procs) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts, ipr.second, true, true, true, true);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

NamedFunc::ScalarType imax_NR(const Baby &b) {
	double max(0.), score(0);
	int imax(-1);
	for(size_t ijet = 0; ijet < b.ak8jets_pt()->size(); ijet++) {
		score = b.ak8jets_nom_raw_top()->at(ijet);
// 		if(score > max && b.ak8jets_pt()->at(ijet) >= 300) {
		if(score > max && b.ak8jets_pt()->at(ijet) >= 30) {
			max = score;
			imax = ijet;
		}
	}
	return imax;
}
