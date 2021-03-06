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
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 19.3;
  // choose processes to include, options are: "ttx", "vjets", "singlet", "qcd", "other", "ttonly"
  // set<string> proc_types = {"ttx", "vjets", "singlet", "qcd", "other"}; // for default data/MC
  set<string> proc_types = {"ttx", "vjets", "singlet"}; // for default data/MC
  // set<string> proc_types = {}; // to make signal plots only
  // set<string> proc_types = {"ttonly"};
  // signal points to include and their colors
  vector<string> sigm = {"250","400"}; 
  vector<int> sig_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
  //for signal plots only
  // vector<string> sigm = {"175","225","350","700","1000"}; 
  // vector<int> sig_colors = {kMagenta+2 , kGreen, kRed, kBlue, kAzure+10}; // need sigm.size() >= sig_colors.size()
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc; //ordered in number of leptons
  foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["ttonly"]  = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TCHiHH/merged_higmc_unskimmed/";
  foldersig += "*TChiHH_mGluino-";

  string c_ps = "pass && stitch";

  vector<shared_ptr<Process>> procs = vector<shared_ptr<Process> >();
  if (proc_types.find("other")!=proc_types.end())       
    procs.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
						    attach_folder(foldermc, mctags["other"]),c_ps));
  if (proc_types.find("singlet")!=proc_types.end())       
    procs.push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, 1,
						    attach_folder(foldermc,mctags["singlet"]),c_ps));
  if (proc_types.find("qcd")!=proc_types.end())       
    procs.push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
						    attach_folder(foldermc, mctags["qcd"]),c_ps)); 
  if (proc_types.find("vjets")!=proc_types.end())       
    procs.push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
						    attach_folder(foldermc,mctags["vjets"]),c_ps));
  if (proc_types.find("ttx")!=proc_types.end()) 
    procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
						    attach_folder(foldermc, mctags["ttx"]),c_ps));
  if (proc_types.find("ttonly")!=proc_types.end()) 
    procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
						    attach_folder(foldermc,mctags["ttonly"]),c_ps));
  if (proc_types.size()==0) { // have to pretend signal is background, otherwise crashes
    for (unsigned isig(0); isig<sigm.size(); isig++)
      procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)",
              Process::Type::background, 2, {foldersig+sigm[isig]+"*.root"}, "1"));
  } else {
    for (unsigned isig(0); isig<sigm.size(); isig++)
      procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
						      Process::Type::signal, 2, {foldersig+sigm[isig]+"*.root"}, "1"));
  }


 

  PlotMaker pm;

  //just to be pretty... already in skim...
  string njcut = "njets>=4 && njets<=5";
  string c_2b = "nbt==2&&nbm==2";
  string c_3b = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b = "nbt>=2&&nbm>=3&&nbl>=4";

  map<string, string> xcuts; // useful additional cut definitions
  xcuts["drmax"] = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40";
  xcuts["fullhig"] = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && (hig_am>100 && hig_am<=140)";
  xcuts["fullsbd"] = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && !(hig_am>100 && hig_am<=140)";

  string baseline = "pass && stitch";// && nvleps==0 && met/met_calo<5";
  string noleps_dphi = "nvleps==0 && ntks==0 && !low_dphi";
  string higsig = "hig_dm<20 && hig_am>100 && hig_am<140";

  //        Cutflow table
  //-------------------------------- 
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;
  pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("No selection", "1",0,0, wgt),
      TableRow("$\\text{4-5 jets}$", baseline+"&&"+njcut,0,0, wgt),
      TableRow("Lepton veto and dphi", baseline+"&&"+njcut+"&&"+noleps_dphi,0,0, wgt),
      TableRow("3b category", baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_3b,0,0, wgt),
      TableRow("3b, $\\Delta R_{max}$", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_3b+"&&hig_drmax<2.2",0,0, wgt),
      TableRow("3b, HIG region", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_3b+"&&hig_drmax<2.2&&"+higsig,0,1, wgt),
      TableRow("3b, HIG region, MET>135", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_3b+"&&hig_drmax<2.2&&"+higsig+"&&met>135",0,1, wgt),
      TableRow("4b category", baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_4b,0,0, wgt),
      TableRow("4b, $\\Delta R_{max}$", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_4b+"&&hig_drmax<2.2",0,0, wgt),
      TableRow("4b, HIG region", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_4b+"&&hig_drmax<2.2&&"+higsig,0,1, wgt),
      TableRow("4b, HIG region, MET>135", 
        baseline+"&&"+njcut+"&&"+noleps_dphi+"&&"+c_4b+"&&hig_drmax<2.2&&"+higsig+"&&met>135",0,1, wgt),

	// TableRow("$E_{T}^{miss} > 150$, 2+ CSVM, $\\text{4-5 jets}$, $0\\ell$", 
	// 	 baseline + " && nbm>=2 &&"+njcut,0,0, wgt),
	// TableRow("2+ CSVT (2b + 3b + 4b)", 
	// 	 baseline + " && nbt>=2 &&"+njcut,0,0, wgt),
	// TableRow("Track veto", 
	// 	 baseline + " && ntks==0 && nbt>=2 &&"+njcut,0,0, wgt),
	// TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$",        
	// 	 baseline + " && ntks==0 && nbt>=2 && nbt>=2 &&"+njcut+"&& !low_dphi",0,0, wgt),
	// TableRow("$|\\Delta m| < 40$",     
	// 	 baseline + " && ntks==0 && nbt>=2 &&"+njcut+"&& !low_dphi && hig_dm<=40",0,0, wgt),
	// TableRow("$\\Delta R_{\\text{max}} < 2.2$",                    
	// 	 baseline +" && ntks==0 && nbt>=2 &&"+njcut+"&& !low_dphi && hig_drmax<=2.2 && hig_dm<=40",0,0,wgt),

	// TableRow("2b 150-200", 
	// 	 baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 && nbt==2 && nbm==2 && met>150 && met<=200 &&"+njcut+"&& !low_dphi",1,0, wgt),
 //  TableRow("2b 200-300", 
 //     baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 && nbt==2 && nbm==2 && met>200 && met<=300 &&"+njcut+"&& !low_dphi",1,0, wgt),
 //  TableRow("2b 300", 
 //     baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 && nbt==2 && nbm==2 && met>300 &&"+njcut+"&& !low_dphi",1,0, wgt),
	// TableRow("3b + 4b", 
	// 	 baseline +" && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 &&("+c_3b+"||"+c_4b+")&&"+njcut+"&& !low_dphi",0,0,wgt),
	// TableRow("4b", 
	// 	 baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 &&"+c_4b+"&&"+njcut+"&& !low_dphi",0,0, wgt),
	// TableRow("$E_{T}^{miss}>200$", 
	// 	 baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 && met>200 &&"+c_4b+"&&"+njcut+"&& !low_dphi",0,0,wgt),
	// TableRow("$E_{T}^{miss}>300$", 
	// 	 baseline + " && ntks==0 && hig_dm <= 40 && hig_am<=200 && hig_drmax<=2.2 && met>300 &&"+c_4b+"&&"+njcut+"&& !low_dphi",0,0,wgt),

  // keep for cross-checking pie charts
  //-------------------------------------
  // TableRow("2b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_2b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("2b, $200<E_{T}^{miss}\\leq$300", 
  //  baseline + " && met>200 && met<=300 &&"+c_2b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("2b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_2b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("3b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_3b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("3b, $200<E_{T}^{miss}\\leq$300", 
  //  baseline + " && met>200 && met<=300 &&"+c_3b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("3b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_3b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("4b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_4b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("4b, $200<E_{T}^{miss}\\leq$300", 
  //  baseline + " && met>200 && met<=300 &&"+c_4b+"&&"+njcut+"&&"+xcuts["drmax"],0,0, wgt),
  // TableRow("4b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_4b+"&&"+njcut+"&&"+xcuts["drmax"],0,2, wgt),

  // TableRow("HIG, 2b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 2b, $200<E_{T}^{miss}\\leq$300", 
  //  baseline + " && met>200 && met<=300 &&"+c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 2b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 3b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 3b, $200<E_{T}^{miss}\\leq$300", 
  //  baseline + " && met>200 && met<=300 &&"+c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 3b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 4b, $150<E_{T}^{miss}\\leq$200", 
  //  baseline + " && met<=200 &&"+c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 4b, $200<E_{T}^{miss}\\leq$300", 
  //   baseline + " && met>200 && met<=300 &&"+c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  // TableRow("HIG, 4b, $E_{T}^{miss}>$300",         
  //  baseline + " && met>300 &&"+c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),

	},procs,0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
