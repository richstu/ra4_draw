#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  bool quick = false;
  bool regions = false;
  double lumi = 135;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  // string foldermc = "/afs/cern.ch/work/a/ana/cms2r0/babymaker/babies/2018_07_25/mc/merged_mcbase_standard/";
  string foldermc = "/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/";
  string foldersig = "/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
 
  Palette colors("txt/colors.txt", "default");

  // Background samples
  string bkg_pass = "pass && stitch_met";
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*.root", foldermc+"*_ZJet*.root", foldermc+"*_ttHTobb_M125_*.root",
        foldermc+"*_TTTT_*.root", foldermc+"*_WH_HToBB*.root", foldermc+"*_ZH_HToBB*.root", 
        foldermc+"*_WWTo*.root", foldermc+"*_WZ*.root", foldermc+"*_ZZ_*.root"}, bkg_pass);
  auto qcd = Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("qcd"),
    {foldermc+"*QCD_HT*0_Tune*.root", foldermc+"*QCD_HT*Inf_Tune*.root"}, bkg_pass);
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTGJets*.root", foldermc+"*_TTWJets*.root", foldermc+"*_TTZ*.root"}, bkg_pass);
  auto st = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, bkg_pass);
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsTo*.root"}, bkg_pass);  
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root"}, "ntruleps==1 &&"+bkg_pass);
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root"}, "ntruleps==2 &&"+bkg_pass);

  auto non_tt = Process::MakeShared<Baby_full>("Non-t#bar{t}", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*.root", foldermc+"*_ZJet*.root", foldermc+"*_ttHTobb_M125_*.root",
        foldermc+"*_TTTT_*.root", foldermc+"*_WH_HToBB*.root", foldermc+"*_ZH_HToBB*.root", 
        foldermc+"*_WWTo*.root", foldermc+"*_WZ*.root", foldermc+"*_ZZ_*.root",
        foldermc+"*QCD_HT*0_Tune*.root", foldermc+"*QCD_HT*Inf_Tune*.root",
        foldermc+"*_TTGJets*.root", foldermc+"*_TTWJets*.root", foldermc+"*_TTZ*.root",
        foldermc+"*_ST_*.root", foldermc+"*_WJetsTo*.root"});

  auto all_bkg = Process::MakeShared<Baby_full>("Non-t#bar{t}", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*.root", foldermc+"*_ZJet*.root", foldermc+"*_ttHTobb_M125_*.root",
        foldermc+"*_TTTT_*.root", foldermc+"*_WH_HToBB*.root", foldermc+"*_ZH_HToBB*.root", 
        foldermc+"*_WWTo*.root", foldermc+"*_WZ*.root", foldermc+"*_ZZ_*.root",
        foldermc+"*QCD_HT*0_Tune*.root", foldermc+"*QCD_HT*Inf_Tune*.root",
        foldermc+"*_TTGJets*.root", foldermc+"*_TTWJets*.root", foldermc+"*_TTZ*.root",
        foldermc+"*_ST_*.root", foldermc+"*_WJetsTo*.root"});

  auto sig_nc = Process::MakeShared<Baby_full>("(2200,100)", Process::Type::signal, colors("t1tttt"),
    {foldersig+"*SMS-T1tttt_mGluino-2200_mLSP-100_*.root"});
  auto sig_c = Process::MakeShared<Baby_full>("(1600,1250)", Process::Type::signal, colors("t1tttt"),
    {foldersig+"*SMS-T1tttt_mGluino-1600_mLSP-1250*.root"});
  sig_c->SetLineStyle(2);

  vector<shared_ptr<Process> > procs_info, procs;
  if (quick){
    procs_info = procs = {sig_nc, sig_c};
  } else {
    procs_info = {other, qcd, ttv, st, wjets, tt1l, tt2l, sig_nc, sig_c};
    procs = {non_tt, tt1l, tt2l, sig_nc, sig_c};
  }

  // --------------------------------------------
  //          Cutflow table
  //---------------------------------------------
  
  string pass = "pass_ra2_badmu && met/met_calo<5";
  string base = pass+"&& nleps==1 && nveto==0 && mj14>250 && st>500 && met>200"; // keep in mind the skim also has njets >=5!!
  vector<TableRow> cutflow_rows = vector<TableRow>();
  cutflow_rows.push_back(TableRow("No selection", pass,0,0,"weight"));
  cutflow_rows.push_back(TableRow("$1\\ell$", pass+" && nleps==1",0,0,"weight"));
  cutflow_rows.push_back(TableRow("Track veto", pass+" && nleps==1 && nveto==0",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$M_J>250$ GeV", pass+" && nleps==1 && nveto==0 && mj14>250",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$S_{T}>$500 GeV", pass+" && nleps==1 && nveto==0 && mj14>250 && st>500",0,0,"weight"));
  cutflow_rows.push_back(TableRow("MET$>$150 GeV", pass+" && nleps==1 && nveto==0 && mj14>250 && st>500 && met>150",0,0,"weight"));
  cutflow_rows.push_back(TableRow("MET$>$200 GeV", pass+" && nleps==1 && nveto==0 && mj14>250 && st>500 && met>200",0,0,"weight"));

  cutflow_rows.push_back(TableRow("$N_{\\rm jets}\\geq6$", base+"&& njets>=6",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$N_{\\rm b}\\geq1$", base+"&& nbm>=1 && njets>=6",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$m_T>140$ GeV", base+"&& mt>140 && nbm>=1 && njets>=6",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$M_J>400$ GeV", base+"&& mt>140 && mj14>400 && nbm>=1 && njets>=6",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$N_{\\rm b}\\geq2$", base+"&& mt>140 && mj14>400 && nbm>=2 && njets>=6",0,0,"weight"));
  cutflow_rows.push_back(TableRow("MET$>350$ GeV", base+"&& mt>140 && mj14>400 && nbm>=2 && njets>=6 && met>350",0,0,"weight"));
  cutflow_rows.push_back(TableRow("MET$>500$ GeV", base+"&& mt>140 && mj14>400 && nbm>=2 && njets>=6 && met>500",0,0,"weight"));
  cutflow_rows.push_back(TableRow("$N_{\\rm jets}\\geq9$", base+"&& mt>140 && mj14>400 && nbm>=2 && njets>=9 && met>500",0,0,"weight"));


  PlotMaker pm;
  pm.Push<Table>("cutflow_mj", cutflow_rows, procs, false);  

  // --------------------------------------------
  //          Regions table
  //---------------------------------------------

  base += "&& njets>=6 && nbm>=1";

  vector<string> abcdcuts; 
  abcdcuts.push_back("mt<=140 && mj14<=400");
  abcdcuts.push_back("mt<=140 && mj14>400");
  abcdcuts.push_back("mt>140  && mj14<=400");
  abcdcuts.push_back("mt>140  && mj14>400");

  vector<string> nj, njcuts;
  nj.push_back("6-8j"); njcuts.push_back("njets>=6 && njets<=8");
  nj.push_back("$\\geq$9j"); njcuts.push_back("njets>=9");

  vector<string> nb, nbcuts;
  nb.push_back("1b"); nbcuts.push_back("nbm==1");
  nb.push_back("2b"); nbcuts.push_back("nbm==2");
  nb.push_back("$\\geq$3b"); nbcuts.push_back("nbm>=3");
  // nb.push_back("3b"); nbcuts.push_back("nbm==3");
  // nb.push_back("$\\geq$4b"); nbcuts.push_back("nbm>=4");

  vector<string> metcuts;
  // metcuts.push_back("met>100 && met<=150");
  // metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");

  vector<TableRow> regions_rows = vector<TableRow>();
  for (auto &imet: metcuts) {
    regions_rows.push_back(TableRow("$"+CodeToLatex(imet)+"$"));
    regions_rows.push_back(TableRow("R1: all $N_{jets}$, $N_b$", 
      base+"&&"+imet+"&&"+abcdcuts[0],0,0,"weight"));
    for (unsigned i(0); i<nb.size(); i++) 
        for (unsigned j(0); j<nj.size(); j++)
          regions_rows.push_back(TableRow("R2: "+nb[i]+", "+nj[j], 
            base+"&&"+abcdcuts[1]+"&&"+imet+"&&"+nbcuts[i]+"&&"+njcuts[j], 0, 0, "weight"));

    regions_rows.push_back(TableRow("R3: all $N_{jets}$, $N_b$", 
      base+"&&"+imet+"&&"+abcdcuts[2],0,1,"weight"));
    for (unsigned i(0); i<nb.size(); i++) 
        for (unsigned j(0); j<nj.size(); j++)
          regions_rows.push_back(TableRow("R4: "+nb[i]+", "+nj[j], 
            base+"&&"+abcdcuts[3]+"&&"+imet+"&&"+nbcuts[i]+"&&"+njcuts[j], 0, 0, "weight"));

    ReplaceAll(imet, ">","_");
    ReplaceAll(imet, " && met<=","to");
  }
  if (regions) pm.Push<Table>("regions_mj", regions_rows, procs, true);  
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},
      {"quick", no_argument, 0, 'q'},
      {"regions", no_argument, 0, 'r'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:qr", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 'q':
      quick = true;
      break;
    case 'r':
      regions = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      printf("Bad option! Found option name %s\n", optname.c_str());
      exit(1);      
      break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}