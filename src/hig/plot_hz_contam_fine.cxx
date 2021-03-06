///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 35.9;
  float bf = 1.;
  bool incl_nonbb = true;
  bool incl_nonhh = true;
  bool hz_only = false;
  bool zz_only = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higtight/");
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higmc_higtight/");
  string folderHZ(bfolder+"/cms2r0/babymaker/babies/2017_06_01/TChiHZ/merged_higmc_higtight/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  string c_ps = "pass && stitch_met";

  vector<shared_ptr<Process> > procs;
  // procs.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
		// 					    attach_folder(foldermc, mctags["other"]),c_ps));
  // procs.push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
		// 					    attach_folder(foldermc, mctags["qcd"]),c_ps)); 
  // procs.push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
		// 					    attach_folder(foldermc,mctags["vjets"]),c_ps));
  // procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                 attach_folder(foldermc, mctags["ttx"]),c_ps));

  vector<float> sigd({225., 400., 700.});
  vector<string> sigm({"225","400", "700"});
  vector<unsigned> sigcol({kRed+1,kGreen+2, kAzure-6});
  for (unsigned isig(0); isig<sigm.size(); isig++) {
    procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet"));
    procs.push_back(Process::MakeShared<Baby_full>("TChiHZ("+sigm[isig]+",1)", 
      Process::Type::signal, 1, {folderHZ+"*TChiHZ_*"+sigm[isig]+"*.root"}, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet"));
  } 
 

  string filters = "pass_ra2_badmu && met/met_calo<5";
  string baseline = filters+"&& njets>=4 && njets<=5 && nvleps==0 && ntks==0 && !low_dphi";

  string c_2b = "nbdt==2&&nbdm==2";
  string c_3b = "nbdt>=2&&nbdm>=3";
  // string c_4b = "nbdt>=2&&nbdm>=3&&nbdl>=4";
  string hig = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && (higd_am>100 && higd_am<=140)";
  string sbd = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && !(higd_am>100 && higd_am<=140)";
  NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  NamedFunc bf_wgt("bf_wgt", [&](const Baby &b){
    float wgt_ = 1;
    if (b.type()==-999999){
      int nh(0), nh_nonbb(0);
      for (unsigned i(0); i<b.mc_id()->size(); i++) {
        if (b.mc_id()->at(i)==25) nh++;
        if (b.mc_mom()->at(i)==25 && abs(b.mc_id()->at(i))!=5) nh_nonbb++;
      }
      nh_nonbb /=2;
      if (hz_only) {
        wgt_ *= 2*(nh==1);
      } else if (zz_only) {
        wgt_ *= 4*(nh==0);
      } else {
        if (incl_nonhh) wgt_ *= (bf*bf/.25*(nh==2) + 2*bf*(1-bf)/.5*(nh==1) + (1-bf)*(1-bf)/.25*(nh==0));
        else wgt_ *= (bf*bf/.25*(nh==2));
        if (!incl_nonbb) wgt_ *= (nh_nonbb==0);
      }
    }
    return wgt_;
  });

  wgt *= bf_wgt;
  
  //        Cutflow table
  //-------------------------------- 
  PlotMaker pm;
  string tabname = "regions";
  pm.Push<Table>(tabname, vector<TableRow>{
  TableRow("SBD, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("SBD, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  TableRow("SBD, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("SBD, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  TableRow("SBD, 2b", baseline + " && met>300 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>300 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>300 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>300 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("SBD, 4b", baseline + " && met>300 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 4b", baseline + " && met>300 &&" +c_4b+"&&"+hig,0,1, wgt)
	},procs,0);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  Table * yield_table = static_cast<Table*>(pm.Figures()[0].get());
  vector<vector<GammaParams> > allyields;
  for (auto &iproc: procs) 
    allyields.push_back(yield_table->Yield(iproc.get(), lumi));

  unsigned nbins = allyields[0].size();

  vector<vector<float> > ratio, err_up, err_dn, xcoord, xerr;
  vector<float> pow_ratio({ 1, -1});
  for (unsigned isig(0); isig<sigm.size(); isig++){
    for (unsigned ireg(0); ireg<2; ireg++){
      ratio.push_back(vector<float>());
      err_up.push_back(vector<float>());
      err_dn.push_back(vector<float>());
      xcoord.push_back(vector<float>());
      xerr.push_back(vector<float>());
      for (unsigned ibin(ireg); ibin<nbins; ibin+=2){ //store HIG and SBD separately
        float val(1.), valup(1.), valdn(1.);
        vector<vector<float> > entries;
        vector<vector<float> > weights;
        for(size_t iobs=1; iobs<2; iobs--){
          size_t sam = 2*isig + iobs;
          entries.push_back(vector<float>());
          weights.push_back(vector<float>());
          entries.back().push_back(allyields[sam][ibin].NEffective());
          weights.back().push_back(allyields[sam][ibin].Weight());
        }
        val = calcKappa(entries, weights, pow_ratio, valdn, valup);
        if (val<0.001) {val=0.001; valup=.001; valdn=0.001;}
        if (val>3) {val=2.99; valup=0; valdn=0;}
        ratio.back().push_back(val);
        err_up.back().push_back(valup);
        err_dn.back().push_back(valdn);
        xcoord.back().push_back(ibin + sigd[isig]/1000.);
        xerr.back().push_back(0.);
      }
    }
  }

  // for (unsigned i=0; i<ratio.size(); i++){
  //   for (unsigned j=0; j<ratio[i].size(); j++){
  //     cout<<"ratio["<<i<<"]["<<j<<"]"<<ratio[i][j]<<endl;
  //   }
  // }

  PlotOpt opts("txt/plot_styles.txt", "Kappa");

  TCanvas can("can","", 1600, 700);
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
  TLatex klab; klab.SetTextFont(42); klab.SetTextAlign(23);

  float minx = 0., maxx = nbins, miny = 0;
  float maxy = 3;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.SetTitleOffset(0.57,"y");
  histo.SetTitleSize(0.07,"y");

  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  if(lumi-floor(lumi)==0) digits_lumi = 0;
  TString lumi_s = RoundNumber(lumi, digits_lumi);
  double legX(opts.LeftMargin()+0.005), legY(1-0.12), legSingle = 0.04;
  legX = 0.595;
  double legW = 0.30, legH = legSingle*sigm.size();

  TString ytitle = (hz_only ? "HZ":"ZZ");
  ytitle += "/HH ";
  // ytitle += (ibin==0 ? "SBD":"HIG"); 
  ytitle += " Yield";
  histo.SetYTitle(ytitle);
  histo.Draw();

  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.15); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; 
  
  for (unsigned isig(0); isig<sigm.size(); isig++){
    for (unsigned ibin(0); ibin<2; ibin++){
      graph[2*isig+ibin] = TGraphAsymmErrors(ratio[0].size(), 
        &(xcoord[2*isig+ibin][0]), &(ratio[2*isig+ibin][0]),
        &(xerr[2*isig+ibin][0]), &(xerr[2*isig+ibin][0]), 
        &(err_dn[2*isig+ibin][0]), &(err_up[2*isig+ibin][0]));
      graph[2*isig+ibin].SetMarkerStyle(20+4*ibin+isig); 
      graph[2*isig+ibin].SetMarkerSize(1.5);
      graph[2*isig+ibin].SetMarkerColor(sigcol[isig]);
      graph[2*isig+ibin].SetLineColor(sigcol[isig]); 
      graph[2*isig+ibin].SetLineWidth(2);
      graph[2*isig+ibin].Draw("p0 same");
      TString pname = sigm[isig];
      pname += (ibin==0 ? ", SBD":", HIG"); 
      leg.AddEntry(&graph[2*isig+ibin], pname, "p");
    }
  }
  leg.Draw();

  // for (unsigned ifit(0); ifit<nbins/3; ifit++){
  //   TF1 f1("f1","gaus",1,3);
  //   graph[].Fit("f1","R");
  // }

  TString fname="plots/contam_";
  // fname += (ibin==0 ? "sbd":"hig");
  if (hz_only) fname += "_hz";
  else fname += "_zz";
  fname += "_fine.pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl; 
  

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"excl_nonbb", no_argument, 0, 0},
      {"excl_nonhh", no_argument, 0, 0}, 
      {"hz", no_argument, 0, 0}, 
      {"zz", no_argument, 0, 0}, 
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 0:
      optname = long_options[option_index].name;
      if(optname == "excl_nonbb"){
        incl_nonbb = false;
      }else if(optname == "excl_nonhh"){
        incl_nonhh = false;
      }else if(optname == "hz"){
        hz_only = true;
      }else if(optname == "zz"){
        zz_only = true;
      }else if(optname == "bf"){
        bf = atof(optarg);
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
