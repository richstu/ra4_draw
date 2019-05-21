///// plot_ratios: plots rMJ and rmT, and combinations of these

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace Functions;

namespace{
  bool debug = false;
  float lumi=1.;
  int year=0;
  TString tag = "";
  enum Regions {r1, r2, r3, r4};

  struct oneplot{
    TString name;
    TString baseline;
    vector<TString> bincuts;
  };
}


void plotRatio(vector<vector<vector<GammaParams> > > &allyields, oneplot &plotdef,
	       vector<vector<vector<int> > > &indices, vector<TString> &leglabels,
	       vector<shared_ptr<Process> > &procs);
void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline, vector<shared_ptr<Process> > &procs);
void GetOptions(int argc, char *argv[]);


int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<int> years;
  if (year==0) years = {2016, 2017, 2018};
  else years = {year};

  map<int, string> foldermc;
  foldermc[2016] = bfolder+"/cms2r0/babymaker/babies/2019_01_11/mc/merged_mcbase_stdnj5/";
  foldermc[2017] = bfolder+"/cms2r0/babymaker/babies/2018_12_17/mc/merged_mcbase_stdnj5/";
  foldermc[2018] = bfolder+"/cms2r0/babymaker/babies/2019_03_30/mc/merged_mcbase_stdnj5/"; 

  set<string> vnames_other = {
    "_WJetsToLNu_HT","_ST_","_TTW","_TTZ", "_DYJetsToLL_M-50_HT","_ZJet","_ttH",
     "_TTGJets","_TTTT","_WH_HToBB","_ZH_HToBB","_WWTo","_WZ","_ZZ_","QCD_HT*0_Tune","QCD_HT*Inf_Tune"
  };

  set<string> tt_files, all_files;
  for (auto &yr: years) {
    // all_files.insert(foldermc[yr]+"*_TTJets_DiLept_T*.root");
    tt_files.insert(foldermc[yr]+"*_TTJets*Lept*.root");
    all_files.insert(tt_files.begin(), tt_files.end());
    for(auto name : vnames_other)
      all_files.insert(foldermc[yr] + "*" + name + "*.root");
  }

  string baseline = "nleps==1 && nveto==0 && mj14>250 && st>500 && met>100 && njets>=5 && nbdm>=1";
  NamedFunc baselinef = baseline && Functions::hem_veto && "st<10000 && pass_ra2_badmu && met/met_calo<5";

  vector<string> nb;
  nb.push_back("nbdm==1");
  nb.push_back("nbdm==2");
  nb.push_back("nbdm>=3");
  vector<int>cols = {kBlue,kRed,kGreen+3};

  vector<shared_ptr<Process> > procs;
  for (size_t inb(0); inb<nb.size(); inb++)
    procs.push_back(Process::MakeShared<Baby_full>(CodeToRootTex(nb[inb]), Process::Type::background, cols[inb], 
      all_files, baselinef && ("pass && stitch_met &&"+nb[inb])));
  
  vector<string> met, mjl, mjh;
  met.push_back("met>200 && met<=350"); mjl.push_back("400"); mjh.push_back("500");
  met.push_back("met>350 && met<=500"); mjl.push_back("450"); mjh.push_back("650");
  met.push_back("met>500");             mjl.push_back("500"); mjh.push_back("800");

  vector<TString> abcdcuts;
  // if (tag.Contains("lowmj")) abcdcuts  = {"mt<=140 && mj14<=MJ1X",
  //                               "mt<=140 && mj14> MJ1X && mj14<=MJ2X",
  //                               "mt>140  && mj14<=MJ1X",
  //                               "mt>140  && mj14> MJ1X && mj14<=MJ2X"};
  // else              abcdcuts  = {"mt<=140 && mj14<=MJ1X",
  //                               "mt<=140 && mj14> MJ2X",
  //                               "mt>140  && mj14<=MJ1X",
  //                               "mt>140  && mj14> MJ2X"};
  abcdcuts  = {"mt<=140",
                                "mt<=140",
                                "mt>140",
                                "mt>140"};

  NamedFunc wgt = Functions::wgt_run2 * Functions::eff_trig_run2;// * wnb;
  size_t Nabcd = abcdcuts.size();

  // Makes a plot for each vector in plotcuts
  vector<oneplot> plotcuts;
  plotcuts.push_back({tag+"_njets", "met>200 && met<=350", {"njets==5 && mj14<=400",             "njets==6 && mj14<=400",             "njets==7 && mj14<=400",             "njets>=8 && mj14<=400",
                                                            "njets==5 && mj14>400 && mj14<=500", "njets==6 && mj14>400 && mj14<=500", "njets==7 && mj14>400 && mj14<=500", "njets>=8 && mj14>400 && mj14<=500", 
                                                            "njets==5 && mj14>500",              "njets==6 && mj14>500",              "njets==7 && mj14>500",              "njets>=8 && mj14>500"}});
  plotcuts.push_back({tag+"_njets", "met>350 && met<=500", {"njets==5 && mj14<=450",             "njets==6 && mj14<=450",             "njets==7 && mj14<=450",             "njets==8 && mj14<=450",
                                                            "njets==5 && mj14>450 && mj14<=650", "njets==6 && mj14>450 && mj14<=650", "njets==7 && mj14>450 && mj14<=650", "njets==8 && mj14>450 && mj14<=650",
                                                            "njets==5 && mj14>650",              "njets==6 && mj14>650",              "njets==7 && mj14>650",              "njets==8 && mj14>650"}});
  plotcuts.push_back({tag+"_njets", "met>500"            , {"njets==5 && mj14<=500",             "njets==6 && mj14<=500",             "njets==7 && mj14<=500",             "njets==8 && mj14<=500",
                                                            "njets==5 && mj14>650 && mj14<=800", "njets==6 && mj14>650 && mj14<=800", "njets==7 && mj14>650 && mj14<=800", "njets==8 && mj14>650 && mj14<=800",
                                                            "njets==5 && mj14>800",              "njets==6 && mj14>800",              "njets==7 && mj14>800",              "njets==8 && mj14>800"}});
  

  PlotMaker pm;
  vector<vector<vector<TString> > > allcuts(plotcuts.size(), vector<vector<TString> > (Nabcd));
  for(size_t iplot=0; iplot<plotcuts.size(); iplot++){
    for(size_t iabcd=0; iabcd<abcdcuts.size(); iabcd++){
      vector<TableRow> table_cuts;
      for(size_t ibin=0; ibin<plotcuts[iplot].bincuts.size(); ibin++){
        TString _abcdcuts = abcdcuts[iabcd];
        // _abcdcuts.ReplaceAll("MJ1X", mjl[iplot]).ReplaceAll("MJ2X", mjh[iplot]); //iplot==imet
        TString totcut=plotcuts[iplot].baseline+" && "+plotcuts[iplot].bincuts[ibin]+" && "+_abcdcuts;
        if (debug) cout<<totcut<<endl;
        table_cuts.push_back(TableRow("", totcut.Data(),0,0, wgt));
        allcuts[iplot][iabcd].push_back(totcut);
      } // Loop over bins
      TString tname = "rmt_"+tag+"_"; tname += iplot; tname += iabcd;
      pm.Push<Table>(tname.Data(),  table_cuts, procs, false, false);
    } // Loop over abcdcuts
  } // Loop over plots

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////

  for(size_t iplot=0; iplot<plotcuts.size(); iplot++){
    // allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
    vector<vector<vector<GammaParams> > > allyields(procs.size(), vector<vector<GammaParams> >(Nabcd));
    for(size_t iabcd=0; iabcd<abcdcuts.size(); iabcd++){
      Table * yield_table = static_cast<Table*>(pm.Figures()[iplot*Nabcd+iabcd].get());
      for(size_t ibkg=0; ibkg<procs.size(); ibkg++) {
        allyields[ibkg][iabcd] = yield_table->Yield(procs[ibkg].get(), lumi);
      }
    } // Loop over ABCD cuts

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(allcuts[iplot], allyields, baseline, procs);

    vector<vector<vector<int> > > indices;
    vector<TString> leglabels;
    for(int ibkg=0; ibkg<static_cast<int>(procs.size()); ibkg++) {
      indices.push_back(vector<vector<int> >({{ibkg, r3, 1}, {ibkg, r1, -1}}));
      leglabels.push_back(procs[ibkg]->name_);
    }

    plotRatio(allyields, plotcuts[iplot], indices, leglabels, procs);


  } // Loop over plots


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<plotcuts.size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////




void plotRatio(vector<vector<vector<GammaParams> > > &allyields, oneplot &plotdef,
	       vector<vector<vector<int> > > &indices, vector<TString> &leglabels, 
	       vector<shared_ptr<Process> > &procs){

  size_t ngraphs = indices.size();
  size_t nbins = allyields[0][0].size();

  //// Finding all ratios for all graphs
  float val(1.), valup(1.), valdown(1.);
  vector<vector<vector<float> > > ratios(ngraphs);
  float maxr=-1., minr=1e6;
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    // Finding powers to calculate ratio
    vector<float> powers;
    for(size_t ipow=0; ipow<indices[igraph].size(); ipow++) powers.push_back(indices[igraph][ipow][2]);

    // Finding ratios for each bin
    for(size_t ibin=0; ibin<nbins; ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      for(size_t ind=0; ind<indices[igraph].size(); ind++) {
        size_t ibkg = indices[igraph][ind][0];
        size_t iabcd = indices[igraph][ind][1];
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[ibkg][iabcd][ibin].NEffective());
        weights.back().push_back(allyields[ibkg][iabcd][ibin].Weight());
      } // Loop over indices

      // Throwing toys to find ratios and uncertainties
      val = calcKappa(entries, weights, powers, valdown, valup);
      if(valdown<0) valdown = 0;
      ratios[igraph].push_back(vector<float>({val, valdown, valup}));
      if(maxr < val+valup) maxr = val+valup;
      if(minr > val-valdown) minr = val-valdown;
    } // Loop over bins
  } // Loop over graphs

  //// Finding ytitle
  TString ytitle="Ratio";
  if(indices[0].size()==2){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    if((ind0==r3&&ind1==r1) || (ind0==r4&&ind1==r2)) ytitle = "R(m_{T})";
    if((ind0==r4&&ind1==r3) || (ind0==r2&&ind1==r1)) ytitle = "R(M_{J})";
  }
  if(indices[0].size()==4){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    size_t ind2=indices[0][2][1], ind3=indices[0][3][1];
    if((ind0==r4&&ind1==r3&&ind2==r2&&ind3==r1)) 
      ytitle = "R(M_{J}^{high}) / R[M_{J}^{low}("+TString(procs[0]->name_)+")]";
  }
  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Ratio");
  setPlotStyle(opts);

  //// Plotting kappas
  TCanvas can("can","",1400,500);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);

  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = 0.26;
  // if(maxy<0.26) maxy = 0.26;
  // if(maxy>6) maxy = 6;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetYaxis()->SetTitleOffset(0.9);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.GetXaxis()->SetLabelSize(0.07);
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  vector<vector<double> > vx(ngraphs), vexh(ngraphs), vexl(ngraphs);
  vector<vector<double> > vy(ngraphs), veyh(ngraphs), veyl(ngraphs);
  for(size_t ibin=0; ibin<nbins; ibin++){
    TString tmp = plotdef.bincuts[ibin];
    tmp.ReplaceAll(" ","");
    tmp.ReplaceAll("&&mj14<=400","").ReplaceAll("&&mj14>400&&mj14<=500","").ReplaceAll("&&mj14>500","");
    tmp.ReplaceAll("&&mj14<=450","").ReplaceAll("&&mj14>450&&mj14<=650","").ReplaceAll("&&mj14>650","");
    tmp.ReplaceAll("&&mj14<=500","").ReplaceAll("&&mj14>500&&mj14<=800","").ReplaceAll("&&mj14>800","");
    histo.GetXaxis()->SetBinLabel(ibin+1, CodeToRootTex(tmp.Data()).c_str());
    // xval is the x position of the first marker in the group
    double xval = ibin+1, minxb = 0.15, binw = 0;
    // If there is more than one point in the group, it starts minxb to the left of the center of the bin
    // binw is the distance between points in the njets group
    if(ngraphs>1) {
      xval -= minxb;
      binw = 2*minxb/(ngraphs-1);
    }
    for(size_t igraph=0; igraph<ngraphs; igraph++){
      vx[igraph].push_back(xval);
      xval += binw;
      vexl[igraph].push_back(0);
      vexh[igraph].push_back(0);
      vy[igraph]  .push_back(ratios[igraph][ibin][0]);
      veyl[igraph].push_back(ratios[igraph][ibin][1]);
      veyh[igraph].push_back(ratios[igraph][ibin][2]);
    } // Loop over TGraphs
  } // Loop over bin cuts

  //// Drawing legend and TGraphs
  double legX(opts.LeftMargin()+0.023), legY(1-opts.TopMargin()-0.03), legSingle = 0.05;
  double legW = 0.19*4, legH = legSingle;
  int Ncol = ngraphs;
  if(ngraphs>3) {
    legH *= 2;
    Ncol = (ngraphs+1)/2;
    legW = 0.25*Ncol;
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(Ncol);

  Palette colors("txt/colors.txt", "default");
  vector<int> mcolors, styles;
  for(size_t ibkg=0; ibkg<procs.size(); ibkg++) {
    mcolors.push_back(procs[ibkg]->color_);
    styles.push_back(19+ibkg);
  }
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    graph[igraph] = TGraphAsymmErrors(vx[igraph].size(), &(vx[igraph][0]), &(vy[igraph][0]),
                                    &(vexl[igraph][0]), &(vexh[igraph][0]), &(veyl[igraph][0]), &(veyh[igraph][0]));
    graph[igraph].SetMarkerStyle(styles[igraph]); 
    if(leglabels[igraph].Contains("All")) graph[igraph].SetMarkerStyle(21);
    graph[igraph].SetMarkerSize(1.4);
    graph[igraph].SetMarkerColor(mcolors[igraph]);
    graph[igraph].SetLineColor(mcolors[igraph]); graph[igraph].SetLineWidth(2);
    graph[igraph].Draw("p0 same");
    leg.AddEntry(&graph[igraph], leglabels[igraph], "p");
  } // Loop over TGraphs
  leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(4.5, miny, 4.5, maxy);
  line.DrawLine(8.5, miny, 8.5, maxy);

  line.SetLineColor(kMagenta+1); line.SetLineStyle(2); line.SetLineWidth(1);
  line.DrawLine(2.5, miny, 2.5, maxy-0.1);
  line.DrawLine(6.5, miny, 6.5, maxy-0.1);
  line.DrawLine(10.5, miny, 10.5, maxy-0.1);

  TLatex mjlabel;
  mjlabel.SetTextSize(0.06);
  if (plotdef.baseline.Contains("met>200")) {
    mjlabel.DrawLatex(1.65, maxy-0.07, "250 #leq M#lower[-0.1]{_{J}} #leq 400");
    mjlabel.DrawLatex(5.65, maxy-0.07, "400 #leq M#lower[-0.1]{_{J}} #leq 500");
    mjlabel.DrawLatex(10, maxy-0.07, "M#lower[-0.1]{_{J}} #geq 500");
  } else if (plotdef.baseline.Contains("met>350")) {
    mjlabel.DrawLatex(1.65, maxy-0.07, "250 #leq M#lower[-0.1]{_{J}} #leq 450");
    mjlabel.DrawLatex(5.65, maxy-0.07, "450 #leq M#lower[-0.1]{_{J}} #leq 650");
    mjlabel.DrawLatex(10, maxy-0.07, "M#lower[-0.1]{_{J}} #geq 650");
  } else {
    mjlabel.DrawLatex(1.65, maxy-0.07, "250 #leq M#lower[-0.1]{_{J}} #leq 500");
    mjlabel.DrawLatex(5.65, maxy-0.07, "500 #leq M#lower[-0.1]{_{J}} #leq 800");
    mjlabel.DrawLatex(10, maxy-0.07, "M#lower[-0.1]{_{J}} #geq 800");
  }

  TString fname = "plots/ratio_"+CodeToPlainText(ytitle.Data())+"_"+plotdef.name+"_"
    +CodeToPlainText(plotdef.baseline.Data())+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

} // plotRatio



// allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline, vector<shared_ptr<Process> > &procs){
  int digits = 3;
  cout<<endl<<endl<<"============================ Printing cuts  ============================"<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl<<endl;
  for(size_t ibin=0; ibin<allcuts[0].size(); ibin++){
    for(size_t iabcd=0; iabcd<allcuts.size(); iabcd++){
      for(size_t ibkg=0; procs.size(); ibkg++)
        cout<<procs[ibkg]->name_<<": "    <<setw(9)<<RoundNumber(allyields[ibkg][iabcd][ibin].Yield(), digits)<<", ";
      cout<<"  - "<< allcuts[iabcd][ibin]<<endl;
    } // Loop over ABCD cuts
    cout<<endl;
  } // Loop over bin cuts

} // printDebug


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"debug", no_argument, 0, 'd'},
      {"year", required_argument, 0, 'y'},
      {"tag", required_argument, 0, 't'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "y:t:d", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'y':
      year = atoi(optarg);
      break;
    case 't':
      tag = optarg;
      break;
    case 'd':
      debug = true;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
