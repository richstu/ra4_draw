//// plot_t5tttt_prl: Plots T1tttt and T5tttt limits in one canvas

// System includes
#include <fstream>
#include <iostream>

// ROOT includes
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TPad.h"
#include "TLine.h"
#include "TH2D.h"
#include "TColor.h"
#include "TPaletteAxis.h"
#include "TPolyLine.h"
#include "TError.h" // Controls error level reporting

// User includes
#include "ra4/plot_t5tttt_prl.hpp"
#include "core/utilities.hpp"

using namespace std;
namespace{
  //double cmsH = 0.075;
  bool do_tev = false;
  double cmsH = 0.03;
  float legLineH = 0.055;
  float legTextSize = 0.04;
  float fillTransparency = 0.5;

  TString lsp = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString stop = "#tilde{t}_{1}";
  TString mlsp_s = "m("+lsp+")";
  TString mglu_s = "m(#tilde{g})";
  TString mstop_s = "m("+stop+")";
  TString t1tttt_s = "#tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #kern[-0.2]{#rightarrow} #kern[-0.2]{t}#kern[0.4]{#bar{t}}#kern[0.4]{"+lsp+"}";
  TString t5tttt_s = "#tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #kern[-0.2]{#rightarrow} #kern[-0.2]{"+stop+"}#bar{t},  "+stop+" #rightarrow t#kern[0.4]{"+lsp+"},  "+mstop_s+" #kern[-0.5]{-} #kern[-0.1]{"+mlsp_s+"} = #kern[-0.1]{175} GeV";
}


void SetupColors();

int main(){
  gErrorIgnoreLevel=kWarning; // Turns off ROOT INFO messages

  // Label definitions
  //TString lsp("#tilde{#chi}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}");
  TString pp_gluglu("pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}");
  TString mj("M#lower[-.1]{_{J}}");
  //TString mt2("M#lower[-.1]{_{T2}}"), mht("#slash{H}#lower[-.1]{_{T}}"), aT("#alpha#lower[-.1]{_{T}}");

  // Folder with root files containing the TGraphs
  TString folder("");
  vector<model_limits> models;

  ///////////////////////////////    Defining T1tttt plot    /////////////////////////////////
  models.push_back(model_limits("T1tttt", pp_gluglu));
  models.back().add(t1tttt_s, folder+"T1tttt_limit_scan_nom.root", 1, "T1ttttObservedLimit", "T1ttttExpectedLimit");
  models.push_back(model_limits("T5tttt", pp_gluglu));
  models.back().add(t5tttt_s, folder+"T5tttt_limit_scan_nom.root", kAzure+7, "T5ttttObservedLimit", "T5ttttExpectedLimit");
  models.back().add(t1tttt_s, folder+"T1tttt_limit_scan_nom.root", 1, "T1ttttObservedLimit", "T1ttttExpectedLimit");

  //////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////    Making plots    //////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////// 
  
  // Creating canvas
  gStyle->SetOptStat(0);  
  SetupColors();

  // Creating base legend for observed/expected
  int wobs(4), wexp(4);
  TH1D hobs("hobs","",1,0,1), hexp("hexp","",1,0,1);
  hobs.SetLineColor(1); hobs.SetLineWidth(wobs);
  hexp.SetLineColor(1); hexp.SetLineStyle(2); hexp.SetLineWidth(wexp);


  // Looping over each model
  cout<<endl;
  for(size_t imodel(0); imodel < models.size(); imodel++){
    model_limits mod(models[imodel]);

    float lMargin(0.13), tMargin(0.08), rMargin(0.16), bMargin(0.11);
    int canW = 800, canH = 700;
    TCanvas can("canvas","", canW, canH);
    if (mod.model == "T5tttt") {
      rMargin = 0.05;
    }
    setCanvas(can, lMargin, tMargin, rMargin, bMargin);

    double legX(1-rMargin-0.04), legY(1-tMargin-0.13);
    double legW = 0.8, legH = 0.07;
    TLegend baseleg(legX-legW, legY-legH, legX, legY);
    baseleg.SetTextSize(0.034); baseleg.SetFillColor(0); 
    baseleg.SetFillStyle(0); baseleg.SetBorderSize(0);
    //baseleg.AddEntry(&hobs, "Observed");
    baseleg.AddEntry(&hexp, "Expected");
    legX = 0.75;
    TLegend obsleg(legX-legW, legY-legH, legX, legY);
    obsleg.SetTextSize(0.034); obsleg.SetFillColor(0); 
    obsleg.SetFillStyle(0); obsleg.SetBorderSize(0);
    obsleg.AddEntry(&hobs, "Observed");

   // Creating base histogram and drawing lumi labels
    float Xmin(700), Xmax(1750), Ymin(0), Ymax(1800), glu_lsp;
    getModelParams(mod.model, Xmin, Xmax, Ymin, Ymax, glu_lsp);

    TH2D hbase = baseHistogram(Xmin, Xmax, Ymin, Ymax);
    hbase.GetYaxis()->SetTitleOffset(1.45);
    hbase.GetYaxis()->SetTitleSize(0.045);
    hbase.GetYaxis()->SetLabelSize(0.039);

    hbase.GetXaxis()->SetTitleOffset(1.1);
    hbase.GetXaxis()->SetTitleSize(0.045);
    hbase.GetXaxis()->SetLabelSize(0.039);
    hbase.Draw();
    TH2D *hxsec_ori = nullptr;

    // Plotting limits
    int widthCentral = 4, widthErr = 2;
    int styleObs = 1, styleObsErr = 3, styleExp = 7;
    int colorExp = 2;
    size_t ncurves(mod.files.size());
    vector<TGraph*> obs(ncurves, 0), exp(ncurves, 0);
    vector<TGraph*> obsUp(ncurves,0), obsDown(ncurves,0), expUp(ncurves,0), expDown(ncurves,0);
    vector<TGraph*> expArea(ncurves,0);
    // Getting all graphs first because the ones that come from TCanvas mess up the colors
    vector<TFile*> flimit(ncurves);
    for(size_t file(0); file < ncurves; file++){
      flimit[file] = new TFile(mod.files[file]);
      exp[file] = getGraph(*flimit[file], mod.expnames[file]);
      obs[file] = getGraph(*flimit[file], mod.obsnames[file]);
      expUp[file] = getGraph(*flimit[file], mod.expnames[file]+"Up");
      expDown[file] = getGraph(*flimit[file], mod.expnames[file]+"Down");
      obsUp[file] = getGraph(*flimit[file], mod.obsnames[file]+"Up");
      obsDown[file] = getGraph(*flimit[file], mod.obsnames[file]+"Down");
      reverseGraph(expDown[file]);
      expArea[file] = joinGraphs(expUp[file], expDown[file]);
      if(file==0 && mod.model=="T1tttt") {
         hxsec_ori = getHist2D(*flimit[file], mod.model+"ObservedExcludedXsec");
         hxsec_ori->SetDirectory(0);
      }
    }

    int yext = 600.;
    if (mod.model=="T1tttt") yext = 400.;

    if (mod.model=="T1tttt") {
      TH2D dummy("","",hxsec_ori->GetNbinsX(),
               hxsec_ori->GetXaxis()->GetXmin(),hxsec_ori->GetXaxis()->GetXmax(), 
               hxsec_ori->GetNbinsY()+yext/hxsec_ori->GetYaxis()->GetBinWidth(1),
               hxsec_ori->GetYaxis()->GetXmin(),hxsec_ori->GetYaxis()->GetXmax()+yext);
      for (int ix(0); ix<hxsec_ori->GetNbinsX(); ix++) {
        for (int iy(0); iy<hxsec_ori->GetNbinsY(); iy++) {
          dummy.SetBinContent(ix+1,iy+1,hxsec_ori->GetBinContent(ix+1, iy+1)*1000.); //switch to fb
          dummy.SetBinError(ix+1,iy+1,hxsec_ori->GetBinError(ix+1, iy+1)*1000.);
          dummy.GetXaxis()->SetTitle(mglu_s + " [GeV]");
          dummy.GetYaxis()->SetTitle(mlsp_s + " [GeV]");
        }
      }

      dummy.SetMinimum(1e-1);
      dummy.GetYaxis()->SetTitleOffset(1.45);
      dummy.GetYaxis()->SetTitleSize(0.045);
      dummy.GetYaxis()->SetLabelSize(0.039);

      dummy.GetXaxis()->SetTitleOffset(1.05);
      dummy.GetXaxis()->SetTitleSize(0.045);
      dummy.GetXaxis()->SetLabelSize(0.039);

      dummy.GetZaxis()->SetLabelSize(0.04);
      dummy.GetZaxis()->SetTitleSize(0.05);
      dummy.GetZaxis()->SetTitleOffset(1);
      dummy.GetZaxis()->SetTitle("Upper limit (95% CL) on #sigma("+t1tttt_s+") [fb]");
      dummy.DrawCopy("colz");
    }
    addLabelsTitle(lMargin, tMargin, rMargin, mod.title);

    gPad->Modified();
    gPad->Update();
    if (mod.model=="T1tttt") {
      TPaletteAxis *palette = static_cast<TPaletteAxis*>(hxsec_ori->GetListOfFunctions()->FindObject("palette"));
      palette->SetX1NDC(1.-rMargin+0.006);
      palette->SetX2NDC(1.-rMargin+0.035);
      palette->SetY1NDC(bMargin);
      palette->SetY2NDC(1.-tMargin);
    }
    for(size_t file(0); file < ncurves; file++){
      TString model_name = "T1tttt";
      if(mod.labels[file].Contains("175") || mod.labels[file].Contains("T5tttt")) {
        model_name = "T5tttt";
        glu_lsp += 40;
      }
      setGraphStyle(obs[file], mod.colors[file], styleObs, widthCentral, glu_lsp, model_name);
      setGraphStyle(obsUp[file], mod.colors[file], styleObsErr, widthErr, glu_lsp, model_name);
      setGraphStyle(obsDown[file], mod.colors[file], styleObsErr, widthErr, glu_lsp, model_name);

      if(mod.labels[file].Contains("T2tt")) colorExp = kOrange+1;
      if(mod.model=="T5tttt") colorExp = mod.colors[file];
      if(mod.labels[file].Contains("175") && imodel==0) setGraphStyle(exp[file], mod.colors[file], styleExp, widthCentral-1, glu_lsp, model_name);
      else setGraphStyle(exp[file], colorExp, styleExp, widthCentral, glu_lsp, model_name);
      setGraphStyle(expUp[file], colorExp, styleExp, widthErr, glu_lsp, model_name);
      setGraphStyle(expDown[file], colorExp, styleExp, widthErr, glu_lsp, model_name);

      TString obsname("obs"); obsname += imodel; obsname += file;
      obs[file]->SetName(obsname);
      TString expAreaname("expArea"); expAreaname += imodel; expAreaname += file;
      expArea[file]->SetName(expAreaname);
      TString expname("exp"); expname += imodel; expname += file;
      exp[file]->SetName(expname);

      // Setting the area style for expected limit
      int areaColor = kOrange;
      expArea[file]->SetFillColor(areaColor);
      expArea[file]->SetFillColorAlpha(areaColor, fillTransparency);
      expArea[file]->SetFillStyle(1001);
      expArea[file]->SetLineColor(colorExp);
      expArea[file]->SetLineStyle(styleExp);
      expArea[file]->SetLineWidth(widthCentral);
    } // Loop over curves in each model

    hbase.Draw("axis same");

    // Drawing legends
    float legEntries = 4;
    if (mod.model=="T1tttt") legEntries = 2.5;
    legX = mod.model=="T1tttt" ? lMargin+0.01 : lMargin+0.04; 
    legY = mod.model=="T1tttt" ? 1-tMargin-0.01: 1-tMargin-0.04;
    legW = mod.model=="T1tttt" ? 0.23 : 0.5; 
    legH = legLineH * legEntries;
    TLegend* limleg[2];

    if (mod.model=="T1tttt") {
      limleg[0] = new TLegend(legX, legY-legH, legX+legW, legY);
      limleg[0]->SetX1NDC(legX); limleg[0]->SetX2NDC(legX+legW); // So that GetX1NDC works in getLegendBoxes
      limleg[0]->SetY1NDC(legY-legH); limleg[0]->SetY2NDC(legY); // So that GetX1NDC works in getLegendBoxes
      limleg[0]->SetTextSize(legTextSize); limleg[0]->SetFillColor(0); 
      limleg[0]->SetFillStyle(0); limleg[0]->SetBorderSize(0);
    } else {
      limleg[0] = new TLegend(legX, legY-legH/2, legX+legW, legY);
      limleg[1] = new TLegend(legX, legY-legH, legX+legW, legY-legH/2);
      for (unsigned ileg(0); ileg<2; ileg++) {
        limleg[ileg]->SetNColumns(2);
        limleg[ileg]->SetTextSize(legTextSize); limleg[ileg]->SetFillColor(0); 
        limleg[ileg]->SetFillStyle(0); limleg[ileg]->SetBorderSize(0);
      }
    }

    Ymax = Ymax+yext-400;
    if (mod.model=="T1tttt") {
      float bheight = (Ymax-Ymin)*legH/(1-tMargin-bMargin)*1.17;
      TBox box;
      Xmin += (Xmax-Xmin)*0.001;
      Xmax -= (Xmax-Xmin)*0.001;
      box.SetFillColor(0); box.SetFillStyle(1001);
      box.SetLineColor(1); box.SetLineWidth(2); box.SetLineStyle(1);
      box.DrawBox(Xmin, Ymax-bheight, Xmax, Ymax);
      box.SetFillColor(0); box.SetFillStyle(0);
      box.SetLineColor(1); box.SetLineWidth(2); box.SetLineStyle(1);
      box.DrawBox(Xmin, Ymax-bheight, Xmax, Ymax);
    }
    
    // Plotting the lines on top of the fills
    for(size_t file(0); file < ncurves; file++){
      if(mod.model!="T5tttt" && (!mod.labels[file].Contains("175") || imodel!=0)) {
	// expArea[file]->Draw("f same");
       expUp[file]->Draw("same");
       expDown[file]->Draw("same");
       obsUp[file]->Draw("same");
       obsDown[file]->Draw("same");
      }
      exp[file]->Draw("same");
      obs[file]->Draw("same");
      obs[0]->Draw("same");
    }// Loop over curves in each model

    //limleg.AddEntry(expArea[1]->GetName(), "95% CL upper limits", "n");
    
    for(size_t file(0); file < ncurves; file++){
      if(!mod.labels[file].Contains("175") && mod.model!="T5tttt") {
      	limleg[0]->AddEntry(exp[file]->GetName(), "Expected ("+mod.labels[file]+") #pm s.d._{experiment}", "l");
      	limleg[0]->AddEntry(obs[file]->GetName(), "Observed ("+mod.labels[file]+") #pm s.d._{theory}", "l");
      } else {
        limleg[file]->SetHeader("Model: "+ mod.labels[file]);
        limleg[file]->AddEntry(exp[file]->GetName(), "Expected", "l");
        limleg[file]->AddEntry(obs[file]->GetName(), "Observed", "l");
      }
    }
    limleg[0]->Draw();
    if(mod.model!="T1tttt") limleg[1]->Draw();


    vector<vector<float> > boxes;
    getLegendBoxes(*limleg[0], boxes);
    int ibox = 0;
    TLine line; TLatex label; 
    // T2tt exclusion
    float x[5] = {800, 2600, 2600, 800};
    float y[5] = {100, 100, 550, 550};
    TPolyLine pline = TPolyLine(4,x,y);
    float x2[5] = {800, 2600, 2600, 800};
    float y2[5] = {0, 0, 33, 33};
    TPolyLine pline2 = TPolyLine(4,x2,y2);

    if (mod.model!="T5tttt") {
      // Drawing expected error lines on legend
      ibox = 0;
      line.SetLineColor(colorExp);line.SetLineWidth(widthErr);line.SetLineStyle(styleExp);
      line.DrawLineNDC(boxes[ibox][0], boxes[ibox][1], boxes[ibox][2], boxes[ibox][1]);
      line.DrawLineNDC(boxes[ibox][0], boxes[ibox][3], boxes[ibox][2], boxes[ibox][3]);

      // Drawing observed error lines on legend
      line.SetLineColor(mod.colors[0]);line.SetLineWidth(widthErr);line.SetLineStyle(styleObsErr);
      ibox = 1;
      line.DrawLineNDC(boxes[ibox][0], boxes[ibox][1], boxes[ibox][2], boxes[ibox][1]);
      line.DrawLineNDC(boxes[ibox][0], boxes[ibox][3], boxes[ibox][2], boxes[ibox][3]);
    } else {


      pline.SetFillColor(kGray+1);
      pline.SetFillStyle(3353);
      pline.SetFillColorAlpha(kGray+1, 0.6);
      pline.SetLineColor(kGray+1);
      pline.SetLineWidth(1);
      pline.Draw("f");
      pline.Draw();

      pline2.SetFillColor(kGray+1);
      pline2.SetFillStyle(3353);
      pline2.SetFillColorAlpha(kGray+1, 0.6);
      pline2.SetLineColor(kGray+1);
      pline2.SetLineWidth(1);
      pline2.Draw("f");
      pline2.Draw();

      // label.SetNDC();  
      label.SetTextAlign(21); label.SetTextFont(72); label.SetTextColor(kGray+2); label.SetTextSize(0.043);
      label.DrawLatex(1400, 400, "Direct stop pair production");
      label.DrawLatex(1400, 250, "excluded");
    }

    legY = 1-legY-legH-0.02-0.1; legH = 0.07;
    obsleg.SetY1NDC(legY-legH); obsleg.SetY2NDC(legY);
    baseleg.SetY1NDC(legY-legH); baseleg.SetY2NDC(legY);
    //baseleg.Draw();
    //obsleg.Draw();

    TString plotname("plots/"+mod.model+"_limit.pdf");
    can.SaveAs(plotname);
    cout<<" open "<<plotname<<endl;
  } // Loop over models
  cout<<endl<<endl;
}

TH2D* getHist2D(TFile &flimit, TString hname){
  TH2D *hist = static_cast<TH2D*>(flimit.Get(hname));
  // If the TH2D is not directly provided in the root file, try to extract it from a TCanvas
  if(hist==nullptr) {
    TPad *current_pad = static_cast<TPad*>(gPad);
    TCanvas *c1 = static_cast<TCanvas*>(flimit.Get("c1"));
    current_pad->cd();
    if(c1!=nullptr){
      hist = static_cast<TH2D*>(c1->GetListOfPrimitives()->FindObject(hname));
    }
  }
  return hist;
}

TGraph* getGraph(TFile &flimit, TString gname){
  TGraph *graph = static_cast<TGraph*>(flimit.Get(gname));
  // If the TGraph is not directly provided in the root file, try to extract it from a TCanvas
  if(graph==nullptr) {
    TPad *current_pad = static_cast<TPad*>(gPad);
    TCanvas *c1 = static_cast<TCanvas*>(flimit.Get("c1"));
    current_pad->cd();
    if(c1!=nullptr){
      graph = static_cast<TGraph*>(c1->GetListOfPrimitives()->FindObject(gname));
    }
  }
  return graph;
}

void setGraphStyle(TGraph* graph, int color, int style, int width, double glu_lsp, TString model_name){
  if(graph==0) return;
  // Setting graph style
  graph->SetLineColor(color);
  graph->SetLineStyle(style);
  int fillcolor(color);
  // graph->SetFillColor(fillcolor);
  // graph->SetFillColorAlpha(fillcolor, fillTransparency);
  // graph->SetFillStyle(1001);
  graph->SetFillColor(0);
  graph->SetFillColorAlpha(fillcolor, fillTransparency);
  graph->SetFillStyle(0);
  graph->SetLineWidth(width); 

  int np(graph->GetN());
  double mglu, iniglu, endglu, mlsp;
  if(do_tev){
    for(int point(0); point < np; point++){
      graph->GetPoint(point, mglu, mlsp);
      graph->SetPoint(point, mglu/1000., mlsp/1000.);
    }
  }

  graph->GetPoint(0, iniglu, mlsp);
  graph->GetPoint(np-1, endglu, mlsp);
  // Reversing graph if printed towards decreasing mgluino
  if(iniglu > endglu) reverseGraph(graph);
  // Adding a point so that it goes down to mLSP = 0
  // graph->SetPoint(graph->GetN(), max(iniglu,endglu), 0);
  // np++;

  reverseGraph(graph);

  // Adding a point at LSP = 0, and removing points beyond the diagonal
  for(int point(0); point < np; point++){
    graph->GetPoint(point, mglu, mlsp);
    if(mlsp > mglu-glu_lsp){
      while(point <= graph->GetN()) {
       graph->RemovePoint(graph->GetN()-1);
       np--;
     }
      break;
    }
  }

  if (model_name=="T5tttt" || model_name=="T2tt") {
    graph->RemovePoint(graph->GetN()-1);
    np--;
  }
  // Finding intersection of line between last 2 points and mlsp = mglu - glu_lsp
  double x1, y1, x2, y2;
  graph->GetPoint(np-1, x1, y1);
  graph->GetPoint(np-2, x2, y2);
  double slope((y1-y2)/(x1-x2)), offset(y1-slope*x1);
  double intersection((offset+glu_lsp)/(1-slope));
  //cout<<"intersection "<<intersection<<", slope "<<slope<<", glu_lsp "<<glu_lsp<<endl;

  // Adding extrapolation into the diagonal, and point for mglu = 0
  if(slope!=1) graph->SetPoint(graph->GetN(), intersection, intersection-glu_lsp);
  graph->SetPoint(graph->GetN(), 0, -glu_lsp);
  if(x1 == x2 || y1 == y2 || slope == 1){
    for(int point(0); point < graph->GetN(); point++){
      graph->GetPoint(point, mglu, mlsp);
      //cout<<point<<": "<<mglu<<", "<<mlsp<<endl;
    }
  }
}

TGraph* joinGraphs(TGraph *graph1, TGraph *graph2){
  TGraph *graph = new TGraph;
  double mglu, mlsp;
  for(int point(0); point < graph1->GetN(); point++) {
    graph1->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  for(int point(0); point < graph2->GetN(); point++) {
    graph2->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  graph1->GetPoint(0, mglu, mlsp);
  graph->SetPoint(graph->GetN(), mglu, mlsp);
  TString gname = graph1->GetName(); gname += graph2->GetName();
  graph->SetName(gname);

  return graph;
}

void reverseGraph(TGraph *graph){
  int np(graph->GetN());
  double mglu, mlsp;
  vector<double> mglus, mlsps;
  for(int point(np-1); point >= 0; point--){
    graph->GetPoint(point, mglu, mlsp);
    mglus.push_back(mglu);
    mlsps.push_back(mlsp);
  }
  for(int point(0); point < np; point++)
    graph->SetPoint(point, mglus[point], mlsps[point]);
}

void getModelParams(TString model, float &Xmin, float &Xmax, float &Ymin, float &Ymax, float &glu_lsp){
  if(model == "T1tttt" || model == "T5tttt" || model == "T2tt"){
    Xmin = 800; Xmax = 2600.;
    Ymin = 0;   Ymax = 2000;
    glu_lsp = 225;
  }
  if(model == "T1bbbb"){
    Xmin = 600; Xmax = 1950;
    Ymin = 0;   Ymax = 1885;
    glu_lsp = 25;
  }    
  if(model == "T1qqqq"){
    Xmin = 600; Xmax = 1950;
    Ymin = 0;   Ymax = 1750;
    glu_lsp = 25;
  }    
  if(do_tev){
    Xmin *= 0.001;
    Xmax *= 0.001;
    Ymin *= 0.001;
    Ymax *= 0.001;
    glu_lsp *= 0.001;
  }
}


void setCanvas(TCanvas &can, float lMargin, float tMargin, float rMargin, float bMargin){
  can.SetLogz();
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(lMargin);
  can.SetTopMargin(tMargin);
  can.SetRightMargin(rMargin);
  can.SetBottomMargin(bMargin);
}

void addLabelsTitle(float lMargin, float tMargin, float rMargin, TString title){
  TLatex label; label.SetNDC();  
  // Printing CMS Preliminary
  double offsetx(0.0), ycms(1-tMargin-cmsH);
  label.SetTextAlign(12); label.SetTextFont(61); label.SetTextSize(0.75*tMargin);
  label.DrawLatex(lMargin+offsetx, 1-tMargin/2., "CMS");
  label.SetTextAlign(12); label.SetTextFont(52); label.SetTextSize(0.038);
  // label.DrawLatex(lMargin+offsetx+0.1, 1-tMargin/2.-0.013, "Preliminary");
  // Printing lumi (energy)
  label.SetTextAlign(31); label.SetTextFont(42); label.SetTextSize(0.6*tMargin);
  label.DrawLatex(1-rMargin, 1-tMargin+0.018, "137 fb^{-1} (13 TeV)");
  
  title += " "; ycms += 1;// To avoid non-used warnings
}

void SetupColors(){
  const unsigned num = 4;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.00, 0.5, 0.8, 1.00};
  double red[num] = {0.50, 1.00, 1.00, 1.00};
  double green[num] = {1.00, 1.00, 0.60, 0.50};
  double blue[num] = {1.00, 0.50, 0.40, 0.50};

  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    colors[i] = fi+i;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}


TH2D baseHistogram(float Xmin, float Xmax, float Ymin, float Ymax){
  TH2D hbase("hbase", "", 1, Xmin, Xmax, 1, Ymin, Ymax);
  hbase.GetXaxis()->SetLabelFont(42);
  hbase.GetXaxis()->SetLabelSize(0.045);
  hbase.GetXaxis()->SetTitleFont(42);
  hbase.GetXaxis()->SetTitleSize(0.05);
  hbase.GetXaxis()->SetTitleOffset(1.2);
  hbase.GetXaxis()->SetLabelOffset(0.015);
  hbase.GetXaxis()->SetTitle(mglu_s+" [GeV]");

  hbase.GetYaxis()->SetLabelFont(42);
  hbase.GetYaxis()->SetLabelSize(0.045);
  hbase.GetYaxis()->SetTitleFont(42);
  hbase.GetYaxis()->SetTitleSize(0.05);
  hbase.GetYaxis()->SetTitleOffset(0.9);
  hbase.GetYaxis()->SetTitle(mlsp_s+" [GeV]");

  return hbase;
}

void model_limits::add(TString label, TString file, int color, TString obsname, TString expname){
  labels.push_back(label);
  files.push_back(file);
  obsnames.push_back(obsname);
  expnames.push_back(expname);
  colors.push_back(color);
}

model_limits::model_limits(TString imodel, TString ititle, float ilegScale):
  model(imodel),
  title(ititle),
  legScale(ilegScale){
  }
