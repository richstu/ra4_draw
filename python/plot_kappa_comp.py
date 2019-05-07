#! /usr/bin/env python
import os, sys
import math
from ROOT import *

yrs = ['0','2016','2017','2018']
mjs = ['lowmj','highmj']
nbs = ['1b','nbdm2','nbdmge3']

# yrs = ['2016','2017']
# mjs = ['lowmj']
# nbs = ['1b']

# gROOT.Reset()
#         Bin names
#-------------------------------------
bins = []
for imet in ['lmet','mmet','hmet']:
    for inj in ['lnj','hnj']:
            bins.append('_'.join([imet,inj]))

nbins = len(bins)

#         Reading the data
#-------------------------------------
# Idx 0: year, 1: mj, 2: nb, 3, 
kap_x = [[[[0. for i in range(nbins)] for i in range(3)] for j in range(2)] for k in range(4)]
kap_y = [[[[0. for i in range(nbins)] for i in range(3)] for j in range(2)] for k in range(4)]
ekup = [[[[0. for i in range(nbins)] for i in range(3)] for j in range(2)] for k in range(4)]
ekdown = [[[[0. for i in range(nbins)] for i in range(3)] for j in range(2)] for k in range(4)]

for iyr,yr in enumerate(yrs):
    for imj,mj in enumerate(mjs):
        file = TFile("root/kappa_signal_"+mj+"_abcd_scen_data_"+yr+".root")
        for inb,nb in enumerate(nbs):
            # gr_in = TGraphAsymmErrors()
            gr_in = file.Get("can").GetListOfPrimitives().FindObject(nb)
            for ibin in range(nbins):
                x, y  = Double(0), Double(0)
                gr_in.GetPoint(ibin, x, y)
                kap_x[iyr][imj][inb][ibin] = x
                kap_y[iyr][imj][inb][ibin] = y
                ekup[iyr][imj][inb][ibin] = gr_in.GetErrorYhigh(ibin)
                ekdown[iyr][imj][inb][ibin] = gr_in.GetErrorYlow(ibin)
        file.Close()

file = TFile("root/kappa_signal_lowmj_abcd_scen_data_0.root")
histo = TH1D()
histo = file.Get("can").GetListOfPrimitives().FindObject("histo")
histo.SetDirectory(0)
file.Close()

#        Draw graphs
#-------------------------------------
gr = [[[TGraphAsymmErrors(nbins) for i in range(3)] for j in range(2)] for k in range(4)]

leg = TLegend(0.4, 0.9, 0.7, 0.98)
# leg.SetTextSize(30)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetNColumns(2)

color = [kRed,kBlue,kGreen+1]
for imj,mj in enumerate(mjs):
    can = TCanvas('c'+mj,'c'+mj,1200,500)
    can.cd()
    if mj=='highmj': histo.GetYaxis().SetTitle("#kappa_{B}")
    histo.DrawCopy()

    leg = TLegend(0.4, 0.9, 0.7, 0.98)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetNColumns(3)

    for iyr,yr in enumerate(yrs):
        for inb,nb in enumerate(nbs):
            gr[iyr][imj][inb].SetName(yr+mj+nb)            
            for ibin in range(nbins):
                offset_yr = 0
                if (yr=='2016'): offset_yr = -0.05
                elif (yr=='2018'): offset_yr = 0.05
                offset_nb = 0
                if inb==0: offset_nb = -0.2
                elif inb==2: offset_nb = 0.2
                gr[iyr][imj][inb].SetPoint(ibin, kap_x[iyr][imj][inb][ibin]+offset_yr+offset_nb, kap_y[iyr][imj][inb][ibin])
                gr[iyr][imj][inb].SetPointEYhigh(ibin, ekup[iyr][imj][inb][ibin])
                gr[iyr][imj][inb].SetPointEYlow(ibin, ekdown[iyr][imj][inb][ibin])
                xerr = 0.
                if (yr=='0'): xerr = 0.15
                gr[iyr][imj][inb].SetPointEXhigh(ibin, xerr)
                gr[iyr][imj][inb].SetPointEXlow(ibin, xerr)
            gr[iyr][imj][inb].SetMarkerStyle(20+inb)
            gr[iyr][imj][inb].SetMarkerSize(1.)
            gr[iyr][imj][inb].SetMarkerColor(color[inb])
            gr[iyr][imj][inb].SetLineColor(color[inb])
            gr[iyr][imj][inb].SetLineWidth(1)
            if (yr=='0'): 
                # gr[iyr][imj][inb].SetFillStyle(3002)
                gr[iyr][imj][inb].SetFillColor(color[inb]-10)
                if (inb==2): gr[iyr][imj][inb].SetFillColor(color[inb]-11)
                gr[iyr][imj][inb].Draw("2")
            else: 
                gr[iyr][imj][inb].Draw("p0")
            if iyr==0: 
                leg_label = str(inb+1)+'b'                
                leg.AddEntry(gr[iyr][imj][inb],leg_label,"p")

    leg.Draw()

    binlabel = TLatex()
    binlabel.SetTextAlign(21)
    binlabel.SetTextSize(0.045)
    ptmiss = "p#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.15]{#lower[0.2]{^{miss}}}}";
    binlabel.DrawLatex(1.5, 2.5,"#font[52]{200 < "+ptmiss+"#leq 350}")
    binlabel.DrawLatex(3.5, 2.5,"#font[52]{350 < "+ptmiss+"#leq 500}")
    binlabel.DrawLatex(5.5, 2.5,"#font[52]{"+ptmiss+"#geq 500}")

    a = TLine()
    a.SetLineWidth(1)
    a.SetLineStyle(2)
    a.DrawLine(2.5, 0,2.5,3)
    a.DrawLine(4.5, 0,4.5,3)
    b = TLine()
    b.SetLineStyle(3)
    # b.SetLineColor(kGray+1)
    b.DrawLine(1.5, 0,1.5,2)
    b.DrawLine(3.5, 0,3.5,2)
    b.DrawLine(5.5, 0,5.5,2)

    can.Print('kappa_comp_'+mj+'.pdf')

sys.exit(0)

#         Plotting the data
#-------------------------------------

# for ibin in range(nbins):
#     print bkg_pre[ibin], ebkg_pre[ibin]

can = TCanvas('c','c',1000,500)
can.cd()

gStyle.SetOptStat(0)
top = TPad("top_pad", "top_pad", 0., 0.3, 1., 1.)
bottom = TPad("bottom_pad", "bottom_pad", 0., 0., 1., 0.3)
can.SetMargin(0., 0., 0., 0.);
can.SetFillStyle(4000);

top.SetTopMargin(0.1)
top.SetBottomMargin(0.)
top.SetLeftMargin(0.1)
top.SetRightMargin(0.05)
top.SetFillStyle(4000);
top.Draw()

bottom.SetTopMargin(0.)
bottom.SetBottomMargin(0.3)
bottom.SetLeftMargin(0.1)
bottom.SetRightMargin(0.05)
bottom.SetFillStyle(4000);
bottom.Draw()

#     Top pad
# -----------------------------

top.cd()
top.SetLogy()
miny, maxy = 0.101, 2000
htopdummy.GetYaxis().SetLabelSize(0.06)
htopdummy.GetYaxis().SetRangeUser(miny,maxy)
htopdummy.GetYaxis().SetTitle("Events / Bin")
htopdummy.GetYaxis().CenterTitle()
htopdummy.GetYaxis().SetTitleSize(0.075)
htopdummy.GetYaxis().SetTitleOffset(0.5)
htopdummy.Draw()

leg = TLegend(0.4, 0.9, 0.7, 0.98)
# leg.SetTextSize(30)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetNColumns(2)

grbkg_pre.SetLineColor(tag1_color)
grbkg_pre.SetLineWidth(1)
grbkg_pre.SetFillColor(0)
grbkg_pre.SetFillStyle(0)
grbkg_pre.Draw('2')
leg.AddEntry(grbkg_pre, tag1_lbl.replace("--","-"), "F")

grbkg_post.SetFillColor(tag2_color)
grbkg_post.SetFillStyle(3144)
grbkg_post.SetLineWidth(0)
grbkg_post.Draw('2')
leg.AddEntry(grbkg_post, tag2_lbl.replace("--","-"), "F")

grdata.SetMarkerStyle(20)
grdata.Draw('P')

leg.Draw()

cmslabel = TLatex()
cmslabel.SetTextSize(0.06)
cmslabel.SetNDC(kTRUE)
cmslabel.SetTextAlign(11)
cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS}") # #scale[0.8]{#font[42]{Preliminary}}")
cmslabel.SetTextAlign(31)
cmslabel.DrawLatex(1-top.GetRightMargin()-0.005, 0.92,"#font[42]{137 fb^{-1} (13 TeV)}")

binlabel = TLatex()
binlabel.SetTextSize(0.05)
# binlabel.SetNDC(kTRUE)
binlabel.SetTextAlign(21)
binlabel.DrawLatex(10, 1000,"Low M#lower[-0.1]{_{J}}")
binlabel.DrawLatex(28, 1000,"High M#lower[-0.1]{_{J}}")
binlabel.SetTextSize(0.045)
ptmiss = "p#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.15]{#lower[0.2]{^{miss}}}}";

for i in range(2):
    binlabel.DrawLatex(4+i*18, 350,"#font[52]{200 < "+ptmiss+"#leq 350}")
    binlabel.DrawLatex(10+i*18, 350,"#font[52]{350 < "+ptmiss+"#leq 500}")
    binlabel.DrawLatex(16+i*18, 350,"#font[52]{"+ptmiss+"#geq 500}")

binlabel.SetTextSize(0.045)
for i in range(6):
    binlabel.DrawLatex(2+i*6, 150,"#font[52]{1b}")
    binlabel.DrawLatex(4+i*6, 150,"#font[52]{2b}")
    binlabel.DrawLatex(6+i*6, 150,"#font[52]{#geq3b}")

binlabel.SetTextAlign(11)
binlabel.SetTextSize(0.045)
binlabel.DrawLatex(31.3, 50,"#font[52]{Low N#lower[-0.1]{_{jets}} (odd bin #)}")
binlabel.DrawLatex(31.3, 25,"#font[52]{High N#lower[-0.1]{_{jets}} (even bin #)}")

a = TLine()
a.SetLineWidth(1)
a.SetLineStyle(3)
a.SetLineColor(kBlack)
a.DrawLine(nhbins/2+1,miny,nhbins/2+1,maxy)
for i in range(0,2):
    a.DrawLine((i+1)*6+1,miny,(i+1)*6+1,0.4*maxy)
    a.DrawLine(nhbins/2+(i+1)*6+1,miny,nhbins/2+(i+1)*6+1,0.4*maxy)


#     Bottom pad
# -----------------------------
bottom.cd()
hbotdummy = TH1D("","",nhbins,0.5,nhbins+1.5)
hbotdummy.GetYaxis().SetRangeUser(-2.9,2.9)
hbotdummy.GetYaxis().SetLabelSize(0.12)
hbotdummy.GetYaxis().SetTitle("Pull")
hbotdummy.GetYaxis().CenterTitle()
hbotdummy.GetYaxis().SetTitleSize(0.15)
hbotdummy.GetYaxis().SetTitleOffset(0.2)

hbotdummy.GetXaxis().SetLabelSize(0.12)
hbotdummy.GetXaxis().SetLabelOffset(0.02)
hbotdummy.GetXaxis().SetTitle("Bin #")
hbotdummy.GetXaxis().CenterTitle()
hbotdummy.GetXaxis().SetTitleSize(0.15)
hbotdummy.GetXaxis().SetTitleOffset(0.9)

hbotdummy.Draw()

grpull_pre.SetLineColor(tag1_color)
grpull_pre.SetLineWidth(1)
grpull_pre.SetFillColor(0)
grpull_pre.SetFillStyle(0)
grpull_pre.Draw('B')

grpull_post.SetFillColor(tag2_color)
grpull_post.SetFillStyle(3144)
grpull_post.Draw('B')

a = TLine()
a.SetLineWidth(1)
a.SetLineStyle(3)
a.DrawLine(nhbins/2+1, -2.9,nhbins/2+1,2.9)
a.SetLineColor(kGray+1)
for i in range(0,2):
    a.DrawLine((i+1)*6+1, -2.9,(i+1)*6+1,2.9)
    a.DrawLine(nhbins/2+(i+1)*6+1, -2.9,nhbins/2+(i+1)*6+1,2.9)


b = TLine()
b.SetLineWidth(1)
b.SetLineColor(kGray+1)
b.SetLineStyle(1)
b.DrawLine(0.5,0, nhbins+1.5,0)
b.SetLineStyle(2)
for i in [-1.,1.]:
    b.DrawLine(0.5,i, nhbins+1.5,i)
b.SetLineStyle(3)
for i in [-2.,2.]:
    b.DrawLine(0.5,i, nhbins+1.5,i)

can.Print('results'+tag1+'_vs'+tag2+'.pdf')


