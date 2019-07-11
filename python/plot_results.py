#! /usr/bin/env python
import os
import math
from ROOT import *

do_sig = True
compile_table = True


# --- Official plot
tag1 = '_nor4_2100'
tag1_lbl = 'R1--R3 fit'
tag1_color = kPink+2

tag2 = '_r4_2100' 
tag2_lbl = 'R1--R4 fit'
tag2_color = kAzure+1

tag_sig_c = '_r4_1900'

# --- Effect of systematics
# tag1 = '_nor4'
# tag1_lbl = 'R1--R3 fit'
# tag1_color = kBlack

# tag2 = '_nor4_nosys' 
# tag2_lbl = 'R1--R3 fit (no syst)'
# tag2_color = kGray+1

#         Bin names
#-------------------------------------
bins = []
for imj in ['lmj','hmj']:
    for imet in ['lmet','mmet','hmet']:
        for ir in ['r1','r2','r3','r4']:
            if ir=='r1' or ir=='r3':
                if imj=='lmj':
                    bins.append('_'.join([ir,imet]))
                else:
                    continue
            else:
                for inb in ['lnb','mnb','hnb']:
                    for inj in ['lnj','hnj']:
                        bins.append('_'.join([ir,imet,inb, inj,imj]))

nbins = len(bins)
# in the combine output they are ordered alphabetically
bins_alpha = sorted(bins)

#         Reading the data
#-------------------------------------
bkg_pre, ebkg_pre, pull_pre = [None]*nbins, [None]*nbins, [None]*nbins
data, edata_up, edata_dn = [None]*nbins, [None]*nbins, [None]*nbins
sig_nc, esig_nc = [None]*nbins, [None]*nbins

file_pre = TFile("root/fitDiagnostics"+tag1+".root")
bkg_in = file_pre.Get('shapes_fit_b/total_background')
data_in = file_pre.Get('shapes_fit_b/total_data')
sig_nc_in = file_pre.Get('shapes_prefit/total_signal')
for ibin in range(0,nbins):
    ibin_alpha = bins_alpha.index(bins[ibin])
    x, y  = Double(0), Double(0)
    data_in.GetPoint(ibin_alpha, x, y)
    data[ibin] = y
    edata_up[ibin] = data_in.GetErrorYhigh(ibin_alpha)
    edata_dn[ibin] = data_in.GetErrorYlow(ibin_alpha)
    bkg_pre[ibin] = bkg_in.GetBinContent(ibin_alpha+1)
    ebkg_pre[ibin] = bkg_in.GetBinError(ibin_alpha+1)
    pull_pre[ibin] = (data[ibin]-bkg_pre[ibin])/math.sqrt(bkg_pre[ibin]+ebkg_pre[ibin]*ebkg_pre[ibin])
    sig_nc[ibin] = sig_nc_in.GetBinContent(ibin_alpha+1)
    esig_nc[ibin] = sig_nc_in.GetBinError(ibin_alpha+1)
file_pre.Close()

bkg_post, ebkg_post, pull_post = [None]*nbins, [None]*nbins, [None]*nbins
file_post = TFile("root/fitDiagnostics"+tag2+".root")
bkg_in = file_post.Get('shapes_fit_b/total_background')
for ibin in range(0,nbins):
    ibin_alpha = bins_alpha.index(bins[ibin])
    bkg_post[ibin] = bkg_in.GetBinContent(ibin_alpha+1)
    ebkg_post[ibin] = bkg_in.GetBinError(ibin_alpha+1)
    pull_post[ibin] = (data[ibin]-bkg_post[ibin])/math.sqrt(bkg_post[ibin]+ebkg_post[ibin]*ebkg_post[ibin])
file_post.Close()

sig_c, esig_c = [None]*nbins, [None]*nbins
file_sig_c = TFile("root/fitDiagnostics"+tag_sig_c+".root")
sig_c_in = file_sig_c.Get('shapes_prefit/total_signal')
for ibin in range(0,nbins):
    ibin_alpha = bins_alpha.index(bins[ibin])
    sig_c[ibin] = sig_c_in.GetBinContent(ibin_alpha+1)
    esig_c[ibin] = sig_c_in.GetBinError(ibin_alpha+1)
file_sig_c.Close()

#         Making table
#-------------------------------------
tab = []
tab.append(open("tables/table_lmj"+tag1+"_vs"+tag2+".tex","w"))
tab.append(open("tables/table_hmj"+tag1+"_vs"+tag2+".tex","w"))
ncols = 6
if (do_sig): ncols +=2

tab_head = "\\begin{tabular}[tbp!]{ l cc"
for i in range(ncols-3): tab_head += " r"
tab_head += "} \n\\hline\\hline\n"
tab_head += "${\\cal L}=137$ fb$^{-1}$ &"
if do_sig:
    tab_head += " SUS-NC & SUS-C &"
tab_head += tag1_lbl+" & Pull & "+tag2_lbl+" & Pull & Obs. \\\\ \\hline\n"
for i in range(2): tab[i].write(tab_head)

save_rows = { 'r1_lmet':'', 'r3_lmet':'', 'r1_mmet':'', 'r3_mmet':'', 'r1_hmet':'', 'r3_hmet':''}
irow = 0
for ibin in range(nbins):
    tmp = bins[ibin].split("_")
    ireg = tmp[0].replace("r","R")
    imet = tmp[1].replace("lmet","$200<p_{T}^{\text{miss}}\\leq350$ GeV")
    imet = imet.replace("mmet","$350<p_{T}^{\text{miss}}\\leq500$ GeV")
    imet = imet.replace("hmet","$p_{T}^{\text{miss}}> 500$ GeV")
    inb,inj = '',''
    if ireg=='R2' or ireg=='R4':
        inb = tmp[2].replace("lnb","1b").replace("mnb","2b").replace("hnb","$\\geq$ 3b")
        if "hmet" in bins[ibin]:
            inj = tmp[3].replace("lnj"," 6--7j").replace("hnj","$\\geq$ 8j")
        else:
            inj = tmp[3].replace("lnj"," 7j").replace("hnj","$\\geq$ 8j")
    
    itab = ibin/42
    if (irow%14==0):
        tab[itab].write("\\hline\n\\multicolumn{"+str(ncols)+"}{c}{"+imet+"}  \\\\ \\hline\n");

    if ibin>41 and (ibin-42)%6==0: 
        irow +=1
        tmp_ = bins[ibin].split("_")
        insert_bin = tmp_[0].replace("r2","r1").replace("r4","r3")+'_'+tmp_[1]
        tab[itab].write(save_rows[insert_bin])
        if '3' in insert_bin:
            tab[itab].write("\\hline\n")


    cols = []
    if ireg=='R2' or ireg=='R4':
        cols.append('{0:<30}'.format(ireg+": "+inb+", "+inj))
    else:
        cols.append('{0:<30}'.format(ireg))
    if do_sig:
        cols.append('{0:>10.1f}'.format(sig_nc[ibin], esig_nc[ibin]))
        cols.append('{0:>10.1f}'.format(sig_c[ibin], esig_c[ibin]))
    cols.append('{0:>20}'.format('${0:.1f} \\pm {1:.1f}$'.format(bkg_pre[ibin], ebkg_pre[ibin])))
    if ireg=="R4": 
        cols.append('{0:>7.1f}'.format(pull_pre[ibin]))
    else: 
        cols.append('')
    cols.append('{0:>20}'.format('${0:.1f} \\pm {1:.1f}$'.format(bkg_post[ibin], ebkg_post[ibin])))
    if ireg=="R4": 
        cols.append('{0:>7.1f}'.format(pull_post[ibin]))
    else: 
        cols.append('')
    cols.append('{0:>10.0f}'.format(data[ibin]))
    tab[itab].write('&'.join(cols)+'\\\\\n')

    if ibin<42 and ('1' in ireg or '3' in ireg): 
        save_rows[bins[ibin]] = '&'.join(cols)+'\\\\\n'

    if '3' in ireg: 
        tab[itab].write("\\hline\n")

    irow += 1

for itab in tab:
    itab.write("\\hline\\hline\n \\end{tabular}\n")
    itab.close()


#         Making tables that can compile standalone
#----------------------------------------------------------
fulltab = []

for i in range(2):
    tabname = tab[i].name.split("/")[-1]
    fulltab.append(open("tables/full"+tabname,"w"))
    with open("txt/header.tex") as head:
        for line in head.readlines():
            fulltab[i].write(line)
    fulltab[i].write("\\begin{document}\n\\begin{preview}\n")
    with open(tab[i].name) as body:
        for line in body.readlines():
            fulltab[i].write(line)
    fulltab[i].write("\\end{preview}\n\\end{document}\n")
    fulltab[i].close()
    if (compile_table):
        print "Converting "+fulltab[i].name+" ..."
        os.system("pdflatex "+fulltab[i].name+" > /dev/null")

if compile_table:
    for i in range(2): 
        print "open "+fulltab[i].name.replace(".tex",".pdf")

#         Plotting the data
#-------------------------------------

nhbins = 36
htopdummy = TH1D("","",nhbins,0.5,nhbins+1.5)
grbkg_pre = TGraphErrors(nhbins)
grbkg_post = TGraphErrors(nhbins)
grdata = TGraphAsymmErrors(nhbins)
grpull_pre = TGraph(nhbins)
grpull_post = TGraph(nhbins)
i = 0
for ibin in range(nbins):
    if 'r4' not in bins[ibin]: continue
    i +=1
    grbkg_pre.SetPoint(i, i+0.5,bkg_pre[ibin])
    grbkg_pre.SetPointError(i, 0.5, ebkg_pre[ibin])    
    grpull_pre.SetPoint(i, i+0.5, pull_pre[ibin])

    grbkg_post.SetPoint(i, i+0.5,bkg_post[ibin])
    grbkg_post.SetPointError(i, 0.5, ebkg_post[ibin])
    grpull_post.SetPoint(i, i+0.5, pull_post[ibin])

    grdata.SetPoint(i, i+0.5, data[ibin])
    grdata.SetPointEYhigh(i, edata_up[ibin])
    grdata.SetPointEYlow(i, edata_dn[ibin])
    grdata.SetPointEXhigh(i, 0.)
    grdata.SetPointEXlow(i, 0.)


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
hbotdummy.GetXaxis().SetTitle("Bin number")
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

pname = 'plots/results'+tag1+'_vs'+tag2+'.pdf'
can.Print(pname)

print 'open', pname


