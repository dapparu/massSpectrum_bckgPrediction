#!/usr/bin/python

import sys, getopt, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import numpy as np

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import CMS_lumi, tdrstyle

ROOT.gROOT.SetBatch(True)

tdrstyle.setTDRStyle()

CMS_lumi.lumi_sqrtS = "13 TeV"
CMS_lumi.lumi_13TeV = "2018 - 51.2 fb^{-1}"
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "In progress"
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

def scale(h1):
    h1.Scale(1./h1.Integral())
    return h1

def ratioHisto(h1,h2):
    h3=h1.Clone()
    h3.Divide(h2)
    return h3

def ratioInt(h1,h2):
    h3=h2.Clone()
    h3.Reset()
    for i in range(0,h1.GetNbinsX()+1):
        e1=ROOT.Double(0.0)
        e2=ROOT.Double(0.0)
        a=h1.IntegralAndError(i,h1.GetNbinsX()+1,e1,"")
        b=h2.IntegralAndError(i,h1.GetNbinsX()+1,e2,"")
        if b != 0 and a != 0:
            c=math.sqrt((e1*e1)/(a*a)+(e2*e2)/(b*b))*a/b
            h3.SetBinContent(i,a/b)
            h3.SetBinError(i,c)
        else:
            h3.SetBinContent(i,0)
    return h3

def diffRel(h1,h2):
    h3=h2.Clone()
    h3.Reset()
    for i in range(0,h1.GetNbinsX()+1):
        a=h1.GetBinContent(i)
        b=h2.GetBinContent(i)
        e1=h1.GetBinError(i)
        e2=h2.GetBinError(i)
        if b != 0 and a != 0:
            c=math.sqrt((e1*e1)/(a*a)+(e2*e2)/(b*b))*b/a
            h3.SetBinContent(i,(a-b)/a)
            h3.SetBinError(i,c)
        else:
            h3.SetBinContent(i,0)
    return h3

def overflowInLastBin(h):
    res=h.Clone()
    res.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
    res.SetBinContent(h.GetNbinsX()+1,0)
    return res

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res


def statErr(h1,name):
    statErr=h1.Clone()
    statErr.SetName(name)
    for i in range (0,statErr.GetNbinsX()+1):
        if statErr.GetBinContent(i)>0:
            statErr.SetBinContent(i,statErr.GetBinError(i)/statErr.GetBinContent(i))
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr

def statErrRInt(h1,name):
    statErr=h1.Clone()
    statErr.SetName(name)
    c=ROOT.Double(0.0)
    e=ROOT.Double(0.0)
    for i in range (0,statErr.GetNbinsX()+1):
        c=h1.IntegralAndError(i,h1.GetNbinsX()+1,e,"")
        if e>0:
            statErr.SetBinContent(i,e/c)
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr


def systMass(nominal,down,up,name,typec,binned):
    if (binned==0):
        ra1=ratioInt(nominal,up)
        ra2=ratioInt(nominal,down)
    elif (binned==1):
        ra1=ratioHisto(nominal,up)
        ra2=ratioHisto(nominal,down)
    res=ra1.Clone()
    res.SetName(name)
    for i in range (0,res.GetNbinsX()+1):
        r1=ra1.GetBinContent(i)
        r2=ra2.GetBinContent(i)
        s1=abs(1-r1)
        s2=abs(1-r2)
        if (typec==0):
            m=max(s1,s2)
        elif (typec==1):
            m=(s1+s2)/2.
        res.SetBinContent(i,100*m)
    return res

def systTotal(list_h):
    res=list_h[0].Clone()
    for i in range (0,res.GetNbinsX()+1):
        systotal=0
        for h in list_h:
            systotal+=h.GetBinContent(i)*h.GetBinContent(i)
        res.SetBinContent(i,math.sqrt(systotal))
    return res


def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    return h1

def superposition(h1,h2,h3,name):
    c=TCanvas()
    '''
    h1=scale(h1)
    h2=scale(h2)
    h3=scale(h3)
    h1.SetMinimum(1e-3)
    h1.SetMaximum(1)
    '''
    h1=setColorAndMarker(h1,20,20)
    h2=setColorAndMarker(h2,30,21)
    h3=setColorAndMarker(h3,46,22)
    h1.Draw()
    h2.Draw("same")
    h3.Draw("same")
    c.SaveAs(name)
    return c

def binWidth(h1):
    res=h1.Clone()
    for i in range (0,h1.GetNbinsX()+1):
        res.SetBinContent(i,h1.GetBinContent(i)/h1.GetBinWidth(i))
        res.SetBinError(i,h1.GetBinError(i)/h1.GetBinWidth(i))
    return res

def plotter_(predNominal,predPullU,legNom,leg1,outDir,outTitle):
    c1=TCanvas("c1","c1",800,800)
    t1=TPad("t1","t1", 0.0, 0.40, 0.95, 0.95)

    t1.Draw()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    t2.SetGridy(1)
    t2.SetTopMargin(0.005)
    t2.SetBottomMargin(0.005)

    t3=TPad("t3","t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetTopMargin(0.005)
    t3.SetBottomMargin(0.4)

    t1.cd()
    c1.SetLogy(1)

    max_mass=2000
    min_entries=5e-5
    max_entries=5e6

    predNominal=setColorAndMarker(predNominal,1,20)
    predNominal.GetXaxis().SetRangeUser(0,max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries,max_entries)
    predNominal.SetTitle(";Mass (GeV);Tracks / bin width")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.05)
    predNominal.Draw()

    predPullU=setColorAndMarker(predPullU,46,21)
    predPullU.Draw("same")

    
    leg=TLegend(0.82,0.9,0.4,0.6);
    leg.AddEntry(predNominal,legNom,"PE1");
    leg.AddEntry(predPullU,leg1,"PE1");

    leg.Draw("same")
    
    LineLastBin=TLine(predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),0,predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(3)
    LineLastBin.SetLineColor(1)
    LineLastBin.Draw("same")


    c1.cd()
    t2.cd()

    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,0, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetXaxis().SetRangeUser(0,max_mass)
    frameR2.GetYaxis().SetTitle("Ratio #int_{m}^{#infty}")
    frameR2.SetMaximum(1.5)
    frameR2.SetMinimum(0.5)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(20) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(20) #font size
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetLabelSize(16) #font size
    frameR2.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetTitleSize(24) #font size
    frameR2.GetXaxis().SetTitleOffset(3.75)
    frameR2.Draw("AXIS")

    LineAtOne=TLine(0,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    LineAtOne.Draw("same")
    
    LineAt1p2=TLine(0,1.2,max_mass,1.2)
    LineAt1p2.SetLineStyle(4)
    LineAt1p2.SetLineColor(1)
    LineAt1p2.Draw("same")
    
    LineAt0p8=TLine(0,0.8,max_mass,0.8)
    LineAt0p8.SetLineStyle(4)
    LineAt0p8.SetLineColor(1)
    LineAt0p8.Draw("same")

    ratioInt2=ratioInt(predNominal,predPullU)
    ratioInt2.Draw("E0 same")

    ofile.cd()
    ratioInt2.Write()

    c1.cd()
    t3.cd()

    frameR3=frameR2.Clone()
    frameR3.GetYaxis().SetTitle("#frac{nominal}{pull}")
    frameR3.GetXaxis().SetRangeUser(0,max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(5)

    LineAtOne.Draw("same")
    LineAt1p2.Draw("same")
    LineAt0p8.Draw("same")
    
    predNominalCl2=predNominal.Clone()
    predNominalCl2.Divide(predPullU)

    predNominalCl2=setColorAndMarker(predNominalCl2,46,21)
    predNominalCl2.Draw("E0 same")

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)

    cmd='mkdir -p '+outDir
    os.system(cmd)

    c1.SaveAs(outDir+"/"+outTitle+".pdf")
    c1.SaveAs(outDir+"/"+outTitle+".root")
    c1.SaveAs(outDir+"/"+outTitle+".C")



def plotter(predNominal,predPullD,predPullU,legNom,leg1,leg2,outDir,outTitle):
    c1=TCanvas("c1","c1",800,800)
    t1=TPad("t1","t1", 0.0, 0.40, 0.95, 0.95)

    t1.Draw()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    t2.SetGridy(1)
    t2.SetTopMargin(0.005)
    t2.SetBottomMargin(0.005)

    t3=TPad("t3","t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetTopMargin(0.005)
    t3.SetBottomMargin(0.4)

    t1.cd()
    c1.SetLogy(1)

    max_mass=2000
    min_entries=predNominal.GetMinimum()/10
    max_entries=predNominal.GetMaximum()*10

    predNominal=setColorAndMarker(predNominal,1,20)
    predNominal.GetXaxis().SetRangeUser(0,max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries,max_entries)
    predNominal.SetTitle(";Mass (GeV);Tracks / bin width")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.05)
    predNominal.Draw()

    predPullD=setColorAndMarker(predPullD,38,21)
    predPullD.Draw("same")

    predPullU=setColorAndMarker(predPullU,46,21)
    predPullU.Draw("same")

    
    leg=TLegend(0.82,0.9,0.4,0.6);
    leg.AddEntry(predNominal,legNom,"PE1");
    leg.AddEntry(predPullD,leg1,"PE1");
    leg.AddEntry(predPullU,leg2,"PE1");

    
    LineLastBin=TLine(predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),0,predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(1)
    LineLastBin.SetLineColor(1)
    LineLastBin.Draw("same")
   
    t=ROOT.TText(0.95,0.7,"+overflow")
    t.SetNDC(True)
    t.SetTextColor(1)
    t.SetTextFont(43)
    t.SetTextSize(24)
    t.SetTextAngle(90)
    t.Draw("same")

    leg.Draw("same")
    
    c1.cd()
    t2.cd()

    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,0, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("Ratio #int_{m}^{#infty}")
    frameR2.GetXaxis().SetRangeUser(0,max_mass)
    frameR2.SetMaximum(1.5)
    frameR2.SetMinimum(0.5)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(20) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(20) #font size
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetLabelSize(16) #font size
    frameR2.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetTitleSize(24) #font size
    frameR2.GetXaxis().SetTitleOffset(3.75)
    frameR2.Draw("AXIS")

    LineAtOne=TLine(0,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    LineAtOne.Draw("same")
    
    LineAt1p2=TLine(0,1.2,max_mass,1.2)
    LineAt1p2.SetLineStyle(4)
    LineAt1p2.SetLineColor(1)
    LineAt1p2.Draw("same")
    
    LineAt0p8=TLine(0,0.8,max_mass,0.8)
    LineAt0p8.SetLineStyle(4)
    LineAt0p8.SetLineColor(1)
    LineAt0p8.Draw("same")

    ratioInt1=ratioInt(predNominal,predPullD)
    ratioInt2=ratioInt(predNominal,predPullU)
    ratioInt1.Draw("E0 same")
    ratioInt2.Draw("E0 same")

    ofile.cd()
    ratioInt1.Write()
    ratioInt2.Write()

    c1.cd()
    t3.cd()

    frameR3=frameR2.Clone()
    frameR3.GetYaxis().SetTitle("#frac{nominal}{pull}")
    frameR3.GetXaxis().SetRangeUser(0,max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(5)

    LineAtOne.Draw("same")
    LineAt1p2.Draw("same")
    LineAt0p8.Draw("same")
    
    predNominalCl1=predNominal.Clone()
    predNominalCl2=predNominal.Clone()
    predNominalCl1.Divide(predPullD)
    predNominalCl2.Divide(predPullU)

    predNominalCl1=setColorAndMarker(predNominalCl1,38,21)
    predNominalCl1.Draw("E0 same")

    predNominalCl2=setColorAndMarker(predNominalCl2,46,21)
    predNominalCl2.Draw("E0 same")

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)

    cmd='mkdir -p '+outDir
    os.system(cmd)

    c1.SaveAs(outDir+"/"+outTitle+".pdf")
    c1.SaveAs(outDir+"/"+outTitle+".root")
    c1.SaveAs(outDir+"/"+outTitle+".C")

def plotSummary(syst_stat,syst_eta,syst_ih,syst_p,sysTot,xtitle,outTitle):
    syst_stat=lowEdge(syst_stat)
    syst_eta=lowEdge(syst_eta)
    syst_ih=lowEdge(syst_ih)
    syst_p=lowEdge(syst_p)
    sysTot=lowEdge(sysTot)

    c2=TCanvas()

    #c2.SetLogy()
    c2.SetGrid()
    syst_stat.SetMinimum(0)
    syst_stat.SetMaximum(1e2)
    syst_stat.GetXaxis().SetTitle(xtitle)
    syst_stat.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    syst_stat.GetXaxis().SetNdivisions(510)
    syst_stat.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetXaxis().SetLabelSize(20) #font size
    syst_stat.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetYaxis().SetLabelSize(20) #font size

    syst_stat.GetXaxis().SetRangeUser(0,2000)

    syst_stat=setColorAndMarker(syst_stat,1,20)
    syst_eta=setColorAndMarker(syst_eta,30,21)
    syst_ih=setColorAndMarker(syst_ih,38,22)
    syst_p=setColorAndMarker(syst_p,46,23)
    sysTot=setColorAndMarker(sysTot,28,34)

    leg2=TLegend(0.2,0.7,0.4,0.9)
    leg2.AddEntry(sysTot,"Total","PE1");
    leg2.AddEntry(syst_stat,"Stat.","PE1");
    leg2.AddEntry(syst_eta,"#eta binning","PE1");
    leg2.AddEntry(syst_ih,"I_{h} binning","PE1");
    leg2.AddEntry(syst_p,"p binning","PE1");

    syst_stat.Draw("AP")
    syst_eta.Draw("P")
    syst_ih.Draw("P")
    syst_p.Draw("P")
    sysTot.Draw("P")
    leg2.Draw("same")

    CMS_lumi.CMS_lumi(c2, iPeriod, iPos)
    
    commandMkdir='mkdir -p '+oDir+'pdf '+oDir+'Cfile '+oDir+'rootfile'
    os.system(commandMkdir)

    c2.SaveAs(oDir+"pdf/summary_"+outTitle+".pdf")
    c2.SaveAs(oDir+"Cfile/summary_"+outTitle+".root")
    c2.SaveAs(oDir+"rootfile/summary_"+outTitle+".C")

#directory="/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/resultSingleMuon_2018D_41p2/"
directory="/opt/sbg/cms/ui14_data2/dapparu/CMSSW_10_6_27/src/BackgroundStudies/tmpOut/"

inputNominal=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta4_rebinIh2_rebinP4_rebinMass1_analysed.root"
inputEtaD=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta8_rebinIh2_rebinP4_rebinMass1_analysed.root"
inputEtaU=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta2_rebinIh2_rebinP4_rebinMass1_analysed.root"
inputIhD=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_analysed.root"
inputIhU=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta4_rebinIh1_rebinP4_rebinMass1_analysed.root"
inputPD=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta4_rebinIh2_rebinP8_rebinMass1_analysed.root"
inputPU=directory+"outfile_SingleMuon_TkOnly_All2018_1overP_noIhCut_v41p2_v4_370092360_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_analysed.root"


ifileNominal=TFile(inputNominal)
ifileEtaD=TFile(inputEtaD)
ifileEtaU=TFile(inputEtaU)
ifileIhD=TFile(inputIhD)
ifileIhU=TFile(inputIhU)
ifilePD=TFile(inputPD)
ifilePU=TFile(inputPU)
   
plotType="mass_predBC_"
region="90ias100"

oDir="systematics_rebin_dir/"
cmd='mkdir -p '+oDir
os.system(cmd)
outTitle="syst"

ofile=TFile(oDir+outTitle+".root","RECREATE")

predNominal_def=ifileNominal.Get(plotType+region)
predEtaD=ifileEtaD.Get(plotType+region)
predEtaU=ifileEtaU.Get(plotType+region)
predIhD=ifileIhD.Get(plotType+region)
predIhU=ifileIhU.Get(plotType+region)
predPD=ifilePD.Get(plotType+region)
predPU=ifilePU.Get(plotType+region)
predTrueTemplates=ifileNominal.Get(plotType+"True_BC")


rebinning=array.array('d',[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,750,800,850,900,1000,1100,1200,1400,2000])
rebinning1=np.arange(0.,500.,50.)
rebinning2=np.arange(500.,900.,100.)
rebinning3=np.arange(900.,2100.,1100.)

rebinning=np.concatenate((rebinning1,rebinning2))
rebinning=np.concatenate((rebinning,rebinning3))
rebinning=array.array('d',rebinning)
#print rebinning

sizeRebinning=len(rebinning)-1

predNominal_def=predNominal_def.Rebin(sizeRebinning,"nominal_new",rebinning)
predEtaD=predEtaD.Rebin(sizeRebinning,"etaD_new",rebinning)
predEtaU=predEtaU.Rebin(sizeRebinning,"etaU_new",rebinning)
predIhD=predIhD.Rebin(sizeRebinning,"ihD_new",rebinning)
predIhU=predIhU.Rebin(sizeRebinning,"ihU_new",rebinning)
predPD=predPD.Rebin(sizeRebinning,"pD_new",rebinning)
predPU=predPU.Rebin(sizeRebinning,"pU_new",rebinning)
predTrueTemplates=predTrueTemplates.Rebin(sizeRebinning,"corr_new",rebinning)

predNominal_def=overflowInLastBin(predNominal_def)
predEtaD=overflowInLastBin(predEtaD)
predEtaU=overflowInLastBin(predEtaU)
predIhD=overflowInLastBin(predIhD)
predIhU=overflowInLastBin(predIhU)
predPD=overflowInLastBin(predPD)
predPU=overflowInLastBin(predPU)
predTrueTemplates=overflowInLastBin(predTrueTemplates)

predNominal_def=binWidth(predNominal_def)
predEtaD=binWidth(predEtaD)
predEtaU=binWidth(predEtaU)
predIhD=binWidth(predIhD)
predIhU=binWidth(predIhU)
predPD=binWidth(predPD)
predPU=binWidth(predPU)
predTrueTemplates=binWidth(predTrueTemplates)

syst_stat=statErrRInt(predNominal_def,"Stat")
syst_stat_binned=statErr(predNominal_def,"Stat")

syst_eta=systMass(predNominal_def,predEtaD,predEtaU,"Eta",0,0)
syst_ih=systMass(predNominal_def,predIhD,predIhU,"Ih",0,0)
syst_p=systMass(predNominal_def,predPD,predPU,"P",0,0)
syst_corr=systMass(predNominal_def,predTrueTemplates,predTrueTemplates,"Correlation",0,0)

syst_eta_mean=systMass(predNominal_def,predEtaD,predEtaU,"Eta_mean",1,0)
syst_ih_mean=systMass(predNominal_def,predIhD,predIhU,"Ih_mean",1,0)
syst_p_mean=systMass(predNominal_def,predPD,predPU,"P_mean",1,0)

syst_eta_binned=systMass(predNominal_def,predEtaD,predEtaU,"Eta",0,1)
syst_ih_binned=systMass(predNominal_def,predIhD,predIhU,"Ih",0,1)
syst_p_binned=systMass(predNominal_def,predPD,predPU,"P",0,1)

syst_eta_binned_mean=systMass(predNominal_def,predEtaD,predEtaU,"Eta",1,1)
syst_ih_binned_mean=systMass(predNominal_def,predIhD,predIhU,"Ih",1,1)
syst_p_binned_mean=systMass(predNominal_def,predPD,predPU,"P",1,1)

syst_stat.Write()
syst_eta.Write()
syst_ih.Write()
syst_p.Write()
syst_eta_mean.Write()
syst_ih_mean.Write()
syst_p_mean.Write()


listOfSyst=[syst_stat,syst_eta,syst_ih,syst_p]
listOfSyst_mean=[syst_stat,syst_eta_mean,syst_ih_mean,syst_p_mean]

sysTot=systTotal(listOfSyst)
sysTot_mean=systTotal(listOfSyst_mean)
sysTot.Write()
sysTot_mean.Write()


listOfSyst_binned=[syst_stat_binned,syst_eta_binned,syst_ih_binned,syst_p_binned]
listOfSyst_binned_mean=[syst_stat_binned,syst_eta_binned_mean,syst_ih_binned_mean,syst_p_binned_mean]

sysTot_binned=systTotal(listOfSyst_binned)
sysTot_binned_mean=systTotal(listOfSyst_binned_mean)

ihB=ifileNominal.Get("ih_eta_regionB_90_py")
ihB_Down=ifileIhD.Get("ih_eta_regionB_90_py")
ihB_Up=ifileIhU.Get("ih_eta_regionB_90_py")
pC=ifileNominal.Get("eta_p_regionC_med_px")
pC_Down=ifilePD.Get("eta_p_regionC_med_px")
pC_Up=ifilePU.Get("eta_p_regionC_med_px")
etaB=ifileNominal.Get("ih_eta_regionB_90_px")
etaB_Down=ifileEtaD.Get("ih_eta_regionB_90_px")
etaB_Up=ifileEtaU.Get("ih_eta_regionB_90_px")

plotter(predNominal_def,predEtaD,predEtaU,"Nominal (50 bins)","Down (25 bins)","Up (100 bins)",oDir+"syst_Eta","plot_Eta")
plotter(predNominal_def,predIhD,predIhU,"Nominal (1000 bins)","Down (500 bins)","Up (2000 bins)",oDir+"syst_Ih","plot_Ih")
plotter(predNominal_def,predPD,predPU,"Nominal (500 bins)","Down (250 bins)","Up (1000 bins)",oDir+"syst_P","plot_P")
plotter_(predNominal_def,predTrueTemplates,"Templates from B and C","Templates from D",oDir+"syst_Corr","plot_Corr")

xtitle="Mass cut [GeV]"
plotSummary(syst_stat,syst_eta,syst_ih,syst_p,sysTot,xtitle,outTitle)
plotSummary(syst_stat,syst_eta_mean,syst_ih_mean,syst_p_mean,sysTot_mean,xtitle,"mean_"+outTitle)

xtitle="Mass bin"
plotSummary(syst_stat_binned,syst_eta_binned,syst_ih_binned,syst_p_binned,sysTot_binned,xtitle,"binned_"+outTitle)
plotSummary(syst_stat_binned,syst_eta_binned_mean,syst_ih_binned_mean,syst_p_binned_mean,sysTot_binned_mean,xtitle,"binned_mean_"+outTitle)

