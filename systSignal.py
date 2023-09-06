import ROOT, sys, os, time, re, numpy
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
from tqdm import tqdm


import sys, getopt, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import numpy as np

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import CMS_lumi, tdrstyle

tdrstyle.setTDRStyle()

CMS_lumi.lumi_sqrtS = "13 TeV"
#CMS_lumi.lumi_13TeV = "t#bar{t} fully leptonic - 26 fb^{-1}"
#CMS_lumi.lumi_13TeV = "W+jets - 26 fb^{-1}"
#CMS_lumi.lumi_13TeV = "2018 - 51.2 fb^{-1}"
#CMS_lumi.lumi_13TeV = "59.7 fb^{-1}"
CMS_lumi.lumi_13TeV = "101 fb^{-1}"
CMS_lumi.writeExtraText = True
CMS_lumi.cmsText = ""
CMS_lumi.extraText = "Private work"
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

codeVersion = sys.argv[1]

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1)
ROOT.gStyle.SetPadBottomMargin(0.14)
ROOT.gStyle.SetPadLeftMargin(0.15)

def scale(h1):
    h1.Scale(1./h1.Integral())
    return h1

def ratioHisto(h2,h1):
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

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    return h1

def overflowInLastBin(h):
    res=h.Clone()
    res.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
    res.SetBinContent(h.GetNbinsX()+1,0)
    return res

def underflowInFirstBin(h):
    res=h.Clone()
    res.SetBinContent(1,h.GetBinContent(0)+h.GetBinContent(1))
    res.SetBinContent(0,0)
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

def superposition(h1,h2,h3,name):
    c=TCanvas()
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

def allSet(h,sizeRebinning,rebinning,st):
    h=h.Rebin(sizeRebinning,st,rebinning)
    h=overflowInLastBin(h)
    h=underflowInFirstBin(h)
    return h

def systMass(nominal,down,up,name,typec,binned,mini=0):
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
            if (mini==1):
                m=min(s1,s2)
        elif (typec==1):
            m=(s1+s2)/2.
        res.SetBinContent(i,100*m)
        if(mini==1):res.SetBinContent(i,0)
    return res

def systTotal(list_h):
    res=list_h[0].Clone()
    for i in range (0,res.GetNbinsX()+1):
        systotal=0
        for h in list_h:
            systotal+=h.GetBinContent(i)*h.GetBinContent(i)
        res.SetBinContent(i,math.sqrt(systotal))
    return res


def plotter(predNominal,predPullD,predPullU,legNom,leg1,leg2,outDir,outTitle):
    c1=TCanvas("c1","c1",800,800)
    t1=TPad("t1","t1", 0.0, 0.40, 0.95, 0.9)

    t1.Draw()
    #t1.cd()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    #t2.cd()
    t2.SetGridy(1)
    t2.SetTopMargin(0.005)
    t2.SetBottomMargin(0.005)

    t3=TPad("t3","t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    #t3.cd()
    t3.SetGridy(1)
    t3.SetTopMargin(0.005)
    t3.SetBottomMargin(0.4)

    t1.cd()
    c1.SetLogy(1)

    max_mass=2000
    min_entries=5e-5
    #min_entries=predNominal.GetMinimum()/10
    max_entries=predNominal.GetMaximum()*10

    predNominal=setColorAndMarker(predNominal,1,20)
    predNominal.GetXaxis().SetRangeUser(0,max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries,max_entries)
    predNominal.SetMinimum(min_entries)
    predNominal.SetMaximum(max_entries)
    predNominal.SetTitle(";Mass (GeV);Tracks")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.05)
    predNominal.Draw()

    predPullD=setColorAndMarker(predPullD,38,21)
    predPullD.Draw("same")

    predPullU=setColorAndMarker(predPullU,46,21)
    predPullU.Draw("same")

    
    leg=TLegend(0.2,0.65,0.4,0.9)
    leg.AddEntry(predNominal,legNom,"PE1")
    leg.AddEntry(predPullD,leg1,"PE1")
    leg.AddEntry(predPullU,leg2,"PE1")

    
    LineLastBin=TLine(predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),0,predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(1)
    LineLastBin.SetLineColor(1)
    LineLastBin.Draw("same")
   
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
    frameR2.SetMaximum(1.1)
    if(outTitle=="plot_Trigger"):
        frameR2.SetMaximum(1.3)
    frameR2.SetMinimum(0.9)
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
    
    LineAt0p8=TLine(0,0.8,max_mass,0.8)
    LineAt0p8.SetLineStyle(4)
    LineAt0p8.SetLineColor(1)

    ratioInt1=ratioInt(predNominal,predPullD)
    ratioInt2=ratioInt(predNominal,predPullU)
    ratioInt1.Draw("E0 same")
    ratioInt2.Draw("E0 same")
    
    ratioInt1.Write()
    ratioInt2.Write()

    c1.cd()
    t3.cd()

    frameR3=frameR2.Clone()
    frameR3.GetYaxis().SetTitle("#frac{nominal}{var}")
    frameR3.GetXaxis().SetRangeUser(0,max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(5)

    LineAtOne.Draw("same")
   
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

def plotSummary(syst_K,syst_C,syst_PU,syst_probq,syst_ias,syst_pt,syst_trigger,sysTot,xtitle,outTitle,outDir):

    syst_K=lowEdge(syst_K)
    syst_C=lowEdge(syst_C)
    syst_PU=lowEdge(syst_PU)
    syst_probq=lowEdge(syst_probq)
    syst_ias=lowEdge(syst_ias)
    syst_pt=lowEdge(syst_pt)
    syst_trigger=lowEdge(syst_trigger)

    sysTot=lowEdge(sysTot)

    c2=TCanvas()

    #c2.SetLogy()
    c2.SetGrid()
    syst_K.SetMinimum(0)
    syst_K.SetMaximum(100)
    syst_K.GetXaxis().SetTitle(xtitle)
    syst_K.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    syst_K.GetXaxis().SetNdivisions(510)
    syst_K.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_K.GetXaxis().SetLabelSize(20) #font size
    syst_K.GetXaxis().SetTitleSize(0.04) #font size
    syst_K.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_K.GetYaxis().SetLabelSize(20) #font size
    syst_K.GetYaxis().SetLabelSize(20) #font size
    syst_K.GetYaxis().SetTitleSize(0.04) #font size

    syst_K.GetXaxis().SetRangeUser(0,2000)

    syst_K=setColorAndMarker(syst_K,30,21)
    syst_C=setColorAndMarker(syst_C,38,22)
    syst_PU=setColorAndMarker(syst_PU,46,23)
    syst_probq=setColorAndMarker(syst_probq,43,43)
    syst_ias=setColorAndMarker(syst_ias,45,45)
    syst_pt=setColorAndMarker(syst_pt,39,29)
    syst_trigger=setColorAndMarker(syst_trigger,40,39)

    sysTot=setColorAndMarker(sysTot,28,34)

    leg2=TLegend(0.2,0.55,0.5,0.8)
    leg2.AddEntry(sysTot,"Total","PE1")
    leg2.AddEntry(syst_K,"K","PE1")
    leg2.AddEntry(syst_C,"C","PE1")
    leg2.AddEntry(syst_PU,"PU","PE1")
    leg2.AddEntry(syst_probq,"F^{pixel}","PE1")
    leg2.AddEntry(syst_ias,"G^{strip}","PE1")
    leg2.AddEntry(syst_pt,"p_{T}","PE1")
    leg2.AddEntry(syst_trigger,"Trigger","PE1")

    #c2.SetLogy()
    
    syst_K.Draw("AP")
    leg2.Draw("same")
    syst_C.Draw("P")
    syst_PU.Draw("P")
    syst_probq.Draw("P")
    syst_ias.Draw("P")
    syst_pt.Draw("P")
    syst_trigger.Draw("P")
    sysTot.Draw("P")
    

    CMS_lumi.CMS_lumi(c2, iPeriod, iPos)
    
    #commandMkdir='mkdir -p '+oDir+'pdf '+oDir+'Cfile '+oDir+'rootfile'
    #os.system(commandMkdir)

    c2.SaveAs(outDir+"pdf/summary_"+outTitle+".pdf")
    c2.SaveAs(outDir+"Cfile/summary_"+outTitle+".C")
    c2.SaveAs(outDir+"rootfile/summary_"+outTitle+".root")

    c2.SaveAs("summary_"+outTitle+".pdf")

def massDown(h1,h2):
    h1.SetName(h1.GetName()+"_down")
    for i in range (0,h1.GetNbinsX()+1):
        h1.SetBinContent(i,h1.GetBinContent(i)*(1-h2.GetBinContent(i)/100.))
    return h1

def massUp(h1,h2):
    h1.SetName(h1.GetName()+"_up")
    for i in range (0,h1.GetNbinsX()+1):
        h1.SetBinContent(i,h1.GetBinContent(i)*(1+h2.GetBinContent(i)/100.))
    return h1

SignalSamples = [
"crab_Analysis_2018_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-800_CodeV"+codeVersion+"_v1.root",

#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",

#"crab_Analysis_2018_HSCPpairStau_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-247_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-308_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-432_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-651_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-745_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-871_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1029_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1218_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1409_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1599_CodeV"+codeVersion+"_v1.root",

#"crab_Analysis_2018_HSCPgmsbStau_M-1029_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-1218_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-1409_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-1599_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-247_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-308_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-432_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-557_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-651_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-745_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgmsbStau_M-871_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstopOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",

"crab_Analysis_2018_HSCPstop_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPstop_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-800_CodeV"+codeVersion+"_v1.root",

#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-500_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-500_CodeV"+codeVersion+"_v1.root",

]


sampleName = {
"crab_Analysis_2018_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root" : "Gluino_M-1000",
"crab_Analysis_2018_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root" : "Gluino_M-1400",
"crab_Analysis_2018_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root" : "Gluino_M-1600",
"crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root" : "Gluino_M-1800",
"crab_Analysis_2018_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root" : "Gluino_M-2000",
"crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root" : "Gluino_M-2200",
"crab_Analysis_2018_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root" : "Gluino_M-2400",
"crab_Analysis_2018_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root" : "Gluino_M-2600",
"crab_Analysis_2018_HSCPgluino_M-500_CodeV"+codeVersion+"_v1.root" : "Gluino_M-500",
"crab_Analysis_2018_HSCPgluino_M-800_CodeV"+codeVersion+"_v1.root" : "Gluino_M-800",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-500",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-800",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-1000",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-1200",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-1400",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-1600",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-1800",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-2000",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-2200",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-2400",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root" : "GluinoOnlyNeutral_M-2600",
"crab_Analysis_2018_HSCPpairStau_M-200_CodeV"+codeVersion+"_v1.root" : "ppStau_M-200",
"crab_Analysis_2018_HSCPpairStau_M-247_CodeV"+codeVersion+"_v1.root" : "ppStau_M-247",
"crab_Analysis_2018_HSCPpairStau_M-308_CodeV"+codeVersion+"_v1.root" : "ppStau_M-308",
"crab_Analysis_2018_HSCPpairStau_M-432_CodeV"+codeVersion+"_v1.root" : "ppStau_M-432",
"crab_Analysis_2018_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root" : "ppStau_M-557",
"crab_Analysis_2018_HSCPpairStau_M-651_CodeV"+codeVersion+"_v1.root" : "ppStau_M-651",
"crab_Analysis_2018_HSCPpairStau_M-745_CodeV"+codeVersion+"_v1.root" : "ppStau_M-745",
"crab_Analysis_2018_HSCPpairStau_M-871_CodeV"+codeVersion+"_v1.root" : "ppStau_M-871",
"crab_Analysis_2018_HSCPpairStau_M-1029_CodeV"+codeVersion+"_v1.root" : "ppStau_M-1029",
"crab_Analysis_2018_HSCPpairStau_M-1218_CodeV"+codeVersion+"_v1.root" : "ppStau_M-1218",
"crab_Analysis_2018_HSCPpairStau_M-1409_CodeV"+codeVersion+"_v1.root" : "ppStau_M-1409",
"crab_Analysis_2018_HSCPpairStau_M-1599_CodeV"+codeVersion+"_v1.root" : "ppStau_M-1599",
"crab_Analysis_2018_HSCPgmsbStau_M-1029_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-1029",
"crab_Analysis_2018_HSCPgmsbStau_M-1218_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-1218",
"crab_Analysis_2018_HSCPgmsbStau_M-1409_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-1409",
"crab_Analysis_2018_HSCPgmsbStau_M-1599_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-1599",
"crab_Analysis_2018_HSCPgmsbStau_M-200_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-200",
"crab_Analysis_2018_HSCPgmsbStau_M-247_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-247",
"crab_Analysis_2018_HSCPgmsbStau_M-308_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-308",
"crab_Analysis_2018_HSCPgmsbStau_M-432_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-432",
"crab_Analysis_2018_HSCPgmsbStau_M-557_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-557",
"crab_Analysis_2018_HSCPgmsbStau_M-651_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-651",
"crab_Analysis_2018_HSCPgmsbStau_M-745_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-745",
"crab_Analysis_2018_HSCPgmsbStau_M-871_CodeV"+codeVersion+"_v1.root" : "gmsbStau_M-871",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-1000",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-1200",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-1400",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-1600",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-1800",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-2000",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-2200",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-2400",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-2600",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-500",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root" : "StopOnlyNeutral_M-800",
"crab_Analysis_2018_HSCPstop_M-1000_CodeV"+codeVersion+"_v1.root" : "Stop_M-1000",
"crab_Analysis_2018_HSCPstop_M-1200_CodeV"+codeVersion+"_v1.root" : "Stop_M-1200",
"crab_Analysis_2018_HSCPstop_M-1400_CodeV"+codeVersion+"_v1.root" : "Stop_M-1400",
"crab_Analysis_2018_HSCPstop_M-1600_CodeV"+codeVersion+"_v1.root" : "Stop_M-1600",
"crab_Analysis_2018_HSCPstop_M-1800_CodeV"+codeVersion+"_v1.root" : "Stop_M-1800",
"crab_Analysis_2018_HSCPstop_M-2000_CodeV"+codeVersion+"_v1.root" : "Stop_M-2000",
"crab_Analysis_2018_HSCPstop_M-2200_CodeV"+codeVersion+"_v1.root" : "Stop_M-2200",
"crab_Analysis_2018_HSCPstop_M-2400_CodeV"+codeVersion+"_v1.root" : "Stop_M-2400",
"crab_Analysis_2018_HSCPstop_M-2600_CodeV"+codeVersion+"_v1.root" : "Stop_M-2600",
"crab_Analysis_2018_HSCPstop_M-500_CodeV"+codeVersion+"_v1.root" : "Stop_M-500",
"crab_Analysis_2018_HSCPstop_M-800_CodeV"+codeVersion+"_v1.root" : "Stop_M-800",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1000_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-1000",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1400_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-1400",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1800_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-1800",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-200_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-200",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2200_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-2200",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2600_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-2600",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-400_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-400",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-500_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-500",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-800_CodeV"+codeVersion+"_v1.root" : "TauPrime1e_M-800",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1000_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-1000",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1400_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-1400",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1800_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-1800",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-200_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-200",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2200_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-2200",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2600_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-2600",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-400_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-400",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-500_CodeV"+codeVersion+"_v1.root" : "TauPrime2e_M-500",
}

RegS="SR3"

idir = "/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/"
region = "HSCParticleAnalyzer/BaseName/PostS_"+RegS+"_Mass"

ofileAllSignal = TFile("massShapeSignal_"+RegS+".root","RECREATE")

for sample in SignalSamples:
    print sample
    directory = "systSignal_25april/"+sampleName.get(sample)+"/"+RegS+"_73p0/"
    commandMkdir='mkdir -p '+directory+'pdf '+directory+'Cfile '+directory+'rootfile'
    os.system(commandMkdir)

    ifile = ROOT.TFile(idir+sample)
    mass_nominal = ifile.Get(region)
    mass_PU_up = ifile.Get(region+"_Pileup_up")
    mass_PU_down = ifile.Get(region+"_Pileup_down")
    mass_ProbQ_up = ifile.Get(region+"_ProbQNoL1_up")
    mass_ProbQ_down = ifile.Get(region+"_ProbQNoL1_down")
    mass_Ias_up = ifile.Get(region+"_Ias_up")
    mass_Ias_down = ifile.Get(region+"_Ias_down")
    mass_Pt_up = ifile.Get(region+"_Pt_up")
    mass_Pt_down = ifile.Get(region+"_Pt_down")
    mass_Trigger_up = ifile.Get(region+"_Trigger_up")
    mass_Trigger_down = ifile.Get(region+"_Trigger_down")
    mass_K_up1 = ifile.Get(region+"_K_up1")
    mass_K_down1 = ifile.Get(region+"_K_down1")
    mass_C_up1 = ifile.Get(region+"_C_up1")
    mass_C_down1 = ifile.Get(region+"_C_down1")
    mass_K_up2 = ifile.Get(region+"_K_up2")
    mass_K_down2 = ifile.Get(region+"_K_down2")
    mass_C_up2 = ifile.Get(region+"_C_up2")
    mass_C_down2 = ifile.Get(region+"_C_down2")

    rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,405.,435.,475.,525.,585.,660.,755.,875.,1025.,1210.,1440.,1730.,2000.])
    
    sizeRebinning=len(rebinning)-1

    mass_nominal = allSet(mass_nominal, sizeRebinning, rebinning, "nominal")
    mass_PU_up = allSet(mass_PU_up, sizeRebinning, rebinning, "PU_up")
    mass_PU_down = allSet(mass_PU_down, sizeRebinning, rebinning, "PU_down")
    mass_ProbQ_up = allSet(mass_ProbQ_up, sizeRebinning, rebinning, "ProbQ_up")
    mass_ProbQ_down = allSet(mass_ProbQ_down, sizeRebinning, rebinning, "ProbQ_down")
    mass_Ias_up = allSet(mass_Ias_up, sizeRebinning, rebinning, "Ias_up")
    mass_Ias_down = allSet(mass_Ias_down, sizeRebinning, rebinning, "Ias_down")
    mass_Pt_up = allSet(mass_Pt_up, sizeRebinning, rebinning, "Pt_up")
    mass_Pt_down = allSet(mass_Pt_down, sizeRebinning, rebinning, "Pt_down")
    mass_Trigger_up = allSet(mass_Trigger_up, sizeRebinning, rebinning, "Trigger_up")
    mass_Trigger_down = allSet(mass_Trigger_down, sizeRebinning, rebinning, "Trigger_down")
    mass_K_up1 = allSet(mass_K_up1, sizeRebinning, rebinning, "K_up1")
    mass_K_down1 = allSet(mass_K_down1, sizeRebinning, rebinning, "K_down1")
    mass_C_up1 = allSet(mass_C_up1, sizeRebinning, rebinning, "C_up1")
    mass_C_down1 = allSet(mass_C_down1, sizeRebinning, rebinning, "C_down1")
    mass_K_up2 = allSet(mass_K_up2, sizeRebinning, rebinning, "K_up2")
    mass_K_down2 = allSet(mass_K_down2, sizeRebinning, rebinning, "K_down2")
    mass_C_up2 = allSet(mass_C_up2, sizeRebinning, rebinning, "C_up2")
    mass_C_down2 = allSet(mass_C_down2, sizeRebinning, rebinning, "C_down2")

    syst_PU = systMass(mass_nominal,mass_PU_down,mass_PU_up,"",0,1)
    syst_ProbQ = systMass(mass_nominal,mass_ProbQ_down,mass_ProbQ_up,"",0,1)
    syst_Ias = systMass(mass_nominal,mass_Ias_down,mass_Ias_up,"",0,1)
    syst_Pt = systMass(mass_nominal,mass_Pt_down,mass_Pt_up,"",0,1)
    syst_Trigger = systMass(mass_nominal,mass_Trigger_down,mass_Trigger_up,"",0,1)
    syst_K1 = systMass(mass_nominal,mass_K_down1,mass_K_up1,"",0,1)
    syst_C1 = systMass(mass_nominal,mass_C_down1,mass_C_up1,"",0,1)
    syst_K2 = systMass(mass_nominal,mass_K_down2,mass_K_up2,"",0,1)
    syst_C2 = systMass(mass_nominal,mass_C_down2,mass_C_up2,"",0,1)
 
    listOfSyst = [syst_PU, syst_ProbQ, syst_Ias, syst_Pt, syst_Trigger, syst_K1, syst_C1]
    
    ofileSample=TFile(directory+"ofile_"+sampleName.get(sample)+".root","RECREATE")
    ofileSample.cd()
    sysTot=systTotal(listOfSyst)
    c1=TCanvas()
    sysTot.Draw()
    sysTot.Write()
    c1.SaveAs(directory+"_sysTot.root")

    plotter(mass_nominal, mass_PU_down, mass_PU_up, "Nominal", "Down", "Up", directory+"_syst_PU", "plot_PU")
    plotter(mass_nominal, mass_ProbQ_down, mass_ProbQ_up, "Nominal", "Down", "Up", directory+"_syst_ProbQ", "plot_ProbQ")
    plotter(mass_nominal, mass_Ias_down, mass_Ias_up, "Nominal", "Down", "Up", directory+"_syst_Ias", "plot_Ias")
    plotter(mass_nominal, mass_Pt_down, mass_Pt_up, "Nominal", "Down", "Up", directory+"_syst_Pt", "plot_Pt")
    plotter(mass_nominal, mass_Trigger_down, mass_Trigger_up, "Nominal", "Down", "Up", directory+"_syst_Trigger", "plot_Trigger")
    plotter(mass_nominal, mass_K_down1, mass_K_up1, "Nominal", "Down", "Up", directory+"_syst_K1", "plot_K1")
    plotter(mass_nominal, mass_C_down1, mass_C_up1, "Nominal", "Down", "Up", directory+"_syst_C1", "plot_C1")
    plotter(mass_nominal, mass_K_down2, mass_K_up2, "Nominal", "Down", "Up", directory+"_syst_K2", "plot_K2")
    plotter(mass_nominal, mass_C_down2, mass_C_up2, "Nominal", "Down", "Up", directory+"_syst_C2", "plot_C2")

    outTitle=sampleName.get(sample)
    plotSummary(syst_K1,syst_C1,syst_PU,syst_ProbQ,syst_Ias,syst_Pt,syst_Trigger,sysTot,"Mass bin",outTitle,directory)

    ofileAllSignal.cd()
    mass_nominal.SetName(sampleName.get(sample))
    mass_nominal.Write()
    massd=massDown(mass_nominal,sysTot)
    massd.Write()
    massu=massUp(mass_nominal,sysTot)
    massu.Write()
    
ofileAllSignal.Close()
