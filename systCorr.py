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
CMS_lumi.lumi_13TeV = "2018 - 59.7 fb^{-1}"
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "In progress"
#CMS_lumi.cmsText=""
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

def pull2D(h1,h2):
    res=h1.Clone()
    res.Reset()
    res.SetName("pull2D")
    for i in range (0,h1.GetNbinsX()+1):
        for j in range (0,h1.GetNbinsY()+1):
            a1=h1.GetBinContent(i,j)
            a2=h2.GetBinContent(i,j)
            err1=h1.GetBinError(i,j)
            if(a1<=0 or a2<=0 or err1<=0):
                continue
            res.SetBinContent(i,j,(a1-a2)/err1)
    return res

def diffRel2D(h1,h2):
    res=h1.Clone()
    res.Reset()
    res.SetName("diffRel2D")
    for i in range (0,h1.GetNbinsX()+1):
        for j in range (0,h1.GetNbinsY()+1):
            a1=h1.GetBinContent(i,j)
            if(a1<=0):
                continue
            a2=h2.GetBinContent(i,j)
            res.SetBinContent(i,j,(a1-a2)/a1)
    return res

def removeDeDxBand(h):
    for i in range (0,h.GetNbinsX()+1):
        for j in range (0,h.GetNbinsY()+1):
            if (h.GetYaxis().GetBinCenter(j) < 3.1) :
                h.SetBinContent(i,j,0)
    return h

inputNominal="/opt/sbg/cms/safe1/cms/dapparu/HSCP_annexe/CMSSW_10_6_30/src/crab_Analysis_SingleMuon_Run2018_CodeV55p0_v1_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE1.root"
ifileNominal=TFile(inputNominal)

oDir="systematic_correlation_dir/"
cmd='mkdir -p '+oDir
os.system(cmd)

c1=TCanvas("c1","c1",1100,700)

cross1D=ifileNominal.Get("cross1D_True_BC")
histo2D=ifileNominal.Get("ih_p_regionD_50ias90")

cross1D=removeDeDxBand(cross1D)
histo2D=removeDeDxBand(histo2D)

cross1D.Scale(1./cross1D.Integral())
histo2D.Scale(1./histo2D.Integral())


histo2DCL=histo2D.Clone()
histo2DCL.Rebin2D(25,2)
histo2DCL.GetYaxis().SetRange(3,5)
histo2DCL.GetYaxis().SetRangeUser(2.9,5)
histo2DCL.GetXaxis().SetTitle("10^{-4}/p [GeV^{-1}]")
histo2DCL.GetXaxis().SetRangeUser(0,180)
histo2DCL.Draw("colz,text")
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
c1.SaveAs(oDir+"histo2D.pdf")
c1.SaveAs(oDir+"histo2D.root")

cross1DCL=cross1D.Clone()
cross1DCL.Rebin2D(25,2)
cross1DCL.GetYaxis().SetRange(3,5)
cross1DCL.GetYaxis().SetRangeUser(2.9,5)
cross1DCL.GetXaxis().SetTitle("10^{-4}/p [GeV^{-1}]")
cross1DCL.GetXaxis().SetRangeUser(0,180)
cross1DCL.Draw("colz,text")
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
c1.SaveAs(oDir+"cross1D.pdf")
c1.SaveAs(oDir+"cross1D.root")


cross1D.Rebin2D(25,2)
histo2D.Rebin2D(25,2)

ofile=TFile(oDir+"systCorrelation_noCorr_regD.root","RECREATE")
h_pull2D=pull2D(histo2D,cross1D)
h_diffRel2D=diffRel2D(histo2D,cross1D)


h_pull2D.GetYaxis().SetRangeUser(2.9,4.5)
h_pull2D.GetXaxis().SetTitle("10^{-4}/p [GeV^{-1}]")
h_diffRel2D.GetYaxis().SetRangeUser(2.9,5)
h_diffRel2D.GetXaxis().SetTitle("10^{-4}/p [GeV^{-1}]")

h_pull2D.GetXaxis().SetRangeUser(0,180)
h_diffRel2D.GetXaxis().SetRangeUser(0,180)

h_pull2D.Draw("colz,text")
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
c1.SaveAs(oDir+"pull2D.pdf")
c1.SaveAs(oDir+"pull2D.root")
h_diffRel2D.Draw("colz,text")
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
c1.SaveAs(oDir+"diffRel2D.pdf")
c1.SaveAs(oDir+"diffRel2D.root")

histo2DCL.Write()
h_pull2D.Write()
h_diffRel2D.Write()
