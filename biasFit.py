#!/usr/bin/python

import sys, getopt, os
sys.argv.append( '-b-' )
import ROOT
import math
from array import array
import numpy as np

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors, TF1
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

def fitParameter(h1,i):
    f1=TF1("f1","gaus",0,8)
    f1.SetParameter(0,1)
    f1.SetParameter(1,3.3)
    f1.SetParameter(2,0.15)
    h1.Reset()
    h1.FillRandom("f1",i)
    h1.Fit("f1","Q","",0,8)
    f1=h1.GetFunction("f1")
    return f1.GetParameter(1), f1.GetParameter(2)

def fitParameterErf(h2,i):
    f2=TF1("f2","[0]*([1]+erf((log(x)-[2])/[3]))",1,70)
    f2.SetParameter(0,1)
    f2.SetParameter(1,1.00006e+00)
    f2.SetParameter(2,5.10116e+00)
    f2.SetParameter(3,1.10152e+00)
    h2.Reset()
    h2.FillRandom("f2",i)
    h2.Fit("f2","Q","",0,70)
    f2=h2.GetFunction("f2")
    return f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3)

h1=ROOT.TH1F("h1","h1",1000,0,8)
h2=ROOT.TH1F("h2","h2",1000,0,70)

#print fitParameter(h1,1000)
n, x, y, x2, y2, z2 = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

for i in range (1,5000,50):
    (x_, y_) = fitParameter(h1,i)
    (x2_, y2_, z2_) = fitParameterErf(h2,i)
    n.append(i)
    x.append((x_-3.3)/3.3)
    y.append((y_-0.15)/0.15)
    x2.append((x2_-1.00006)/1.00006)
    y2.append((y2_-5.10116)/5.10116)
    z2.append((z2_-1.10152)/1.10152)

gr1=TGraph(len(n),n,x)
gr2=TGraph(len(n),n,y)

gr3=TGraph(len(n),n,x2)
gr4=TGraph(len(n),n,y2)
gr5=TGraph(len(n),n,z2)

c1=TCanvas("","Gaussian bias test",700,700)
c1.Divide(1,2)

c1.cd(1)
gr1.Draw("AP")
gr1.GetXaxis().SetTitle("PE")
gr1.GetYaxis().SetTitle("#frac{#mu_{fit}-#mu}{#mu}")

c1.cd(2)
gr2.Draw("AP")
gr2.GetXaxis().SetTitle("PE")
gr2.GetYaxis().SetTitle("#frac{#sigma_{fit}-#sigma}{#sigma}")

c1.SaveAs("biasStudyGaus.pdf")

c2=TCanvas("","",700,700)
c2.Divide(1,3)
c2.cd(1)
gr3.Draw("AP")
gr3.GetXaxis().SetTitle("PE")
gr3.GetYaxis().SetTitle("#frac{par1_{fit}-par1}{par1}")
c2.cd(2)
gr4.Draw("AP")
gr4.GetXaxis().SetTitle("PE")
gr4.GetYaxis().SetTitle("#frac{par2_{fit}-par2}{par2}")
c2.cd(3)
gr5.Draw("AP")
gr5.GetXaxis().SetTitle("PE")
gr5.GetYaxis().SetTitle("#frac{par3_{fit}-par3}{par3}")
c2.SaveAs("biasStudyErf.pdf")