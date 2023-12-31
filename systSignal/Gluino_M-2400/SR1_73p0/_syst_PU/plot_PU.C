void plot_PU()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:02 2023) by ROOT version 6.14/09
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,800);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.15);
   c1->SetBottomMargin(0.14);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: t1
   TPad *t1 = new TPad("t1", "t1",0,0.4,0.95,0.9);
   t1->Draw();
   t1->cd();
   t1->Range(-428.5714,-4.326304,2428.571,0.7285576);
   t1->SetFillColor(0);
   t1->SetBorderMode(0);
   t1->SetBorderSize(2);
   t1->SetLogy();
   t1->SetGridx();
   t1->SetGridy();
   t1->SetTickx(1);
   t1->SetTicky(1);
   t1->SetLeftMargin(0.15);
   t1->SetRightMargin(0.15);
   t1->SetTopMargin(0.005);
   t1->SetBottomMargin(0.005);
   t1->SetFrameFillStyle(0);
   t1->SetFrameBorderMode(0);
   t1->SetFrameFillStyle(0);
   t1->SetFrameBorderMode(0);
   Double_t xAxis379[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__487 = new TH1F("nominal__487","",32, xAxis379);
   nominal__487->SetBinContent(12,0.0001509822);
   nominal__487->SetBinContent(22,0.0001332393);
   nominal__487->SetBinContent(23,0.0001202936);
   nominal__487->SetBinContent(24,0.0001385778);
   nominal__487->SetBinContent(26,0.0001376264);
   nominal__487->SetBinContent(27,0.0002742345);
   nominal__487->SetBinContent(28,0.001101506);
   nominal__487->SetBinContent(29,0.003481095);
   nominal__487->SetBinContent(30,0.01318569);
   nominal__487->SetBinContent(31,0.0642131);
   nominal__487->SetBinContent(32,0.5049906);
   nominal__487->SetBinError(12,0.0001509822);
   nominal__487->SetBinError(22,0.0001332393);
   nominal__487->SetBinError(23,0.0001202936);
   nominal__487->SetBinError(24,0.0001385778);
   nominal__487->SetBinError(26,0.0001376264);
   nominal__487->SetBinError(27,0.0001939144);
   nominal__487->SetBinError(28,0.0003903977);
   nominal__487->SetBinError(29,0.0006969875);
   nominal__487->SetBinError(30,0.001347568);
   nominal__487->SetBinError(31,0.002977733);
   nominal__487->SetBinError(32,0.003673154);
   nominal__487->SetBinError(33,0.0075031);
   nominal__487->SetMinimum(5e-05);
   nominal__487->SetMaximum(5.049906);
   nominal__487->SetEntries(4276);
   nominal__487->SetFillColor(1);
   nominal__487->SetMarkerStyle(20);
   nominal__487->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__487->GetXaxis()->SetRange(1,32);
   nominal__487->GetXaxis()->SetLabelFont(42);
   nominal__487->GetXaxis()->SetLabelSize(0.035);
   nominal__487->GetXaxis()->SetTitleSize(0.035);
   nominal__487->GetXaxis()->SetTitleFont(42);
   nominal__487->GetYaxis()->SetTitle("Tracks");
   nominal__487->GetYaxis()->SetLabelFont(42);
   nominal__487->GetYaxis()->SetLabelSize(0.05);
   nominal__487->GetYaxis()->SetTitleSize(0.07);
   nominal__487->GetYaxis()->SetTitleOffset(0);
   nominal__487->GetYaxis()->SetTitleFont(42);
   nominal__487->GetZaxis()->SetLabelFont(42);
   nominal__487->GetZaxis()->SetLabelSize(0.035);
   nominal__487->GetZaxis()->SetTitleSize(0.035);
   nominal__487->GetZaxis()->SetTitleFont(42);
   nominal__487->Draw("");
   Double_t xAxis380[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__488 = new TH1F("PU_down__488","",32, xAxis380);
   PU_down__488->SetBinContent(12,0.0002016105);
   PU_down__488->SetBinContent(22,0.000114438);
   PU_down__488->SetBinContent(23,9.621504e-05);
   PU_down__488->SetBinContent(24,0.0001303931);
   PU_down__488->SetBinContent(26,0.0001276535);
   PU_down__488->SetBinContent(27,0.000249856);
   PU_down__488->SetBinContent(28,0.001115194);
   PU_down__488->SetBinContent(29,0.003515749);
   PU_down__488->SetBinContent(30,0.01287071);
   PU_down__488->SetBinContent(31,0.06371346);
   PU_down__488->SetBinContent(32,0.5022117);
   PU_down__488->SetBinError(12,0.0002016105);
   PU_down__488->SetBinError(22,0.000114438);
   PU_down__488->SetBinError(23,9.621504e-05);
   PU_down__488->SetBinError(24,0.0001303931);
   PU_down__488->SetBinError(26,0.0001276535);
   PU_down__488->SetBinError(27,0.0001767169);
   PU_down__488->SetBinError(28,0.0004083958);
   PU_down__488->SetBinError(29,0.0007153473);
   PU_down__488->SetBinError(30,0.001344883);
   PU_down__488->SetBinError(31,0.003041555);
   PU_down__488->SetBinError(32,0.00370506);
   PU_down__488->SetBinError(33,0.007710356);
   PU_down__488->SetEntries(4276);
   PU_down__488->SetFillColor(38);
   PU_down__488->SetLineColor(38);
   PU_down__488->SetMarkerColor(38);
   PU_down__488->SetMarkerStyle(21);
   PU_down__488->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__488->GetXaxis()->SetRange(1,400);
   PU_down__488->GetXaxis()->SetLabelFont(42);
   PU_down__488->GetXaxis()->SetLabelSize(0.035);
   PU_down__488->GetXaxis()->SetTitleSize(0.035);
   PU_down__488->GetXaxis()->SetTitleFont(42);
   PU_down__488->GetYaxis()->SetTitle("Events / bin");
   PU_down__488->GetYaxis()->SetLabelFont(42);
   PU_down__488->GetYaxis()->SetLabelSize(0.035);
   PU_down__488->GetYaxis()->SetTitleSize(0.035);
   PU_down__488->GetYaxis()->SetTitleOffset(0);
   PU_down__488->GetYaxis()->SetTitleFont(42);
   PU_down__488->GetZaxis()->SetLabelFont(42);
   PU_down__488->GetZaxis()->SetLabelSize(0.035);
   PU_down__488->GetZaxis()->SetTitleSize(0.035);
   PU_down__488->GetZaxis()->SetTitleFont(42);
   PU_down__488->Draw("same");
   Double_t xAxis381[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__489 = new TH1F("PU_up__489","",32, xAxis381);
   PU_up__489->SetBinContent(12,9.919715e-05);
   PU_up__489->SetBinContent(22,0.0001527162);
   PU_up__489->SetBinContent(23,0.0001512105);
   PU_up__489->SetBinContent(24,0.0001453477);
   PU_up__489->SetBinContent(26,0.0001481521);
   PU_up__489->SetBinContent(27,0.0002990724);
   PU_up__489->SetBinContent(28,0.001088385);
   PU_up__489->SetBinContent(29,0.003394963);
   PU_up__489->SetBinContent(30,0.01344951);
   PU_up__489->SetBinContent(31,0.06458578);
   PU_up__489->SetBinContent(32,0.507418);
   PU_up__489->SetBinError(12,9.919716e-05);
   PU_up__489->SetBinError(22,0.0001527162);
   PU_up__489->SetBinError(23,0.0001512105);
   PU_up__489->SetBinError(24,0.0001453477);
   PU_up__489->SetBinError(26,0.0001481521);
   PU_up__489->SetBinError(27,0.0002114852);
   PU_up__489->SetBinError(28,0.0003885389);
   PU_up__489->SetBinError(29,0.0006843836);
   PU_up__489->SetBinError(30,0.001380919);
   PU_up__489->SetBinError(31,0.003012379);
   PU_up__489->SetBinError(32,0.003741702);
   PU_up__489->SetBinError(33,0.007571149);
   PU_up__489->SetEntries(4276);
   PU_up__489->SetFillColor(46);
   PU_up__489->SetLineColor(46);
   PU_up__489->SetMarkerColor(46);
   PU_up__489->SetMarkerStyle(21);
   PU_up__489->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__489->GetXaxis()->SetRange(1,400);
   PU_up__489->GetXaxis()->SetLabelFont(42);
   PU_up__489->GetXaxis()->SetLabelSize(0.035);
   PU_up__489->GetXaxis()->SetTitleSize(0.035);
   PU_up__489->GetXaxis()->SetTitleFont(42);
   PU_up__489->GetYaxis()->SetTitle("Events / bin");
   PU_up__489->GetYaxis()->SetLabelFont(42);
   PU_up__489->GetYaxis()->SetLabelSize(0.035);
   PU_up__489->GetYaxis()->SetTitleSize(0.035);
   PU_up__489->GetYaxis()->SetTitleOffset(0);
   PU_up__489->GetYaxis()->SetTitleFont(42);
   PU_up__489->GetZaxis()->SetLabelFont(42);
   PU_up__489->GetZaxis()->SetLabelSize(0.035);
   PU_up__489->GetZaxis()->SetTitleSize(0.035);
   PU_up__489->GetZaxis()->SetTitleFont(42);
   PU_up__489->Draw("same");
   TLine *line = new TLine(1730,0,1730,5.049906);
   line->Draw();
   
   TLegend *leg = new TLegend(0.2,0.65,0.4,0.9,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("nominal","Nominal","PE1");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("PU_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("PU_up","Up","PE1");
   entry->SetLineColor(46);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(46);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   t1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: t2
   TPad *t2 = new TPad("t2", "t2",0,0.225,0.95,0.375);
   t2->Draw();
   t2->cd();
   t2->Range(-428.5714,0.4949495,2428.571,1.505051);
   t2->SetFillColor(0);
   t2->SetBorderMode(0);
   t2->SetBorderSize(2);
   t2->SetGridy();
   t2->SetTickx(1);
   t2->SetTicky(1);
   t2->SetLeftMargin(0.15);
   t2->SetRightMargin(0.15);
   t2->SetTopMargin(0.005);
   t2->SetBottomMargin(0.005);
   t2->SetFrameFillStyle(0);
   t2->SetFrameBorderMode(0);
   t2->SetFrameFillStyle(0);
   t2->SetFrameBorderMode(0);
   
   TH1D *frameR2__490 = new TH1D("frameR2__490","",1,0,2000);
   frameR2__490->SetMinimum(0.5);
   frameR2__490->SetMaximum(1.5);
   frameR2__490->SetStats(0);
   frameR2__490->SetLineStyle(0);
   frameR2__490->SetMarkerStyle(20);
   frameR2__490->GetXaxis()->SetRange(1,1);
   frameR2__490->GetXaxis()->SetLabelFont(43);
   frameR2__490->GetXaxis()->SetLabelOffset(0.007);
   frameR2__490->GetXaxis()->SetLabelSize(16);
   frameR2__490->GetXaxis()->SetTitleSize(24);
   frameR2__490->GetXaxis()->SetTitleOffset(3.75);
   frameR2__490->GetXaxis()->SetTitleFont(43);
   frameR2__490->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__490->GetYaxis()->SetNdivisions(205);
   frameR2__490->GetYaxis()->SetLabelFont(43);
   frameR2__490->GetYaxis()->SetLabelOffset(0.007);
   frameR2__490->GetYaxis()->SetLabelSize(20);
   frameR2__490->GetYaxis()->SetTitleSize(20);
   frameR2__490->GetYaxis()->SetTitleOffset(2);
   frameR2__490->GetYaxis()->SetTitleFont(43);
   frameR2__490->GetZaxis()->SetLabelFont(42);
   frameR2__490->GetZaxis()->SetLabelOffset(0.007);
   frameR2__490->GetZaxis()->SetLabelSize(0.05);
   frameR2__490->GetZaxis()->SetTitleSize(0.06);
   frameR2__490->GetZaxis()->SetTitleFont(42);
   frameR2__490->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis382[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__491 = new TH1F("PU_down__491","",32, xAxis382);
   PU_down__491->SetBinContent(0,1.006126);
   PU_down__491->SetBinContent(1,1.006126);
   PU_down__491->SetBinContent(2,1.006126);
   PU_down__491->SetBinContent(3,1.006126);
   PU_down__491->SetBinContent(4,1.006126);
   PU_down__491->SetBinContent(5,1.006126);
   PU_down__491->SetBinContent(6,1.006126);
   PU_down__491->SetBinContent(7,1.006126);
   PU_down__491->SetBinContent(8,1.006126);
   PU_down__491->SetBinContent(9,1.006126);
   PU_down__491->SetBinContent(10,1.006126);
   PU_down__491->SetBinContent(11,1.006126);
   PU_down__491->SetBinContent(12,1.006126);
   PU_down__491->SetBinContent(13,1.006215);
   PU_down__491->SetBinContent(14,1.006215);
   PU_down__491->SetBinContent(15,1.006215);
   PU_down__491->SetBinContent(16,1.006215);
   PU_down__491->SetBinContent(17,1.006215);
   PU_down__491->SetBinContent(18,1.006215);
   PU_down__491->SetBinContent(19,1.006215);
   PU_down__491->SetBinContent(20,1.006215);
   PU_down__491->SetBinContent(21,1.006215);
   PU_down__491->SetBinContent(22,1.006215);
   PU_down__491->SetBinContent(23,1.006184);
   PU_down__491->SetBinContent(24,1.006144);
   PU_down__491->SetBinContent(25,1.006131);
   PU_down__491->SetBinContent(26,1.006131);
   PU_down__491->SetBinContent(27,1.006116);
   PU_down__491->SetBinContent(28,1.006076);
   PU_down__491->SetBinContent(29,1.006112);
   PU_down__491->SetBinContent(30,1.006209);
   PU_down__491->SetBinContent(31,1.005793);
   PU_down__491->SetBinContent(32,1.005533);
   PU_down__491->SetBinError(0,0.02213676);
   PU_down__491->SetBinError(1,0.02213676);
   PU_down__491->SetBinError(2,0.02213676);
   PU_down__491->SetBinError(3,0.02213676);
   PU_down__491->SetBinError(4,0.02213676);
   PU_down__491->SetBinError(5,0.02213676);
   PU_down__491->SetBinError(6,0.02213676);
   PU_down__491->SetBinError(7,0.02213676);
   PU_down__491->SetBinError(8,0.02213676);
   PU_down__491->SetBinError(9,0.02213676);
   PU_down__491->SetBinError(10,0.02213676);
   PU_down__491->SetBinError(11,0.02213676);
   PU_down__491->SetBinError(12,0.02213676);
   PU_down__491->SetBinError(13,0.02214117);
   PU_down__491->SetBinError(14,0.02214117);
   PU_down__491->SetBinError(15,0.02214117);
   PU_down__491->SetBinError(16,0.02214117);
   PU_down__491->SetBinError(17,0.02214117);
   PU_down__491->SetBinError(18,0.02214117);
   PU_down__491->SetBinError(19,0.02214117);
   PU_down__491->SetBinError(20,0.02214117);
   PU_down__491->SetBinError(21,0.02214117);
   PU_down__491->SetBinError(22,0.02214117);
   PU_down__491->SetBinError(23,0.02214311);
   PU_down__491->SetBinError(24,0.02214472);
   PU_down__491->SetBinError(25,0.02214711);
   PU_down__491->SetBinError(26,0.02214711);
   PU_down__491->SetBinError(27,0.02214943);
   PU_down__491->SetBinError(28,0.02215388);
   PU_down__491->SetBinError(29,0.02217536);
   PU_down__491->SetBinError(30,0.02224391);
   PU_down__491->SetBinError(31,0.02249334);
   PU_down__491->SetBinError(32,0.02387585);
   PU_down__491->SetEntries(33);
   PU_down__491->SetFillColor(38);
   PU_down__491->SetLineColor(38);
   PU_down__491->SetMarkerColor(38);
   PU_down__491->SetMarkerStyle(21);
   PU_down__491->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__491->GetXaxis()->SetRange(1,400);
   PU_down__491->GetXaxis()->SetLabelFont(42);
   PU_down__491->GetXaxis()->SetLabelSize(0.035);
   PU_down__491->GetXaxis()->SetTitleSize(0.035);
   PU_down__491->GetXaxis()->SetTitleFont(42);
   PU_down__491->GetYaxis()->SetTitle("Events / bin");
   PU_down__491->GetYaxis()->SetLabelFont(42);
   PU_down__491->GetYaxis()->SetLabelSize(0.035);
   PU_down__491->GetYaxis()->SetTitleSize(0.035);
   PU_down__491->GetYaxis()->SetTitleOffset(0);
   PU_down__491->GetYaxis()->SetTitleFont(42);
   PU_down__491->GetZaxis()->SetLabelFont(42);
   PU_down__491->GetZaxis()->SetLabelSize(0.035);
   PU_down__491->GetZaxis()->SetTitleSize(0.035);
   PU_down__491->GetZaxis()->SetTitleFont(42);
   PU_down__491->Draw("E0 same");
   Double_t xAxis383[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__492 = new TH1F("PU_up__492","",32, xAxis383);
   PU_up__492->SetBinContent(0,0.9949142);
   PU_up__492->SetBinContent(1,0.9949142);
   PU_up__492->SetBinContent(2,0.9949142);
   PU_up__492->SetBinContent(3,0.9949142);
   PU_up__492->SetBinContent(4,0.9949142);
   PU_up__492->SetBinContent(5,0.9949142);
   PU_up__492->SetBinContent(6,0.9949142);
   PU_up__492->SetBinContent(7,0.9949142);
   PU_up__492->SetBinContent(8,0.9949142);
   PU_up__492->SetBinContent(9,0.9949142);
   PU_up__492->SetBinContent(10,0.9949142);
   PU_up__492->SetBinContent(11,0.9949142);
   PU_up__492->SetBinContent(12,0.9949142);
   PU_up__492->SetBinContent(13,0.9948257);
   PU_up__492->SetBinContent(14,0.9948257);
   PU_up__492->SetBinContent(15,0.9948257);
   PU_up__492->SetBinContent(16,0.9948257);
   PU_up__492->SetBinContent(17,0.9948257);
   PU_up__492->SetBinContent(18,0.9948257);
   PU_up__492->SetBinContent(19,0.9948257);
   PU_up__492->SetBinContent(20,0.9948257);
   PU_up__492->SetBinContent(21,0.9948257);
   PU_up__492->SetBinContent(22,0.9948257);
   PU_up__492->SetBinContent(23,0.9948574);
   PU_up__492->SetBinContent(24,0.9949084);
   PU_up__492->SetBinContent(25,0.9949186);
   PU_up__492->SetBinContent(26,0.9949186);
   PU_up__492->SetBinContent(27,0.9949352);
   PU_up__492->SetBinContent(28,0.9949747);
   PU_up__492->SetBinContent(29,0.9949431);
   PU_up__492->SetBinContent(30,0.9947667);
   PU_up__492->SetBinContent(31,0.9951049);
   PU_up__492->SetBinContent(32,0.9952163);
   PU_up__492->SetBinError(0,0.02163587);
   PU_up__492->SetBinError(1,0.02163587);
   PU_up__492->SetBinError(2,0.02163587);
   PU_up__492->SetBinError(3,0.02163587);
   PU_up__492->SetBinError(4,0.02163587);
   PU_up__492->SetBinError(5,0.02163587);
   PU_up__492->SetBinError(6,0.02163587);
   PU_up__492->SetBinError(7,0.02163587);
   PU_up__492->SetBinError(8,0.02163587);
   PU_up__492->SetBinError(9,0.02163587);
   PU_up__492->SetBinError(10,0.02163587);
   PU_up__492->SetBinError(11,0.02163587);
   PU_up__492->SetBinError(12,0.02163587);
   PU_up__492->SetBinError(13,0.02163638);
   PU_up__492->SetBinError(14,0.02163638);
   PU_up__492->SetBinError(15,0.02163638);
   PU_up__492->SetBinError(16,0.02163638);
   PU_up__492->SetBinError(17,0.02163638);
   PU_up__492->SetBinError(18,0.02163638);
   PU_up__492->SetBinError(19,0.02163638);
   PU_up__492->SetBinError(20,0.02163638);
   PU_up__492->SetBinError(21,0.02163638);
   PU_up__492->SetBinError(22,0.02163638);
   PU_up__492->SetBinError(23,0.02163961);
   PU_up__492->SetBinError(24,0.02164326);
   PU_up__492->SetBinError(25,0.02164604);
   PU_up__492->SetBinError(26,0.02164604);
   PU_up__492->SetBinError(27,0.02164896);
   PU_up__492->SetBinError(28,0.02165493);
   PU_up__492->SetBinError(29,0.02167454);
   PU_up__492->SetBinError(30,0.02173466);
   PU_up__492->SetBinError(31,0.02199362);
   PU_up__492->SetBinError(32,0.02335419);
   PU_up__492->SetEntries(33);
   PU_up__492->SetFillColor(46);
   PU_up__492->SetLineColor(46);
   PU_up__492->SetMarkerColor(46);
   PU_up__492->SetMarkerStyle(21);
   PU_up__492->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__492->GetXaxis()->SetRange(1,400);
   PU_up__492->GetXaxis()->SetLabelFont(42);
   PU_up__492->GetXaxis()->SetLabelSize(0.035);
   PU_up__492->GetXaxis()->SetTitleSize(0.035);
   PU_up__492->GetXaxis()->SetTitleFont(42);
   PU_up__492->GetYaxis()->SetTitle("Events / bin");
   PU_up__492->GetYaxis()->SetLabelFont(42);
   PU_up__492->GetYaxis()->SetLabelSize(0.035);
   PU_up__492->GetYaxis()->SetTitleSize(0.035);
   PU_up__492->GetYaxis()->SetTitleOffset(0);
   PU_up__492->GetYaxis()->SetTitleFont(42);
   PU_up__492->GetZaxis()->SetLabelFont(42);
   PU_up__492->GetZaxis()->SetLabelSize(0.035);
   PU_up__492->GetZaxis()->SetTitleSize(0.035);
   PU_up__492->GetZaxis()->SetTitleFont(42);
   PU_up__492->Draw("E0 same");
   t2->Modified();
   c1->cd();
  
// ------------>Primitives in pad: t3
   TPad *t3 = new TPad("t3", "t3",0,0,0.95,0.2);
   t3->Draw();
   t3->cd();
   t3->Range(-428.5714,-0.1722689,2428.571,1.508403);
   t3->SetFillColor(0);
   t3->SetBorderMode(0);
   t3->SetBorderSize(2);
   t3->SetGridy();
   t3->SetTickx(1);
   t3->SetTicky(1);
   t3->SetLeftMargin(0.15);
   t3->SetRightMargin(0.15);
   t3->SetTopMargin(0.005);
   t3->SetBottomMargin(0.4);
   t3->SetFrameFillStyle(0);
   t3->SetFrameBorderMode(0);
   t3->SetFrameFillStyle(0);
   t3->SetFrameBorderMode(0);
   
   TH1D *frameR2__493 = new TH1D("frameR2__493","",1,0,2000);
   frameR2__493->SetMinimum(0.5);
   frameR2__493->SetMaximum(1.5);
   frameR2__493->SetStats(0);
   frameR2__493->SetLineStyle(0);
   frameR2__493->SetMarkerStyle(20);
   frameR2__493->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__493->GetXaxis()->SetRange(1,1);
   frameR2__493->GetXaxis()->SetLabelFont(43);
   frameR2__493->GetXaxis()->SetLabelOffset(0.007);
   frameR2__493->GetXaxis()->SetLabelSize(16);
   frameR2__493->GetXaxis()->SetTitleSize(24);
   frameR2__493->GetXaxis()->SetTitleOffset(5);
   frameR2__493->GetXaxis()->SetTitleFont(43);
   frameR2__493->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__493->GetYaxis()->SetNdivisions(205);
   frameR2__493->GetYaxis()->SetLabelFont(43);
   frameR2__493->GetYaxis()->SetLabelOffset(0.007);
   frameR2__493->GetYaxis()->SetLabelSize(20);
   frameR2__493->GetYaxis()->SetTitleSize(20);
   frameR2__493->GetYaxis()->SetTitleOffset(2);
   frameR2__493->GetYaxis()->SetTitleFont(43);
   frameR2__493->GetZaxis()->SetLabelFont(42);
   frameR2__493->GetZaxis()->SetLabelOffset(0.007);
   frameR2__493->GetZaxis()->SetLabelSize(0.05);
   frameR2__493->GetZaxis()->SetTitleSize(0.06);
   frameR2__493->GetZaxis()->SetTitleFont(42);
   frameR2__493->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis384[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__494 = new TH1F("nominal__494","",32, xAxis384);
   nominal__494->SetBinContent(12,0.7488809);
   nominal__494->SetBinContent(22,1.164293);
   nominal__494->SetBinContent(23,1.250258);
   nominal__494->SetBinContent(24,1.06277);
   nominal__494->SetBinContent(26,1.078125);
   nominal__494->SetBinContent(27,1.09757);
   nominal__494->SetBinContent(28,0.9877263);
   nominal__494->SetBinContent(29,0.9901431);
   nominal__494->SetBinContent(30,1.024472);
   nominal__494->SetBinContent(31,1.007842);
   nominal__494->SetBinContent(32,1.005533);
   nominal__494->SetBinError(12,1.059077);
   nominal__494->SetBinError(22,1.646559);
   nominal__494->SetBinError(23,1.768132);
   nominal__494->SetBinError(24,1.502983);
   nominal__494->SetBinError(26,1.524698);
   nominal__494->SetBinError(27,1.097705);
   nominal__494->SetBinError(28,0.5033771);
   nominal__494->SetBinError(29,0.2826475);
   nominal__494->SetBinError(30,0.1497386);
   nominal__494->SetBinError(31,0.06707523);
   nominal__494->SetBinError(32,0.01041754);
   nominal__494->SetMinimum(5e-05);
   nominal__494->SetMaximum(5.049906);
   nominal__494->SetEntries(9.945202);
   nominal__494->SetFillColor(38);
   nominal__494->SetLineColor(38);
   nominal__494->SetMarkerColor(38);
   nominal__494->SetMarkerStyle(21);
   nominal__494->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__494->GetXaxis()->SetRange(1,32);
   nominal__494->GetXaxis()->SetLabelFont(42);
   nominal__494->GetXaxis()->SetLabelSize(0.035);
   nominal__494->GetXaxis()->SetTitleSize(0.035);
   nominal__494->GetXaxis()->SetTitleFont(42);
   nominal__494->GetYaxis()->SetTitle("Tracks");
   nominal__494->GetYaxis()->SetLabelFont(42);
   nominal__494->GetYaxis()->SetLabelSize(0.05);
   nominal__494->GetYaxis()->SetTitleSize(0.07);
   nominal__494->GetYaxis()->SetTitleOffset(0);
   nominal__494->GetYaxis()->SetTitleFont(42);
   nominal__494->GetZaxis()->SetLabelFont(42);
   nominal__494->GetZaxis()->SetLabelSize(0.035);
   nominal__494->GetZaxis()->SetTitleSize(0.035);
   nominal__494->GetZaxis()->SetTitleFont(42);
   nominal__494->Draw("E0 same");
   Double_t xAxis385[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__495 = new TH1F("nominal__495","",32, xAxis385);
   nominal__495->SetBinContent(12,1.522042);
   nominal__495->SetBinContent(22,0.8724636);
   nominal__495->SetBinContent(23,0.7955373);
   nominal__495->SetBinContent(24,0.953423);
   nominal__495->SetBinContent(26,0.9289529);
   nominal__495->SetBinContent(27,0.9169502);
   nominal__495->SetBinContent(28,1.012056);
   nominal__495->SetBinContent(29,1.02537);
   nominal__495->SetBinContent(30,0.9803844);
   nominal__495->SetBinContent(31,0.9942296);
   nominal__495->SetBinContent(32,0.9952163);
   nominal__495->SetBinError(12,2.152492);
   nominal__495->SetBinError(22,1.23385);
   nominal__495->SetBinError(23,1.12506);
   nominal__495->SetBinError(24,1.348344);
   nominal__495->SetBinError(26,1.313738);
   nominal__495->SetBinError(27,0.916973);
   nominal__495->SetBinError(28,0.5091093);
   nominal__495->SetBinError(29,0.2913317);
   nominal__495->SetBinError(30,0.1420263);
   nominal__495->SetBinError(31,0.06539171);
   nominal__495->SetBinError(32,0.01030819);
   nominal__495->SetMinimum(5e-05);
   nominal__495->SetMaximum(5.049906);
   nominal__495->SetEntries(9.932493);
   nominal__495->SetFillColor(46);
   nominal__495->SetLineColor(46);
   nominal__495->SetMarkerColor(46);
   nominal__495->SetMarkerStyle(21);
   nominal__495->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__495->GetXaxis()->SetRange(1,32);
   nominal__495->GetXaxis()->SetLabelFont(42);
   nominal__495->GetXaxis()->SetLabelSize(0.035);
   nominal__495->GetXaxis()->SetTitleSize(0.035);
   nominal__495->GetXaxis()->SetTitleFont(42);
   nominal__495->GetYaxis()->SetTitle("Tracks");
   nominal__495->GetYaxis()->SetLabelFont(42);
   nominal__495->GetYaxis()->SetLabelSize(0.05);
   nominal__495->GetYaxis()->SetTitleSize(0.07);
   nominal__495->GetYaxis()->SetTitleOffset(0);
   nominal__495->GetYaxis()->SetTitleFont(42);
   nominal__495->GetZaxis()->SetLabelFont(42);
   nominal__495->GetZaxis()->SetLabelSize(0.035);
   nominal__495->GetZaxis()->SetTitleSize(0.035);
   nominal__495->GetZaxis()->SetTitleFont(42);
   nominal__495->Draw("E0 same");
   t3->Modified();
   c1->cd();
   TLatex *   tex = new TLatex(0.85,0.92,"#scale[0.85]{101 fb^{-1} (13 TeV)}");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.92,"Private work");
tex->SetNDC();
   tex->SetTextFont(52);
   tex->SetTextSize(0.057);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
