void plot_Pt()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:05 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.347128,2428.571,4.872406);
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
   Double_t xAxis715[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__919 = new TH1F("nominal__919","",32, xAxis715);
   nominal__919->SetBinContent(2,0.0004271438);
   nominal__919->SetBinContent(5,0.0004657558);
   nominal__919->SetBinContent(6,0.0004665128);
   nominal__919->SetBinContent(7,0.0004416279);
   nominal__919->SetBinContent(8,0.0004562392);
   nominal__919->SetBinContent(9,0.0004493202);
   nominal__919->SetBinContent(10,0.0004450588);
   nominal__919->SetBinContent(11,0.0004493202);
   nominal__919->SetBinContent(12,0.001827087);
   nominal__919->SetBinContent(13,0.002300727);
   nominal__919->SetBinContent(14,0.003799015);
   nominal__919->SetBinContent(15,0.004216676);
   nominal__919->SetBinContent(16,0.01508476);
   nominal__919->SetBinContent(17,0.02368568);
   nominal__919->SetBinContent(18,0.05945345);
   nominal__919->SetBinContent(19,0.1197385);
   nominal__919->SetBinContent(20,0.3419138);
   nominal__919->SetBinContent(21,0.7384515);
   nominal__919->SetBinContent(22,2.091074);
   nominal__919->SetBinContent(23,5.199684);
   nominal__919->SetBinContent(24,6.7036);
   nominal__919->SetBinContent(25,3.288148);
   nominal__919->SetBinContent(26,0.906163);
   nominal__919->SetBinContent(27,0.1376991);
   nominal__919->SetBinContent(28,0.03025302);
   nominal__919->SetBinContent(29,0.008260727);
   nominal__919->SetBinContent(30,0.002395275);
   nominal__919->SetBinContent(31,0.0009228621);
   nominal__919->SetBinContent(32,0.0009801134);
   nominal__919->SetBinError(2,0.0004271438);
   nominal__919->SetBinError(5,0.0004657558);
   nominal__919->SetBinError(6,0.0004665128);
   nominal__919->SetBinError(7,0.0004416279);
   nominal__919->SetBinError(8,0.0004562392);
   nominal__919->SetBinError(9,0.0004493202);
   nominal__919->SetBinError(10,0.0004450588);
   nominal__919->SetBinError(11,0.0004493202);
   nominal__919->SetBinError(12,0.0009148069);
   nominal__919->SetBinError(13,0.001029015);
   nominal__919->SetBinError(14,0.001344761);
   nominal__919->SetBinError(15,0.001411037);
   nominal__919->SetBinError(16,0.002629427);
   nominal__919->SetBinError(17,0.003321459);
   nominal__919->SetBinError(18,0.005243704);
   nominal__919->SetBinError(19,0.007457969);
   nominal__919->SetBinError(20,0.01259405);
   nominal__919->SetBinError(21,0.01853888);
   nominal__919->SetBinError(22,0.03117085);
   nominal__919->SetBinError(23,0.0492001);
   nominal__919->SetBinError(24,0.05589285);
   nominal__919->SetBinError(25,0.03917436);
   nominal__919->SetBinError(26,0.02059337);
   nominal__919->SetBinError(27,0.008005042);
   nominal__919->SetBinError(28,0.00378786);
   nominal__919->SetBinError(29,0.001949863);
   nominal__919->SetBinError(30,0.001073186);
   nominal__919->SetBinError(31,0.0006526152);
   nominal__919->SetBinError(32,0.0005159987);
   nominal__919->SetBinError(33,0.0004641146);
   nominal__919->SetMinimum(5e-05);
   nominal__919->SetMaximum(67036);
   nominal__919->SetEntries(42418);
   nominal__919->SetFillColor(1);
   nominal__919->SetMarkerStyle(20);
   nominal__919->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__919->GetXaxis()->SetRange(1,32);
   nominal__919->GetXaxis()->SetLabelFont(42);
   nominal__919->GetXaxis()->SetLabelSize(0.035);
   nominal__919->GetXaxis()->SetTitleSize(0.035);
   nominal__919->GetXaxis()->SetTitleFont(42);
   nominal__919->GetYaxis()->SetTitle("Tracks");
   nominal__919->GetYaxis()->SetLabelFont(42);
   nominal__919->GetYaxis()->SetLabelSize(0.05);
   nominal__919->GetYaxis()->SetTitleSize(0.07);
   nominal__919->GetYaxis()->SetTitleOffset(0);
   nominal__919->GetYaxis()->SetTitleFont(42);
   nominal__919->GetZaxis()->SetLabelFont(42);
   nominal__919->GetZaxis()->SetLabelSize(0.035);
   nominal__919->GetZaxis()->SetTitleSize(0.035);
   nominal__919->GetZaxis()->SetTitleFont(42);
   nominal__919->Draw("");
   Double_t xAxis716[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_down__920 = new TH1F("Pt_down__920","",32, xAxis716);
   Pt_down__920->SetBinContent(2,0.0004271438);
   Pt_down__920->SetBinContent(5,0.0004657558);
   Pt_down__920->SetBinContent(6,0.0004665128);
   Pt_down__920->SetBinContent(7,0.0004416279);
   Pt_down__920->SetBinContent(8,0.0004562392);
   Pt_down__920->SetBinContent(9,0.0004493202);
   Pt_down__920->SetBinContent(10,0.0004450588);
   Pt_down__920->SetBinContent(11,0.0004493202);
   Pt_down__920->SetBinContent(12,0.001827087);
   Pt_down__920->SetBinContent(13,0.002300727);
   Pt_down__920->SetBinContent(14,0.003799015);
   Pt_down__920->SetBinContent(15,0.004216676);
   Pt_down__920->SetBinContent(16,0.01508476);
   Pt_down__920->SetBinContent(17,0.02368568);
   Pt_down__920->SetBinContent(18,0.05945345);
   Pt_down__920->SetBinContent(19,0.1197385);
   Pt_down__920->SetBinContent(20,0.3419138);
   Pt_down__920->SetBinContent(21,0.7384515);
   Pt_down__920->SetBinContent(22,2.091074);
   Pt_down__920->SetBinContent(23,5.199684);
   Pt_down__920->SetBinContent(24,6.7036);
   Pt_down__920->SetBinContent(25,3.288148);
   Pt_down__920->SetBinContent(26,0.906163);
   Pt_down__920->SetBinContent(27,0.1376991);
   Pt_down__920->SetBinContent(28,0.03025302);
   Pt_down__920->SetBinContent(29,0.008260727);
   Pt_down__920->SetBinContent(30,0.002395275);
   Pt_down__920->SetBinContent(31,0.0009228621);
   Pt_down__920->SetBinContent(32,0.0009801134);
   Pt_down__920->SetBinError(2,0.0004271438);
   Pt_down__920->SetBinError(5,0.0004657558);
   Pt_down__920->SetBinError(6,0.0004665128);
   Pt_down__920->SetBinError(7,0.0004416279);
   Pt_down__920->SetBinError(8,0.0004562392);
   Pt_down__920->SetBinError(9,0.0004493202);
   Pt_down__920->SetBinError(10,0.0004450588);
   Pt_down__920->SetBinError(11,0.0004493202);
   Pt_down__920->SetBinError(12,0.0009148069);
   Pt_down__920->SetBinError(13,0.001029015);
   Pt_down__920->SetBinError(14,0.001344761);
   Pt_down__920->SetBinError(15,0.001411037);
   Pt_down__920->SetBinError(16,0.002629427);
   Pt_down__920->SetBinError(17,0.003321459);
   Pt_down__920->SetBinError(18,0.005243704);
   Pt_down__920->SetBinError(19,0.007457969);
   Pt_down__920->SetBinError(20,0.01259405);
   Pt_down__920->SetBinError(21,0.01853888);
   Pt_down__920->SetBinError(22,0.03117085);
   Pt_down__920->SetBinError(23,0.0492001);
   Pt_down__920->SetBinError(24,0.05589285);
   Pt_down__920->SetBinError(25,0.03917436);
   Pt_down__920->SetBinError(26,0.02059337);
   Pt_down__920->SetBinError(27,0.008005042);
   Pt_down__920->SetBinError(28,0.00378786);
   Pt_down__920->SetBinError(29,0.001949863);
   Pt_down__920->SetBinError(30,0.001073186);
   Pt_down__920->SetBinError(31,0.0006526152);
   Pt_down__920->SetBinError(32,0.0005159987);
   Pt_down__920->SetBinError(33,0.0004641146);
   Pt_down__920->SetEntries(42418);
   Pt_down__920->SetFillColor(38);
   Pt_down__920->SetLineColor(38);
   Pt_down__920->SetMarkerColor(38);
   Pt_down__920->SetMarkerStyle(21);
   Pt_down__920->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_down__920->GetXaxis()->SetRange(1,400);
   Pt_down__920->GetXaxis()->SetLabelFont(42);
   Pt_down__920->GetXaxis()->SetLabelSize(0.035);
   Pt_down__920->GetXaxis()->SetTitleSize(0.035);
   Pt_down__920->GetXaxis()->SetTitleFont(42);
   Pt_down__920->GetYaxis()->SetTitle("Events / bin");
   Pt_down__920->GetYaxis()->SetLabelFont(42);
   Pt_down__920->GetYaxis()->SetLabelSize(0.035);
   Pt_down__920->GetYaxis()->SetTitleSize(0.035);
   Pt_down__920->GetYaxis()->SetTitleOffset(0);
   Pt_down__920->GetYaxis()->SetTitleFont(42);
   Pt_down__920->GetZaxis()->SetLabelFont(42);
   Pt_down__920->GetZaxis()->SetLabelSize(0.035);
   Pt_down__920->GetZaxis()->SetTitleSize(0.035);
   Pt_down__920->GetZaxis()->SetTitleFont(42);
   Pt_down__920->Draw("same");
   Double_t xAxis717[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_up__921 = new TH1F("Pt_up__921","",32, xAxis717);
   Pt_up__921->SetBinContent(2,0.0004271438);
   Pt_up__921->SetBinContent(5,0.0004657558);
   Pt_up__921->SetBinContent(6,0.0004665128);
   Pt_up__921->SetBinContent(7,0.0004416279);
   Pt_up__921->SetBinContent(8,0.0004562392);
   Pt_up__921->SetBinContent(9,0.0004493202);
   Pt_up__921->SetBinContent(10,0.0004450588);
   Pt_up__921->SetBinContent(11,0.0004493202);
   Pt_up__921->SetBinContent(12,0.001827087);
   Pt_up__921->SetBinContent(13,0.002300727);
   Pt_up__921->SetBinContent(14,0.003799015);
   Pt_up__921->SetBinContent(15,0.004216676);
   Pt_up__921->SetBinContent(16,0.01508476);
   Pt_up__921->SetBinContent(17,0.02368568);
   Pt_up__921->SetBinContent(18,0.05945345);
   Pt_up__921->SetBinContent(19,0.1197385);
   Pt_up__921->SetBinContent(20,0.3419138);
   Pt_up__921->SetBinContent(21,0.7384515);
   Pt_up__921->SetBinContent(22,2.091074);
   Pt_up__921->SetBinContent(23,5.199684);
   Pt_up__921->SetBinContent(24,6.7036);
   Pt_up__921->SetBinContent(25,3.288148);
   Pt_up__921->SetBinContent(26,0.906163);
   Pt_up__921->SetBinContent(27,0.1376991);
   Pt_up__921->SetBinContent(28,0.03025302);
   Pt_up__921->SetBinContent(29,0.008260727);
   Pt_up__921->SetBinContent(30,0.002395275);
   Pt_up__921->SetBinContent(31,0.0009228621);
   Pt_up__921->SetBinContent(32,0.0009801134);
   Pt_up__921->SetBinError(2,0.0004271438);
   Pt_up__921->SetBinError(5,0.0004657558);
   Pt_up__921->SetBinError(6,0.0004665128);
   Pt_up__921->SetBinError(7,0.0004416279);
   Pt_up__921->SetBinError(8,0.0004562392);
   Pt_up__921->SetBinError(9,0.0004493202);
   Pt_up__921->SetBinError(10,0.0004450588);
   Pt_up__921->SetBinError(11,0.0004493202);
   Pt_up__921->SetBinError(12,0.0009148069);
   Pt_up__921->SetBinError(13,0.001029015);
   Pt_up__921->SetBinError(14,0.001344761);
   Pt_up__921->SetBinError(15,0.001411037);
   Pt_up__921->SetBinError(16,0.002629427);
   Pt_up__921->SetBinError(17,0.003321459);
   Pt_up__921->SetBinError(18,0.005243704);
   Pt_up__921->SetBinError(19,0.007457969);
   Pt_up__921->SetBinError(20,0.01259405);
   Pt_up__921->SetBinError(21,0.01853888);
   Pt_up__921->SetBinError(22,0.03117085);
   Pt_up__921->SetBinError(23,0.0492001);
   Pt_up__921->SetBinError(24,0.05589285);
   Pt_up__921->SetBinError(25,0.03917436);
   Pt_up__921->SetBinError(26,0.02059337);
   Pt_up__921->SetBinError(27,0.008005042);
   Pt_up__921->SetBinError(28,0.00378786);
   Pt_up__921->SetBinError(29,0.001949863);
   Pt_up__921->SetBinError(30,0.001073186);
   Pt_up__921->SetBinError(31,0.0006526152);
   Pt_up__921->SetBinError(32,0.0005159987);
   Pt_up__921->SetBinError(33,0.0004641146);
   Pt_up__921->SetEntries(42418);
   Pt_up__921->SetFillColor(46);
   Pt_up__921->SetLineColor(46);
   Pt_up__921->SetMarkerColor(46);
   Pt_up__921->SetMarkerStyle(21);
   Pt_up__921->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_up__921->GetXaxis()->SetRange(1,400);
   Pt_up__921->GetXaxis()->SetLabelFont(42);
   Pt_up__921->GetXaxis()->SetLabelSize(0.035);
   Pt_up__921->GetXaxis()->SetTitleSize(0.035);
   Pt_up__921->GetXaxis()->SetTitleFont(42);
   Pt_up__921->GetYaxis()->SetTitle("Events / bin");
   Pt_up__921->GetYaxis()->SetLabelFont(42);
   Pt_up__921->GetYaxis()->SetLabelSize(0.035);
   Pt_up__921->GetYaxis()->SetTitleSize(0.035);
   Pt_up__921->GetYaxis()->SetTitleOffset(0);
   Pt_up__921->GetYaxis()->SetTitleFont(42);
   Pt_up__921->GetZaxis()->SetLabelFont(42);
   Pt_up__921->GetZaxis()->SetLabelSize(0.035);
   Pt_up__921->GetZaxis()->SetTitleSize(0.035);
   Pt_up__921->GetZaxis()->SetTitleFont(42);
   Pt_up__921->Draw("same");
   TLine *line = new TLine(1730,0,1730,67036);
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
   entry=leg->AddEntry("Pt_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Pt_up","Up","PE1");
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
   
   TH1D *frameR2__922 = new TH1D("frameR2__922","",1,0,2000);
   frameR2__922->SetMinimum(0.5);
   frameR2__922->SetMaximum(1.5);
   frameR2__922->SetStats(0);
   frameR2__922->SetLineStyle(0);
   frameR2__922->SetMarkerStyle(20);
   frameR2__922->GetXaxis()->SetRange(1,1);
   frameR2__922->GetXaxis()->SetLabelFont(43);
   frameR2__922->GetXaxis()->SetLabelOffset(0.007);
   frameR2__922->GetXaxis()->SetLabelSize(16);
   frameR2__922->GetXaxis()->SetTitleSize(24);
   frameR2__922->GetXaxis()->SetTitleOffset(3.75);
   frameR2__922->GetXaxis()->SetTitleFont(43);
   frameR2__922->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__922->GetYaxis()->SetNdivisions(205);
   frameR2__922->GetYaxis()->SetLabelFont(43);
   frameR2__922->GetYaxis()->SetLabelOffset(0.007);
   frameR2__922->GetYaxis()->SetLabelSize(20);
   frameR2__922->GetYaxis()->SetTitleSize(20);
   frameR2__922->GetYaxis()->SetTitleOffset(2);
   frameR2__922->GetYaxis()->SetTitleFont(43);
   frameR2__922->GetZaxis()->SetLabelFont(42);
   frameR2__922->GetZaxis()->SetLabelOffset(0.007);
   frameR2__922->GetZaxis()->SetLabelSize(0.05);
   frameR2__922->GetZaxis()->SetTitleSize(0.06);
   frameR2__922->GetZaxis()->SetTitleFont(42);
   frameR2__922->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis718[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_down__923 = new TH1F("Pt_down__923","",32, xAxis718);
   Pt_down__923->SetBinContent(0,1);
   Pt_down__923->SetBinContent(1,1);
   Pt_down__923->SetBinContent(2,1);
   Pt_down__923->SetBinContent(3,1);
   Pt_down__923->SetBinContent(4,1);
   Pt_down__923->SetBinContent(5,1);
   Pt_down__923->SetBinContent(6,1);
   Pt_down__923->SetBinContent(7,1);
   Pt_down__923->SetBinContent(8,1);
   Pt_down__923->SetBinContent(9,1);
   Pt_down__923->SetBinContent(10,1);
   Pt_down__923->SetBinContent(11,1);
   Pt_down__923->SetBinContent(12,1);
   Pt_down__923->SetBinContent(13,1);
   Pt_down__923->SetBinContent(14,1);
   Pt_down__923->SetBinContent(15,1);
   Pt_down__923->SetBinContent(16,1);
   Pt_down__923->SetBinContent(17,1);
   Pt_down__923->SetBinContent(18,1);
   Pt_down__923->SetBinContent(19,1);
   Pt_down__923->SetBinContent(20,1);
   Pt_down__923->SetBinContent(21,1);
   Pt_down__923->SetBinContent(22,1);
   Pt_down__923->SetBinContent(23,1);
   Pt_down__923->SetBinContent(24,1);
   Pt_down__923->SetBinContent(25,1);
   Pt_down__923->SetBinContent(26,1);
   Pt_down__923->SetBinContent(27,1);
   Pt_down__923->SetBinContent(28,1);
   Pt_down__923->SetBinContent(29,1);
   Pt_down__923->SetBinContent(30,1);
   Pt_down__923->SetBinContent(31,1);
   Pt_down__923->SetBinContent(32,1);
   Pt_down__923->SetBinError(0,0.006880216);
   Pt_down__923->SetBinError(1,0.006880216);
   Pt_down__923->SetBinError(2,0.006880216);
   Pt_down__923->SetBinError(3,0.006880297);
   Pt_down__923->SetBinError(4,0.006880297);
   Pt_down__923->SetBinError(5,0.006880297);
   Pt_down__923->SetBinError(6,0.006880379);
   Pt_down__923->SetBinError(7,0.00688046);
   Pt_down__923->SetBinError(8,0.006880541);
   Pt_down__923->SetBinError(9,0.006880623);
   Pt_down__923->SetBinError(10,0.006880704);
   Pt_down__923->SetBinError(11,0.006880785);
   Pt_down__923->SetBinError(12,0.006880867);
   Pt_down__923->SetBinError(13,0.006881191);
   Pt_down__923->SetBinError(14,0.006881599);
   Pt_down__923->SetBinError(15,0.006882249);
   Pt_down__923->SetBinError(16,0.006882976);
   Pt_down__923->SetBinError(17,0.00688566);
   Pt_down__923->SetBinError(18,0.006889814);
   Pt_down__923->SetBinError(19,0.006900349);
   Pt_down__923->SetBinError(20,0.006921619);
   Pt_down__923->SetBinError(21,0.006983543);
   Pt_down__923->SetBinError(22,0.00712269);
   Pt_down__923->SetBinError(23,0.007567724);
   Pt_down__923->SetBinError(24,0.009176078);
   Pt_down__923->SetBinError(25,0.01461132);
   Pt_down__923->SetBinError(26,0.02934327);
   Pt_down__923->SetBinError(27,0.07192811);
   Pt_down__923->SetBinError(28,0.1484989);
   Pt_down__923->SetBinError(29,0.2726189);
   Pt_down__923->SetBinError(30,0.4721524);
   Pt_down__923->SetBinError(31,0.707979);
   Pt_down__923->SetBinError(32,1.0014);
   Pt_down__923->SetEntries(33);
   Pt_down__923->SetFillColor(38);
   Pt_down__923->SetLineColor(38);
   Pt_down__923->SetMarkerColor(38);
   Pt_down__923->SetMarkerStyle(21);
   Pt_down__923->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_down__923->GetXaxis()->SetRange(1,400);
   Pt_down__923->GetXaxis()->SetLabelFont(42);
   Pt_down__923->GetXaxis()->SetLabelSize(0.035);
   Pt_down__923->GetXaxis()->SetTitleSize(0.035);
   Pt_down__923->GetXaxis()->SetTitleFont(42);
   Pt_down__923->GetYaxis()->SetTitle("Events / bin");
   Pt_down__923->GetYaxis()->SetLabelFont(42);
   Pt_down__923->GetYaxis()->SetLabelSize(0.035);
   Pt_down__923->GetYaxis()->SetTitleSize(0.035);
   Pt_down__923->GetYaxis()->SetTitleOffset(0);
   Pt_down__923->GetYaxis()->SetTitleFont(42);
   Pt_down__923->GetZaxis()->SetLabelFont(42);
   Pt_down__923->GetZaxis()->SetLabelSize(0.035);
   Pt_down__923->GetZaxis()->SetTitleSize(0.035);
   Pt_down__923->GetZaxis()->SetTitleFont(42);
   Pt_down__923->Draw("E0 same");
   Double_t xAxis719[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_up__924 = new TH1F("Pt_up__924","",32, xAxis719);
   Pt_up__924->SetBinContent(0,1);
   Pt_up__924->SetBinContent(1,1);
   Pt_up__924->SetBinContent(2,1);
   Pt_up__924->SetBinContent(3,1);
   Pt_up__924->SetBinContent(4,1);
   Pt_up__924->SetBinContent(5,1);
   Pt_up__924->SetBinContent(6,1);
   Pt_up__924->SetBinContent(7,1);
   Pt_up__924->SetBinContent(8,1);
   Pt_up__924->SetBinContent(9,1);
   Pt_up__924->SetBinContent(10,1);
   Pt_up__924->SetBinContent(11,1);
   Pt_up__924->SetBinContent(12,1);
   Pt_up__924->SetBinContent(13,1);
   Pt_up__924->SetBinContent(14,1);
   Pt_up__924->SetBinContent(15,1);
   Pt_up__924->SetBinContent(16,1);
   Pt_up__924->SetBinContent(17,1);
   Pt_up__924->SetBinContent(18,1);
   Pt_up__924->SetBinContent(19,1);
   Pt_up__924->SetBinContent(20,1);
   Pt_up__924->SetBinContent(21,1);
   Pt_up__924->SetBinContent(22,1);
   Pt_up__924->SetBinContent(23,1);
   Pt_up__924->SetBinContent(24,1);
   Pt_up__924->SetBinContent(25,1);
   Pt_up__924->SetBinContent(26,1);
   Pt_up__924->SetBinContent(27,1);
   Pt_up__924->SetBinContent(28,1);
   Pt_up__924->SetBinContent(29,1);
   Pt_up__924->SetBinContent(30,1);
   Pt_up__924->SetBinContent(31,1);
   Pt_up__924->SetBinContent(32,1);
   Pt_up__924->SetBinError(0,0.006880216);
   Pt_up__924->SetBinError(1,0.006880216);
   Pt_up__924->SetBinError(2,0.006880216);
   Pt_up__924->SetBinError(3,0.006880297);
   Pt_up__924->SetBinError(4,0.006880297);
   Pt_up__924->SetBinError(5,0.006880297);
   Pt_up__924->SetBinError(6,0.006880379);
   Pt_up__924->SetBinError(7,0.00688046);
   Pt_up__924->SetBinError(8,0.006880541);
   Pt_up__924->SetBinError(9,0.006880623);
   Pt_up__924->SetBinError(10,0.006880704);
   Pt_up__924->SetBinError(11,0.006880785);
   Pt_up__924->SetBinError(12,0.006880867);
   Pt_up__924->SetBinError(13,0.006881191);
   Pt_up__924->SetBinError(14,0.006881599);
   Pt_up__924->SetBinError(15,0.006882249);
   Pt_up__924->SetBinError(16,0.006882976);
   Pt_up__924->SetBinError(17,0.00688566);
   Pt_up__924->SetBinError(18,0.006889814);
   Pt_up__924->SetBinError(19,0.006900349);
   Pt_up__924->SetBinError(20,0.006921619);
   Pt_up__924->SetBinError(21,0.006983543);
   Pt_up__924->SetBinError(22,0.00712269);
   Pt_up__924->SetBinError(23,0.007567724);
   Pt_up__924->SetBinError(24,0.009176078);
   Pt_up__924->SetBinError(25,0.01461132);
   Pt_up__924->SetBinError(26,0.02934327);
   Pt_up__924->SetBinError(27,0.07192811);
   Pt_up__924->SetBinError(28,0.1484989);
   Pt_up__924->SetBinError(29,0.2726189);
   Pt_up__924->SetBinError(30,0.4721524);
   Pt_up__924->SetBinError(31,0.707979);
   Pt_up__924->SetBinError(32,1.0014);
   Pt_up__924->SetEntries(33);
   Pt_up__924->SetFillColor(46);
   Pt_up__924->SetLineColor(46);
   Pt_up__924->SetMarkerColor(46);
   Pt_up__924->SetMarkerStyle(21);
   Pt_up__924->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_up__924->GetXaxis()->SetRange(1,400);
   Pt_up__924->GetXaxis()->SetLabelFont(42);
   Pt_up__924->GetXaxis()->SetLabelSize(0.035);
   Pt_up__924->GetXaxis()->SetTitleSize(0.035);
   Pt_up__924->GetXaxis()->SetTitleFont(42);
   Pt_up__924->GetYaxis()->SetTitle("Events / bin");
   Pt_up__924->GetYaxis()->SetLabelFont(42);
   Pt_up__924->GetYaxis()->SetLabelSize(0.035);
   Pt_up__924->GetYaxis()->SetTitleSize(0.035);
   Pt_up__924->GetYaxis()->SetTitleOffset(0);
   Pt_up__924->GetYaxis()->SetTitleFont(42);
   Pt_up__924->GetZaxis()->SetLabelFont(42);
   Pt_up__924->GetZaxis()->SetLabelSize(0.035);
   Pt_up__924->GetZaxis()->SetTitleSize(0.035);
   Pt_up__924->GetZaxis()->SetTitleFont(42);
   Pt_up__924->Draw("E0 same");
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
   
   TH1D *frameR2__925 = new TH1D("frameR2__925","",1,0,2000);
   frameR2__925->SetMinimum(0.5);
   frameR2__925->SetMaximum(1.5);
   frameR2__925->SetStats(0);
   frameR2__925->SetLineStyle(0);
   frameR2__925->SetMarkerStyle(20);
   frameR2__925->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__925->GetXaxis()->SetRange(1,1);
   frameR2__925->GetXaxis()->SetLabelFont(43);
   frameR2__925->GetXaxis()->SetLabelOffset(0.007);
   frameR2__925->GetXaxis()->SetLabelSize(16);
   frameR2__925->GetXaxis()->SetTitleSize(24);
   frameR2__925->GetXaxis()->SetTitleOffset(5);
   frameR2__925->GetXaxis()->SetTitleFont(43);
   frameR2__925->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__925->GetYaxis()->SetNdivisions(205);
   frameR2__925->GetYaxis()->SetLabelFont(43);
   frameR2__925->GetYaxis()->SetLabelOffset(0.007);
   frameR2__925->GetYaxis()->SetLabelSize(20);
   frameR2__925->GetYaxis()->SetTitleSize(20);
   frameR2__925->GetYaxis()->SetTitleOffset(2);
   frameR2__925->GetYaxis()->SetTitleFont(43);
   frameR2__925->GetZaxis()->SetLabelFont(42);
   frameR2__925->GetZaxis()->SetLabelOffset(0.007);
   frameR2__925->GetZaxis()->SetLabelSize(0.05);
   frameR2__925->GetZaxis()->SetTitleSize(0.06);
   frameR2__925->GetZaxis()->SetTitleFont(42);
   frameR2__925->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis720[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__926 = new TH1F("nominal__926","",32, xAxis720);
   nominal__926->SetBinContent(2,1);
   nominal__926->SetBinContent(5,1);
   nominal__926->SetBinContent(6,1);
   nominal__926->SetBinContent(7,1);
   nominal__926->SetBinContent(8,1);
   nominal__926->SetBinContent(9,1);
   nominal__926->SetBinContent(10,1);
   nominal__926->SetBinContent(11,1);
   nominal__926->SetBinContent(12,1);
   nominal__926->SetBinContent(13,1);
   nominal__926->SetBinContent(14,1);
   nominal__926->SetBinContent(15,1);
   nominal__926->SetBinContent(16,1);
   nominal__926->SetBinContent(17,1);
   nominal__926->SetBinContent(18,1);
   nominal__926->SetBinContent(19,1);
   nominal__926->SetBinContent(20,1);
   nominal__926->SetBinContent(21,1);
   nominal__926->SetBinContent(22,1);
   nominal__926->SetBinContent(23,1);
   nominal__926->SetBinContent(24,1);
   nominal__926->SetBinContent(25,1);
   nominal__926->SetBinContent(26,1);
   nominal__926->SetBinContent(27,1);
   nominal__926->SetBinContent(28,1);
   nominal__926->SetBinContent(29,1);
   nominal__926->SetBinContent(30,1);
   nominal__926->SetBinContent(31,1);
   nominal__926->SetBinContent(32,1);
   nominal__926->SetBinError(2,1.414214);
   nominal__926->SetBinError(5,1.414214);
   nominal__926->SetBinError(6,1.414214);
   nominal__926->SetBinError(7,1.414214);
   nominal__926->SetBinError(8,1.414214);
   nominal__926->SetBinError(9,1.414214);
   nominal__926->SetBinError(10,1.414214);
   nominal__926->SetBinError(11,1.414214);
   nominal__926->SetBinError(12,0.7080846);
   nominal__926->SetBinError(13,0.6325159);
   nominal__926->SetBinError(14,0.5005981);
   nominal__926->SetBinError(15,0.473242);
   nominal__926->SetBinError(16,0.2465117);
   nominal__926->SetBinError(17,0.1983161);
   nominal__926->SetBinError(18,0.1247315);
   nominal__926->SetBinError(19,0.08808497);
   nominal__926->SetBinError(20,0.05209112);
   nominal__926->SetBinError(21,0.03550394);
   nominal__926->SetBinError(22,0.02108114);
   nominal__926->SetBinError(23,0.01338148);
   nominal__926->SetBinError(24,0.01179134);
   nominal__926->SetBinError(25,0.01684867);
   nominal__926->SetBinError(26,0.03213928);
   nominal__926->SetBinError(27,0.08221432);
   nominal__926->SetBinError(28,0.177068);
   nominal__926->SetBinError(29,0.3338111);
   nominal__926->SetBinError(30,0.6336282);
   nominal__926->SetBinError(31,1.000081);
   nominal__926->SetBinError(32,0.7445387);
   nominal__926->SetMinimum(5e-05);
   nominal__926->SetMaximum(67036);
   nominal__926->SetEntries(42.88407);
   nominal__926->SetFillColor(38);
   nominal__926->SetLineColor(38);
   nominal__926->SetMarkerColor(38);
   nominal__926->SetMarkerStyle(21);
   nominal__926->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__926->GetXaxis()->SetRange(1,32);
   nominal__926->GetXaxis()->SetLabelFont(42);
   nominal__926->GetXaxis()->SetLabelSize(0.035);
   nominal__926->GetXaxis()->SetTitleSize(0.035);
   nominal__926->GetXaxis()->SetTitleFont(42);
   nominal__926->GetYaxis()->SetTitle("Tracks");
   nominal__926->GetYaxis()->SetLabelFont(42);
   nominal__926->GetYaxis()->SetLabelSize(0.05);
   nominal__926->GetYaxis()->SetTitleSize(0.07);
   nominal__926->GetYaxis()->SetTitleOffset(0);
   nominal__926->GetYaxis()->SetTitleFont(42);
   nominal__926->GetZaxis()->SetLabelFont(42);
   nominal__926->GetZaxis()->SetLabelSize(0.035);
   nominal__926->GetZaxis()->SetTitleSize(0.035);
   nominal__926->GetZaxis()->SetTitleFont(42);
   nominal__926->Draw("E0 same");
   Double_t xAxis721[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__927 = new TH1F("nominal__927","",32, xAxis721);
   nominal__927->SetBinContent(2,1);
   nominal__927->SetBinContent(5,1);
   nominal__927->SetBinContent(6,1);
   nominal__927->SetBinContent(7,1);
   nominal__927->SetBinContent(8,1);
   nominal__927->SetBinContent(9,1);
   nominal__927->SetBinContent(10,1);
   nominal__927->SetBinContent(11,1);
   nominal__927->SetBinContent(12,1);
   nominal__927->SetBinContent(13,1);
   nominal__927->SetBinContent(14,1);
   nominal__927->SetBinContent(15,1);
   nominal__927->SetBinContent(16,1);
   nominal__927->SetBinContent(17,1);
   nominal__927->SetBinContent(18,1);
   nominal__927->SetBinContent(19,1);
   nominal__927->SetBinContent(20,1);
   nominal__927->SetBinContent(21,1);
   nominal__927->SetBinContent(22,1);
   nominal__927->SetBinContent(23,1);
   nominal__927->SetBinContent(24,1);
   nominal__927->SetBinContent(25,1);
   nominal__927->SetBinContent(26,1);
   nominal__927->SetBinContent(27,1);
   nominal__927->SetBinContent(28,1);
   nominal__927->SetBinContent(29,1);
   nominal__927->SetBinContent(30,1);
   nominal__927->SetBinContent(31,1);
   nominal__927->SetBinContent(32,1);
   nominal__927->SetBinError(2,1.414214);
   nominal__927->SetBinError(5,1.414214);
   nominal__927->SetBinError(6,1.414214);
   nominal__927->SetBinError(7,1.414214);
   nominal__927->SetBinError(8,1.414214);
   nominal__927->SetBinError(9,1.414214);
   nominal__927->SetBinError(10,1.414214);
   nominal__927->SetBinError(11,1.414214);
   nominal__927->SetBinError(12,0.7080846);
   nominal__927->SetBinError(13,0.6325159);
   nominal__927->SetBinError(14,0.5005981);
   nominal__927->SetBinError(15,0.473242);
   nominal__927->SetBinError(16,0.2465117);
   nominal__927->SetBinError(17,0.1983161);
   nominal__927->SetBinError(18,0.1247315);
   nominal__927->SetBinError(19,0.08808497);
   nominal__927->SetBinError(20,0.05209112);
   nominal__927->SetBinError(21,0.03550394);
   nominal__927->SetBinError(22,0.02108114);
   nominal__927->SetBinError(23,0.01338148);
   nominal__927->SetBinError(24,0.01179134);
   nominal__927->SetBinError(25,0.01684867);
   nominal__927->SetBinError(26,0.03213928);
   nominal__927->SetBinError(27,0.08221432);
   nominal__927->SetBinError(28,0.177068);
   nominal__927->SetBinError(29,0.3338111);
   nominal__927->SetBinError(30,0.6336282);
   nominal__927->SetBinError(31,1.000081);
   nominal__927->SetBinError(32,0.7445387);
   nominal__927->SetMinimum(5e-05);
   nominal__927->SetMaximum(67036);
   nominal__927->SetEntries(42.88407);
   nominal__927->SetFillColor(46);
   nominal__927->SetLineColor(46);
   nominal__927->SetMarkerColor(46);
   nominal__927->SetMarkerStyle(21);
   nominal__927->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__927->GetXaxis()->SetRange(1,32);
   nominal__927->GetXaxis()->SetLabelFont(42);
   nominal__927->GetXaxis()->SetLabelSize(0.035);
   nominal__927->GetXaxis()->SetTitleSize(0.035);
   nominal__927->GetXaxis()->SetTitleFont(42);
   nominal__927->GetYaxis()->SetTitle("Tracks");
   nominal__927->GetYaxis()->SetLabelFont(42);
   nominal__927->GetYaxis()->SetLabelSize(0.05);
   nominal__927->GetYaxis()->SetTitleSize(0.07);
   nominal__927->GetYaxis()->SetTitleOffset(0);
   nominal__927->GetYaxis()->SetTitleFont(42);
   nominal__927->GetZaxis()->SetLabelFont(42);
   nominal__927->GetZaxis()->SetLabelSize(0.035);
   nominal__927->GetZaxis()->SetTitleSize(0.035);
   nominal__927->GetZaxis()->SetTitleFont(42);
   nominal__927->Draw("E0 same");
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
