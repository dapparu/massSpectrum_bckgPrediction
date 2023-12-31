void plot_K1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:06 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.354892,2428.571,6.417608);
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
   Double_t xAxis792[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1018 = new TH1F("nominal__1018","",32, xAxis792);
   nominal__1018->SetBinContent(2,0.0001276471);
   nominal__1018->SetBinContent(4,0.0001254598);
   nominal__1018->SetBinContent(5,0.0005029399);
   nominal__1018->SetBinContent(7,0.000125725);
   nominal__1018->SetBinContent(8,0.0002434331);
   nominal__1018->SetBinContent(9,0.0001226969);
   nominal__1018->SetBinContent(10,0.0001225562);
   nominal__1018->SetBinContent(11,0.0001118656);
   nominal__1018->SetBinContent(12,0.0005197951);
   nominal__1018->SetBinContent(13,0.0002517582);
   nominal__1018->SetBinContent(14,0.0004995628);
   nominal__1018->SetBinContent(15,0.0003769499);
   nominal__1018->SetBinContent(16,0.0008585473);
   nominal__1018->SetBinContent(17,0.002139057);
   nominal__1018->SetBinContent(18,0.00229686);
   nominal__1018->SetBinContent(19,0.004799545);
   nominal__1018->SetBinContent(20,0.02039263);
   nominal__1018->SetBinContent(21,0.04602493);
   nominal__1018->SetBinContent(22,0.1375895);
   nominal__1018->SetBinContent(23,0.4305354);
   nominal__1018->SetBinContent(24,1.283651);
   nominal__1018->SetBinContent(25,2.310712);
   nominal__1018->SetBinContent(26,1.861487);
   nominal__1018->SetBinContent(27,0.4265608);
   nominal__1018->SetBinContent(28,0.0716294);
   nominal__1018->SetBinContent(29,0.01315903);
   nominal__1018->SetBinContent(30,0.002491383);
   nominal__1018->SetBinContent(31,0.0005185502);
   nominal__1018->SetBinContent(32,0.0005040495);
   nominal__1018->SetBinError(2,0.0001276471);
   nominal__1018->SetBinError(4,0.0001254598);
   nominal__1018->SetBinError(5,0.000252398);
   nominal__1018->SetBinError(7,0.000125725);
   nominal__1018->SetBinError(8,0.0001724051);
   nominal__1018->SetBinError(9,0.0001226969);
   nominal__1018->SetBinError(10,0.0001225562);
   nominal__1018->SetBinError(11,0.0001118656);
   nominal__1018->SetBinError(12,0.0002600355);
   nominal__1018->SetBinError(13,0.0001780375);
   nominal__1018->SetBinError(14,0.000249798);
   nominal__1018->SetBinError(15,0.0002186157);
   nominal__1018->SetBinError(16,0.0003249841);
   nominal__1018->SetBinError(17,0.0005198812);
   nominal__1018->SetBinError(18,0.0005419345);
   nominal__1018->SetBinError(19,0.0007798717);
   nominal__1018->SetBinError(20,0.001600887);
   nominal__1018->SetBinError(21,0.002400276);
   nominal__1018->SetBinError(22,0.004149472);
   nominal__1018->SetBinError(23,0.007337162);
   nominal__1018->SetBinError(24,0.01267707);
   nominal__1018->SetBinError(25,0.01701336);
   nominal__1018->SetBinError(26,0.01528545);
   nominal__1018->SetBinError(27,0.007326473);
   nominal__1018->SetBinError(28,0.003003646);
   nominal__1018->SetBinError(29,0.001292828);
   nominal__1018->SetBinError(30,0.000558052);
   nominal__1018->SetBinError(31,0.000259427);
   nominal__1018->SetBinError(32,0.0001270409);
   nominal__1018->SetBinError(33,0.0002177277);
   nominal__1018->SetMinimum(5e-05);
   nominal__1018->SetMaximum(2310712);
   nominal__1018->SetEntries(53005);
   nominal__1018->SetFillColor(1);
   nominal__1018->SetMarkerStyle(20);
   nominal__1018->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1018->GetXaxis()->SetRange(1,32);
   nominal__1018->GetXaxis()->SetLabelFont(42);
   nominal__1018->GetXaxis()->SetLabelSize(0.035);
   nominal__1018->GetXaxis()->SetTitleSize(0.035);
   nominal__1018->GetXaxis()->SetTitleFont(42);
   nominal__1018->GetYaxis()->SetTitle("Tracks");
   nominal__1018->GetYaxis()->SetLabelFont(42);
   nominal__1018->GetYaxis()->SetLabelSize(0.05);
   nominal__1018->GetYaxis()->SetTitleSize(0.07);
   nominal__1018->GetYaxis()->SetTitleOffset(0);
   nominal__1018->GetYaxis()->SetTitleFont(42);
   nominal__1018->GetZaxis()->SetLabelFont(42);
   nominal__1018->GetZaxis()->SetLabelSize(0.035);
   nominal__1018->GetZaxis()->SetTitleSize(0.035);
   nominal__1018->GetZaxis()->SetTitleFont(42);
   nominal__1018->Draw("");
   Double_t xAxis793[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__1019 = new TH1F("K_down1__1019","",32, xAxis793);
   K_down1__1019->SetBinContent(2,0.0001276471);
   K_down1__1019->SetBinContent(4,0.0001254598);
   K_down1__1019->SetBinContent(5,0.0005029399);
   K_down1__1019->SetBinContent(7,0.000125725);
   K_down1__1019->SetBinContent(8,0.0002434331);
   K_down1__1019->SetBinContent(9,0.0001226969);
   K_down1__1019->SetBinContent(10,0.0001225562);
   K_down1__1019->SetBinContent(11,0.0001118656);
   K_down1__1019->SetBinContent(12,0.0005197951);
   K_down1__1019->SetBinContent(13,0.0002517582);
   K_down1__1019->SetBinContent(14,0.0004995628);
   K_down1__1019->SetBinContent(15,0.0003769499);
   K_down1__1019->SetBinContent(16,0.0007351505);
   K_down1__1019->SetBinContent(17,0.002014033);
   K_down1__1019->SetBinContent(18,0.002290235);
   K_down1__1019->SetBinContent(19,0.004574508);
   K_down1__1019->SetBinContent(20,0.02009956);
   K_down1__1019->SetBinContent(21,0.04467551);
   K_down1__1019->SetBinContent(22,0.1355132);
   K_down1__1019->SetBinContent(23,0.4198046);
   K_down1__1019->SetBinContent(24,1.257954);
   K_down1__1019->SetBinContent(25,2.308392);
   K_down1__1019->SetBinContent(26,1.890308);
   K_down1__1019->SetBinContent(27,0.4385396);
   K_down1__1019->SetBinContent(28,0.07315595);
   K_down1__1019->SetBinContent(29,0.01377958);
   K_down1__1019->SetBinContent(30,0.002491383);
   K_down1__1019->SetBinContent(31,0.0005185502);
   K_down1__1019->SetBinContent(32,0.0005040495);
   K_down1__1019->SetBinError(2,0.0001276471);
   K_down1__1019->SetBinError(4,0.0001254598);
   K_down1__1019->SetBinError(5,0.000252398);
   K_down1__1019->SetBinError(7,0.000125725);
   K_down1__1019->SetBinError(8,0.0001724051);
   K_down1__1019->SetBinError(9,0.0001226969);
   K_down1__1019->SetBinError(10,0.0001225562);
   K_down1__1019->SetBinError(11,0.0001118656);
   K_down1__1019->SetBinError(12,0.0002600355);
   K_down1__1019->SetBinError(13,0.0001780375);
   K_down1__1019->SetBinError(14,0.000249798);
   K_down1__1019->SetBinError(15,0.0002186157);
   K_down1__1019->SetBinError(16,0.0003006458);
   K_down1__1019->SetBinError(17,0.0005046199);
   K_down1__1019->SetBinError(18,0.0005403556);
   K_down1__1019->SetBinError(19,0.0007634934);
   K_down1__1019->SetBinError(20,0.001587502);
   K_down1__1019->SetBinError(21,0.002365411);
   K_down1__1019->SetBinError(22,0.004118724);
   K_down1__1019->SetBinError(23,0.0072448);
   K_down1__1019->SetBinError(24,0.01255038);
   K_down1__1019->SetBinError(25,0.01700452);
   K_down1__1019->SetBinError(26,0.01540209);
   K_down1__1019->SetBinError(27,0.007428974);
   K_down1__1019->SetBinError(28,0.003035947);
   K_down1__1019->SetBinError(29,0.001322443);
   K_down1__1019->SetBinError(30,0.000558052);
   K_down1__1019->SetBinError(31,0.000259427);
   K_down1__1019->SetBinError(32,0.0001270409);
   K_down1__1019->SetBinError(33,0.0002177277);
   K_down1__1019->SetEntries(53005);
   K_down1__1019->SetFillColor(38);
   K_down1__1019->SetLineColor(38);
   K_down1__1019->SetMarkerColor(38);
   K_down1__1019->SetMarkerStyle(21);
   K_down1__1019->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__1019->GetXaxis()->SetRange(1,400);
   K_down1__1019->GetXaxis()->SetLabelFont(42);
   K_down1__1019->GetXaxis()->SetLabelSize(0.035);
   K_down1__1019->GetXaxis()->SetTitleSize(0.035);
   K_down1__1019->GetXaxis()->SetTitleFont(42);
   K_down1__1019->GetYaxis()->SetTitle("Events / bin");
   K_down1__1019->GetYaxis()->SetLabelFont(42);
   K_down1__1019->GetYaxis()->SetLabelSize(0.035);
   K_down1__1019->GetYaxis()->SetTitleSize(0.035);
   K_down1__1019->GetYaxis()->SetTitleOffset(0);
   K_down1__1019->GetYaxis()->SetTitleFont(42);
   K_down1__1019->GetZaxis()->SetLabelFont(42);
   K_down1__1019->GetZaxis()->SetLabelSize(0.035);
   K_down1__1019->GetZaxis()->SetTitleSize(0.035);
   K_down1__1019->GetZaxis()->SetTitleFont(42);
   K_down1__1019->Draw("same");
   Double_t xAxis794[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__1020 = new TH1F("K_up1__1020","",32, xAxis794);
   K_up1__1020->SetBinContent(2,0.0001276471);
   K_up1__1020->SetBinContent(4,0.0001254598);
   K_up1__1020->SetBinContent(5,0.0005029399);
   K_up1__1020->SetBinContent(7,0.000125725);
   K_up1__1020->SetBinContent(8,0.0002434331);
   K_up1__1020->SetBinContent(9,0.0001226969);
   K_up1__1020->SetBinContent(10,0.0001225562);
   K_up1__1020->SetBinContent(11,0.0001118656);
   K_up1__1020->SetBinContent(12,0.0005197951);
   K_up1__1020->SetBinContent(13,0.0002517582);
   K_up1__1020->SetBinContent(14,0.0004995628);
   K_up1__1020->SetBinContent(15,0.0003769499);
   K_up1__1020->SetBinContent(16,0.0008585473);
   K_up1__1020->SetBinContent(17,0.002139057);
   K_up1__1020->SetBinContent(18,0.002425421);
   K_up1__1020->SetBinContent(19,0.004926129);
   K_up1__1020->SetBinContent(20,0.02102119);
   K_up1__1020->SetBinContent(21,0.04773964);
   K_up1__1020->SetBinContent(22,0.1388434);
   K_up1__1020->SetBinContent(23,0.4418609);
   K_up1__1020->SetBinContent(24,1.308814);
   K_up1__1020->SetBinContent(25,2.31703);
   K_up1__1020->SetBinContent(26,1.828466);
   K_up1__1020->SetBinContent(27,0.414019);
   K_up1__1020->SetBinContent(28,0.07090919);
   K_up1__1020->SetBinContent(29,0.01278296);
   K_up1__1020->SetBinContent(30,0.002491383);
   K_up1__1020->SetBinContent(31,0.0005185502);
   K_up1__1020->SetBinContent(32,0.0005040495);
   K_up1__1020->SetBinError(2,0.0001276471);
   K_up1__1020->SetBinError(4,0.0001254598);
   K_up1__1020->SetBinError(5,0.000252398);
   K_up1__1020->SetBinError(7,0.000125725);
   K_up1__1020->SetBinError(8,0.0001724051);
   K_up1__1020->SetBinError(9,0.0001226969);
   K_up1__1020->SetBinError(10,0.0001225562);
   K_up1__1020->SetBinError(11,0.0001118656);
   K_up1__1020->SetBinError(12,0.0002600355);
   K_up1__1020->SetBinError(13,0.0001780375);
   K_up1__1020->SetBinError(14,0.000249798);
   K_up1__1020->SetBinError(15,0.0002186157);
   K_up1__1020->SetBinError(16,0.0003249841);
   K_up1__1020->SetBinError(17,0.0005198812);
   K_up1__1020->SetBinError(18,0.0005569748);
   K_up1__1020->SetBinError(19,0.0007900825);
   K_up1__1020->SetBinError(20,0.00162543);
   K_up1__1020->SetBinError(21,0.002443801);
   K_up1__1020->SetBinError(22,0.004168548);
   K_up1__1020->SetBinError(23,0.007433829);
   K_up1__1020->SetBinError(24,0.01280019);
   K_up1__1020->SetBinError(25,0.01703634);
   K_up1__1020->SetBinError(26,0.01515078);
   K_up1__1020->SetBinError(27,0.00721673);
   K_up1__1020->SetBinError(28,0.002989189);
   K_up1__1020->SetBinError(29,0.001274458);
   K_up1__1020->SetBinError(30,0.000558052);
   K_up1__1020->SetBinError(31,0.000259427);
   K_up1__1020->SetBinError(32,0.0001270409);
   K_up1__1020->SetBinError(33,0.0002177277);
   K_up1__1020->SetEntries(53005);
   K_up1__1020->SetFillColor(46);
   K_up1__1020->SetLineColor(46);
   K_up1__1020->SetMarkerColor(46);
   K_up1__1020->SetMarkerStyle(21);
   K_up1__1020->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__1020->GetXaxis()->SetRange(1,400);
   K_up1__1020->GetXaxis()->SetLabelFont(42);
   K_up1__1020->GetXaxis()->SetLabelSize(0.035);
   K_up1__1020->GetXaxis()->SetTitleSize(0.035);
   K_up1__1020->GetXaxis()->SetTitleFont(42);
   K_up1__1020->GetYaxis()->SetTitle("Events / bin");
   K_up1__1020->GetYaxis()->SetLabelFont(42);
   K_up1__1020->GetYaxis()->SetLabelSize(0.035);
   K_up1__1020->GetYaxis()->SetTitleSize(0.035);
   K_up1__1020->GetYaxis()->SetTitleOffset(0);
   K_up1__1020->GetYaxis()->SetTitleFont(42);
   K_up1__1020->GetZaxis()->SetLabelFont(42);
   K_up1__1020->GetZaxis()->SetLabelSize(0.035);
   K_up1__1020->GetZaxis()->SetTitleSize(0.035);
   K_up1__1020->GetZaxis()->SetTitleFont(42);
   K_up1__1020->Draw("same");
   TLine *line = new TLine(1730,0,1730,2310712);
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
   entry=leg->AddEntry("K_down1","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("K_up1","Up","PE1");
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
   
   TH1D *frameR2__1021 = new TH1D("frameR2__1021","",1,0,2000);
   frameR2__1021->SetMinimum(0.5);
   frameR2__1021->SetMaximum(1.5);
   frameR2__1021->SetStats(0);
   frameR2__1021->SetLineStyle(0);
   frameR2__1021->SetMarkerStyle(20);
   frameR2__1021->GetXaxis()->SetRange(1,1);
   frameR2__1021->GetXaxis()->SetLabelFont(43);
   frameR2__1021->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1021->GetXaxis()->SetLabelSize(16);
   frameR2__1021->GetXaxis()->SetTitleSize(24);
   frameR2__1021->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1021->GetXaxis()->SetTitleFont(43);
   frameR2__1021->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1021->GetYaxis()->SetNdivisions(205);
   frameR2__1021->GetYaxis()->SetLabelFont(43);
   frameR2__1021->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1021->GetYaxis()->SetLabelSize(20);
   frameR2__1021->GetYaxis()->SetTitleSize(20);
   frameR2__1021->GetYaxis()->SetTitleOffset(2);
   frameR2__1021->GetYaxis()->SetTitleFont(43);
   frameR2__1021->GetZaxis()->SetLabelFont(42);
   frameR2__1021->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1021->GetZaxis()->SetLabelSize(0.05);
   frameR2__1021->GetZaxis()->SetTitleSize(0.06);
   frameR2__1021->GetZaxis()->SetTitleFont(42);
   frameR2__1021->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis795[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__1022 = new TH1F("K_down1__1022","",32, xAxis795);
   K_down1__1022->SetBinContent(0,0.9999999);
   K_down1__1022->SetBinContent(1,0.9999999);
   K_down1__1022->SetBinContent(2,0.9999999);
   K_down1__1022->SetBinContent(3,0.9999999);
   K_down1__1022->SetBinContent(4,0.9999999);
   K_down1__1022->SetBinContent(5,0.9999999);
   K_down1__1022->SetBinContent(6,0.9999999);
   K_down1__1022->SetBinContent(7,0.9999999);
   K_down1__1022->SetBinContent(8,0.9999999);
   K_down1__1022->SetBinContent(9,0.9999999);
   K_down1__1022->SetBinContent(10,0.9999999);
   K_down1__1022->SetBinContent(11,0.9999999);
   K_down1__1022->SetBinContent(12,0.9999999);
   K_down1__1022->SetBinContent(13,0.9999999);
   K_down1__1022->SetBinContent(14,0.9999999);
   K_down1__1022->SetBinContent(15,0.9999999);
   K_down1__1022->SetBinContent(16,0.9999999);
   K_down1__1022->SetBinContent(17,0.9999813);
   K_down1__1022->SetBinContent(18,0.9999624);
   K_down1__1022->SetBinContent(19,0.9999614);
   K_down1__1022->SetBinContent(20,0.9999273);
   K_down1__1022->SetBinContent(21,0.9998826);
   K_down1__1022->SetBinContent(22,0.9996755);
   K_down1__1022->SetBinContent(23,0.9993445);
   K_down1__1022->SetBinContent(24,0.9975057);
   K_down1__1022->SetBinContent(25,0.9914066);
   K_down1__1022->SetBinContent(26,0.9822478);
   K_down1__1022->SetBinContent(27,0.9732964);
   K_down1__1022->SetBinContent(28,0.9762619);
   K_down1__1022->SetBinContent(29,0.9641169);
   K_down1__1022->SetBinContent(30,1);
   K_down1__1022->SetBinContent(31,1);
   K_down1__1022->SetBinContent(32,1);
   K_down1__1022->SetBinError(0,0.006154726);
   K_down1__1022->SetBinError(1,0.006154726);
   K_down1__1022->SetBinError(2,0.006154726);
   K_down1__1022->SetBinError(3,0.006154784);
   K_down1__1022->SetBinError(4,0.006154784);
   K_down1__1022->SetBinError(5,0.006154842);
   K_down1__1022->SetBinError(6,0.006155074);
   K_down1__1022->SetBinError(7,0.006155074);
   K_down1__1022->SetBinError(8,0.006155132);
   K_down1__1022->SetBinError(9,0.006155248);
   K_down1__1022->SetBinError(10,0.006155306);
   K_down1__1022->SetBinError(11,0.006155365);
   K_down1__1022->SetBinError(12,0.006155422);
   K_down1__1022->SetBinError(13,0.006155655);
   K_down1__1022->SetBinError(14,0.006155772);
   K_down1__1022->SetBinError(15,0.006156005);
   K_down1__1022->SetBinError(16,0.006156178);
   K_down1__1022->SetBinError(17,0.006156441);
   K_down1__1022->SetBinError(18,0.006157283);
   K_down1__1022->SetBinError(19,0.006158325);
   K_down1__1022->SetBinError(20,0.00616027);
   K_down1__1022->SetBinError(21,0.006169443);
   K_down1__1022->SetBinError(22,0.006189539);
   K_down1__1022->SetBinError(23,0.006253219);
   K_down1__1022->SetBinError(24,0.006460458);
   K_down1__1022->SetBinError(25,0.007237379);
   K_down1__1022->SetBinError(26,0.01005358);
   K_down1__1022->SetBinError(27,0.02137916);
   K_down1__1022->SetBinError(28,0.05186672);
   K_down1__1022->SetBinError(29,0.1178037);
   K_down1__1022->SetBinError(30,0.2676448);
   K_down1__1022->SetBinError(31,0.500255);
   K_down1__1022->SetBinError(32,0.7072641);
   K_down1__1022->SetEntries(33);
   K_down1__1022->SetFillColor(38);
   K_down1__1022->SetLineColor(38);
   K_down1__1022->SetMarkerColor(38);
   K_down1__1022->SetMarkerStyle(21);
   K_down1__1022->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__1022->GetXaxis()->SetRange(1,400);
   K_down1__1022->GetXaxis()->SetLabelFont(42);
   K_down1__1022->GetXaxis()->SetLabelSize(0.035);
   K_down1__1022->GetXaxis()->SetTitleSize(0.035);
   K_down1__1022->GetXaxis()->SetTitleFont(42);
   K_down1__1022->GetYaxis()->SetTitle("Events / bin");
   K_down1__1022->GetYaxis()->SetLabelFont(42);
   K_down1__1022->GetYaxis()->SetLabelSize(0.035);
   K_down1__1022->GetYaxis()->SetTitleSize(0.035);
   K_down1__1022->GetYaxis()->SetTitleOffset(0);
   K_down1__1022->GetYaxis()->SetTitleFont(42);
   K_down1__1022->GetZaxis()->SetLabelFont(42);
   K_down1__1022->GetZaxis()->SetLabelSize(0.035);
   K_down1__1022->GetZaxis()->SetTitleSize(0.035);
   K_down1__1022->GetZaxis()->SetTitleFont(42);
   K_down1__1022->Draw("E0 same");
   Double_t xAxis796[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__1023 = new TH1F("K_up1__1023","",32, xAxis796);
   K_up1__1023->SetBinContent(0,0.9999999);
   K_up1__1023->SetBinContent(1,0.9999999);
   K_up1__1023->SetBinContent(2,0.9999999);
   K_up1__1023->SetBinContent(3,0.9999999);
   K_up1__1023->SetBinContent(4,0.9999999);
   K_up1__1023->SetBinContent(5,0.9999999);
   K_up1__1023->SetBinContent(6,0.9999999);
   K_up1__1023->SetBinContent(7,0.9999999);
   K_up1__1023->SetBinContent(8,0.9999999);
   K_up1__1023->SetBinContent(9,0.9999999);
   K_up1__1023->SetBinContent(10,0.9999999);
   K_up1__1023->SetBinContent(11,0.9999999);
   K_up1__1023->SetBinContent(12,0.9999999);
   K_up1__1023->SetBinContent(13,0.9999999);
   K_up1__1023->SetBinContent(14,0.9999999);
   K_up1__1023->SetBinContent(15,0.9999999);
   K_up1__1023->SetBinContent(16,0.9999999);
   K_up1__1023->SetBinContent(17,0.9999999);
   K_up1__1023->SetBinContent(18,0.9999999);
   K_up1__1023->SetBinContent(19,1.000019);
   K_up1__1023->SetBinContent(20,1.000039);
   K_up1__1023->SetBinContent(21,1.000134);
   K_up1__1023->SetBinContent(22,1.000398);
   K_up1__1023->SetBinContent(23,1.000602);
   K_up1__1023->SetBinContent(24,1.002548);
   K_up1__1023->SetBinContent(25,1.008682);
   K_up1__1023->SetBinContent(26,1.020028);
   K_up1__1023->SetBinContent(27,1.027209);
   K_up1__1023->SetBinContent(28,1.012571);
   K_up1__1023->SetBinContent(29,1.023076);
   K_up1__1023->SetBinContent(30,0.9999999);
   K_up1__1023->SetBinContent(31,1);
   K_up1__1023->SetBinContent(32,1);
   K_up1__1023->SetBinError(0,0.006154726);
   K_up1__1023->SetBinError(1,0.006154726);
   K_up1__1023->SetBinError(2,0.006154726);
   K_up1__1023->SetBinError(3,0.006154784);
   K_up1__1023->SetBinError(4,0.006154784);
   K_up1__1023->SetBinError(5,0.006154842);
   K_up1__1023->SetBinError(6,0.006155074);
   K_up1__1023->SetBinError(7,0.006155074);
   K_up1__1023->SetBinError(8,0.006155132);
   K_up1__1023->SetBinError(9,0.006155248);
   K_up1__1023->SetBinError(10,0.006155306);
   K_up1__1023->SetBinError(11,0.006155365);
   K_up1__1023->SetBinError(12,0.006155422);
   K_up1__1023->SetBinError(13,0.006155655);
   K_up1__1023->SetBinError(14,0.006155772);
   K_up1__1023->SetBinError(15,0.006156005);
   K_up1__1023->SetBinError(16,0.006156178);
   K_up1__1023->SetBinError(17,0.006156585);
   K_up1__1023->SetBinError(18,0.006157573);
   K_up1__1023->SetBinError(19,0.00615877);
   K_up1__1023->SetBinError(20,0.00616113);
   K_up1__1023->SetBinError(21,0.006171375);
   K_up1__1023->SetBinError(22,0.006195133);
   K_up1__1023->SetBinError(23,0.006263071);
   K_up1__1023->SetBinError(24,0.006501304);
   K_up1__1023->SetBinError(25,0.007395502);
   K_up1__1023->SetBinError(26,0.01053983);
   K_up1__1023->SetBinError(27,0.02286778);
   K_up1__1023->SetBinError(28,0.05429113);
   K_up1__1023->SetBinError(29,0.1268955);
   K_up1__1023->SetBinError(30,0.2676447);
   K_up1__1023->SetBinError(31,0.500255);
   K_up1__1023->SetBinError(32,0.7072641);
   K_up1__1023->SetEntries(33);
   K_up1__1023->SetFillColor(46);
   K_up1__1023->SetLineColor(46);
   K_up1__1023->SetMarkerColor(46);
   K_up1__1023->SetMarkerStyle(21);
   K_up1__1023->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__1023->GetXaxis()->SetRange(1,400);
   K_up1__1023->GetXaxis()->SetLabelFont(42);
   K_up1__1023->GetXaxis()->SetLabelSize(0.035);
   K_up1__1023->GetXaxis()->SetTitleSize(0.035);
   K_up1__1023->GetXaxis()->SetTitleFont(42);
   K_up1__1023->GetYaxis()->SetTitle("Events / bin");
   K_up1__1023->GetYaxis()->SetLabelFont(42);
   K_up1__1023->GetYaxis()->SetLabelSize(0.035);
   K_up1__1023->GetYaxis()->SetTitleSize(0.035);
   K_up1__1023->GetYaxis()->SetTitleOffset(0);
   K_up1__1023->GetYaxis()->SetTitleFont(42);
   K_up1__1023->GetZaxis()->SetLabelFont(42);
   K_up1__1023->GetZaxis()->SetLabelSize(0.035);
   K_up1__1023->GetZaxis()->SetTitleSize(0.035);
   K_up1__1023->GetZaxis()->SetTitleFont(42);
   K_up1__1023->Draw("E0 same");
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
   
   TH1D *frameR2__1024 = new TH1D("frameR2__1024","",1,0,2000);
   frameR2__1024->SetMinimum(0.5);
   frameR2__1024->SetMaximum(1.5);
   frameR2__1024->SetStats(0);
   frameR2__1024->SetLineStyle(0);
   frameR2__1024->SetMarkerStyle(20);
   frameR2__1024->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1024->GetXaxis()->SetRange(1,1);
   frameR2__1024->GetXaxis()->SetLabelFont(43);
   frameR2__1024->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1024->GetXaxis()->SetLabelSize(16);
   frameR2__1024->GetXaxis()->SetTitleSize(24);
   frameR2__1024->GetXaxis()->SetTitleOffset(5);
   frameR2__1024->GetXaxis()->SetTitleFont(43);
   frameR2__1024->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1024->GetYaxis()->SetNdivisions(205);
   frameR2__1024->GetYaxis()->SetLabelFont(43);
   frameR2__1024->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1024->GetYaxis()->SetLabelSize(20);
   frameR2__1024->GetYaxis()->SetTitleSize(20);
   frameR2__1024->GetYaxis()->SetTitleOffset(2);
   frameR2__1024->GetYaxis()->SetTitleFont(43);
   frameR2__1024->GetZaxis()->SetLabelFont(42);
   frameR2__1024->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1024->GetZaxis()->SetLabelSize(0.05);
   frameR2__1024->GetZaxis()->SetTitleSize(0.06);
   frameR2__1024->GetZaxis()->SetTitleFont(42);
   frameR2__1024->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis797[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1025 = new TH1F("nominal__1025","",32, xAxis797);
   nominal__1025->SetBinContent(2,1);
   nominal__1025->SetBinContent(4,1);
   nominal__1025->SetBinContent(5,1);
   nominal__1025->SetBinContent(7,1);
   nominal__1025->SetBinContent(8,1);
   nominal__1025->SetBinContent(9,1);
   nominal__1025->SetBinContent(10,1);
   nominal__1025->SetBinContent(11,1);
   nominal__1025->SetBinContent(12,1);
   nominal__1025->SetBinContent(13,1);
   nominal__1025->SetBinContent(14,1);
   nominal__1025->SetBinContent(15,1);
   nominal__1025->SetBinContent(16,1.167852);
   nominal__1025->SetBinContent(17,1.062076);
   nominal__1025->SetBinContent(18,1.002893);
   nominal__1025->SetBinContent(19,1.049194);
   nominal__1025->SetBinContent(20,1.014581);
   nominal__1025->SetBinContent(21,1.030205);
   nominal__1025->SetBinContent(22,1.015322);
   nominal__1025->SetBinContent(23,1.025561);
   nominal__1025->SetBinContent(24,1.020428);
   nominal__1025->SetBinContent(25,1.001005);
   nominal__1025->SetBinContent(26,0.9847528);
   nominal__1025->SetBinContent(27,0.9726848);
   nominal__1025->SetBinContent(28,0.979133);
   nominal__1025->SetBinContent(29,0.9549662);
   nominal__1025->SetBinContent(30,1);
   nominal__1025->SetBinContent(31,1);
   nominal__1025->SetBinContent(32,1);
   nominal__1025->SetBinError(2,1.414214);
   nominal__1025->SetBinError(4,1.414214);
   nominal__1025->SetBinError(5,0.7097162);
   nominal__1025->SetBinError(7,1.414214);
   nominal__1025->SetBinError(8,1.00158);
   nominal__1025->SetBinError(9,1.414214);
   nominal__1025->SetBinError(10,1.414213);
   nominal__1025->SetBinError(11,1.414214);
   nominal__1025->SetBinError(12,0.7074821);
   nominal__1025->SetBinError(13,1.000099);
   nominal__1025->SetBinError(14,0.7071538);
   nominal__1025->SetBinError(15,0.8201866);
   nominal__1025->SetBinError(16,0.6507884);
   nominal__1025->SetBinError(17,0.3707329);
   nominal__1025->SetBinError(18,0.3346381);
   nominal__1025->SetBinError(19,0.244394);
   nominal__1025->SetBinError(20,0.1129831);
   nominal__1025->SetBinError(21,0.07656249);
   nominal__1025->SetBinError(22,0.043473);
   nominal__1025->SetBinError(23,0.02487385);
   nominal__1025->SetBinError(24,0.01432486);
   nominal__1025->SetBinError(25,0.0104256);
   nominal__1025->SetBinError(26,0.01139151);
   nominal__1025->SetBinError(27,0.02346523);
   nominal__1025->SetBinError(28,0.05776562);
   nominal__1025->SetBinError(29,0.131157);
   nominal__1025->SetBinError(30,0.3167738);
   nominal__1025->SetBinError(31,0.7075212);
   nominal__1025->SetBinError(32,0.3564392);
   nominal__1025->SetMinimum(5e-05);
   nominal__1025->SetMaximum(2310712);
   nominal__1025->SetEntries(48.48291);
   nominal__1025->SetFillColor(38);
   nominal__1025->SetLineColor(38);
   nominal__1025->SetMarkerColor(38);
   nominal__1025->SetMarkerStyle(21);
   nominal__1025->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1025->GetXaxis()->SetRange(1,32);
   nominal__1025->GetXaxis()->SetLabelFont(42);
   nominal__1025->GetXaxis()->SetLabelSize(0.035);
   nominal__1025->GetXaxis()->SetTitleSize(0.035);
   nominal__1025->GetXaxis()->SetTitleFont(42);
   nominal__1025->GetYaxis()->SetTitle("Tracks");
   nominal__1025->GetYaxis()->SetLabelFont(42);
   nominal__1025->GetYaxis()->SetLabelSize(0.05);
   nominal__1025->GetYaxis()->SetTitleSize(0.07);
   nominal__1025->GetYaxis()->SetTitleOffset(0);
   nominal__1025->GetYaxis()->SetTitleFont(42);
   nominal__1025->GetZaxis()->SetLabelFont(42);
   nominal__1025->GetZaxis()->SetLabelSize(0.035);
   nominal__1025->GetZaxis()->SetTitleSize(0.035);
   nominal__1025->GetZaxis()->SetTitleFont(42);
   nominal__1025->Draw("E0 same");
   Double_t xAxis798[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1026 = new TH1F("nominal__1026","",32, xAxis798);
   nominal__1026->SetBinContent(2,1);
   nominal__1026->SetBinContent(4,1);
   nominal__1026->SetBinContent(5,1);
   nominal__1026->SetBinContent(7,1);
   nominal__1026->SetBinContent(8,1);
   nominal__1026->SetBinContent(9,1);
   nominal__1026->SetBinContent(10,1);
   nominal__1026->SetBinContent(11,1);
   nominal__1026->SetBinContent(12,1);
   nominal__1026->SetBinContent(13,1);
   nominal__1026->SetBinContent(14,1);
   nominal__1026->SetBinContent(15,1);
   nominal__1026->SetBinContent(16,1);
   nominal__1026->SetBinContent(17,1);
   nominal__1026->SetBinContent(18,0.9469944);
   nominal__1026->SetBinContent(19,0.9743035);
   nominal__1026->SetBinContent(20,0.9700989);
   nominal__1026->SetBinContent(21,0.9640819);
   nominal__1026->SetBinContent(22,0.9909688);
   nominal__1026->SetBinContent(23,0.9743687);
   nominal__1026->SetBinContent(24,0.9807742);
   nominal__1026->SetBinContent(25,0.9972732);
   nominal__1026->SetBinContent(26,1.018059);
   nominal__1026->SetBinContent(27,1.030293);
   nominal__1026->SetBinContent(28,1.010157);
   nominal__1026->SetBinContent(29,1.02942);
   nominal__1026->SetBinContent(30,0.9999999);
   nominal__1026->SetBinContent(31,1);
   nominal__1026->SetBinContent(32,1);
   nominal__1026->SetBinError(2,1.414214);
   nominal__1026->SetBinError(4,1.414214);
   nominal__1026->SetBinError(5,0.7097162);
   nominal__1026->SetBinError(7,1.414214);
   nominal__1026->SetBinError(8,1.00158);
   nominal__1026->SetBinError(9,1.414214);
   nominal__1026->SetBinError(10,1.414213);
   nominal__1026->SetBinError(11,1.414214);
   nominal__1026->SetBinError(12,0.7074821);
   nominal__1026->SetBinError(13,1.000099);
   nominal__1026->SetBinError(14,0.7071538);
   nominal__1026->SetBinError(15,0.8201866);
   nominal__1026->SetBinError(16,0.5353192);
   nominal__1026->SetBinError(17,0.3437137);
   nominal__1026->SetBinError(18,0.3117973);
   nominal__1026->SetBinError(19,0.2224449);
   nominal__1026->SetBinError(20,0.1068944);
   nominal__1026->SetBinError(21,0.0704521);
   nominal__1026->SetBinError(22,0.0421707);
   nominal__1026->SetBinError(23,0.02333348);
   nominal__1026->SetBinError(24,0.01363168);
   nominal__1026->SetBinError(25,0.01037704);
   nominal__1026->SetBinError(26,0.01187627);
   nominal__1026->SetBinError(27,0.02521252);
   nominal__1026->SetBinError(28,0.06006354);
   nominal__1026->SetBinError(29,0.1440909);
   nominal__1026->SetBinError(30,0.3167737);
   nominal__1026->SetBinError(31,0.7075212);
   nominal__1026->SetBinError(32,0.3564392);
   nominal__1026->SetMinimum(5e-05);
   nominal__1026->SetMaximum(2310712);
   nominal__1026->SetEntries(47.67205);
   nominal__1026->SetFillColor(46);
   nominal__1026->SetLineColor(46);
   nominal__1026->SetMarkerColor(46);
   nominal__1026->SetMarkerStyle(21);
   nominal__1026->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1026->GetXaxis()->SetRange(1,32);
   nominal__1026->GetXaxis()->SetLabelFont(42);
   nominal__1026->GetXaxis()->SetLabelSize(0.035);
   nominal__1026->GetXaxis()->SetTitleSize(0.035);
   nominal__1026->GetXaxis()->SetTitleFont(42);
   nominal__1026->GetYaxis()->SetTitle("Tracks");
   nominal__1026->GetYaxis()->SetLabelFont(42);
   nominal__1026->GetYaxis()->SetLabelSize(0.05);
   nominal__1026->GetYaxis()->SetTitleSize(0.07);
   nominal__1026->GetYaxis()->SetTitleOffset(0);
   nominal__1026->GetYaxis()->SetTitleFont(42);
   nominal__1026->GetZaxis()->SetLabelFont(42);
   nominal__1026->GetZaxis()->SetLabelSize(0.035);
   nominal__1026->GetZaxis()->SetTitleSize(0.035);
   nominal__1026->GetZaxis()->SetTitleFont(42);
   nominal__1026->Draw("E0 same");
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
