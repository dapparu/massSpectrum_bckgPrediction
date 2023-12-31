void plot_PU()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:01 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.328416,2428.571,1.148829);
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
   Double_t xAxis316[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__406 = new TH1F("nominal__406","",32, xAxis316);
   nominal__406->SetBinContent(18,0.0003669347);
   nominal__406->SetBinContent(21,0.0003689668);
   nominal__406->SetBinContent(24,0.0003709604);
   nominal__406->SetBinContent(26,0.0003891264);
   nominal__406->SetBinContent(27,0.002262764);
   nominal__406->SetBinContent(28,0.003053271);
   nominal__406->SetBinContent(29,0.01591948);
   nominal__406->SetBinContent(30,0.06855007);
   nominal__406->SetBinContent(31,0.2704401);
   nominal__406->SetBinContent(32,1.322644);
   nominal__406->SetBinError(18,0.0003669347);
   nominal__406->SetBinError(21,0.0003689668);
   nominal__406->SetBinError(24,0.0003709604);
   nominal__406->SetBinError(26,0.0003891264);
   nominal__406->SetBinError(27,0.0009241219);
   nominal__406->SetBinError(28,0.001082269);
   nominal__406->SetBinError(29,0.002434662);
   nominal__406->SetBinError(30,0.005025268);
   nominal__406->SetBinError(31,0.01001539);
   nominal__406->SetBinError(32,0.01178473);
   nominal__406->SetBinError(33,0.01879447);
   nominal__406->SetMinimum(5e-05);
   nominal__406->SetMaximum(13.22644);
   nominal__406->SetEntries(4553);
   nominal__406->SetFillColor(1);
   nominal__406->SetMarkerStyle(20);
   nominal__406->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__406->GetXaxis()->SetRange(1,32);
   nominal__406->GetXaxis()->SetLabelFont(42);
   nominal__406->GetXaxis()->SetLabelSize(0.035);
   nominal__406->GetXaxis()->SetTitleSize(0.035);
   nominal__406->GetXaxis()->SetTitleFont(42);
   nominal__406->GetYaxis()->SetTitle("Tracks");
   nominal__406->GetYaxis()->SetLabelFont(42);
   nominal__406->GetYaxis()->SetLabelSize(0.05);
   nominal__406->GetYaxis()->SetTitleSize(0.07);
   nominal__406->GetYaxis()->SetTitleOffset(0);
   nominal__406->GetYaxis()->SetTitleFont(42);
   nominal__406->GetZaxis()->SetLabelFont(42);
   nominal__406->GetZaxis()->SetLabelSize(0.035);
   nominal__406->GetZaxis()->SetTitleSize(0.035);
   nominal__406->GetZaxis()->SetTitleFont(42);
   nominal__406->Draw("");
   Double_t xAxis317[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__407 = new TH1F("PU_down__407","",32, xAxis317);
   PU_down__407->SetBinContent(18,0.0003365517);
   PU_down__407->SetBinContent(21,0.0003417139);
   PU_down__407->SetBinContent(24,0.0003440794);
   PU_down__407->SetBinContent(26,0.0004343394);
   PU_down__407->SetBinContent(27,0.002249045);
   PU_down__407->SetBinContent(28,0.003392647);
   PU_down__407->SetBinContent(29,0.01617095);
   PU_down__407->SetBinContent(30,0.06678736);
   PU_down__407->SetBinContent(31,0.2651918);
   PU_down__407->SetBinContent(32,1.315619);
   PU_down__407->SetBinError(18,0.0003365518);
   PU_down__407->SetBinError(21,0.0003417139);
   PU_down__407->SetBinError(24,0.0003440794);
   PU_down__407->SetBinError(26,0.0004343394);
   PU_down__407->SetBinError(27,0.0009270628);
   PU_down__407->SetBinError(28,0.001292322);
   PU_down__407->SetBinError(29,0.002557165);
   PU_down__407->SetBinError(30,0.005030164);
   PU_down__407->SetBinError(31,0.01009504);
   PU_down__407->SetBinError(32,0.01210306);
   PU_down__407->SetBinError(33,0.01929053);
   PU_down__407->SetEntries(4553);
   PU_down__407->SetFillColor(38);
   PU_down__407->SetLineColor(38);
   PU_down__407->SetMarkerColor(38);
   PU_down__407->SetMarkerStyle(21);
   PU_down__407->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__407->GetXaxis()->SetRange(1,400);
   PU_down__407->GetXaxis()->SetLabelFont(42);
   PU_down__407->GetXaxis()->SetLabelSize(0.035);
   PU_down__407->GetXaxis()->SetTitleSize(0.035);
   PU_down__407->GetXaxis()->SetTitleFont(42);
   PU_down__407->GetYaxis()->SetTitle("Events / bin");
   PU_down__407->GetYaxis()->SetLabelFont(42);
   PU_down__407->GetYaxis()->SetLabelSize(0.035);
   PU_down__407->GetYaxis()->SetTitleSize(0.035);
   PU_down__407->GetYaxis()->SetTitleOffset(0);
   PU_down__407->GetYaxis()->SetTitleFont(42);
   PU_down__407->GetZaxis()->SetLabelFont(42);
   PU_down__407->GetZaxis()->SetLabelSize(0.035);
   PU_down__407->GetZaxis()->SetTitleSize(0.035);
   PU_down__407->GetZaxis()->SetTitleFont(42);
   PU_down__407->Draw("same");
   Double_t xAxis318[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__408 = new TH1F("PU_up__408","",32, xAxis318);
   PU_up__408->SetBinContent(18,0.0003990454);
   PU_up__408->SetBinContent(21,0.0003987929);
   PU_up__408->SetBinContent(24,0.0003993318);
   PU_up__408->SetBinContent(26,0.0003183201);
   PU_up__408->SetBinContent(27,0.002221487);
   PU_up__408->SetBinContent(28,0.002875072);
   PU_up__408->SetBinContent(29,0.01571578);
   PU_up__408->SetBinContent(30,0.07041905);
   PU_up__408->SetBinContent(31,0.2758451);
   PU_up__408->SetBinContent(32,1.330649);
   PU_up__408->SetBinError(18,0.0003990454);
   PU_up__408->SetBinError(21,0.0003987928);
   PU_up__408->SetBinError(24,0.0003993318);
   PU_up__408->SetBinError(26,0.0003183202);
   PU_up__408->SetBinError(27,0.0009119571);
   PU_up__408->SetBinError(28,0.00103969);
   PU_up__408->SetBinError(29,0.002422587);
   PU_up__408->SetBinError(30,0.005182612);
   PU_up__408->SetBinError(31,0.0102679);
   PU_up__408->SetBinError(32,0.01198389);
   PU_up__408->SetBinError(33,0.01898579);
   PU_up__408->SetEntries(4553);
   PU_up__408->SetFillColor(46);
   PU_up__408->SetLineColor(46);
   PU_up__408->SetMarkerColor(46);
   PU_up__408->SetMarkerStyle(21);
   PU_up__408->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__408->GetXaxis()->SetRange(1,400);
   PU_up__408->GetXaxis()->SetLabelFont(42);
   PU_up__408->GetXaxis()->SetLabelSize(0.035);
   PU_up__408->GetXaxis()->SetTitleSize(0.035);
   PU_up__408->GetXaxis()->SetTitleFont(42);
   PU_up__408->GetYaxis()->SetTitle("Events / bin");
   PU_up__408->GetYaxis()->SetLabelFont(42);
   PU_up__408->GetYaxis()->SetLabelSize(0.035);
   PU_up__408->GetYaxis()->SetTitleSize(0.035);
   PU_up__408->GetYaxis()->SetTitleOffset(0);
   PU_up__408->GetYaxis()->SetTitleFont(42);
   PU_up__408->GetZaxis()->SetLabelFont(42);
   PU_up__408->GetZaxis()->SetLabelSize(0.035);
   PU_up__408->GetZaxis()->SetTitleSize(0.035);
   PU_up__408->GetZaxis()->SetTitleFont(42);
   PU_up__408->Draw("same");
   TLine *line = new TLine(1730,0,1730,13.22644);
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
   
   TH1D *frameR2__409 = new TH1D("frameR2__409","",1,0,2000);
   frameR2__409->SetMinimum(0.5);
   frameR2__409->SetMaximum(1.5);
   frameR2__409->SetStats(0);
   frameR2__409->SetLineStyle(0);
   frameR2__409->SetMarkerStyle(20);
   frameR2__409->GetXaxis()->SetRange(1,1);
   frameR2__409->GetXaxis()->SetLabelFont(43);
   frameR2__409->GetXaxis()->SetLabelOffset(0.007);
   frameR2__409->GetXaxis()->SetLabelSize(16);
   frameR2__409->GetXaxis()->SetTitleSize(24);
   frameR2__409->GetXaxis()->SetTitleOffset(3.75);
   frameR2__409->GetXaxis()->SetTitleFont(43);
   frameR2__409->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__409->GetYaxis()->SetNdivisions(205);
   frameR2__409->GetYaxis()->SetLabelFont(43);
   frameR2__409->GetYaxis()->SetLabelOffset(0.007);
   frameR2__409->GetYaxis()->SetLabelSize(20);
   frameR2__409->GetYaxis()->SetTitleSize(20);
   frameR2__409->GetYaxis()->SetTitleOffset(2);
   frameR2__409->GetYaxis()->SetTitleFont(43);
   frameR2__409->GetZaxis()->SetLabelFont(42);
   frameR2__409->GetZaxis()->SetLabelOffset(0.007);
   frameR2__409->GetZaxis()->SetLabelSize(0.05);
   frameR2__409->GetZaxis()->SetTitleSize(0.06);
   frameR2__409->GetZaxis()->SetTitleFont(42);
   frameR2__409->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis319[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__410 = new TH1F("PU_down__410","",32, xAxis319);
   PU_down__410->SetBinContent(0,1.008078);
   PU_down__410->SetBinContent(1,1.008078);
   PU_down__410->SetBinContent(2,1.008078);
   PU_down__410->SetBinContent(3,1.008078);
   PU_down__410->SetBinContent(4,1.008078);
   PU_down__410->SetBinContent(5,1.008078);
   PU_down__410->SetBinContent(6,1.008078);
   PU_down__410->SetBinContent(7,1.008078);
   PU_down__410->SetBinContent(8,1.008078);
   PU_down__410->SetBinContent(9,1.008078);
   PU_down__410->SetBinContent(10,1.008078);
   PU_down__410->SetBinContent(11,1.008078);
   PU_down__410->SetBinContent(12,1.008078);
   PU_down__410->SetBinContent(13,1.008078);
   PU_down__410->SetBinContent(14,1.008078);
   PU_down__410->SetBinContent(15,1.008078);
   PU_down__410->SetBinContent(16,1.008078);
   PU_down__410->SetBinContent(17,1.008078);
   PU_down__410->SetBinContent(18,1.008078);
   PU_down__410->SetBinContent(19,1.008062);
   PU_down__410->SetBinContent(20,1.008062);
   PU_down__410->SetBinContent(21,1.008062);
   PU_down__410->SetBinContent(22,1.008047);
   PU_down__410->SetBinContent(23,1.008047);
   PU_down__410->SetBinContent(24,1.008047);
   PU_down__410->SetBinContent(25,1.008033);
   PU_down__410->SetBinContent(26,1.008033);
   PU_down__410->SetBinContent(27,1.008062);
   PU_down__410->SetBinContent(28,1.008064);
   PU_down__410->SetBinContent(29,1.008285);
   PU_down__410->SetBinContent(30,1.008519);
   PU_down__410->SetBinContent(31,1.007764);
   PU_down__410->SetBinContent(32,1.005339);
   PU_down__410->SetBinError(0,0.02151445);
   PU_down__410->SetBinError(1,0.02151445);
   PU_down__410->SetBinError(2,0.02151445);
   PU_down__410->SetBinError(3,0.02151445);
   PU_down__410->SetBinError(4,0.02151445);
   PU_down__410->SetBinError(5,0.02151445);
   PU_down__410->SetBinError(6,0.02151445);
   PU_down__410->SetBinError(7,0.02151445);
   PU_down__410->SetBinError(8,0.02151445);
   PU_down__410->SetBinError(9,0.02151445);
   PU_down__410->SetBinError(10,0.02151445);
   PU_down__410->SetBinError(11,0.02151445);
   PU_down__410->SetBinError(12,0.02151445);
   PU_down__410->SetBinError(13,0.02151445);
   PU_down__410->SetBinError(14,0.02151445);
   PU_down__410->SetBinError(15,0.02151445);
   PU_down__410->SetBinError(16,0.02151445);
   PU_down__410->SetBinError(17,0.02151445);
   PU_down__410->SetBinError(18,0.02151445);
   PU_down__410->SetBinError(19,0.02151652);
   PU_down__410->SetBinError(20,0.02151652);
   PU_down__410->SetBinError(21,0.02151652);
   PU_down__410->SetBinError(22,0.02151864);
   PU_down__410->SetBinError(23,0.02151864);
   PU_down__410->SetBinError(24,0.02151864);
   PU_down__410->SetBinError(25,0.02152077);
   PU_down__410->SetBinError(26,0.02152077);
   PU_down__410->SetBinError(27,0.02152383);
   PU_down__410->SetBinError(28,0.02153847);
   PU_down__410->SetBinError(29,0.0215608);
   PU_down__410->SetBinError(30,0.02166836);
   PU_down__410->SetBinError(31,0.02211893);
   PU_down__410->SetBinError(32,0.0242312);
   PU_down__410->SetEntries(33);
   PU_down__410->SetFillColor(38);
   PU_down__410->SetLineColor(38);
   PU_down__410->SetMarkerColor(38);
   PU_down__410->SetMarkerStyle(21);
   PU_down__410->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__410->GetXaxis()->SetRange(1,400);
   PU_down__410->GetXaxis()->SetLabelFont(42);
   PU_down__410->GetXaxis()->SetLabelSize(0.035);
   PU_down__410->GetXaxis()->SetTitleSize(0.035);
   PU_down__410->GetXaxis()->SetTitleFont(42);
   PU_down__410->GetYaxis()->SetTitle("Events / bin");
   PU_down__410->GetYaxis()->SetLabelFont(42);
   PU_down__410->GetYaxis()->SetLabelSize(0.035);
   PU_down__410->GetYaxis()->SetTitleSize(0.035);
   PU_down__410->GetYaxis()->SetTitleOffset(0);
   PU_down__410->GetYaxis()->SetTitleFont(42);
   PU_down__410->GetZaxis()->SetLabelFont(42);
   PU_down__410->GetZaxis()->SetLabelSize(0.035);
   PU_down__410->GetZaxis()->SetTitleSize(0.035);
   PU_down__410->GetZaxis()->SetTitleFont(42);
   PU_down__410->Draw("E0 same");
   Double_t xAxis320[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__411 = new TH1F("PU_up__411","",32, xAxis320);
   PU_up__411->SetBinContent(0,0.991246);
   PU_up__411->SetBinContent(1,0.991246);
   PU_up__411->SetBinContent(2,0.991246);
   PU_up__411->SetBinContent(3,0.991246);
   PU_up__411->SetBinContent(4,0.991246);
   PU_up__411->SetBinContent(5,0.991246);
   PU_up__411->SetBinContent(6,0.991246);
   PU_up__411->SetBinContent(7,0.991246);
   PU_up__411->SetBinContent(8,0.991246);
   PU_up__411->SetBinContent(9,0.991246);
   PU_up__411->SetBinContent(10,0.991246);
   PU_up__411->SetBinContent(11,0.991246);
   PU_up__411->SetBinContent(12,0.991246);
   PU_up__411->SetBinContent(13,0.991246);
   PU_up__411->SetBinContent(14,0.991246);
   PU_up__411->SetBinContent(15,0.991246);
   PU_up__411->SetBinContent(16,0.991246);
   PU_up__411->SetBinContent(17,0.991246);
   PU_up__411->SetBinContent(18,0.991246);
   PU_up__411->SetBinContent(19,0.9912629);
   PU_up__411->SetBinContent(20,0.9912629);
   PU_up__411->SetBinContent(21,0.9912629);
   PU_up__411->SetBinContent(22,0.9912784);
   PU_up__411->SetBinContent(23,0.9912784);
   PU_up__411->SetBinContent(24,0.9912784);
   PU_up__411->SetBinContent(25,0.9912931);
   PU_up__411->SetBinContent(26,0.9912931);
   PU_up__411->SetBinContent(27,0.9912497);
   PU_up__411->SetBinContent(28,0.9912139);
   PU_up__411->SetBinContent(29,0.9910937);
   PU_up__411->SetBinContent(30,0.9908888);
   PU_up__411->SetBinContent(31,0.9916528);
   PU_up__411->SetBinContent(32,0.9939843);
   PU_up__411->SetBinError(0,0.02088689);
   PU_up__411->SetBinError(1,0.02088689);
   PU_up__411->SetBinError(2,0.02088689);
   PU_up__411->SetBinError(3,0.02088689);
   PU_up__411->SetBinError(4,0.02088689);
   PU_up__411->SetBinError(5,0.02088689);
   PU_up__411->SetBinError(6,0.02088689);
   PU_up__411->SetBinError(7,0.02088689);
   PU_up__411->SetBinError(8,0.02088689);
   PU_up__411->SetBinError(9,0.02088689);
   PU_up__411->SetBinError(10,0.02088689);
   PU_up__411->SetBinError(11,0.02088689);
   PU_up__411->SetBinError(12,0.02088689);
   PU_up__411->SetBinError(13,0.02088689);
   PU_up__411->SetBinError(14,0.02088689);
   PU_up__411->SetBinError(15,0.02088689);
   PU_up__411->SetBinError(16,0.02088689);
   PU_up__411->SetBinError(17,0.02088689);
   PU_up__411->SetBinError(18,0.02088689);
   PU_up__411->SetBinError(19,0.02088956);
   PU_up__411->SetBinError(20,0.02088956);
   PU_up__411->SetBinError(21,0.02088956);
   PU_up__411->SetBinError(22,0.0208922);
   PU_up__411->SetBinError(23,0.0208922);
   PU_up__411->SetBinError(24,0.0208922);
   PU_up__411->SetBinError(25,0.02089483);
   PU_up__411->SetBinError(26,0.02089483);
   PU_up__411->SetBinError(27,0.0208962);
   PU_up__411->SetBinError(28,0.0209093);
   PU_up__411->SetBinError(29,0.02092493);
   PU_up__411->SetBinError(30,0.02102018);
   PU_up__411->SetBinError(31,0.02148923);
   PU_up__411->SetBinError(32,0.02364747);
   PU_up__411->SetEntries(33);
   PU_up__411->SetFillColor(46);
   PU_up__411->SetLineColor(46);
   PU_up__411->SetMarkerColor(46);
   PU_up__411->SetMarkerStyle(21);
   PU_up__411->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__411->GetXaxis()->SetRange(1,400);
   PU_up__411->GetXaxis()->SetLabelFont(42);
   PU_up__411->GetXaxis()->SetLabelSize(0.035);
   PU_up__411->GetXaxis()->SetTitleSize(0.035);
   PU_up__411->GetXaxis()->SetTitleFont(42);
   PU_up__411->GetYaxis()->SetTitle("Events / bin");
   PU_up__411->GetYaxis()->SetLabelFont(42);
   PU_up__411->GetYaxis()->SetLabelSize(0.035);
   PU_up__411->GetYaxis()->SetTitleSize(0.035);
   PU_up__411->GetYaxis()->SetTitleOffset(0);
   PU_up__411->GetYaxis()->SetTitleFont(42);
   PU_up__411->GetZaxis()->SetLabelFont(42);
   PU_up__411->GetZaxis()->SetLabelSize(0.035);
   PU_up__411->GetZaxis()->SetTitleSize(0.035);
   PU_up__411->GetZaxis()->SetTitleFont(42);
   PU_up__411->Draw("E0 same");
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
   
   TH1D *frameR2__412 = new TH1D("frameR2__412","",1,0,2000);
   frameR2__412->SetMinimum(0.5);
   frameR2__412->SetMaximum(1.5);
   frameR2__412->SetStats(0);
   frameR2__412->SetLineStyle(0);
   frameR2__412->SetMarkerStyle(20);
   frameR2__412->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__412->GetXaxis()->SetRange(1,1);
   frameR2__412->GetXaxis()->SetLabelFont(43);
   frameR2__412->GetXaxis()->SetLabelOffset(0.007);
   frameR2__412->GetXaxis()->SetLabelSize(16);
   frameR2__412->GetXaxis()->SetTitleSize(24);
   frameR2__412->GetXaxis()->SetTitleOffset(5);
   frameR2__412->GetXaxis()->SetTitleFont(43);
   frameR2__412->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__412->GetYaxis()->SetNdivisions(205);
   frameR2__412->GetYaxis()->SetLabelFont(43);
   frameR2__412->GetYaxis()->SetLabelOffset(0.007);
   frameR2__412->GetYaxis()->SetLabelSize(20);
   frameR2__412->GetYaxis()->SetTitleSize(20);
   frameR2__412->GetYaxis()->SetTitleOffset(2);
   frameR2__412->GetYaxis()->SetTitleFont(43);
   frameR2__412->GetZaxis()->SetLabelFont(42);
   frameR2__412->GetZaxis()->SetLabelOffset(0.007);
   frameR2__412->GetZaxis()->SetLabelSize(0.05);
   frameR2__412->GetZaxis()->SetTitleSize(0.06);
   frameR2__412->GetZaxis()->SetTitleFont(42);
   frameR2__412->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis321[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__413 = new TH1F("nominal__413","",32, xAxis321);
   nominal__413->SetBinContent(18,1.090277);
   nominal__413->SetBinContent(21,1.079754);
   nominal__413->SetBinContent(24,1.078125);
   nominal__413->SetBinContent(26,0.8959039);
   nominal__413->SetBinContent(27,1.0061);
   nominal__413->SetBinContent(28,0.8999673);
   nominal__413->SetBinContent(29,0.9844495);
   nominal__413->SetBinContent(30,1.026393);
   nominal__413->SetBinContent(31,1.019791);
   nominal__413->SetBinContent(32,1.005339);
   nominal__413->SetBinError(18,1.541885);
   nominal__413->SetBinError(21,1.527002);
   nominal__413->SetBinError(24,1.524698);
   nominal__413->SetBinError(26,1.266999);
   nominal__413->SetBinError(27,0.5838027);
   nominal__413->SetBinError(28,0.468279);
   nominal__413->SetBinError(29,0.2165689);
   nominal__413->SetBinError(30,0.1078767);
   nominal__413->SetBinError(31,0.05416024);
   nominal__413->SetBinError(32,0.01287536);
   nominal__413->SetMinimum(5e-05);
   nominal__413->SetMaximum(13.22644);
   nominal__413->SetEntries(10.98484);
   nominal__413->SetFillColor(38);
   nominal__413->SetLineColor(38);
   nominal__413->SetMarkerColor(38);
   nominal__413->SetMarkerStyle(21);
   nominal__413->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__413->GetXaxis()->SetRange(1,32);
   nominal__413->GetXaxis()->SetLabelFont(42);
   nominal__413->GetXaxis()->SetLabelSize(0.035);
   nominal__413->GetXaxis()->SetTitleSize(0.035);
   nominal__413->GetXaxis()->SetTitleFont(42);
   nominal__413->GetYaxis()->SetTitle("Tracks");
   nominal__413->GetYaxis()->SetLabelFont(42);
   nominal__413->GetYaxis()->SetLabelSize(0.05);
   nominal__413->GetYaxis()->SetTitleSize(0.07);
   nominal__413->GetYaxis()->SetTitleOffset(0);
   nominal__413->GetYaxis()->SetTitleFont(42);
   nominal__413->GetZaxis()->SetLabelFont(42);
   nominal__413->GetZaxis()->SetLabelSize(0.035);
   nominal__413->GetZaxis()->SetTitleSize(0.035);
   nominal__413->GetZaxis()->SetTitleFont(42);
   nominal__413->Draw("E0 same");
   Double_t xAxis322[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__414 = new TH1F("nominal__414","",32, xAxis322);
   nominal__414->SetBinContent(18,0.919531);
   nominal__414->SetBinContent(21,0.9252093);
   nominal__414->SetBinContent(24,0.9289529);
   nominal__414->SetBinContent(26,1.222437);
   nominal__414->SetBinContent(27,1.018581);
   nominal__414->SetBinContent(28,1.061981);
   nominal__414->SetBinContent(29,1.012962);
   nominal__414->SetBinContent(30,0.9734592);
   nominal__414->SetBinContent(31,0.9804057);
   nominal__414->SetBinContent(32,0.9939843);
   nominal__414->SetBinError(18,1.300413);
   nominal__414->SetBinError(21,1.308443);
   nominal__414->SetBinError(24,1.313738);
   nominal__414->SetBinError(26,1.728787);
   nominal__414->SetBinError(27,0.5898259);
   nominal__414->SetBinError(28,0.5377587);
   nominal__414->SetBinError(29,0.2199588);
   nominal__414->SetBinError(30,0.1011205);
   nominal__414->SetBinError(31,0.05147901);
   nominal__414->SetBinError(32,0.01259251);
   nominal__414->SetMinimum(5e-05);
   nominal__414->SetMaximum(13.22644);
   nominal__414->SetEntries(11.42798);
   nominal__414->SetFillColor(46);
   nominal__414->SetLineColor(46);
   nominal__414->SetMarkerColor(46);
   nominal__414->SetMarkerStyle(21);
   nominal__414->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__414->GetXaxis()->SetRange(1,32);
   nominal__414->GetXaxis()->SetLabelFont(42);
   nominal__414->GetXaxis()->SetLabelSize(0.035);
   nominal__414->GetXaxis()->SetTitleSize(0.035);
   nominal__414->GetXaxis()->SetTitleFont(42);
   nominal__414->GetYaxis()->SetTitle("Tracks");
   nominal__414->GetYaxis()->SetLabelFont(42);
   nominal__414->GetYaxis()->SetLabelSize(0.05);
   nominal__414->GetYaxis()->SetTitleSize(0.07);
   nominal__414->GetYaxis()->SetTitleOffset(0);
   nominal__414->GetYaxis()->SetTitleFont(42);
   nominal__414->GetZaxis()->SetLabelFont(42);
   nominal__414->GetZaxis()->SetLabelSize(0.035);
   nominal__414->GetZaxis()->SetTitleSize(0.035);
   nominal__414->GetZaxis()->SetTitleFont(42);
   nominal__414->Draw("E0 same");
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
