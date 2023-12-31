void plot_Trigger()
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
   t1->Range(-428.5714,-4.350845,2428.571,5.612124);
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
   Double_t xAxis281[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__361 = new TH1F("nominal__361","",32, xAxis281);
   nominal__361->SetBinContent(14,0.001057806);
   nominal__361->SetBinContent(20,0.00106468);
   nominal__361->SetBinContent(26,0.00816225);
   nominal__361->SetBinContent(27,0.008581624);
   nominal__361->SetBinContent(28,0.02230647);
   nominal__361->SetBinContent(29,0.1038955);
   nominal__361->SetBinContent(30,0.4463632);
   nominal__361->SetBinContent(31,1.27958);
   nominal__361->SetBinContent(32,3.650138);
   nominal__361->SetBinError(14,0.001057806);
   nominal__361->SetBinError(20,0.00106468);
   nominal__361->SetBinError(26,0.002889091);
   nominal__361->SetBinError(27,0.003038866);
   nominal__361->SetBinError(28,0.004882737);
   nominal__361->SetBinError(29,0.01045514);
   nominal__361->SetBinError(30,0.02174857);
   nominal__361->SetBinError(31,0.03673089);
   nominal__361->SetBinError(32,0.03803299);
   nominal__361->SetBinError(33,0.04904298);
   nominal__361->SetMinimum(5e-05);
   nominal__361->SetMaximum(365013.8);
   nominal__361->SetEntries(5256);
   nominal__361->SetFillColor(1);
   nominal__361->SetMarkerStyle(20);
   nominal__361->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__361->GetXaxis()->SetRange(1,32);
   nominal__361->GetXaxis()->SetLabelFont(42);
   nominal__361->GetXaxis()->SetLabelSize(0.035);
   nominal__361->GetXaxis()->SetTitleSize(0.035);
   nominal__361->GetXaxis()->SetTitleFont(42);
   nominal__361->GetYaxis()->SetTitle("Tracks");
   nominal__361->GetYaxis()->SetLabelFont(42);
   nominal__361->GetYaxis()->SetLabelSize(0.05);
   nominal__361->GetYaxis()->SetTitleSize(0.07);
   nominal__361->GetYaxis()->SetTitleOffset(0);
   nominal__361->GetYaxis()->SetTitleFont(42);
   nominal__361->GetZaxis()->SetLabelFont(42);
   nominal__361->GetZaxis()->SetLabelSize(0.035);
   nominal__361->GetZaxis()->SetTitleSize(0.035);
   nominal__361->GetZaxis()->SetTitleFont(42);
   nominal__361->Draw("");
   Double_t xAxis282[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__362 = new TH1F("Trigger_down__362","",32, xAxis282);
   Trigger_down__362->SetBinContent(14,0.001131853);
   Trigger_down__362->SetBinContent(20,0.001011446);
   Trigger_down__362->SetBinContent(26,0.007447215);
   Trigger_down__362->SetBinContent(27,0.008015553);
   Trigger_down__362->SetBinContent(28,0.02081024);
   Trigger_down__362->SetBinContent(29,0.095031);
   Trigger_down__362->SetBinContent(30,0.3996986);
   Trigger_down__362->SetBinContent(31,1.115402);
   Trigger_down__362->SetBinContent(32,3.063723);
   Trigger_down__362->SetBinError(14,0.001131853);
   Trigger_down__362->SetBinError(20,0.001011446);
   Trigger_down__362->SetBinError(26,0.002640547);
   Trigger_down__362->SetBinError(27,0.002875978);
   Trigger_down__362->SetBinError(28,0.004578572);
   Trigger_down__362->SetBinError(29,0.009660917);
   Trigger_down__362->SetBinError(30,0.01964589);
   Trigger_down__362->SetBinError(31,0.03230201);
   Trigger_down__362->SetBinError(32,0.03239663);
   Trigger_down__362->SetBinError(33,0.0412493);
   Trigger_down__362->SetEntries(5256);
   Trigger_down__362->SetFillColor(38);
   Trigger_down__362->SetLineColor(38);
   Trigger_down__362->SetMarkerColor(38);
   Trigger_down__362->SetMarkerStyle(21);
   Trigger_down__362->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__362->GetXaxis()->SetRange(1,400);
   Trigger_down__362->GetXaxis()->SetLabelFont(42);
   Trigger_down__362->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__362->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__362->GetXaxis()->SetTitleFont(42);
   Trigger_down__362->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__362->GetYaxis()->SetLabelFont(42);
   Trigger_down__362->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__362->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__362->GetYaxis()->SetTitleOffset(0);
   Trigger_down__362->GetYaxis()->SetTitleFont(42);
   Trigger_down__362->GetZaxis()->SetLabelFont(42);
   Trigger_down__362->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__362->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__362->GetZaxis()->SetTitleFont(42);
   Trigger_down__362->Draw("same");
   Double_t xAxis283[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__363 = new TH1F("Trigger_up__363","",32, xAxis283);
   Trigger_up__363->SetBinContent(14,0.001068384);
   Trigger_up__363->SetBinContent(20,0.000958212);
   Trigger_up__363->SetBinContent(26,0.008177425);
   Trigger_up__363->SetBinContent(27,0.0084982);
   Trigger_up__363->SetBinContent(28,0.02205073);
   Trigger_up__363->SetBinContent(29,0.103631);
   Trigger_up__363->SetBinContent(30,0.4470586);
   Trigger_up__363->SetBinContent(31,1.287152);
   Trigger_up__363->SetBinContent(32,3.655972);
   Trigger_up__363->SetBinError(14,0.001068384);
   Trigger_up__363->SetBinError(20,0.000958212);
   Trigger_up__363->SetBinError(26,0.002893958);
   Trigger_up__363->SetBinError(27,0.003010424);
   Trigger_up__363->SetBinError(28,0.004828414);
   Trigger_up__363->SetBinError(29,0.01043612);
   Trigger_up__363->SetBinError(30,0.02180272);
   Trigger_up__363->SetBinError(31,0.03699088);
   Trigger_up__363->SetBinError(32,0.03821124);
   Trigger_up__363->SetBinError(33,0.04911632);
   Trigger_up__363->SetEntries(5256);
   Trigger_up__363->SetFillColor(46);
   Trigger_up__363->SetLineColor(46);
   Trigger_up__363->SetMarkerColor(46);
   Trigger_up__363->SetMarkerStyle(21);
   Trigger_up__363->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__363->GetXaxis()->SetRange(1,400);
   Trigger_up__363->GetXaxis()->SetLabelFont(42);
   Trigger_up__363->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__363->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__363->GetXaxis()->SetTitleFont(42);
   Trigger_up__363->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__363->GetYaxis()->SetLabelFont(42);
   Trigger_up__363->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__363->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__363->GetYaxis()->SetTitleOffset(0);
   Trigger_up__363->GetYaxis()->SetTitleFont(42);
   Trigger_up__363->GetZaxis()->SetLabelFont(42);
   Trigger_up__363->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__363->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__363->GetZaxis()->SetTitleFont(42);
   Trigger_up__363->Draw("same");
   TLine *line = new TLine(1730,0,1730,365013.8);
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
   entry=leg->AddEntry("Trigger_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Trigger_up","Up","PE1");
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
   
   TH1D *frameR2__364 = new TH1D("frameR2__364","",1,0,2000);
   frameR2__364->SetMinimum(0.5);
   frameR2__364->SetMaximum(1.5);
   frameR2__364->SetStats(0);
   frameR2__364->SetLineStyle(0);
   frameR2__364->SetMarkerStyle(20);
   frameR2__364->GetXaxis()->SetRange(1,1);
   frameR2__364->GetXaxis()->SetLabelFont(43);
   frameR2__364->GetXaxis()->SetLabelOffset(0.007);
   frameR2__364->GetXaxis()->SetLabelSize(16);
   frameR2__364->GetXaxis()->SetTitleSize(24);
   frameR2__364->GetXaxis()->SetTitleOffset(3.75);
   frameR2__364->GetXaxis()->SetTitleFont(43);
   frameR2__364->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__364->GetYaxis()->SetNdivisions(205);
   frameR2__364->GetYaxis()->SetLabelFont(43);
   frameR2__364->GetYaxis()->SetLabelOffset(0.007);
   frameR2__364->GetYaxis()->SetLabelSize(20);
   frameR2__364->GetYaxis()->SetTitleSize(20);
   frameR2__364->GetYaxis()->SetTitleOffset(2);
   frameR2__364->GetYaxis()->SetTitleFont(43);
   frameR2__364->GetZaxis()->SetLabelFont(42);
   frameR2__364->GetZaxis()->SetLabelOffset(0.007);
   frameR2__364->GetZaxis()->SetLabelSize(0.05);
   frameR2__364->GetZaxis()->SetTitleSize(0.06);
   frameR2__364->GetZaxis()->SetTitleFont(42);
   frameR2__364->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis284[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__365 = new TH1F("Trigger_down__365","",32, xAxis284);
   Trigger_down__365->SetBinContent(0,1.171654);
   Trigger_down__365->SetBinContent(1,1.171654);
   Trigger_down__365->SetBinContent(2,1.171654);
   Trigger_down__365->SetBinContent(3,1.171654);
   Trigger_down__365->SetBinContent(4,1.171654);
   Trigger_down__365->SetBinContent(5,1.171654);
   Trigger_down__365->SetBinContent(6,1.171654);
   Trigger_down__365->SetBinContent(7,1.171654);
   Trigger_down__365->SetBinContent(8,1.171654);
   Trigger_down__365->SetBinContent(9,1.171654);
   Trigger_down__365->SetBinContent(10,1.171654);
   Trigger_down__365->SetBinContent(11,1.171654);
   Trigger_down__365->SetBinContent(12,1.171654);
   Trigger_down__365->SetBinContent(13,1.171654);
   Trigger_down__365->SetBinContent(14,1.171654);
   Trigger_down__365->SetBinContent(15,1.171711);
   Trigger_down__365->SetBinContent(16,1.171711);
   Trigger_down__365->SetBinContent(17,1.171711);
   Trigger_down__365->SetBinContent(18,1.171711);
   Trigger_down__365->SetBinContent(19,1.171711);
   Trigger_down__365->SetBinContent(20,1.171711);
   Trigger_down__365->SetBinContent(21,1.171736);
   Trigger_down__365->SetBinContent(22,1.171736);
   Trigger_down__365->SetBinContent(23,1.171736);
   Trigger_down__365->SetBinContent(24,1.171736);
   Trigger_down__365->SetBinContent(25,1.171736);
   Trigger_down__365->SetBinContent(26,1.171736);
   Trigger_down__365->SetBinContent(27,1.171856);
   Trigger_down__365->SetBinContent(28,1.172029);
   Trigger_down__365->SetBinContent(29,1.172475);
   Trigger_down__365->SetBinContent(30,1.174119);
   Trigger_down__365->SetBinContent(31,1.179605);
   Trigger_down__365->SetBinContent(32,1.191406);
   Trigger_down__365->SetBinError(0,0.02299993);
   Trigger_down__365->SetBinError(1,0.02299993);
   Trigger_down__365->SetBinError(2,0.02299993);
   Trigger_down__365->SetBinError(3,0.02299993);
   Trigger_down__365->SetBinError(4,0.02299993);
   Trigger_down__365->SetBinError(5,0.02299993);
   Trigger_down__365->SetBinError(6,0.02299993);
   Trigger_down__365->SetBinError(7,0.02299993);
   Trigger_down__365->SetBinError(8,0.02299993);
   Trigger_down__365->SetBinError(9,0.02299993);
   Trigger_down__365->SetBinError(10,0.02299993);
   Trigger_down__365->SetBinError(11,0.02299993);
   Trigger_down__365->SetBinError(12,0.02299993);
   Trigger_down__365->SetBinError(13,0.02299993);
   Trigger_down__365->SetBinError(14,0.02299993);
   Trigger_down__365->SetBinError(15,0.0230032);
   Trigger_down__365->SetBinError(16,0.0230032);
   Trigger_down__365->SetBinError(17,0.0230032);
   Trigger_down__365->SetBinError(18,0.0230032);
   Trigger_down__365->SetBinError(19,0.0230032);
   Trigger_down__365->SetBinError(20,0.0230032);
   Trigger_down__365->SetBinError(21,0.02300591);
   Trigger_down__365->SetBinError(22,0.02300591);
   Trigger_down__365->SetBinError(23,0.02300591);
   Trigger_down__365->SetBinError(24,0.02300591);
   Trigger_down__365->SetBinError(25,0.02300591);
   Trigger_down__365->SetBinError(26,0.02300591);
   Trigger_down__365->SetBinError(27,0.02302594);
   Trigger_down__365->SetBinError(28,0.0230467);
   Trigger_down__365->SetBinError(29,0.02310164);
   Trigger_down__365->SetBinError(30,0.02335602);
   Trigger_down__365->SetBinError(31,0.0244969);
   Trigger_down__365->SetBinError(32,0.02874675);
   Trigger_down__365->SetEntries(33);
   Trigger_down__365->SetFillColor(38);
   Trigger_down__365->SetLineColor(38);
   Trigger_down__365->SetMarkerColor(38);
   Trigger_down__365->SetMarkerStyle(21);
   Trigger_down__365->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__365->GetXaxis()->SetRange(1,400);
   Trigger_down__365->GetXaxis()->SetLabelFont(42);
   Trigger_down__365->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__365->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__365->GetXaxis()->SetTitleFont(42);
   Trigger_down__365->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__365->GetYaxis()->SetLabelFont(42);
   Trigger_down__365->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__365->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__365->GetYaxis()->SetTitleOffset(0);
   Trigger_down__365->GetYaxis()->SetTitleFont(42);
   Trigger_down__365->GetZaxis()->SetLabelFont(42);
   Trigger_down__365->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__365->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__365->GetZaxis()->SetTitleFont(42);
   Trigger_down__365->Draw("E0 same");
   Double_t xAxis285[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__366 = new TH1F("Trigger_up__366","",32, xAxis285);
   Trigger_up__366->SetBinContent(0,0.9975758);
   Trigger_up__366->SetBinContent(1,0.9975758);
   Trigger_up__366->SetBinContent(2,0.9975758);
   Trigger_up__366->SetBinContent(3,0.9975758);
   Trigger_up__366->SetBinContent(4,0.9975758);
   Trigger_up__366->SetBinContent(5,0.9975758);
   Trigger_up__366->SetBinContent(6,0.9975758);
   Trigger_up__366->SetBinContent(7,0.9975758);
   Trigger_up__366->SetBinContent(8,0.9975758);
   Trigger_up__366->SetBinContent(9,0.9975758);
   Trigger_up__366->SetBinContent(10,0.9975758);
   Trigger_up__366->SetBinContent(11,0.9975758);
   Trigger_up__366->SetBinContent(12,0.9975758);
   Trigger_up__366->SetBinContent(13,0.9975758);
   Trigger_up__366->SetBinContent(14,0.9975758);
   Trigger_up__366->SetBinContent(15,0.9975773);
   Trigger_up__366->SetBinContent(16,0.9975773);
   Trigger_up__366->SetBinContent(17,0.9975773);
   Trigger_up__366->SetBinContent(18,0.9975773);
   Trigger_up__366->SetBinContent(19,0.9975773);
   Trigger_up__366->SetBinContent(20,0.9975773);
   Trigger_up__366->SetBinContent(21,0.9975576);
   Trigger_up__366->SetBinContent(22,0.9975576);
   Trigger_up__366->SetBinContent(23,0.9975576);
   Trigger_up__366->SetBinContent(24,0.9975576);
   Trigger_up__366->SetBinContent(25,0.9975576);
   Trigger_up__366->SetBinContent(26,0.9975576);
   Trigger_up__366->SetBinContent(27,0.9975567);
   Trigger_up__366->SetBinContent(28,0.9975378);
   Trigger_up__366->SetBinContent(29,0.9974813);
   Trigger_up__366->SetBinContent(30,0.9973839);
   Trigger_up__366->SetBinContent(31,0.997288);
   Trigger_up__366->SetBinContent(32,0.9984042);
   Trigger_up__366->SetBinError(0,0.0195157);
   Trigger_up__366->SetBinError(1,0.0195157);
   Trigger_up__366->SetBinError(2,0.0195157);
   Trigger_up__366->SetBinError(3,0.0195157);
   Trigger_up__366->SetBinError(4,0.0195157);
   Trigger_up__366->SetBinError(5,0.0195157);
   Trigger_up__366->SetBinError(6,0.0195157);
   Trigger_up__366->SetBinError(7,0.0195157);
   Trigger_up__366->SetBinError(8,0.0195157);
   Trigger_up__366->SetBinError(9,0.0195157);
   Trigger_up__366->SetBinError(10,0.0195157);
   Trigger_up__366->SetBinError(11,0.0195157);
   Trigger_up__366->SetBinError(12,0.0195157);
   Trigger_up__366->SetBinError(13,0.0195157);
   Trigger_up__366->SetBinError(14,0.0195157);
   Trigger_up__366->SetBinError(15,0.0195176);
   Trigger_up__366->SetBinError(16,0.0195176);
   Trigger_up__366->SetBinError(17,0.0195176);
   Trigger_up__366->SetBinError(18,0.0195176);
   Trigger_up__366->SetBinError(19,0.0195176);
   Trigger_up__366->SetBinError(20,0.0195176);
   Trigger_up__366->SetBinError(21,0.01951907);
   Trigger_up__366->SetBinError(22,0.01951907);
   Trigger_up__366->SetBinError(23,0.01951907);
   Trigger_up__366->SetBinError(24,0.01951907);
   Trigger_up__366->SetBinError(25,0.01951907);
   Trigger_up__366->SetBinError(26,0.01951907);
   Trigger_up__366->SetBinError(27,0.01953397);
   Trigger_up__366->SetBinError(28,0.01954854);
   Trigger_up__366->SetBinError(29,0.01958671);
   Trigger_up__366->SetBinError(30,0.01977377);
   Trigger_up__366->SetBinError(31,0.02064397);
   Trigger_up__366->SetBinError(32,0.02402024);
   Trigger_up__366->SetEntries(33);
   Trigger_up__366->SetFillColor(46);
   Trigger_up__366->SetLineColor(46);
   Trigger_up__366->SetMarkerColor(46);
   Trigger_up__366->SetMarkerStyle(21);
   Trigger_up__366->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__366->GetXaxis()->SetRange(1,400);
   Trigger_up__366->GetXaxis()->SetLabelFont(42);
   Trigger_up__366->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__366->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__366->GetXaxis()->SetTitleFont(42);
   Trigger_up__366->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__366->GetYaxis()->SetLabelFont(42);
   Trigger_up__366->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__366->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__366->GetYaxis()->SetTitleOffset(0);
   Trigger_up__366->GetYaxis()->SetTitleFont(42);
   Trigger_up__366->GetZaxis()->SetLabelFont(42);
   Trigger_up__366->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__366->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__366->GetZaxis()->SetTitleFont(42);
   Trigger_up__366->Draw("E0 same");
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
   
   TH1D *frameR2__367 = new TH1D("frameR2__367","",1,0,2000);
   frameR2__367->SetMinimum(0.5);
   frameR2__367->SetMaximum(1.5);
   frameR2__367->SetStats(0);
   frameR2__367->SetLineStyle(0);
   frameR2__367->SetMarkerStyle(20);
   frameR2__367->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__367->GetXaxis()->SetRange(1,1);
   frameR2__367->GetXaxis()->SetLabelFont(43);
   frameR2__367->GetXaxis()->SetLabelOffset(0.007);
   frameR2__367->GetXaxis()->SetLabelSize(16);
   frameR2__367->GetXaxis()->SetTitleSize(24);
   frameR2__367->GetXaxis()->SetTitleOffset(5);
   frameR2__367->GetXaxis()->SetTitleFont(43);
   frameR2__367->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__367->GetYaxis()->SetNdivisions(205);
   frameR2__367->GetYaxis()->SetLabelFont(43);
   frameR2__367->GetYaxis()->SetLabelOffset(0.007);
   frameR2__367->GetYaxis()->SetLabelSize(20);
   frameR2__367->GetYaxis()->SetTitleSize(20);
   frameR2__367->GetYaxis()->SetTitleOffset(2);
   frameR2__367->GetYaxis()->SetTitleFont(43);
   frameR2__367->GetZaxis()->SetLabelFont(42);
   frameR2__367->GetZaxis()->SetLabelOffset(0.007);
   frameR2__367->GetZaxis()->SetLabelSize(0.05);
   frameR2__367->GetZaxis()->SetTitleSize(0.06);
   frameR2__367->GetZaxis()->SetTitleFont(42);
   frameR2__367->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis286[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__368 = new TH1F("nominal__368","",32, xAxis286);
   nominal__368->SetBinContent(14,0.9345794);
   nominal__368->SetBinContent(20,1.052632);
   nominal__368->SetBinContent(26,1.096014);
   nominal__368->SetBinContent(27,1.070622);
   nominal__368->SetBinContent(28,1.071899);
   nominal__368->SetBinContent(29,1.09328);
   nominal__368->SetBinContent(30,1.11675);
   nominal__368->SetBinContent(31,1.147192);
   nominal__368->SetBinContent(32,1.191406);
   nominal__368->SetBinError(14,1.321695);
   nominal__368->SetBinError(20,1.488646);
   nominal__368->SetBinError(26,0.549107);
   nominal__368->SetBinError(27,0.5397178);
   nominal__368->SetBinError(28,0.3326706);
   nominal__368->SetBinError(29,0.156387);
   nominal__368->SetBinError(30,0.07728936);
   nominal__368->SetBinError(31,0.04677788);
   nominal__368->SetBinError(32,0.0176868);
   nominal__368->SetMinimum(5e-05);
   nominal__368->SetMaximum(365013.8);
   nominal__368->SetEntries(20.33008);
   nominal__368->SetFillColor(38);
   nominal__368->SetLineColor(38);
   nominal__368->SetMarkerColor(38);
   nominal__368->SetMarkerStyle(21);
   nominal__368->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__368->GetXaxis()->SetRange(1,32);
   nominal__368->GetXaxis()->SetLabelFont(42);
   nominal__368->GetXaxis()->SetLabelSize(0.035);
   nominal__368->GetXaxis()->SetTitleSize(0.035);
   nominal__368->GetXaxis()->SetTitleFont(42);
   nominal__368->GetYaxis()->SetTitle("Tracks");
   nominal__368->GetYaxis()->SetLabelFont(42);
   nominal__368->GetYaxis()->SetLabelSize(0.05);
   nominal__368->GetYaxis()->SetTitleSize(0.07);
   nominal__368->GetYaxis()->SetTitleOffset(0);
   nominal__368->GetYaxis()->SetTitleFont(42);
   nominal__368->GetZaxis()->SetLabelFont(42);
   nominal__368->GetZaxis()->SetLabelSize(0.035);
   nominal__368->GetZaxis()->SetTitleSize(0.035);
   nominal__368->GetZaxis()->SetTitleFont(42);
   nominal__368->Draw("E0 same");
   Double_t xAxis287[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__369 = new TH1F("nominal__369","",32, xAxis287);
   nominal__369->SetBinContent(14,0.990099);
   nominal__369->SetBinContent(20,1.111111);
   nominal__369->SetBinContent(26,0.9981443);
   nominal__369->SetBinContent(27,1.009817);
   nominal__369->SetBinContent(28,1.011598);
   nominal__369->SetBinContent(29,1.002552);
   nominal__369->SetBinContent(30,0.9984444);
   nominal__369->SetBinContent(31,0.9941174);
   nominal__369->SetBinContent(32,0.9984042);
   nominal__369->SetBinError(14,1.400211);
   nominal__369->SetBinError(20,1.571348);
   nominal__369->SetBinError(26,0.4995993);
   nominal__369->SetBinError(27,0.5058001);
   nominal__369->SetBinError(28,0.3132059);
   nominal__369->SetBinError(29,0.1427294);
   nominal__369->SetBinError(30,0.06883087);
   nominal__369->SetBinError(31,0.04038007);
   nominal__369->SetBinError(32,0.01473473);
   nominal__369->SetMinimum(5e-05);
   nominal__369->SetMaximum(365013.8);
   nominal__369->SetEntries(16.41634);
   nominal__369->SetFillColor(46);
   nominal__369->SetLineColor(46);
   nominal__369->SetMarkerColor(46);
   nominal__369->SetMarkerStyle(21);
   nominal__369->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__369->GetXaxis()->SetRange(1,32);
   nominal__369->GetXaxis()->SetLabelFont(42);
   nominal__369->GetXaxis()->SetLabelSize(0.035);
   nominal__369->GetXaxis()->SetTitleSize(0.035);
   nominal__369->GetXaxis()->SetTitleFont(42);
   nominal__369->GetYaxis()->SetTitle("Tracks");
   nominal__369->GetYaxis()->SetLabelFont(42);
   nominal__369->GetYaxis()->SetLabelSize(0.05);
   nominal__369->GetYaxis()->SetTitleSize(0.07);
   nominal__369->GetYaxis()->SetTitleOffset(0);
   nominal__369->GetYaxis()->SetTitleFont(42);
   nominal__369->GetZaxis()->SetLabelFont(42);
   nominal__369->GetZaxis()->SetLabelSize(0.035);
   nominal__369->GetZaxis()->SetTitleSize(0.035);
   nominal__369->GetZaxis()->SetTitleFont(42);
   nominal__369->Draw("E0 same");
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
