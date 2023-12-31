void plot_Trigger()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:17:58 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.364001,2428.571,8.230275);
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
   Double_t xAxis29[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__37 = new TH1F("nominal__37","",32, xAxis29);
   nominal__37->SetBinContent(8,0.9293805);
   nominal__37->SetBinContent(11,0.4490763);
   nominal__37->SetBinContent(14,0.4442573);
   nominal__37->SetBinContent(17,0.4678394);
   nominal__37->SetBinContent(19,0.7901867);
   nominal__37->SetBinContent(20,0.4734108);
   nominal__37->SetBinContent(21,1.312151);
   nominal__37->SetBinContent(22,2.795414);
   nominal__37->SetBinContent(23,3.999052);
   nominal__37->SetBinContent(24,11.92961);
   nominal__37->SetBinContent(25,42.3964);
   nominal__37->SetBinContent(26,159.6213);
   nominal__37->SetBinContent(27,718.7883);
   nominal__37->SetBinContent(28,1469.955);
   nominal__37->SetBinContent(29,910.9312);
   nominal__37->SetBinContent(30,268.57);
   nominal__37->SetBinContent(31,60.78858);
   nominal__37->SetBinContent(32,38.59114);
   nominal__37->SetBinError(8,0.6574474);
   nominal__37->SetBinError(11,0.4490763);
   nominal__37->SetBinError(14,0.4442573);
   nominal__37->SetBinError(17,0.4678393);
   nominal__37->SetBinError(19,0.5617946);
   nominal__37->SetBinError(20,0.4734108);
   nominal__37->SetBinError(21,0.7587229);
   nominal__37->SetBinError(22,1.145493);
   nominal__37->SetBinError(23,1.333819);
   nominal__37->SetBinError(24,2.303587);
   nominal__37->SetBinError(25,4.35733);
   nominal__37->SetBinError(26,8.40769);
   nominal__37->SetBinError(27,17.84672);
   nominal__37->SetBinError(28,25.53434);
   nominal__37->SetBinError(29,20.10793);
   nominal__37->SetBinError(30,10.91164);
   nominal__37->SetBinError(31,5.203216);
   nominal__37->SetBinError(32,2.595201);
   nominal__37->SetBinError(33,3.232536);
   nominal__37->SetMinimum(5e-05);
   nominal__37->SetMaximum(1.469955e+08);
   nominal__37->SetEntries(8361);
   nominal__37->SetFillColor(1);
   nominal__37->SetMarkerStyle(20);
   nominal__37->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__37->GetXaxis()->SetRange(1,32);
   nominal__37->GetXaxis()->SetLabelFont(42);
   nominal__37->GetXaxis()->SetLabelSize(0.035);
   nominal__37->GetXaxis()->SetTitleSize(0.035);
   nominal__37->GetXaxis()->SetTitleFont(42);
   nominal__37->GetYaxis()->SetTitle("Tracks");
   nominal__37->GetYaxis()->SetLabelFont(42);
   nominal__37->GetYaxis()->SetLabelSize(0.05);
   nominal__37->GetYaxis()->SetTitleSize(0.07);
   nominal__37->GetYaxis()->SetTitleOffset(0);
   nominal__37->GetYaxis()->SetTitleFont(42);
   nominal__37->GetZaxis()->SetLabelFont(42);
   nominal__37->GetZaxis()->SetLabelSize(0.035);
   nominal__37->GetZaxis()->SetTitleSize(0.035);
   nominal__37->GetZaxis()->SetTitleFont(42);
   nominal__37->Draw("");
   Double_t xAxis30[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__38 = new TH1F("Trigger_down__38","",32, xAxis30);
   Trigger_down__38->SetBinContent(8,0.9436573);
   Trigger_down__38->SetBinContent(11,0.3457887);
   Trigger_down__38->SetBinContent(14,0.3598484);
   Trigger_down__38->SetBinContent(17,0.3274876);
   Trigger_down__38->SetBinContent(19,0.711168);
   Trigger_down__38->SetBinContent(20,0.4639426);
   Trigger_down__38->SetBinContent(21,1.299092);
   Trigger_down__38->SetBinContent(22,2.797081);
   Trigger_down__38->SetBinContent(23,4.000702);
   Trigger_down__38->SetBinContent(24,12.05633);
   Trigger_down__38->SetBinContent(25,41.05981);
   Trigger_down__38->SetBinContent(26,151.1552);
   Trigger_down__38->SetBinContent(27,641.7568);
   Trigger_down__38->SetBinContent(28,1260.632);
   Trigger_down__38->SetBinContent(29,778.6368);
   Trigger_down__38->SetBinContent(30,234.8779);
   Trigger_down__38->SetBinContent(31,54.40232);
   Trigger_down__38->SetBinContent(32,34.23172);
   Trigger_down__38->SetBinError(8,0.6727493);
   Trigger_down__38->SetBinError(11,0.3457887);
   Trigger_down__38->SetBinError(14,0.3598484);
   Trigger_down__38->SetBinError(17,0.3274875);
   Trigger_down__38->SetBinError(19,0.5056151);
   Trigger_down__38->SetBinError(20,0.4639425);
   Trigger_down__38->SetBinError(21,0.7562689);
   Trigger_down__38->SetBinError(22,1.150235);
   Trigger_down__38->SetBinError(23,1.342615);
   Trigger_down__38->SetBinError(24,2.34214);
   Trigger_down__38->SetBinError(25,4.246832);
   Trigger_down__38->SetBinError(26,8.004856);
   Trigger_down__38->SetBinError(27,16.05038);
   Trigger_down__38->SetBinError(28,22.05224);
   Trigger_down__38->SetBinError(29,17.28984);
   Trigger_down__38->SetBinError(30,9.608925);
   Trigger_down__38->SetBinError(31,4.675332);
   Trigger_down__38->SetBinError(32,2.370798);
   Trigger_down__38->SetBinError(33,2.88063);
   Trigger_down__38->SetEntries(8361);
   Trigger_down__38->SetFillColor(38);
   Trigger_down__38->SetLineColor(38);
   Trigger_down__38->SetMarkerColor(38);
   Trigger_down__38->SetMarkerStyle(21);
   Trigger_down__38->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__38->GetXaxis()->SetRange(1,400);
   Trigger_down__38->GetXaxis()->SetLabelFont(42);
   Trigger_down__38->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__38->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__38->GetXaxis()->SetTitleFont(42);
   Trigger_down__38->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__38->GetYaxis()->SetLabelFont(42);
   Trigger_down__38->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__38->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__38->GetYaxis()->SetTitleOffset(0);
   Trigger_down__38->GetYaxis()->SetTitleFont(42);
   Trigger_down__38->GetZaxis()->SetLabelFont(42);
   Trigger_down__38->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__38->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__38->GetZaxis()->SetTitleFont(42);
   Trigger_down__38->Draw("same");
   Double_t xAxis31[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__39 = new TH1F("Trigger_up__39","",32, xAxis31);
   Trigger_up__39->SetBinContent(8,0.9009603);
   Trigger_up__39->SetBinContent(11,0.4400948);
   Trigger_up__39->SetBinContent(14,0.4620275);
   Trigger_up__39->SetBinContent(17,0.5005881);
   Trigger_up__39->SetBinContent(19,0.7901867);
   Trigger_up__39->SetBinContent(20,0.4923472);
   Trigger_up__39->SetBinContent(21,1.307606);
   Trigger_up__39->SetBinContent(22,2.752937);
   Trigger_up__39->SetBinContent(23,3.982635);
   Trigger_up__39->SetBinContent(24,11.89789);
   Trigger_up__39->SetBinContent(25,42.02271);
   Trigger_up__39->SetBinContent(26,158.6495);
   Trigger_up__39->SetBinContent(27,715.9306);
   Trigger_up__39->SetBinContent(28,1468.938);
   Trigger_up__39->SetBinContent(29,904.6055);
   Trigger_up__39->SetBinContent(30,265.2089);
   Trigger_up__39->SetBinContent(31,59.8304);
   Trigger_up__39->SetBinContent(32,38.34474);
   Trigger_up__39->SetBinError(8,0.6370975);
   Trigger_up__39->SetBinError(11,0.4400948);
   Trigger_up__39->SetBinError(14,0.4620276);
   Trigger_up__39->SetBinError(17,0.5005881);
   Trigger_up__39->SetBinError(19,0.5617946);
   Trigger_up__39->SetBinError(20,0.4923472);
   Trigger_up__39->SetBinError(21,0.7560643);
   Trigger_up__39->SetBinError(22,1.127994);
   Trigger_up__39->SetBinError(23,1.329532);
   Trigger_up__39->SetBinError(24,2.297996);
   Trigger_up__39->SetBinError(25,4.320062);
   Trigger_up__39->SetBinError(26,8.361396);
   Trigger_up__39->SetBinError(27,17.79123);
   Trigger_up__39->SetBinError(28,25.54237);
   Trigger_up__39->SetBinError(29,19.98495);
   Trigger_up__39->SetBinError(30,10.78173);
   Trigger_up__39->SetBinError(31,5.123754);
   Trigger_up__39->SetBinError(32,2.57872);
   Trigger_up__39->SetBinError(33,3.213556);
   Trigger_up__39->SetEntries(8361);
   Trigger_up__39->SetFillColor(46);
   Trigger_up__39->SetLineColor(46);
   Trigger_up__39->SetMarkerColor(46);
   Trigger_up__39->SetMarkerStyle(21);
   Trigger_up__39->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__39->GetXaxis()->SetRange(1,400);
   Trigger_up__39->GetXaxis()->SetLabelFont(42);
   Trigger_up__39->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__39->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__39->GetXaxis()->SetTitleFont(42);
   Trigger_up__39->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__39->GetYaxis()->SetLabelFont(42);
   Trigger_up__39->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__39->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__39->GetYaxis()->SetTitleOffset(0);
   Trigger_up__39->GetYaxis()->SetTitleFont(42);
   Trigger_up__39->GetZaxis()->SetLabelFont(42);
   Trigger_up__39->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__39->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__39->GetZaxis()->SetTitleFont(42);
   Trigger_up__39->Draw("same");
   TLine *line = new TLine(1730,0,1730,1.469955e+08);
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
   
   TH1D *frameR2__40 = new TH1D("frameR2__40","",1,0,2000);
   frameR2__40->SetMinimum(0.5);
   frameR2__40->SetMaximum(1.5);
   frameR2__40->SetStats(0);
   frameR2__40->SetLineStyle(0);
   frameR2__40->SetMarkerStyle(20);
   frameR2__40->GetXaxis()->SetRange(1,1);
   frameR2__40->GetXaxis()->SetLabelFont(43);
   frameR2__40->GetXaxis()->SetLabelOffset(0.007);
   frameR2__40->GetXaxis()->SetLabelSize(16);
   frameR2__40->GetXaxis()->SetTitleSize(24);
   frameR2__40->GetXaxis()->SetTitleOffset(3.75);
   frameR2__40->GetXaxis()->SetTitleFont(43);
   frameR2__40->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__40->GetYaxis()->SetNdivisions(205);
   frameR2__40->GetYaxis()->SetLabelFont(43);
   frameR2__40->GetYaxis()->SetLabelOffset(0.007);
   frameR2__40->GetYaxis()->SetLabelSize(20);
   frameR2__40->GetYaxis()->SetTitleSize(20);
   frameR2__40->GetYaxis()->SetTitleOffset(2);
   frameR2__40->GetYaxis()->SetTitleFont(43);
   frameR2__40->GetZaxis()->SetLabelFont(42);
   frameR2__40->GetZaxis()->SetLabelOffset(0.007);
   frameR2__40->GetZaxis()->SetLabelSize(0.05);
   frameR2__40->GetZaxis()->SetTitleSize(0.06);
   frameR2__40->GetZaxis()->SetTitleFont(42);
   frameR2__40->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis32[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__41 = new TH1F("Trigger_down__41","",32, xAxis32);
   Trigger_down__41->SetBinContent(0,1.146946);
   Trigger_down__41->SetBinContent(1,1.146946);
   Trigger_down__41->SetBinContent(2,1.146946);
   Trigger_down__41->SetBinContent(3,1.146946);
   Trigger_down__41->SetBinContent(4,1.146946);
   Trigger_down__41->SetBinContent(5,1.146946);
   Trigger_down__41->SetBinContent(6,1.146946);
   Trigger_down__41->SetBinContent(7,1.146946);
   Trigger_down__41->SetBinContent(8,1.146946);
   Trigger_down__41->SetBinContent(9,1.146994);
   Trigger_down__41->SetBinContent(10,1.146994);
   Trigger_down__41->SetBinContent(11,1.146994);
   Trigger_down__41->SetBinContent(12,1.146977);
   Trigger_down__41->SetBinContent(13,1.146977);
   Trigger_down__41->SetBinContent(14,1.146977);
   Trigger_down__41->SetBinContent(15,1.146968);
   Trigger_down__41->SetBinContent(16,1.146968);
   Trigger_down__41->SetBinContent(17,1.146968);
   Trigger_down__41->SetBinContent(18,1.146939);
   Trigger_down__41->SetBinContent(19,1.146939);
   Trigger_down__41->SetBinContent(20,1.146947);
   Trigger_down__41->SetBinContent(21,1.146965);
   Trigger_down__41->SetBinContent(22,1.14702);
   Trigger_down__41->SetBinContent(23,1.147149);
   Trigger_down__41->SetBinContent(24,1.147333);
   Trigger_down__41->SetBinContent(25,1.147928);
   Trigger_down__41->SetBinContent(26,1.149429);
   Trigger_down__41->SetBinContent(27,1.154129);
   Trigger_down__41->SetBinContent(28,1.16339);
   Trigger_down__41->SetBinContent(29,1.160352);
   Trigger_down__41->SetBinContent(30,1.137361);
   Trigger_down__41->SetBinContent(31,1.121237);
   Trigger_down__41->SetBinContent(32,1.12735);
   Trigger_down__41->SetBinError(0,0.01784147);
   Trigger_down__41->SetBinError(1,0.01784147);
   Trigger_down__41->SetBinError(2,0.01784147);
   Trigger_down__41->SetBinError(3,0.01784147);
   Trigger_down__41->SetBinError(4,0.01784147);
   Trigger_down__41->SetBinError(5,0.01784147);
   Trigger_down__41->SetBinError(6,0.01784147);
   Trigger_down__41->SetBinError(7,0.01784147);
   Trigger_down__41->SetBinError(8,0.01784147);
   Trigger_down__41->SetBinError(9,0.01784429);
   Trigger_down__41->SetBinError(10,0.01784429);
   Trigger_down__41->SetBinError(11,0.01784429);
   Trigger_down__41->SetBinError(12,0.01784511);
   Trigger_down__41->SetBinError(13,0.01784511);
   Trigger_down__41->SetBinError(14,0.01784511);
   Trigger_down__41->SetBinError(15,0.01784603);
   Trigger_down__41->SetBinError(16,0.01784603);
   Trigger_down__41->SetBinError(17,0.01784603);
   Trigger_down__41->SetBinError(18,0.01784665);
   Trigger_down__41->SetBinError(19,0.01784665);
   Trigger_down__41->SetBinError(20,0.0178489);
   Trigger_down__41->SetBinError(21,0.01785024);
   Trigger_down__41->SetBinError(22,0.01785429);
   Trigger_down__41->SetBinError(23,0.01786256);
   Trigger_down__41->SetBinError(24,0.017875);
   Trigger_down__41->SetBinError(25,0.0179129);
   Trigger_down__41->SetBinError(26,0.01803915);
   Trigger_down__41->SetBinError(27,0.01852491);
   Trigger_down__41->SetBinError(28,0.02097382);
   Trigger_down__41->SetBinError(29,0.03066998);
   Trigger_down__41->SetBinError(30,0.05607062);
   Trigger_down__41->SetBinError(31,0.1065786);
   Trigger_down__41->SetBinError(32,0.1725129);
   Trigger_down__41->SetEntries(33);
   Trigger_down__41->SetFillColor(38);
   Trigger_down__41->SetLineColor(38);
   Trigger_down__41->SetMarkerColor(38);
   Trigger_down__41->SetMarkerStyle(21);
   Trigger_down__41->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__41->GetXaxis()->SetRange(1,400);
   Trigger_down__41->GetXaxis()->SetLabelFont(42);
   Trigger_down__41->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__41->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__41->GetXaxis()->SetTitleFont(42);
   Trigger_down__41->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__41->GetYaxis()->SetLabelFont(42);
   Trigger_down__41->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__41->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__41->GetYaxis()->SetTitleOffset(0);
   Trigger_down__41->GetYaxis()->SetTitleFont(42);
   Trigger_down__41->GetZaxis()->SetLabelFont(42);
   Trigger_down__41->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__41->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__41->GetZaxis()->SetTitleFont(42);
   Trigger_down__41->Draw("E0 same");
   Double_t xAxis33[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__42 = new TH1F("Trigger_up__42","",32, xAxis33);
   Trigger_up__42->SetBinContent(0,1.004399);
   Trigger_up__42->SetBinContent(1,1.004399);
   Trigger_up__42->SetBinContent(2,1.004399);
   Trigger_up__42->SetBinContent(3,1.004399);
   Trigger_up__42->SetBinContent(4,1.004399);
   Trigger_up__42->SetBinContent(5,1.004399);
   Trigger_up__42->SetBinContent(6,1.004399);
   Trigger_up__42->SetBinContent(7,1.004399);
   Trigger_up__42->SetBinContent(8,1.004399);
   Trigger_up__42->SetBinContent(9,1.004392);
   Trigger_up__42->SetBinContent(10,1.004392);
   Trigger_up__42->SetBinContent(11,1.004392);
   Trigger_up__42->SetBinContent(12,1.00439);
   Trigger_up__42->SetBinContent(13,1.00439);
   Trigger_up__42->SetBinContent(14,1.00439);
   Trigger_up__42->SetBinContent(15,1.004396);
   Trigger_up__42->SetBinContent(16,1.004396);
   Trigger_up__42->SetBinContent(17,1.004396);
   Trigger_up__42->SetBinContent(18,1.004405);
   Trigger_up__42->SetBinContent(19,1.004405);
   Trigger_up__42->SetBinContent(20,1.004406);
   Trigger_up__42->SetBinContent(21,1.004412);
   Trigger_up__42->SetBinContent(22,1.004412);
   Trigger_up__42->SetBinContent(23,1.004404);
   Trigger_up__42->SetBinContent(24,1.004404);
   Trigger_up__42->SetBinContent(25,1.00441);
   Trigger_up__42->SetBinContent(26,1.004358);
   Trigger_up__42->SetBinContent(27,1.004276);
   Trigger_up__42->SetBinContent(28,1.004351);
   Trigger_up__42->SetBinContent(29,1.00859);
   Trigger_up__42->SetBinContent(30,1.012565);
   Trigger_up__42->SetBinContent(31,1.01227);
   Trigger_up__42->SetBinContent(32,1.006426);
   Trigger_up__42->SetBinError(0,0.01557446);
   Trigger_up__42->SetBinError(1,0.01557446);
   Trigger_up__42->SetBinError(2,0.01557446);
   Trigger_up__42->SetBinError(3,0.01557446);
   Trigger_up__42->SetBinError(4,0.01557446);
   Trigger_up__42->SetBinError(5,0.01557446);
   Trigger_up__42->SetBinError(6,0.01557446);
   Trigger_up__42->SetBinError(7,0.01557446);
   Trigger_up__42->SetBinError(8,0.01557446);
   Trigger_up__42->SetBinError(9,0.01557623);
   Trigger_up__42->SetBinError(10,0.01557623);
   Trigger_up__42->SetBinError(11,0.01557623);
   Trigger_up__42->SetBinError(12,0.01557713);
   Trigger_up__42->SetBinError(13,0.01557713);
   Trigger_up__42->SetBinError(14,0.01557713);
   Trigger_up__42->SetBinError(15,0.01557815);
   Trigger_up__42->SetBinError(16,0.01557815);
   Trigger_up__42->SetBinError(17,0.01557815);
   Trigger_up__42->SetBinError(18,0.01557923);
   Trigger_up__42->SetBinError(19,0.01557923);
   Trigger_up__42->SetBinError(20,0.01558108);
   Trigger_up__42->SetBinError(21,0.0155821);
   Trigger_up__42->SetBinError(22,0.01558491);
   Trigger_up__42->SetBinError(23,0.01559035);
   Trigger_up__42->SetBinError(24,0.0155988);
   Trigger_up__42->SetBinError(25,0.01562417);
   Trigger_up__42->SetBinError(26,0.01571358);
   Trigger_up__42->SetBinError(27,0.01607055);
   Trigger_up__42->SetBinError(28,0.01805366);
   Trigger_up__42->SetBinError(29,0.02658308);
   Trigger_up__42->SetBinError(30,0.04975021);
   Trigger_up__42->SetBinError(31,0.09585271);
   Trigger_up__42->SetBinError(32,0.1529139);
   Trigger_up__42->SetEntries(33);
   Trigger_up__42->SetFillColor(46);
   Trigger_up__42->SetLineColor(46);
   Trigger_up__42->SetMarkerColor(46);
   Trigger_up__42->SetMarkerStyle(21);
   Trigger_up__42->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__42->GetXaxis()->SetRange(1,400);
   Trigger_up__42->GetXaxis()->SetLabelFont(42);
   Trigger_up__42->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__42->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__42->GetXaxis()->SetTitleFont(42);
   Trigger_up__42->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__42->GetYaxis()->SetLabelFont(42);
   Trigger_up__42->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__42->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__42->GetYaxis()->SetTitleOffset(0);
   Trigger_up__42->GetYaxis()->SetTitleFont(42);
   Trigger_up__42->GetZaxis()->SetLabelFont(42);
   Trigger_up__42->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__42->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__42->GetZaxis()->SetTitleFont(42);
   Trigger_up__42->Draw("E0 same");
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
   
   TH1D *frameR2__43 = new TH1D("frameR2__43","",1,0,2000);
   frameR2__43->SetMinimum(0.5);
   frameR2__43->SetMaximum(1.5);
   frameR2__43->SetStats(0);
   frameR2__43->SetLineStyle(0);
   frameR2__43->SetMarkerStyle(20);
   frameR2__43->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__43->GetXaxis()->SetRange(1,1);
   frameR2__43->GetXaxis()->SetLabelFont(43);
   frameR2__43->GetXaxis()->SetLabelOffset(0.007);
   frameR2__43->GetXaxis()->SetLabelSize(16);
   frameR2__43->GetXaxis()->SetTitleSize(24);
   frameR2__43->GetXaxis()->SetTitleOffset(5);
   frameR2__43->GetXaxis()->SetTitleFont(43);
   frameR2__43->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__43->GetYaxis()->SetNdivisions(205);
   frameR2__43->GetYaxis()->SetLabelFont(43);
   frameR2__43->GetYaxis()->SetLabelOffset(0.007);
   frameR2__43->GetYaxis()->SetLabelSize(20);
   frameR2__43->GetYaxis()->SetTitleSize(20);
   frameR2__43->GetYaxis()->SetTitleOffset(2);
   frameR2__43->GetYaxis()->SetTitleFont(43);
   frameR2__43->GetZaxis()->SetLabelFont(42);
   frameR2__43->GetZaxis()->SetLabelOffset(0.007);
   frameR2__43->GetZaxis()->SetLabelSize(0.05);
   frameR2__43->GetZaxis()->SetTitleSize(0.06);
   frameR2__43->GetZaxis()->SetTitleFont(42);
   frameR2__43->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis34[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__44 = new TH1F("nominal__44","",32, xAxis34);
   nominal__44->SetBinContent(8,0.9848708);
   nominal__44->SetBinContent(11,1.298701);
   nominal__44->SetBinContent(14,1.234568);
   nominal__44->SetBinContent(17,1.428571);
   nominal__44->SetBinContent(19,1.111111);
   nominal__44->SetBinContent(20,1.020408);
   nominal__44->SetBinContent(21,1.010052);
   nominal__44->SetBinContent(22,0.9994043);
   nominal__44->SetBinContent(23,0.9995875);
   nominal__44->SetBinContent(24,0.9894891);
   nominal__44->SetBinContent(25,1.032552);
   nominal__44->SetBinContent(26,1.056009);
   nominal__44->SetBinContent(27,1.120032);
   nominal__44->SetBinContent(28,1.166046);
   nominal__44->SetBinContent(29,1.169905);
   nominal__44->SetBinContent(30,1.143445);
   nominal__44->SetBinContent(31,1.11739);
   nominal__44->SetBinContent(32,1.12735);
   nominal__44->SetBinError(8,0.9891314);
   nominal__44->SetBinError(11,1.836641);
   nominal__44->SetBinError(14,1.745943);
   nominal__44->SetBinError(17,2.020305);
   nominal__44->SetBinError(19,1.117173);
   nominal__44->SetBinError(20,1.443075);
   nominal__44->SetBinError(21,0.8287652);
   nominal__44->SetBinError(22,0.5801916);
   nominal__44->SetBinError(23,0.4729524);
   nominal__44->SetBinError(24,0.2710305);
   nominal__44->SetBinError(25,0.1505571);
   nominal__44->SetBinError(26,0.07887584);
   nominal__44->SetBinError(27,0.03947185);
   nominal__44->SetBinError(28,0.02874608);
   nominal__44->SetBinError(29,0.03663013);
   nominal__44->SetBinError(30,0.06592774);
   nominal__44->SetBinError(31,0.1355326);
   nominal__44->SetBinError(32,0.1088285);
   nominal__44->SetMinimum(5e-05);
   nominal__44->SetMaximum(1.469955e+08);
   nominal__44->SetEntries(24.71494);
   nominal__44->SetFillColor(38);
   nominal__44->SetLineColor(38);
   nominal__44->SetMarkerColor(38);
   nominal__44->SetMarkerStyle(21);
   nominal__44->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__44->GetXaxis()->SetRange(1,32);
   nominal__44->GetXaxis()->SetLabelFont(42);
   nominal__44->GetXaxis()->SetLabelSize(0.035);
   nominal__44->GetXaxis()->SetTitleSize(0.035);
   nominal__44->GetXaxis()->SetTitleFont(42);
   nominal__44->GetYaxis()->SetTitle("Tracks");
   nominal__44->GetYaxis()->SetLabelFont(42);
   nominal__44->GetYaxis()->SetLabelSize(0.05);
   nominal__44->GetYaxis()->SetTitleSize(0.07);
   nominal__44->GetYaxis()->SetTitleOffset(0);
   nominal__44->GetYaxis()->SetTitleFont(42);
   nominal__44->GetZaxis()->SetLabelFont(42);
   nominal__44->GetZaxis()->SetLabelSize(0.035);
   nominal__44->GetZaxis()->SetTitleSize(0.035);
   nominal__44->GetZaxis()->SetTitleFont(42);
   nominal__44->Draw("E0 same");
   Double_t xAxis35[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__45 = new TH1F("nominal__45","",32, xAxis35);
   nominal__45->SetBinContent(8,1.031544);
   nominal__45->SetBinContent(11,1.020408);
   nominal__45->SetBinContent(14,0.9615384);
   nominal__45->SetBinContent(17,0.9345794);
   nominal__45->SetBinContent(19,1);
   nominal__45->SetBinContent(20,0.9615386);
   nominal__45->SetBinContent(21,1.003475);
   nominal__45->SetBinContent(22,1.01543);
   nominal__45->SetBinContent(23,1.004122);
   nominal__45->SetBinContent(24,1.002666);
   nominal__45->SetBinContent(25,1.008893);
   nominal__45->SetBinContent(26,1.006125);
   nominal__45->SetBinContent(27,1.003992);
   nominal__45->SetBinContent(28,1.000692);
   nominal__45->SetBinContent(29,1.006993);
   nominal__45->SetBinContent(30,1.012674);
   nominal__45->SetBinContent(31,1.016015);
   nominal__45->SetBinContent(32,1.006426);
   nominal__45->SetBinError(8,1.031779);
   nominal__45->SetBinError(11,1.443075);
   nominal__45->SetBinError(14,1.359821);
   nominal__45->SetBinError(17,1.321695);
   nominal__45->SetBinError(19,1.005455);
   nominal__45->SetBinError(20,1.359821);
   nominal__45->SetBinError(21,0.8205638);
   nominal__45->SetBinError(22,0.5884279);
   nominal__45->SetBinError(23,0.4738444);
   nominal__45->SetBinError(24,0.2738422);
   nominal__45->SetBinError(25,0.146659);
   nominal__45->SetBinError(26,0.07496872);
   nominal__45->SetBinError(27,0.03526887);
   nominal__45->SetBinError(28,0.02459545);
   nominal__45->SetBinError(29,0.03144881);
   nominal__45->SetBinError(30,0.05820375);
   nominal__45->SetBinError(31,0.1230193);
   nominal__45->SetBinError(32,0.0957167);
   nominal__45->SetMinimum(5e-05);
   nominal__45->SetMaximum(1.469955e+08);
   nominal__45->SetEntries(29.49922);
   nominal__45->SetFillColor(46);
   nominal__45->SetLineColor(46);
   nominal__45->SetMarkerColor(46);
   nominal__45->SetMarkerStyle(21);
   nominal__45->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__45->GetXaxis()->SetRange(1,32);
   nominal__45->GetXaxis()->SetLabelFont(42);
   nominal__45->GetXaxis()->SetLabelSize(0.035);
   nominal__45->GetXaxis()->SetTitleSize(0.035);
   nominal__45->GetXaxis()->SetTitleFont(42);
   nominal__45->GetYaxis()->SetTitle("Tracks");
   nominal__45->GetYaxis()->SetLabelFont(42);
   nominal__45->GetYaxis()->SetLabelSize(0.05);
   nominal__45->GetYaxis()->SetTitleSize(0.07);
   nominal__45->GetYaxis()->SetTitleOffset(0);
   nominal__45->GetYaxis()->SetTitleFont(42);
   nominal__45->GetZaxis()->SetLabelFont(42);
   nominal__45->GetZaxis()->SetLabelSize(0.035);
   nominal__45->GetZaxis()->SetTitleSize(0.035);
   nominal__45->GetZaxis()->SetTitleFont(42);
   nominal__45->Draw("E0 same");
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
