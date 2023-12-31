void plot_PU()
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
   t1->Range(-428.5714,-4.343799,2428.571,4.210073);
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
   Double_t xAxis1[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1 = new TH1F("nominal__1","",32, xAxis1);
   nominal__1->SetBinContent(8,0.9293805);
   nominal__1->SetBinContent(11,0.4490763);
   nominal__1->SetBinContent(14,0.4442573);
   nominal__1->SetBinContent(17,0.4678394);
   nominal__1->SetBinContent(19,0.7901867);
   nominal__1->SetBinContent(20,0.4734108);
   nominal__1->SetBinContent(21,1.312151);
   nominal__1->SetBinContent(22,2.795414);
   nominal__1->SetBinContent(23,3.999052);
   nominal__1->SetBinContent(24,11.92961);
   nominal__1->SetBinContent(25,42.3964);
   nominal__1->SetBinContent(26,159.6213);
   nominal__1->SetBinContent(27,718.7883);
   nominal__1->SetBinContent(28,1469.955);
   nominal__1->SetBinContent(29,910.9312);
   nominal__1->SetBinContent(30,268.57);
   nominal__1->SetBinContent(31,60.78858);
   nominal__1->SetBinContent(32,38.59114);
   nominal__1->SetBinError(8,0.6574474);
   nominal__1->SetBinError(11,0.4490763);
   nominal__1->SetBinError(14,0.4442573);
   nominal__1->SetBinError(17,0.4678393);
   nominal__1->SetBinError(19,0.5617946);
   nominal__1->SetBinError(20,0.4734108);
   nominal__1->SetBinError(21,0.7587229);
   nominal__1->SetBinError(22,1.145493);
   nominal__1->SetBinError(23,1.333819);
   nominal__1->SetBinError(24,2.303587);
   nominal__1->SetBinError(25,4.35733);
   nominal__1->SetBinError(26,8.40769);
   nominal__1->SetBinError(27,17.84672);
   nominal__1->SetBinError(28,25.53434);
   nominal__1->SetBinError(29,20.10793);
   nominal__1->SetBinError(30,10.91164);
   nominal__1->SetBinError(31,5.203216);
   nominal__1->SetBinError(32,2.595201);
   nominal__1->SetBinError(33,3.232536);
   nominal__1->SetMinimum(5e-05);
   nominal__1->SetMaximum(14699.55);
   nominal__1->SetEntries(8361);
   nominal__1->SetFillColor(1);
   nominal__1->SetMarkerStyle(20);
   nominal__1->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1->GetXaxis()->SetRange(1,32);
   nominal__1->GetXaxis()->SetLabelFont(42);
   nominal__1->GetXaxis()->SetLabelSize(0.035);
   nominal__1->GetXaxis()->SetTitleSize(0.035);
   nominal__1->GetXaxis()->SetTitleFont(42);
   nominal__1->GetYaxis()->SetTitle("Tracks");
   nominal__1->GetYaxis()->SetLabelFont(42);
   nominal__1->GetYaxis()->SetLabelSize(0.05);
   nominal__1->GetYaxis()->SetTitleSize(0.07);
   nominal__1->GetYaxis()->SetTitleOffset(0);
   nominal__1->GetYaxis()->SetTitleFont(42);
   nominal__1->GetZaxis()->SetLabelFont(42);
   nominal__1->GetZaxis()->SetLabelSize(0.035);
   nominal__1->GetZaxis()->SetTitleSize(0.035);
   nominal__1->GetZaxis()->SetTitleFont(42);
   nominal__1->Draw("");
   Double_t xAxis2[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__2 = new TH1F("PU_down__2","",32, xAxis2);
   PU_down__2->SetBinContent(8,1.066027);
   PU_down__2->SetBinContent(11,0.4458178);
   PU_down__2->SetBinContent(14,0.4180184);
   PU_down__2->SetBinContent(17,0.5434421);
   PU_down__2->SetBinContent(19,0.6846303);
   PU_down__2->SetBinContent(20,0.5743546);
   PU_down__2->SetBinContent(21,1.283756);
   PU_down__2->SetBinContent(22,3.554881);
   PU_down__2->SetBinContent(23,3.921008);
   PU_down__2->SetBinContent(24,12.47859);
   PU_down__2->SetBinContent(25,43.58586);
   PU_down__2->SetBinContent(26,161.099);
   PU_down__2->SetBinContent(27,712.6426);
   PU_down__2->SetBinContent(28,1468.743);
   PU_down__2->SetBinContent(29,911.7599);
   PU_down__2->SetBinContent(30,267.3781);
   PU_down__2->SetBinContent(31,61.58892);
   PU_down__2->SetBinContent(32,39.40254);
   PU_down__2->SetBinError(8,0.7611817);
   PU_down__2->SetBinError(11,0.4458178);
   PU_down__2->SetBinError(14,0.4180184);
   PU_down__2->SetBinError(17,0.5434421);
   PU_down__2->SetBinError(19,0.4909987);
   PU_down__2->SetBinError(20,0.5743547);
   PU_down__2->SetBinError(21,0.7497858);
   PU_down__2->SetBinError(22,1.585806);
   PU_down__2->SetBinError(23,1.327214);
   PU_down__2->SetBinError(24,2.615791);
   PU_down__2->SetBinError(25,4.634877);
   PU_down__2->SetBinError(26,8.813367);
   PU_down__2->SetBinError(27,18.18352);
   PU_down__2->SetBinError(28,26.31102);
   PU_down__2->SetBinError(29,20.78295);
   PU_down__2->SetBinError(30,11.1515);
   PU_down__2->SetBinError(31,5.455497);
   PU_down__2->SetBinError(32,2.814612);
   PU_down__2->SetBinError(33,3.363026);
   PU_down__2->SetEntries(8361);
   PU_down__2->SetFillColor(38);
   PU_down__2->SetLineColor(38);
   PU_down__2->SetMarkerColor(38);
   PU_down__2->SetMarkerStyle(21);
   PU_down__2->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__2->GetXaxis()->SetRange(1,400);
   PU_down__2->GetXaxis()->SetLabelFont(42);
   PU_down__2->GetXaxis()->SetLabelSize(0.035);
   PU_down__2->GetXaxis()->SetTitleSize(0.035);
   PU_down__2->GetXaxis()->SetTitleFont(42);
   PU_down__2->GetYaxis()->SetTitle("Events / bin");
   PU_down__2->GetYaxis()->SetLabelFont(42);
   PU_down__2->GetYaxis()->SetLabelSize(0.035);
   PU_down__2->GetYaxis()->SetTitleSize(0.035);
   PU_down__2->GetYaxis()->SetTitleOffset(0);
   PU_down__2->GetYaxis()->SetTitleFont(42);
   PU_down__2->GetZaxis()->SetLabelFont(42);
   PU_down__2->GetZaxis()->SetLabelSize(0.035);
   PU_down__2->GetZaxis()->SetTitleSize(0.035);
   PU_down__2->GetZaxis()->SetTitleFont(42);
   PU_down__2->Draw("same");
   Double_t xAxis3[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__3 = new TH1F("PU_up__3","",32, xAxis3);
   PU_up__3->SetBinContent(8,0.7529135);
   PU_up__3->SetBinContent(11,0.4324154);
   PU_up__3->SetBinContent(14,0.4659604);
   PU_up__3->SetBinContent(17,0.3636434);
   PU_up__3->SetBinContent(19,0.9219347);
   PU_up__3->SetBinContent(20,0.3485974);
   PU_up__3->SetBinContent(21,1.307596);
   PU_up__3->SetBinContent(22,2.317645);
   PU_up__3->SetBinContent(23,4.045685);
   PU_up__3->SetBinContent(24,11.76653);
   PU_up__3->SetBinContent(25,41.26015);
   PU_up__3->SetBinContent(26,158.8047);
   PU_up__3->SetBinContent(27,723.3053);
   PU_up__3->SetBinContent(28,1470.573);
   PU_up__3->SetBinContent(29,910.3069);
   PU_up__3->SetBinContent(30,268.5987);
   PU_up__3->SetBinContent(31,60.00992);
   PU_up__3->SetBinContent(32,37.74822);
   PU_up__3->SetBinError(8,0.5359561);
   PU_up__3->SetBinError(11,0.4324154);
   PU_up__3->SetBinError(14,0.4659603);
   PU_up__3->SetBinError(17,0.3636434);
   PU_up__3->SetBinError(19,0.6521917);
   PU_up__3->SetBinError(20,0.3485974);
   PU_up__3->SetBinError(21,0.7582286);
   PU_up__3->SetBinError(22,0.9768017);
   PU_up__3->SetBinError(23,1.356289);
   PU_up__3->SetBinError(24,2.291312);
   PU_up__3->SetBinError(25,4.277226);
   PU_up__3->SetBinError(26,8.423988);
   PU_up__3->SetBinError(27,18.06624);
   PU_up__3->SetBinError(28,25.7099);
   PU_up__3->SetBinError(29,20.22421);
   PU_up__3->SetBinError(30,10.9766);
   PU_up__3->SetBinError(31,5.175665);
   PU_up__3->SetBinError(32,2.524209);
   PU_up__3->SetBinError(33,3.214886);
   PU_up__3->SetEntries(8361);
   PU_up__3->SetFillColor(46);
   PU_up__3->SetLineColor(46);
   PU_up__3->SetMarkerColor(46);
   PU_up__3->SetMarkerStyle(21);
   PU_up__3->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__3->GetXaxis()->SetRange(1,400);
   PU_up__3->GetXaxis()->SetLabelFont(42);
   PU_up__3->GetXaxis()->SetLabelSize(0.035);
   PU_up__3->GetXaxis()->SetTitleSize(0.035);
   PU_up__3->GetXaxis()->SetTitleFont(42);
   PU_up__3->GetYaxis()->SetTitle("Events / bin");
   PU_up__3->GetYaxis()->SetLabelFont(42);
   PU_up__3->GetYaxis()->SetLabelSize(0.035);
   PU_up__3->GetYaxis()->SetTitleSize(0.035);
   PU_up__3->GetYaxis()->SetTitleOffset(0);
   PU_up__3->GetYaxis()->SetTitleFont(42);
   PU_up__3->GetZaxis()->SetLabelFont(42);
   PU_up__3->GetZaxis()->SetLabelSize(0.035);
   PU_up__3->GetZaxis()->SetTitleSize(0.035);
   PU_up__3->GetZaxis()->SetTitleFont(42);
   PU_up__3->Draw("same");
   TLine *line = new TLine(1730,0,1730,14699.55);
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
   
   TH1D *frameR2__4 = new TH1D("frameR2__4","",1,0,2000);
   frameR2__4->SetMinimum(0.5);
   frameR2__4->SetMaximum(1.5);
   frameR2__4->SetStats(0);
   frameR2__4->SetLineStyle(0);
   frameR2__4->SetMarkerStyle(20);
   frameR2__4->GetXaxis()->SetRange(1,1);
   frameR2__4->GetXaxis()->SetLabelFont(43);
   frameR2__4->GetXaxis()->SetLabelOffset(0.007);
   frameR2__4->GetXaxis()->SetLabelSize(16);
   frameR2__4->GetXaxis()->SetTitleSize(24);
   frameR2__4->GetXaxis()->SetTitleOffset(3.75);
   frameR2__4->GetXaxis()->SetTitleFont(43);
   frameR2__4->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__4->GetYaxis()->SetNdivisions(205);
   frameR2__4->GetYaxis()->SetLabelFont(43);
   frameR2__4->GetYaxis()->SetLabelOffset(0.007);
   frameR2__4->GetYaxis()->SetLabelSize(20);
   frameR2__4->GetYaxis()->SetTitleSize(20);
   frameR2__4->GetYaxis()->SetTitleOffset(2);
   frameR2__4->GetYaxis()->SetTitleFont(43);
   frameR2__4->GetZaxis()->SetLabelFont(42);
   frameR2__4->GetZaxis()->SetLabelOffset(0.007);
   frameR2__4->GetZaxis()->SetLabelSize(0.05);
   frameR2__4->GetZaxis()->SetTitleSize(0.06);
   frameR2__4->GetZaxis()->SetTitleFont(42);
   frameR2__4->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis4[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__5 = new TH1F("PU_down__5","",32, xAxis4);
   PU_down__5->SetBinContent(0,1.000559);
   PU_down__5->SetBinContent(1,1.000559);
   PU_down__5->SetBinContent(2,1.000559);
   PU_down__5->SetBinContent(3,1.000559);
   PU_down__5->SetBinContent(4,1.000559);
   PU_down__5->SetBinContent(5,1.000559);
   PU_down__5->SetBinContent(6,1.000559);
   PU_down__5->SetBinContent(7,1.000559);
   PU_down__5->SetBinContent(8,1.000559);
   PU_down__5->SetBinContent(9,1.000596);
   PU_down__5->SetBinContent(10,1.000596);
   PU_down__5->SetBinContent(11,1.000596);
   PU_down__5->SetBinContent(12,1.000595);
   PU_down__5->SetBinContent(13,1.000595);
   PU_down__5->SetBinContent(14,1.000595);
   PU_down__5->SetBinContent(15,1.000588);
   PU_down__5->SetBinContent(16,1.000588);
   PU_down__5->SetBinContent(17,1.000588);
   PU_down__5->SetBinContent(18,1.000609);
   PU_down__5->SetBinContent(19,1.000609);
   PU_down__5->SetBinContent(20,1.00058);
   PU_down__5->SetBinContent(21,1.000607);
   PU_down__5->SetBinContent(22,1.0006);
   PU_down__5->SetBinContent(23,1.000807);
   PU_down__5->SetBinContent(24,1.000786);
   PU_down__5->SetBinContent(25,1.000939);
   PU_down__5->SetBinContent(26,1.001279);
   PU_down__5->SetBinContent(27,1.001765);
   PU_down__5->SetBinContent(28,0.9999868);
   PU_down__5->SetBinContent(29,0.9990247);
   PU_down__5->SetBinContent(30,0.9988604);
   PU_down__5->SetBinContent(31,0.9840409);
   PU_down__5->SetBinContent(32,0.9794075);
   PU_down__5->SetBinError(0,0.01575327);
   PU_down__5->SetBinError(1,0.01575327);
   PU_down__5->SetBinError(2,0.01575327);
   PU_down__5->SetBinError(3,0.01575327);
   PU_down__5->SetBinError(4,0.01575327);
   PU_down__5->SetBinError(5,0.01575327);
   PU_down__5->SetBinError(6,0.01575327);
   PU_down__5->SetBinError(7,0.01575327);
   PU_down__5->SetBinError(8,0.01575327);
   PU_down__5->SetBinError(9,0.01575577);
   PU_down__5->SetBinError(10,0.01575577);
   PU_down__5->SetBinError(11,0.01575577);
   PU_down__5->SetBinError(12,0.01575673);
   PU_down__5->SetBinError(13,0.01575673);
   PU_down__5->SetBinError(14,0.01575673);
   PU_down__5->SetBinError(15,0.01575759);
   PU_down__5->SetBinError(16,0.01575759);
   PU_down__5->SetBinError(17,0.01575759);
   PU_down__5->SetBinError(18,0.01575888);
   PU_down__5->SetBinError(19,0.01575888);
   PU_down__5->SetBinError(20,0.01576027);
   PU_down__5->SetBinError(21,0.01576166);
   PU_down__5->SetBinError(22,0.01576443);
   PU_down__5->SetBinError(23,0.01577238);
   PU_down__5->SetBinError(24,0.01578073);
   PU_down__5->SetBinError(25,0.01580709);
   PU_down__5->SetBinError(26,0.01590328);
   PU_down__5->SetBinError(27,0.01627074);
   PU_down__5->SetBinError(28,0.01825163);
   PU_down__5->SetBinError(29,0.02674016);
   PU_down__5->SetBinError(30,0.04978116);
   PU_down__5->SetBinError(31,0.09482049);
   PU_down__5->SetBinError(32,0.1514951);
   PU_down__5->SetEntries(33);
   PU_down__5->SetFillColor(38);
   PU_down__5->SetLineColor(38);
   PU_down__5->SetMarkerColor(38);
   PU_down__5->SetMarkerStyle(21);
   PU_down__5->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__5->GetXaxis()->SetRange(1,400);
   PU_down__5->GetXaxis()->SetLabelFont(42);
   PU_down__5->GetXaxis()->SetLabelSize(0.035);
   PU_down__5->GetXaxis()->SetTitleSize(0.035);
   PU_down__5->GetXaxis()->SetTitleFont(42);
   PU_down__5->GetYaxis()->SetTitle("Events / bin");
   PU_down__5->GetYaxis()->SetLabelFont(42);
   PU_down__5->GetYaxis()->SetLabelSize(0.035);
   PU_down__5->GetYaxis()->SetTitleSize(0.035);
   PU_down__5->GetYaxis()->SetTitleOffset(0);
   PU_down__5->GetYaxis()->SetTitleFont(42);
   PU_down__5->GetZaxis()->SetLabelFont(42);
   PU_down__5->GetZaxis()->SetLabelSize(0.035);
   PU_down__5->GetZaxis()->SetTitleSize(0.035);
   PU_down__5->GetZaxis()->SetTitleFont(42);
   PU_down__5->Draw("E0 same");
   Double_t xAxis5[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__6 = new TH1F("PU_up__6","",32, xAxis5);
   PU_up__6->SetBinContent(0,0.9999737);
   PU_up__6->SetBinContent(1,0.9999737);
   PU_up__6->SetBinContent(2,0.9999737);
   PU_up__6->SetBinContent(3,0.9999737);
   PU_up__6->SetBinContent(4,0.9999737);
   PU_up__6->SetBinContent(5,0.9999737);
   PU_up__6->SetBinContent(6,0.9999737);
   PU_up__6->SetBinContent(7,0.9999737);
   PU_up__6->SetBinContent(8,0.9999737);
   PU_up__6->SetBinContent(9,0.9999259);
   PU_up__6->SetBinContent(10,0.9999259);
   PU_up__6->SetBinContent(11,0.9999259);
   PU_up__6->SetBinContent(12,0.9999214);
   PU_up__6->SetBinContent(13,0.9999214);
   PU_up__6->SetBinContent(14,0.9999214);
   PU_up__6->SetBinContent(15,0.9999273);
   PU_up__6->SetBinContent(16,0.9999273);
   PU_up__6->SetBinContent(17,0.9999273);
   PU_up__6->SetBinContent(18,0.999899);
   PU_up__6->SetBinContent(19,0.999899);
   PU_up__6->SetBinContent(20,0.9999347);
   PU_up__6->SetBinContent(21,0.9999009);
   PU_up__6->SetBinContent(22,0.9998996);
   PU_up__6->SetBinContent(23,0.9997699);
   PU_up__6->SetBinContent(24,0.9997823);
   PU_up__6->SetBinContent(25,0.9997372);
   PU_up__6->SetBinContent(26,0.9994212);
   PU_up__6->SetBinContent(27,0.9991594);
   PU_up__6->SetBinContent(28,1.000582);
   PU_up__6->SetBinContent(29,1.001737);
   PU_up__6->SetBinContent(30,1.004348);
   PU_up__6->SetBinContent(31,1.016588);
   PU_up__6->SetBinContent(32,1.02233);
   PU_up__6->SetBinError(0,0.01554904);
   PU_up__6->SetBinError(1,0.01554904);
   PU_up__6->SetBinError(2,0.01554904);
   PU_up__6->SetBinError(3,0.01554904);
   PU_up__6->SetBinError(4,0.01554904);
   PU_up__6->SetBinError(5,0.01554904);
   PU_up__6->SetBinError(6,0.01554904);
   PU_up__6->SetBinError(7,0.01554904);
   PU_up__6->SetBinError(8,0.01554904);
   PU_up__6->SetBinError(9,0.01555014);
   PU_up__6->SetBinError(10,0.01555014);
   PU_up__6->SetBinError(11,0.01555014);
   PU_up__6->SetBinError(12,0.01555101);
   PU_up__6->SetBinError(13,0.01555101);
   PU_up__6->SetBinError(14,0.01555101);
   PU_up__6->SetBinError(15,0.01555204);
   PU_up__6->SetBinError(16,0.01555204);
   PU_up__6->SetBinError(17,0.01555204);
   PU_up__6->SetBinError(18,0.01555252);
   PU_up__6->SetBinError(19,0.01555252);
   PU_up__6->SetBinError(20,0.01555494);
   PU_up__6->SetBinError(21,0.01555533);
   PU_up__6->SetBinError(22,0.01555812);
   PU_up__6->SetBinError(23,0.01556153);
   PU_up__6->SetBinError(24,0.01557015);
   PU_up__6->SetBinError(25,0.01559462);
   PU_up__6->SetBinError(26,0.01567945);
   PU_up__6->SetBinError(27,0.01603237);
   PU_up__6->SetBinError(28,0.01803589);
   PU_up__6->SetBinError(29,0.02647759);
   PU_up__6->SetBinError(30,0.04948985);
   PU_up__6->SetBinError(31,0.09661564);
   PU_up__6->SetBinError(32,0.15593);
   PU_up__6->SetEntries(33);
   PU_up__6->SetFillColor(46);
   PU_up__6->SetLineColor(46);
   PU_up__6->SetMarkerColor(46);
   PU_up__6->SetMarkerStyle(21);
   PU_up__6->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__6->GetXaxis()->SetRange(1,400);
   PU_up__6->GetXaxis()->SetLabelFont(42);
   PU_up__6->GetXaxis()->SetLabelSize(0.035);
   PU_up__6->GetXaxis()->SetTitleSize(0.035);
   PU_up__6->GetXaxis()->SetTitleFont(42);
   PU_up__6->GetYaxis()->SetTitle("Events / bin");
   PU_up__6->GetYaxis()->SetLabelFont(42);
   PU_up__6->GetYaxis()->SetLabelSize(0.035);
   PU_up__6->GetYaxis()->SetTitleSize(0.035);
   PU_up__6->GetYaxis()->SetTitleOffset(0);
   PU_up__6->GetYaxis()->SetTitleFont(42);
   PU_up__6->GetZaxis()->SetLabelFont(42);
   PU_up__6->GetZaxis()->SetLabelSize(0.035);
   PU_up__6->GetZaxis()->SetTitleSize(0.035);
   PU_up__6->GetZaxis()->SetTitleFont(42);
   PU_up__6->Draw("E0 same");
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
   
   TH1D *frameR2__7 = new TH1D("frameR2__7","",1,0,2000);
   frameR2__7->SetMinimum(0.5);
   frameR2__7->SetMaximum(1.5);
   frameR2__7->SetStats(0);
   frameR2__7->SetLineStyle(0);
   frameR2__7->SetMarkerStyle(20);
   frameR2__7->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__7->GetXaxis()->SetRange(1,1);
   frameR2__7->GetXaxis()->SetLabelFont(43);
   frameR2__7->GetXaxis()->SetLabelOffset(0.007);
   frameR2__7->GetXaxis()->SetLabelSize(16);
   frameR2__7->GetXaxis()->SetTitleSize(24);
   frameR2__7->GetXaxis()->SetTitleOffset(5);
   frameR2__7->GetXaxis()->SetTitleFont(43);
   frameR2__7->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__7->GetYaxis()->SetNdivisions(205);
   frameR2__7->GetYaxis()->SetLabelFont(43);
   frameR2__7->GetYaxis()->SetLabelOffset(0.007);
   frameR2__7->GetYaxis()->SetLabelSize(20);
   frameR2__7->GetYaxis()->SetTitleSize(20);
   frameR2__7->GetYaxis()->SetTitleOffset(2);
   frameR2__7->GetYaxis()->SetTitleFont(43);
   frameR2__7->GetZaxis()->SetLabelFont(42);
   frameR2__7->GetZaxis()->SetLabelOffset(0.007);
   frameR2__7->GetZaxis()->SetLabelSize(0.05);
   frameR2__7->GetZaxis()->SetTitleSize(0.06);
   frameR2__7->GetZaxis()->SetTitleFont(42);
   frameR2__7->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis6[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__8 = new TH1F("nominal__8","",32, xAxis6);
   nominal__8->SetBinContent(8,0.8718169);
   nominal__8->SetBinContent(11,1.007309);
   nominal__8->SetBinContent(14,1.06277);
   nominal__8->SetBinContent(17,0.8608817);
   nominal__8->SetBinContent(19,1.15418);
   nominal__8->SetBinContent(20,0.8242482);
   nominal__8->SetBinContent(21,1.022118);
   nominal__8->SetBinContent(22,0.7863595);
   nominal__8->SetBinContent(23,1.019904);
   nominal__8->SetBinContent(24,0.9560061);
   nominal__8->SetBinContent(25,0.97271);
   nominal__8->SetBinContent(26,0.9908271);
   nominal__8->SetBinContent(27,1.008624);
   nominal__8->SetBinContent(28,1.000825);
   nominal__8->SetBinContent(29,0.9990911);
   nominal__8->SetBinContent(30,1.004458);
   nominal__8->SetBinContent(31,0.9870051);
   nominal__8->SetBinContent(32,0.9794075);
   nominal__8->SetBinError(8,0.8762812);
   nominal__8->SetBinError(11,1.42455);
   nominal__8->SetBinError(14,1.502983);
   nominal__8->SetBinError(17,1.217471);
   nominal__8->SetBinError(19,1.165555);
   nominal__8->SetBinError(20,1.165663);
   nominal__8->SetBinError(21,0.8400479);
   nominal__8->SetBinError(22,0.4763254);
   nominal__8->SetBinError(23,0.4846626);
   nominal__8->SetBinError(24,0.2724676);
   nominal__8->SetBinError(25,0.1438522);
   nominal__8->SetBinError(26,0.07524649);
   nominal__8->SetBinError(27,0.03590929);
   nominal__8->SetBinError(28,0.02497367);
   nominal__8->SetBinError(29,0.03170197);
   nominal__8->SetBinError(30,0.05848456);
   nominal__8->SetBinError(31,0.1215774);
   nominal__8->SetBinError(32,0.09608654);
   nominal__8->SetMinimum(5e-05);
   nominal__8->SetMaximum(14699.55);
   nominal__8->SetEntries(29.04548);
   nominal__8->SetFillColor(38);
   nominal__8->SetLineColor(38);
   nominal__8->SetMarkerColor(38);
   nominal__8->SetMarkerStyle(21);
   nominal__8->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__8->GetXaxis()->SetRange(1,32);
   nominal__8->GetXaxis()->SetLabelFont(42);
   nominal__8->GetXaxis()->SetLabelSize(0.035);
   nominal__8->GetXaxis()->SetTitleSize(0.035);
   nominal__8->GetXaxis()->SetTitleFont(42);
   nominal__8->GetYaxis()->SetTitle("Tracks");
   nominal__8->GetYaxis()->SetLabelFont(42);
   nominal__8->GetYaxis()->SetLabelSize(0.05);
   nominal__8->GetYaxis()->SetTitleSize(0.07);
   nominal__8->GetYaxis()->SetTitleOffset(0);
   nominal__8->GetYaxis()->SetTitleFont(42);
   nominal__8->GetZaxis()->SetLabelFont(42);
   nominal__8->GetZaxis()->SetLabelSize(0.035);
   nominal__8->GetZaxis()->SetTitleSize(0.035);
   nominal__8->GetZaxis()->SetTitleFont(42);
   nominal__8->Draw("E0 same");
   Double_t xAxis7[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__9 = new TH1F("nominal__9","",32, xAxis7);
   nominal__9->SetBinContent(8,1.234379);
   nominal__9->SetBinContent(11,1.03853);
   nominal__9->SetBinContent(14,0.9534228);
   nominal__9->SetBinContent(17,1.286533);
   nominal__9->SetBinContent(19,0.8570962);
   nominal__9->SetBinContent(20,1.358045);
   nominal__9->SetBinContent(21,1.003483);
   nominal__9->SetBinContent(22,1.206144);
   nominal__9->SetBinContent(23,0.9884734);
   nominal__9->SetBinContent(24,1.01386);
   nominal__9->SetBinContent(25,1.027539);
   nominal__9->SetBinContent(26,1.005142);
   nominal__9->SetBinContent(27,0.9937551);
   nominal__9->SetBinContent(28,0.9995801);
   nominal__9->SetBinContent(29,1.000686);
   nominal__9->SetBinContent(30,0.9998932);
   nominal__9->SetBinContent(31,1.012975);
   nominal__9->SetBinContent(32,1.02233);
   nominal__9->SetBinError(8,1.238778);
   nominal__9->SetBinError(11,1.468703);
   nominal__9->SetBinError(14,1.348344);
   nominal__9->SetBinError(17,1.819433);
   nominal__9->SetBinError(19,0.8596245);
   nominal__9->SetBinError(20,1.920565);
   nominal__9->SetBinError(21,0.8217482);
   nominal__9->SetBinError(22,0.7090108);
   nominal__9->SetBinError(23,0.4674474);
   nominal__9->SetBinError(24,0.2780403);
   nominal__9->SetBinError(25,0.1499971);
   nominal__9->SetBinError(26,0.07513942);
   nominal__9->SetBinError(27,0.03499855);
   nominal__9->SetBinError(28,0.0246351);
   nominal__9->SetBinError(29,0.03134007);
   nominal__9->SetBinError(30,0.05761964);
   nominal__9->SetBinError(31,0.1230883);
   nominal__9->SetBinError(32,0.096954);
   nominal__9->SetMinimum(5e-05);
   nominal__9->SetMaximum(14699.55);
   nominal__9->SetEntries(24.42943);
   nominal__9->SetFillColor(46);
   nominal__9->SetLineColor(46);
   nominal__9->SetMarkerColor(46);
   nominal__9->SetMarkerStyle(21);
   nominal__9->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__9->GetXaxis()->SetRange(1,32);
   nominal__9->GetXaxis()->SetLabelFont(42);
   nominal__9->GetXaxis()->SetLabelSize(0.035);
   nominal__9->GetXaxis()->SetTitleSize(0.035);
   nominal__9->GetXaxis()->SetTitleFont(42);
   nominal__9->GetYaxis()->SetTitle("Tracks");
   nominal__9->GetYaxis()->SetLabelFont(42);
   nominal__9->GetYaxis()->SetLabelSize(0.05);
   nominal__9->GetYaxis()->SetTitleSize(0.07);
   nominal__9->GetYaxis()->SetTitleOffset(0);
   nominal__9->GetYaxis()->SetTitleFont(42);
   nominal__9->GetZaxis()->SetLabelFont(42);
   nominal__9->GetZaxis()->SetLabelSize(0.035);
   nominal__9->GetZaxis()->SetTitleSize(0.035);
   nominal__9->GetZaxis()->SetTitleFont(42);
   nominal__9->Draw("E0 same");
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
