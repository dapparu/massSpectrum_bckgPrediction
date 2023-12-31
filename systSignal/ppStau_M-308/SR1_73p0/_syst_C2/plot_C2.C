void plot_C2()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:04 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.376516,2428.571,10.72063);
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
   Double_t xAxis624[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__802 = new TH1F("nominal__802","",32, xAxis624);
   nominal__802->SetBinContent(2,0.004900225);
   nominal__802->SetBinContent(3,0.03129842);
   nominal__802->SetBinContent(4,0.01028152);
   nominal__802->SetBinContent(5,0.01646051);
   nominal__802->SetBinContent(6,0.005006641);
   nominal__802->SetBinContent(7,0.02525013);
   nominal__802->SetBinContent(8,0.03048938);
   nominal__802->SetBinContent(9,0.1046266);
   nominal__802->SetBinContent(10,0.5669704);
   nominal__802->SetBinContent(11,2.218784);
   nominal__802->SetBinContent(12,6.319841);
   nominal__802->SetBinContent(13,14.93229);
   nominal__802->SetBinContent(14,29.01737);
   nominal__802->SetBinContent(15,44.17146);
   nominal__802->SetBinContent(16,43.6709);
   nominal__802->SetBinContent(17,25.78515);
   nominal__802->SetBinContent(18,11.90388);
   nominal__802->SetBinContent(19,4.887452);
   nominal__802->SetBinContent(20,2.550593);
   nominal__802->SetBinContent(21,1.034574);
   nominal__802->SetBinContent(22,0.4728294);
   nominal__802->SetBinContent(23,0.269127);
   nominal__802->SetBinContent(24,0.110839);
   nominal__802->SetBinContent(25,0.05180309);
   nominal__802->SetBinContent(26,0.0291109);
   nominal__802->SetBinContent(27,0.01644041);
   nominal__802->SetBinContent(28,0.009588117);
   nominal__802->SetBinContent(29,0.004968045);
   nominal__802->SetBinError(2,0.004900225);
   nominal__802->SetBinError(3,0.01278375);
   nominal__802->SetBinError(4,0.007270153);
   nominal__802->SetBinError(5,0.00953037);
   nominal__802->SetBinError(6,0.005006641);
   nominal__802->SetBinError(7,0.01129546);
   nominal__802->SetBinError(8,0.01244841);
   nominal__802->SetBinError(9,0.02291324);
   nominal__802->SetBinError(10,0.05345823);
   nominal__802->SetBinError(11,0.1059821);
   nominal__802->SetBinError(12,0.1788637);
   nominal__802->SetBinError(13,0.2751561);
   nominal__802->SetBinError(14,0.3836934);
   nominal__802->SetBinError(15,0.4731686);
   nominal__802->SetBinError(16,0.4710379);
   nominal__802->SetBinError(17,0.3620423);
   nominal__802->SetBinError(18,0.2463207);
   nominal__802->SetBinError(19,0.1579102);
   nominal__802->SetBinError(20,0.1142782);
   nominal__802->SetBinError(21,0.07273229);
   nominal__802->SetBinError(22,0.04910741);
   nominal__802->SetBinError(23,0.03702299);
   nominal__802->SetBinError(24,0.02369192);
   nominal__802->SetBinError(25,0.01639507);
   nominal__802->SetBinError(26,0.0119664);
   nominal__802->SetBinError(27,0.009510709);
   nominal__802->SetBinError(28,0.006795033);
   nominal__802->SetBinError(29,0.004968045);
   nominal__802->SetMinimum(5e-05);
   nominal__802->SetMaximum(4.417146e+10);
   nominal__802->SetEntries(37224);
   nominal__802->SetFillColor(1);
   nominal__802->SetMarkerStyle(20);
   nominal__802->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__802->GetXaxis()->SetRange(1,32);
   nominal__802->GetXaxis()->SetLabelFont(42);
   nominal__802->GetXaxis()->SetLabelSize(0.035);
   nominal__802->GetXaxis()->SetTitleSize(0.035);
   nominal__802->GetXaxis()->SetTitleFont(42);
   nominal__802->GetYaxis()->SetTitle("Tracks");
   nominal__802->GetYaxis()->SetLabelFont(42);
   nominal__802->GetYaxis()->SetLabelSize(0.05);
   nominal__802->GetYaxis()->SetTitleSize(0.07);
   nominal__802->GetYaxis()->SetTitleOffset(0);
   nominal__802->GetYaxis()->SetTitleFont(42);
   nominal__802->GetZaxis()->SetLabelFont(42);
   nominal__802->GetZaxis()->SetLabelSize(0.035);
   nominal__802->GetZaxis()->SetTitleSize(0.035);
   nominal__802->GetZaxis()->SetTitleFont(42);
   nominal__802->Draw("");
   Double_t xAxis625[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_down2__803 = new TH1F("C_down2__803","",32, xAxis625);
   C_down2__803->SetBinContent(3,0.02016019);
   C_down2__803->SetBinContent(4,0.02631998);
   C_down2__803->SetBinContent(5,0.01038991);
   C_down2__803->SetBinContent(6,0.01107724);
   C_down2__803->SetBinContent(7,0.02525013);
   C_down2__803->SetBinContent(8,0.02028629);
   C_down2__803->SetBinContent(9,0.0945803);
   C_down2__803->SetBinContent(10,0.3572257);
   C_down2__803->SetBinContent(11,1.867678);
   C_down2__803->SetBinContent(12,5.416124);
   C_down2__803->SetBinContent(13,13.81671);
   C_down2__803->SetBinContent(14,27.75295);
   C_down2__803->SetBinContent(15,43.8954);
   C_down2__803->SetBinContent(16,44.53931);
   C_down2__803->SetBinContent(17,27.25611);
   C_down2__803->SetBinContent(18,12.865);
   C_down2__803->SetBinContent(19,5.284717);
   C_down2__803->SetBinContent(20,2.866648);
   C_down2__803->SetBinContent(21,1.099858);
   C_down2__803->SetBinContent(22,0.4887711);
   C_down2__803->SetBinContent(23,0.3045148);
   C_down2__803->SetBinContent(24,0.1066702);
   C_down2__803->SetBinContent(25,0.06640694);
   C_down2__803->SetBinContent(26,0.0291109);
   C_down2__803->SetBinContent(27,0.00525908);
   C_down2__803->SetBinContent(28,0.02076945);
   C_down2__803->SetBinContent(29,0.004968045);
   C_down2__803->SetBinError(3,0.01008213);
   C_down2__803->SetBinError(4,0.01177465);
   C_down2__803->SetBinError(5,0.007346822);
   C_down2__803->SetBinError(6,0.007868838);
   C_down2__803->SetBinError(7,0.01129546);
   C_down2__803->SetBinError(8,0.0101445);
   C_down2__803->SetBinError(9,0.02178399);
   C_down2__803->SetBinError(10,0.04245806);
   C_down2__803->SetBinError(11,0.09717854);
   C_down2__803->SetBinError(12,0.1655394);
   C_down2__803->SetBinError(13,0.2647945);
   C_down2__803->SetBinError(14,0.3750434);
   C_down2__803->SetBinError(15,0.4718119);
   C_down2__803->SetBinError(16,0.4756238);
   C_down2__803->SetBinError(17,0.3721895);
   C_down2__803->SetBinError(18,0.2560999);
   C_down2__803->SetBinError(19,0.1641347);
   C_down2__803->SetBinError(20,0.1211485);
   C_down2__803->SetBinError(21,0.07496062);
   C_down2__803->SetBinError(22,0.04996514);
   C_down2__803->SetBinError(23,0.03937072);
   C_down2__803->SetBinError(24,0.02333641);
   C_down2__803->SetBinError(25,0.01844055);
   C_down2__803->SetBinError(26,0.0119664);
   C_down2__803->SetBinError(27,0.00525908);
   C_down2__803->SetBinError(28,0.01043878);
   C_down2__803->SetBinError(29,0.004968045);
   C_down2__803->SetEntries(37224);
   C_down2__803->SetFillColor(38);
   C_down2__803->SetLineColor(38);
   C_down2__803->SetMarkerColor(38);
   C_down2__803->SetMarkerStyle(21);
   C_down2__803->GetXaxis()->SetTitle("Mass [GeV]");
   C_down2__803->GetXaxis()->SetRange(1,400);
   C_down2__803->GetXaxis()->SetLabelFont(42);
   C_down2__803->GetXaxis()->SetLabelSize(0.035);
   C_down2__803->GetXaxis()->SetTitleSize(0.035);
   C_down2__803->GetXaxis()->SetTitleFont(42);
   C_down2__803->GetYaxis()->SetTitle("Events / bin");
   C_down2__803->GetYaxis()->SetLabelFont(42);
   C_down2__803->GetYaxis()->SetLabelSize(0.035);
   C_down2__803->GetYaxis()->SetTitleSize(0.035);
   C_down2__803->GetYaxis()->SetTitleOffset(0);
   C_down2__803->GetYaxis()->SetTitleFont(42);
   C_down2__803->GetZaxis()->SetLabelFont(42);
   C_down2__803->GetZaxis()->SetLabelSize(0.035);
   C_down2__803->GetZaxis()->SetTitleSize(0.035);
   C_down2__803->GetZaxis()->SetTitleFont(42);
   C_down2__803->Draw("same");
   Double_t xAxis626[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_up2__804 = new TH1F("C_up2__804","",32, xAxis626);
   C_up2__804->SetBinContent(1,0.004900225);
   C_up2__804->SetBinContent(2,0.005176427);
   C_up2__804->SetBinContent(3,0.02612199);
   C_up2__804->SetBinContent(4,0.01635212);
   C_up2__804->SetBinContent(5,0.01038991);
   C_up2__804->SetBinContent(6,0.009824847);
   C_up2__804->SetBinContent(7,0.02043193);
   C_up2__804->SetBinContent(8,0.05940935);
   C_up2__804->SetBinContent(9,0.1308517);
   C_up2__804->SetBinContent(10,0.776821);
   C_up2__804->SetBinContent(11,2.707711);
   C_up2__804->SetBinContent(12,7.140691);
   C_up2__804->SetBinContent(13,15.99357);
   C_up2__804->SetBinContent(14,30.0258);
   C_up2__804->SetBinContent(15,44.88433);
   C_up2__804->SetBinContent(16,42.26543);
   C_up2__804->SetBinContent(17,24.47414);
   C_up2__804->SetBinContent(18,11.14432);
   C_up2__804->SetBinContent(19,4.396082);
   C_up2__804->SetBinContent(20,2.322708);
   C_up2__804->SetBinContent(21,0.9433003);
   C_up2__804->SetBinContent(22,0.4490614);
   C_up2__804->SetBinContent(23,0.2523448);
   C_up2__804->SetBinContent(24,0.1008265);
   C_up2__804->SetBinContent(25,0.03160117);
   C_up2__804->SetBinContent(26,0.0291109);
   C_up2__804->SetBinContent(27,0.01644041);
   C_up2__804->SetBinContent(28,0.01455616);
   C_up2__804->SetBinError(0,0.004900225);
   C_up2__804->SetBinError(2,0.005176427);
   C_up2__804->SetBinError(3,0.01168884);
   C_up2__804->SetBinError(4,0.009471393);
   C_up2__804->SetBinError(5,0.007346822);
   C_up2__804->SetBinError(6,0.006948493);
   C_up2__804->SetBinError(7,0.01021628);
   C_up2__804->SetBinError(8,0.01719879);
   C_up2__804->SetBinError(9,0.02573789);
   C_up2__804->SetBinError(10,0.06253052);
   C_up2__804->SetBinError(11,0.1170805);
   C_up2__804->SetBinError(12,0.1902449);
   C_up2__804->SetBinError(13,0.2846825);
   C_up2__804->SetBinError(14,0.3903037);
   C_up2__804->SetBinError(15,0.4771146);
   C_up2__804->SetBinError(16,0.46334);
   C_up2__804->SetBinError(17,0.3527414);
   C_up2__804->SetBinError(18,0.2383562);
   C_up2__804->SetBinError(19,0.1497906);
   C_up2__804->SetBinError(20,0.1089735);
   C_up2__804->SetBinError(21,0.0694756);
   C_up2__804->SetBinError(22,0.0479363);
   C_up2__804->SetBinError(23,0.03573611);
   C_up2__804->SetBinError(24,0.02260741);
   C_up2__804->SetBinError(25,0.01291351);
   C_up2__804->SetBinError(26,0.0119664);
   C_up2__804->SetBinError(27,0.009510709);
   C_up2__804->SetBinError(28,0.008417478);
   C_up2__804->SetEntries(37224);
   C_up2__804->SetFillColor(46);
   C_up2__804->SetLineColor(46);
   C_up2__804->SetMarkerColor(46);
   C_up2__804->SetMarkerStyle(21);
   C_up2__804->GetXaxis()->SetTitle("Mass [GeV]");
   C_up2__804->GetXaxis()->SetRange(1,400);
   C_up2__804->GetXaxis()->SetLabelFont(42);
   C_up2__804->GetXaxis()->SetLabelSize(0.035);
   C_up2__804->GetXaxis()->SetTitleSize(0.035);
   C_up2__804->GetXaxis()->SetTitleFont(42);
   C_up2__804->GetYaxis()->SetTitle("Events / bin");
   C_up2__804->GetYaxis()->SetLabelFont(42);
   C_up2__804->GetYaxis()->SetLabelSize(0.035);
   C_up2__804->GetYaxis()->SetTitleSize(0.035);
   C_up2__804->GetYaxis()->SetTitleOffset(0);
   C_up2__804->GetYaxis()->SetTitleFont(42);
   C_up2__804->GetZaxis()->SetLabelFont(42);
   C_up2__804->GetZaxis()->SetLabelSize(0.035);
   C_up2__804->GetZaxis()->SetTitleSize(0.035);
   C_up2__804->GetZaxis()->SetTitleFont(42);
   C_up2__804->Draw("same");
   TLine *line = new TLine(1730,0,1730,4.417146e+10);
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
   entry=leg->AddEntry("C_down2","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("C_up2","Up","PE1");
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
   
   TH1D *frameR2__805 = new TH1D("frameR2__805","",1,0,2000);
   frameR2__805->SetMinimum(0.5);
   frameR2__805->SetMaximum(1.5);
   frameR2__805->SetStats(0);
   frameR2__805->SetLineStyle(0);
   frameR2__805->SetMarkerStyle(20);
   frameR2__805->GetXaxis()->SetRange(1,1);
   frameR2__805->GetXaxis()->SetLabelFont(43);
   frameR2__805->GetXaxis()->SetLabelOffset(0.007);
   frameR2__805->GetXaxis()->SetLabelSize(16);
   frameR2__805->GetXaxis()->SetTitleSize(24);
   frameR2__805->GetXaxis()->SetTitleOffset(3.75);
   frameR2__805->GetXaxis()->SetTitleFont(43);
   frameR2__805->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__805->GetYaxis()->SetNdivisions(205);
   frameR2__805->GetYaxis()->SetLabelFont(43);
   frameR2__805->GetYaxis()->SetLabelOffset(0.007);
   frameR2__805->GetYaxis()->SetLabelSize(20);
   frameR2__805->GetYaxis()->SetTitleSize(20);
   frameR2__805->GetYaxis()->SetTitleOffset(2);
   frameR2__805->GetYaxis()->SetTitleFont(43);
   frameR2__805->GetZaxis()->SetLabelFont(42);
   frameR2__805->GetZaxis()->SetLabelOffset(0.007);
   frameR2__805->GetZaxis()->SetLabelSize(0.05);
   frameR2__805->GetZaxis()->SetTitleSize(0.06);
   frameR2__805->GetZaxis()->SetTitleFont(42);
   frameR2__805->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis627[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_down2__806 = new TH1F("C_down2__806","",32, xAxis627);
   C_down2__806->SetBinContent(0,1);
   C_down2__806->SetBinContent(1,1);
   C_down2__806->SetBinContent(2,1);
   C_down2__806->SetBinContent(3,0.999974);
   C_down2__806->SetBinContent(4,0.9999148);
   C_down2__806->SetBinContent(5,1);
   C_down2__806->SetBinContent(6,0.9999678);
   C_down2__806->SetBinContent(7,1);
   C_down2__806->SetBinContent(8,1);
   C_down2__806->SetBinContent(9,0.9999458);
   C_down2__806->SetBinContent(10,0.9998924);
   C_down2__806->SetBinContent(11,0.9987746);
   C_down2__806->SetBinContent(12,0.9968728);
   C_down2__806->SetBinContent(13,0.9917694);
   C_down2__806->SetBinContent(14,0.9843901);
   C_down2__806->SetBinContent(15,0.9721622);
   C_down2__806->SetBinContent(16,0.9563834);
   C_down2__806->SetBinContent(17,0.9350687);
   C_down2__806->SetBinContent(18,0.922157);
   C_down2__806->SetBinContent(19,0.9182338);
   C_down2__806->SetBinContent(20,0.9112546);
   C_down2__806->SetBinContent(21,0.9402496);
   C_down2__806->SetBinContent(22,0.9398282);
   C_down2__806->SetBinContent(23,0.9147797);
   C_down2__806->SetBinContent(24,0.9552497);
   C_down2__806->SetBinContent(25,0.8845677);
   C_down2__806->SetBinContent(26,1);
   C_down2__806->SetBinContent(27,1);
   C_down2__806->SetBinContent(28,0.5655625);
   C_down2__806->SetBinContent(29,1);
   C_down2__806->SetBinError(0,0.007344791);
   C_down2__806->SetBinError(1,0.007344791);
   C_down2__806->SetBinError(2,0.007344791);
   C_down2__806->SetBinError(3,0.007344649);
   C_down2__806->SetBinError(4,0.007344709);
   C_down2__806->SetBinError(5,0.007345681);
   C_down2__806->SetBinError(6,0.00734569);
   C_down2__806->SetBinError(7,0.007346074);
   C_down2__806->SetBinError(8,0.007346569);
   C_down2__806->SetBinError(9,0.007346666);
   C_down2__806->SetBinError(10,0.007348242);
   C_down2__806->SetBinError(11,0.007349127);
   C_down2__806->SetBinError(12,0.007375567);
   C_down2__806->SetBinError(13,0.00745714);
   C_down2__806->SetBinError(14,0.007717376);
   C_down2__806->SetBinError(15,0.008376406);
   C_down2__806->SetBinError(16,0.01001168);
   C_down2__806->SetBinError(17,0.01352001);
   C_down2__806->SetBinError(18,0.01976695);
   C_down2__806->SetBinError(19,0.02958308);
   C_down2__806->SetBinError(20,0.04223398);
   C_down2__806->SetBinError(21,0.06618166);
   C_down2__806->SetBinError(22,0.09516128);
   C_down2__806->SetBinError(23,0.1288158);
   C_down2__806->SetBinError(24,0.2020448);
   C_down2__806->SetBinError(25,0.2594758);
   C_down2__806->SetBinError(26,0.4105649);
   C_down2__806->SetBinError(27,0.5794668);
   C_down2__806->SetBinError(28,0.4141227);
   C_down2__806->SetBinError(29,1.414214);
   C_down2__806->SetEntries(33);
   C_down2__806->SetFillColor(38);
   C_down2__806->SetLineColor(38);
   C_down2__806->SetMarkerColor(38);
   C_down2__806->SetMarkerStyle(21);
   C_down2__806->GetXaxis()->SetTitle("Mass [GeV]");
   C_down2__806->GetXaxis()->SetRange(1,400);
   C_down2__806->GetXaxis()->SetLabelFont(42);
   C_down2__806->GetXaxis()->SetLabelSize(0.035);
   C_down2__806->GetXaxis()->SetTitleSize(0.035);
   C_down2__806->GetXaxis()->SetTitleFont(42);
   C_down2__806->GetYaxis()->SetTitle("Events / bin");
   C_down2__806->GetYaxis()->SetLabelFont(42);
   C_down2__806->GetYaxis()->SetLabelSize(0.035);
   C_down2__806->GetYaxis()->SetTitleSize(0.035);
   C_down2__806->GetYaxis()->SetTitleOffset(0);
   C_down2__806->GetYaxis()->SetTitleFont(42);
   C_down2__806->GetZaxis()->SetLabelFont(42);
   C_down2__806->GetZaxis()->SetLabelSize(0.035);
   C_down2__806->GetZaxis()->SetTitleSize(0.035);
   C_down2__806->GetZaxis()->SetTitleFont(42);
   C_down2__806->Draw("E0 same");
   Double_t xAxis628[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_up2__807 = new TH1F("C_up2__807","",32, xAxis628);
   C_up2__807->SetBinContent(0,1);
   C_up2__807->SetBinContent(1,1);
   C_up2__807->SetBinContent(2,1.000026);
   C_up2__807->SetBinContent(3,1.000028);
   C_up2__807->SetBinContent(4,1);
   C_up2__807->SetBinContent(5,1.000032);
   C_up2__807->SetBinContent(6,1);
   C_up2__807->SetBinContent(7,1.000026);
   C_up2__807->SetBinContent(8,1);
   C_up2__807->SetBinContent(9,1.000154);
   C_up2__807->SetBinContent(10,1.000293);
   C_up2__807->SetBinContent(11,1.001416);
   C_up2__807->SetBinContent(12,1.004087);
   C_up2__807->SetBinContent(13,1.00888);
   C_up2__807->SetBinContent(14,1.016337);
   C_up2__807->SetBinContent(15,1.027752);
   C_up2__807->SetBinContent(16,1.050409);
   C_up2__807->SetBinContent(17,1.066823);
   C_up2__807->SetBinContent(18,1.083291);
   C_up2__807->SetBinContent(19,1.103003);
   C_up2__807->SetBinContent(20,1.093733);
   C_up2__807->SetBinContent(21,1.088196);
   C_up2__807->SetBinContent(22,1.07916);
   C_up2__807->SetBinContent(23,1.105639);
   C_up2__807->SetBinContent(24,1.15693);
   C_up2__807->SetBinContent(25,1.220284);
   C_up2__807->SetBinContent(26,1);
   C_up2__807->SetBinContent(27,1);
   C_up2__807->SetBinContent(28,1);
   C_up2__807->SetBinError(0,0.00734479);
   C_up2__807->SetBinError(1,0.007344744);
   C_up2__807->SetBinError(2,0.007345031);
   C_up2__807->SetBinError(3,0.007345141);
   C_up2__807->SetBinError(4,0.007345483);
   C_up2__807->SetBinError(5,0.007345965);
   C_up2__807->SetBinError(6,0.007345974);
   C_up2__807->SetBinError(7,0.007346311);
   C_up2__807->SetBinError(8,0.007346568);
   C_up2__807->SetBinError(9,0.007348587);
   C_up2__807->SetBinError(10,0.00735193);
   C_up2__807->SetBinError(11,0.007373475);
   C_up2__807->SetBinError(12,0.0074424);
   C_up2__807->SetBinError(13,0.007618297);
   C_up2__807->SetBinError(14,0.008032023);
   C_up2__807->SetBinError(15,0.008979338);
   C_up2__807->SetBinError(16,0.01125731);
   C_up2__807->SetBinError(17,0.01594223);
   C_up2__807->SetBinError(18,0.02417537);
   C_up2__807->SetBinError(19,0.03721082);
   C_up2__807->SetBinError(20,0.05304777);
   C_up2__807->SetBinError(21,0.0794731);
   C_up2__807->SetBinError(22,0.1131114);
   C_up2__807->SetBinError(23,0.1631368);
   C_up2__807->SetBinError(24,0.2570554);
   C_up2__807->SetBinError(25,0.3893621);
   C_up2__807->SetBinError(26,0.4105649);
   C_up2__807->SetBinError(27,0.5794668);
   C_up2__807->SetBinError(28,0.8178057);
   C_up2__807->SetEntries(33);
   C_up2__807->SetFillColor(46);
   C_up2__807->SetLineColor(46);
   C_up2__807->SetMarkerColor(46);
   C_up2__807->SetMarkerStyle(21);
   C_up2__807->GetXaxis()->SetTitle("Mass [GeV]");
   C_up2__807->GetXaxis()->SetRange(1,400);
   C_up2__807->GetXaxis()->SetLabelFont(42);
   C_up2__807->GetXaxis()->SetLabelSize(0.035);
   C_up2__807->GetXaxis()->SetTitleSize(0.035);
   C_up2__807->GetXaxis()->SetTitleFont(42);
   C_up2__807->GetYaxis()->SetTitle("Events / bin");
   C_up2__807->GetYaxis()->SetLabelFont(42);
   C_up2__807->GetYaxis()->SetLabelSize(0.035);
   C_up2__807->GetYaxis()->SetTitleSize(0.035);
   C_up2__807->GetYaxis()->SetTitleOffset(0);
   C_up2__807->GetYaxis()->SetTitleFont(42);
   C_up2__807->GetZaxis()->SetLabelFont(42);
   C_up2__807->GetZaxis()->SetLabelSize(0.035);
   C_up2__807->GetZaxis()->SetTitleSize(0.035);
   C_up2__807->GetZaxis()->SetTitleFont(42);
   C_up2__807->Draw("E0 same");
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
   
   TH1D *frameR2__808 = new TH1D("frameR2__808","",1,0,2000);
   frameR2__808->SetMinimum(0.5);
   frameR2__808->SetMaximum(1.5);
   frameR2__808->SetStats(0);
   frameR2__808->SetLineStyle(0);
   frameR2__808->SetMarkerStyle(20);
   frameR2__808->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__808->GetXaxis()->SetRange(1,1);
   frameR2__808->GetXaxis()->SetLabelFont(43);
   frameR2__808->GetXaxis()->SetLabelOffset(0.007);
   frameR2__808->GetXaxis()->SetLabelSize(16);
   frameR2__808->GetXaxis()->SetTitleSize(24);
   frameR2__808->GetXaxis()->SetTitleOffset(5);
   frameR2__808->GetXaxis()->SetTitleFont(43);
   frameR2__808->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__808->GetYaxis()->SetNdivisions(205);
   frameR2__808->GetYaxis()->SetLabelFont(43);
   frameR2__808->GetYaxis()->SetLabelOffset(0.007);
   frameR2__808->GetYaxis()->SetLabelSize(20);
   frameR2__808->GetYaxis()->SetTitleSize(20);
   frameR2__808->GetYaxis()->SetTitleOffset(2);
   frameR2__808->GetYaxis()->SetTitleFont(43);
   frameR2__808->GetZaxis()->SetLabelFont(42);
   frameR2__808->GetZaxis()->SetLabelOffset(0.007);
   frameR2__808->GetZaxis()->SetLabelSize(0.05);
   frameR2__808->GetZaxis()->SetTitleSize(0.06);
   frameR2__808->GetZaxis()->SetTitleFont(42);
   frameR2__808->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis629[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__809 = new TH1F("nominal__809","",32, xAxis629);
   nominal__809->SetBinContent(3,1.552486);
   nominal__809->SetBinContent(4,0.3906358);
   nominal__809->SetBinContent(5,1.584278);
   nominal__809->SetBinContent(6,0.4519756);
   nominal__809->SetBinContent(7,1);
   nominal__809->SetBinContent(8,1.502955);
   nominal__809->SetBinContent(9,1.10622);
   nominal__809->SetBinContent(10,1.587149);
   nominal__809->SetBinContent(11,1.187991);
   nominal__809->SetBinContent(12,1.166857);
   nominal__809->SetBinContent(13,1.080742);
   nominal__809->SetBinContent(14,1.04556);
   nominal__809->SetBinContent(15,1.006289);
   nominal__809->SetBinContent(16,0.9805022);
   nominal__809->SetBinContent(17,0.9460318);
   nominal__809->SetBinContent(18,0.9252912);
   nominal__809->SetBinContent(19,0.9248276);
   nominal__809->SetBinContent(20,0.8897477);
   nominal__809->SetBinContent(21,0.9406428);
   nominal__809->SetBinContent(22,0.9673842);
   nominal__809->SetBinContent(23,0.8837895);
   nominal__809->SetBinContent(24,1.039081);
   nominal__809->SetBinContent(25,0.7800855);
   nominal__809->SetBinContent(26,1);
   nominal__809->SetBinContent(27,3.126101);
   nominal__809->SetBinContent(28,0.4616452);
   nominal__809->SetBinContent(29,1);
   nominal__809->SetBinError(3,1.002442);
   nominal__809->SetBinError(4,0.3268616);
   nominal__809->SetBinError(5,1.447885);
   nominal__809->SetBinError(6,0.5544052);
   nominal__809->SetBinError(7,0.6326381);
   nominal__809->SetBinError(8,0.9702676);
   nominal__809->SetBinError(9,0.3515789);
   nominal__809->SetBinError(10,0.2407902);
   nominal__809->SetBinError(11,0.08391014);
   nominal__809->SetBinError(12,0.04860584);
   nominal__809->SetBinError(13,0.02873311);
   nominal__809->SetBinError(14,0.01976808);
   nominal__809->SetBinError(15,0.01527042);
   nominal__809->SetBinError(16,0.01488217);
   nominal__809->SetBinError(17,0.01852891);
   nominal__809->SetBinError(18,0.0265682);
   nominal__809->SetBinError(19,0.04144749);
   nominal__809->SetBinError(20,0.05480061);
   nominal__809->SetBinError(21,0.09210333);
   nominal__809->SetBinError(22,0.1409754);
   nominal__809->SetBinError(23,0.1668481);
   nominal__809->SetBinError(24,0.3178134);
   nominal__809->SetBinError(25,0.3284489);
   nominal__809->SetBinError(26,0.5813303);
   nominal__809->SetBinError(27,3.611502);
   nominal__809->SetBinError(28,0.4010885);
   nominal__809->SetBinError(29,1.414214);
   nominal__809->SetMinimum(5e-05);
   nominal__809->SetMaximum(4.417146e+10);
   nominal__809->SetEntries(41.7938);
   nominal__809->SetFillColor(38);
   nominal__809->SetLineColor(38);
   nominal__809->SetMarkerColor(38);
   nominal__809->SetMarkerStyle(21);
   nominal__809->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__809->GetXaxis()->SetRange(1,32);
   nominal__809->GetXaxis()->SetLabelFont(42);
   nominal__809->GetXaxis()->SetLabelSize(0.035);
   nominal__809->GetXaxis()->SetTitleSize(0.035);
   nominal__809->GetXaxis()->SetTitleFont(42);
   nominal__809->GetYaxis()->SetTitle("Tracks");
   nominal__809->GetYaxis()->SetLabelFont(42);
   nominal__809->GetYaxis()->SetLabelSize(0.05);
   nominal__809->GetYaxis()->SetTitleSize(0.07);
   nominal__809->GetYaxis()->SetTitleOffset(0);
   nominal__809->GetYaxis()->SetTitleFont(42);
   nominal__809->GetZaxis()->SetLabelFont(42);
   nominal__809->GetZaxis()->SetLabelSize(0.035);
   nominal__809->GetZaxis()->SetTitleSize(0.035);
   nominal__809->GetZaxis()->SetTitleFont(42);
   nominal__809->Draw("E0 same");
   Double_t xAxis630[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__810 = new TH1F("nominal__810","",32, xAxis630);
   nominal__810->SetBinContent(2,0.9466422);
   nominal__810->SetBinContent(3,1.198164);
   nominal__810->SetBinContent(4,0.6287578);
   nominal__810->SetBinContent(5,1.584278);
   nominal__810->SetBinContent(6,0.5095897);
   nominal__810->SetBinContent(7,1.235817);
   nominal__810->SetBinContent(8,0.5132084);
   nominal__810->SetBinContent(9,0.7995817);
   nominal__810->SetBinContent(10,0.7298598);
   nominal__810->SetBinContent(11,0.8194315);
   nominal__810->SetBinContent(12,0.8850462);
   nominal__810->SetBinContent(13,0.9336438);
   nominal__810->SetBinContent(14,0.9664146);
   nominal__810->SetBinContent(15,0.9841176);
   nominal__810->SetBinContent(16,1.033253);
   nominal__810->SetBinContent(17,1.053567);
   nominal__810->SetBinContent(18,1.068157);
   nominal__810->SetBinContent(19,1.111774);
   nominal__810->SetBinContent(20,1.098112);
   nominal__810->SetBinContent(21,1.096759);
   nominal__810->SetBinContent(22,1.052928);
   nominal__810->SetBinContent(23,1.066505);
   nominal__810->SetBinContent(24,1.099305);
   nominal__810->SetBinContent(25,1.639278);
   nominal__810->SetBinContent(26,1);
   nominal__810->SetBinContent(27,1);
   nominal__810->SetBinContent(28,0.6586981);
   nominal__810->SetBinError(2,1.338754);
   nominal__810->SetBinError(3,0.7259127);
   nominal__810->SetBinError(4,0.5747177);
   nominal__810->SetBinError(5,1.447885);
   nominal__810->SetBinError(6,0.6241557);
   nominal__810->SetBinError(7,0.8291321);
   nominal__810->SetBinError(8,0.2568638);
   nominal__810->SetBinError(9,0.2353679);
   nominal__810->SetBinError(10,0.0904839);
   nominal__810->SetBinError(11,0.05279608);
   nominal__810->SetBinError(12,0.03440104);
   nominal__810->SetBinError(13,0.02391995);
   nominal__810->SetBinError(14,0.01791956);
   nominal__810->SetBinError(15,0.01485147);
   nominal__810->SetBinError(16,0.01589058);
   nominal__810->SetBinError(17,0.02119927);
   nominal__810->SetBinError(18,0.03178787);
   nominal__810->SetBinError(19,0.05220495);
   nominal__810->SetBinError(20,0.07123872);
   nominal__810->SetBinError(21,0.1116698);
   nominal__810->SetBinError(22,0.1568181);
   nominal__810->SetBinError(23,0.2105633);
   nominal__810->SetBinError(24,0.340544);
   nominal__810->SetBinError(25,0.8472887);
   nominal__810->SetBinError(26,0.5813303);
   nominal__810->SetBinError(27,0.8181165);
   nominal__810->SetBinError(28,0.6025014);
   nominal__810->SetMinimum(5e-05);
   nominal__810->SetMaximum(4.417146e+10);
   nominal__810->SetEntries(86.46699);
   nominal__810->SetFillColor(46);
   nominal__810->SetLineColor(46);
   nominal__810->SetMarkerColor(46);
   nominal__810->SetMarkerStyle(21);
   nominal__810->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__810->GetXaxis()->SetRange(1,32);
   nominal__810->GetXaxis()->SetLabelFont(42);
   nominal__810->GetXaxis()->SetLabelSize(0.035);
   nominal__810->GetXaxis()->SetTitleSize(0.035);
   nominal__810->GetXaxis()->SetTitleFont(42);
   nominal__810->GetYaxis()->SetTitle("Tracks");
   nominal__810->GetYaxis()->SetLabelFont(42);
   nominal__810->GetYaxis()->SetLabelSize(0.05);
   nominal__810->GetYaxis()->SetTitleSize(0.07);
   nominal__810->GetYaxis()->SetTitleOffset(0);
   nominal__810->GetYaxis()->SetTitleFont(42);
   nominal__810->GetZaxis()->SetLabelFont(42);
   nominal__810->GetZaxis()->SetLabelSize(0.035);
   nominal__810->GetZaxis()->SetTitleSize(0.035);
   nominal__810->GetZaxis()->SetTitleFont(42);
   nominal__810->Draw("E0 same");
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
