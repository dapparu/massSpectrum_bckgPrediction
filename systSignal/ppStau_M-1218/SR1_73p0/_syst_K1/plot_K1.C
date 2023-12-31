void plot_K1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:09 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.343272,2428.571,4.105082);
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
   Double_t xAxis1044[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1342 = new TH1F("nominal__1342","",32, xAxis1044);
   nominal__1342->SetBinContent(9,6.686583e-07);
   nominal__1342->SetBinContent(14,1.271793e-06);
   nominal__1342->SetBinContent(17,1.662335e-06);
   nominal__1342->SetBinContent(18,6.686583e-07);
   nominal__1342->SetBinContent(19,1.776726e-06);
   nominal__1342->SetBinContent(20,3.525109e-06);
   nominal__1342->SetBinContent(21,1.846568e-06);
   nominal__1342->SetBinContent(22,2.949615e-06);
   nominal__1342->SetBinContent(23,5.817542e-06);
   nominal__1342->SetBinContent(24,8.994208e-06);
   nominal__1342->SetBinContent(25,3.708652e-05);
   nominal__1342->SetBinContent(26,0.000187238);
   nominal__1342->SetBinContent(27,0.0009509593);
   nominal__1342->SetBinContent(28,0.004909313);
   nominal__1342->SetBinContent(29,0.01155687);
   nominal__1342->SetBinContent(30,0.01003278);
   nominal__1342->SetBinContent(31,0.003467502);
   nominal__1342->SetBinContent(32,0.0009908726);
   nominal__1342->SetBinError(9,6.686583e-07);
   nominal__1342->SetBinError(14,9.020513e-07);
   nominal__1342->SetBinError(17,9.653158e-07);
   nominal__1342->SetBinError(18,6.686583e-07);
   nominal__1342->SetBinError(19,1.025849e-06);
   nominal__1342->SetBinError(20,1.459443e-06);
   nominal__1342->SetBinError(21,1.066628e-06);
   nominal__1342->SetBinError(22,1.321852e-06);
   nominal__1342->SetBinError(23,1.844384e-06);
   nominal__1342->SetBinError(24,2.326118e-06);
   nominal__1342->SetBinError(25,4.681799e-06);
   nominal__1342->SetBinError(26,1.057096e-05);
   nominal__1342->SetBinError(27,2.378964e-05);
   nominal__1342->SetBinError(28,5.404915e-05);
   nominal__1342->SetBinError(29,8.293589e-05);
   nominal__1342->SetBinError(30,7.733293e-05);
   nominal__1342->SetBinError(31,4.548308e-05);
   nominal__1342->SetBinError(32,1.963845e-05);
   nominal__1342->SetBinError(33,1.43543e-05);
   nominal__1342->SetMinimum(5e-05);
   nominal__1342->SetMaximum(11556.87);
   nominal__1342->SetEntries(54202);
   nominal__1342->SetFillColor(1);
   nominal__1342->SetMarkerStyle(20);
   nominal__1342->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1342->GetXaxis()->SetRange(1,32);
   nominal__1342->GetXaxis()->SetLabelFont(42);
   nominal__1342->GetXaxis()->SetLabelSize(0.035);
   nominal__1342->GetXaxis()->SetTitleSize(0.035);
   nominal__1342->GetXaxis()->SetTitleFont(42);
   nominal__1342->GetYaxis()->SetTitle("Tracks");
   nominal__1342->GetYaxis()->SetLabelFont(42);
   nominal__1342->GetYaxis()->SetLabelSize(0.05);
   nominal__1342->GetYaxis()->SetTitleSize(0.07);
   nominal__1342->GetYaxis()->SetTitleOffset(0);
   nominal__1342->GetYaxis()->SetTitleFont(42);
   nominal__1342->GetZaxis()->SetLabelFont(42);
   nominal__1342->GetZaxis()->SetLabelSize(0.035);
   nominal__1342->GetZaxis()->SetTitleSize(0.035);
   nominal__1342->GetZaxis()->SetTitleFont(42);
   nominal__1342->Draw("");
   Double_t xAxis1045[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__1343 = new TH1F("K_down1__1343","",32, xAxis1045);
   K_down1__1343->SetBinContent(9,6.686583e-07);
   K_down1__1343->SetBinContent(14,1.271793e-06);
   K_down1__1343->SetBinContent(17,1.662335e-06);
   K_down1__1343->SetBinContent(18,6.686583e-07);
   K_down1__1343->SetBinContent(19,1.187122e-06);
   K_down1__1343->SetBinContent(20,4.114713e-06);
   K_down1__1343->SetBinContent(21,1.846568e-06);
   K_down1__1343->SetBinContent(22,2.949615e-06);
   K_down1__1343->SetBinContent(23,5.817542e-06);
   K_down1__1343->SetBinContent(24,8.994207e-06);
   K_down1__1343->SetBinContent(25,3.532332e-05);
   K_down1__1343->SetBinContent(26,0.0001841485);
   K_down1__1343->SetBinContent(27,0.0009246812);
   K_down1__1343->SetBinContent(28,0.004837002);
   K_down1__1343->SetBinContent(29,0.01149327);
   K_down1__1343->SetBinContent(30,0.01012716);
   K_down1__1343->SetBinContent(31,0.003521028);
   K_down1__1343->SetBinContent(32,0.001010015);
   K_down1__1343->SetBinError(9,6.686583e-07);
   K_down1__1343->SetBinError(14,9.020513e-07);
   K_down1__1343->SetBinError(17,9.653158e-07);
   K_down1__1343->SetBinError(18,6.686583e-07);
   K_down1__1343->SetBinError(19,8.394841e-07);
   K_down1__1343->SetBinError(20,1.574041e-06);
   K_down1__1343->SetBinError(21,1.066628e-06);
   K_down1__1343->SetBinError(22,1.321852e-06);
   K_down1__1343->SetBinError(23,1.844384e-06);
   K_down1__1343->SetBinError(24,2.326118e-06);
   K_down1__1343->SetBinError(25,4.569176e-06);
   K_down1__1343->SetBinError(26,1.048042e-05);
   K_down1__1343->SetBinError(27,2.345662e-05);
   K_down1__1343->SetBinError(28,5.364945e-05);
   K_down1__1343->SetBinError(29,8.270726e-05);
   K_down1__1343->SetBinError(30,7.76931e-05);
   K_down1__1343->SetBinError(31,4.583645e-05);
   K_down1__1343->SetBinError(32,1.979501e-05);
   K_down1__1343->SetBinError(33,1.453837e-05);
   K_down1__1343->SetEntries(54202);
   K_down1__1343->SetFillColor(38);
   K_down1__1343->SetLineColor(38);
   K_down1__1343->SetMarkerColor(38);
   K_down1__1343->SetMarkerStyle(21);
   K_down1__1343->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__1343->GetXaxis()->SetRange(1,400);
   K_down1__1343->GetXaxis()->SetLabelFont(42);
   K_down1__1343->GetXaxis()->SetLabelSize(0.035);
   K_down1__1343->GetXaxis()->SetTitleSize(0.035);
   K_down1__1343->GetXaxis()->SetTitleFont(42);
   K_down1__1343->GetYaxis()->SetTitle("Events / bin");
   K_down1__1343->GetYaxis()->SetLabelFont(42);
   K_down1__1343->GetYaxis()->SetLabelSize(0.035);
   K_down1__1343->GetYaxis()->SetTitleSize(0.035);
   K_down1__1343->GetYaxis()->SetTitleOffset(0);
   K_down1__1343->GetYaxis()->SetTitleFont(42);
   K_down1__1343->GetZaxis()->SetLabelFont(42);
   K_down1__1343->GetZaxis()->SetLabelSize(0.035);
   K_down1__1343->GetZaxis()->SetTitleSize(0.035);
   K_down1__1343->GetZaxis()->SetTitleFont(42);
   K_down1__1343->Draw("same");
   Double_t xAxis1046[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__1344 = new TH1F("K_up1__1344","",32, xAxis1046);
   K_up1__1344->SetBinContent(9,6.686583e-07);
   K_up1__1344->SetBinContent(14,1.271793e-06);
   K_up1__1344->SetBinContent(17,1.662335e-06);
   K_up1__1344->SetBinContent(18,6.686583e-07);
   K_up1__1344->SetBinContent(19,1.776726e-06);
   K_up1__1344->SetBinContent(20,3.525109e-06);
   K_up1__1344->SetBinContent(21,1.846568e-06);
   K_up1__1344->SetBinContent(22,2.949615e-06);
   K_up1__1344->SetBinContent(23,5.817542e-06);
   K_up1__1344->SetBinContent(24,8.994208e-06);
   K_up1__1344->SetBinContent(25,3.830998e-05);
   K_up1__1344->SetBinContent(26,0.0001901636);
   K_up1__1344->SetBinContent(27,0.0009807335);
   K_up1__1344->SetBinContent(28,0.004990528);
   K_up1__1344->SetBinContent(29,0.01162045);
   K_up1__1344->SetBinContent(30,0.009928904);
   K_up1__1344->SetBinContent(31,0.003407067);
   K_up1__1344->SetBinContent(32,0.0009764647);
   K_up1__1344->SetBinError(9,6.686583e-07);
   K_up1__1344->SetBinError(14,9.020513e-07);
   K_up1__1344->SetBinError(17,9.653158e-07);
   K_up1__1344->SetBinError(18,6.686583e-07);
   K_up1__1344->SetBinError(19,1.025849e-06);
   K_up1__1344->SetBinError(20,1.459443e-06);
   K_up1__1344->SetBinError(21,1.066628e-06);
   K_up1__1344->SetBinError(22,1.321852e-06);
   K_up1__1344->SetBinError(23,1.844384e-06);
   K_up1__1344->SetBinError(24,2.326118e-06);
   K_up1__1344->SetBinError(25,4.761117e-06);
   K_up1__1344->SetBinError(26,1.065182e-05);
   K_up1__1344->SetBinError(27,2.416078e-05);
   K_up1__1344->SetBinError(28,5.449727e-05);
   K_up1__1344->SetBinError(29,8.316679e-05);
   K_up1__1344->SetBinError(30,7.692499e-05);
   K_up1__1344->SetBinError(31,4.508928e-05);
   K_up1__1344->SetBinError(32,1.950411e-05);
   K_up1__1344->SetBinError(33,1.423525e-05);
   K_up1__1344->SetEntries(54202);
   K_up1__1344->SetFillColor(46);
   K_up1__1344->SetLineColor(46);
   K_up1__1344->SetMarkerColor(46);
   K_up1__1344->SetMarkerStyle(21);
   K_up1__1344->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__1344->GetXaxis()->SetRange(1,400);
   K_up1__1344->GetXaxis()->SetLabelFont(42);
   K_up1__1344->GetXaxis()->SetLabelSize(0.035);
   K_up1__1344->GetXaxis()->SetTitleSize(0.035);
   K_up1__1344->GetXaxis()->SetTitleFont(42);
   K_up1__1344->GetYaxis()->SetTitle("Events / bin");
   K_up1__1344->GetYaxis()->SetLabelFont(42);
   K_up1__1344->GetYaxis()->SetLabelSize(0.035);
   K_up1__1344->GetYaxis()->SetTitleSize(0.035);
   K_up1__1344->GetYaxis()->SetTitleOffset(0);
   K_up1__1344->GetYaxis()->SetTitleFont(42);
   K_up1__1344->GetZaxis()->SetLabelFont(42);
   K_up1__1344->GetZaxis()->SetLabelSize(0.035);
   K_up1__1344->GetZaxis()->SetTitleSize(0.035);
   K_up1__1344->GetZaxis()->SetTitleFont(42);
   K_up1__1344->Draw("same");
   TLine *line = new TLine(1730,0,1730,11556.87);
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
   
   TH1D *frameR2__1345 = new TH1D("frameR2__1345","",1,0,2000);
   frameR2__1345->SetMinimum(0.5);
   frameR2__1345->SetMaximum(1.5);
   frameR2__1345->SetStats(0);
   frameR2__1345->SetLineStyle(0);
   frameR2__1345->SetMarkerStyle(20);
   frameR2__1345->GetXaxis()->SetRange(1,1);
   frameR2__1345->GetXaxis()->SetLabelFont(43);
   frameR2__1345->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1345->GetXaxis()->SetLabelSize(16);
   frameR2__1345->GetXaxis()->SetTitleSize(24);
   frameR2__1345->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1345->GetXaxis()->SetTitleFont(43);
   frameR2__1345->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1345->GetYaxis()->SetNdivisions(205);
   frameR2__1345->GetYaxis()->SetLabelFont(43);
   frameR2__1345->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1345->GetYaxis()->SetLabelSize(20);
   frameR2__1345->GetYaxis()->SetTitleSize(20);
   frameR2__1345->GetYaxis()->SetTitleOffset(2);
   frameR2__1345->GetYaxis()->SetTitleFont(43);
   frameR2__1345->GetZaxis()->SetLabelFont(42);
   frameR2__1345->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1345->GetZaxis()->SetLabelSize(0.05);
   frameR2__1345->GetZaxis()->SetTitleSize(0.06);
   frameR2__1345->GetZaxis()->SetTitleFont(42);
   frameR2__1345->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1047[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__1346 = new TH1F("K_down1__1346","",32, xAxis1047);
   K_down1__1346->SetBinContent(0,1);
   K_down1__1346->SetBinContent(1,1);
   K_down1__1346->SetBinContent(2,1);
   K_down1__1346->SetBinContent(3,1);
   K_down1__1346->SetBinContent(4,1);
   K_down1__1346->SetBinContent(5,1);
   K_down1__1346->SetBinContent(6,1);
   K_down1__1346->SetBinContent(7,1);
   K_down1__1346->SetBinContent(8,1);
   K_down1__1346->SetBinContent(9,1);
   K_down1__1346->SetBinContent(10,1);
   K_down1__1346->SetBinContent(11,1);
   K_down1__1346->SetBinContent(12,1);
   K_down1__1346->SetBinContent(13,1);
   K_down1__1346->SetBinContent(14,1);
   K_down1__1346->SetBinContent(15,1);
   K_down1__1346->SetBinContent(16,1);
   K_down1__1346->SetBinContent(17,1);
   K_down1__1346->SetBinContent(18,1);
   K_down1__1346->SetBinContent(19,1);
   K_down1__1346->SetBinContent(20,0.9999816);
   K_down1__1346->SetBinContent(21,1);
   K_down1__1346->SetBinContent(22,1);
   K_down1__1346->SetBinContent(23,1);
   K_down1__1346->SetBinContent(24,1);
   K_down1__1346->SetBinContent(25,1);
   K_down1__1346->SetBinContent(26,0.999945);
   K_down1__1346->SetBinContent(27,0.9998479);
   K_down1__1346->SetBinContent(28,0.9989954);
   K_down1__1346->SetBinContent(29,0.9960445);
   K_down1__1346->SetBinContent(30,0.9886039);
   K_down1__1346->SetBinContent(31,0.983962);
   K_down1__1346->SetBinContent(32,0.9810473);
   K_down1__1346->SetBinError(0,0.006086186);
   K_down1__1346->SetBinError(1,0.006086186);
   K_down1__1346->SetBinError(2,0.006086186);
   K_down1__1346->SetBinError(3,0.006086186);
   K_down1__1346->SetBinError(4,0.006086186);
   K_down1__1346->SetBinError(5,0.006086186);
   K_down1__1346->SetBinError(6,0.006086186);
   K_down1__1346->SetBinError(7,0.006086186);
   K_down1__1346->SetBinError(8,0.006086186);
   K_down1__1346->SetBinError(9,0.006086186);
   K_down1__1346->SetBinError(10,0.006086242);
   K_down1__1346->SetBinError(11,0.006086242);
   K_down1__1346->SetBinError(12,0.006086242);
   K_down1__1346->SetBinError(13,0.006086242);
   K_down1__1346->SetBinError(14,0.006086242);
   K_down1__1346->SetBinError(15,0.006086353);
   K_down1__1346->SetBinError(16,0.006086353);
   K_down1__1346->SetBinError(17,0.006086353);
   K_down1__1346->SetBinError(18,0.00608652);
   K_down1__1346->SetBinError(19,0.006086575);
   K_down1__1346->SetBinError(20,0.006086605);
   K_down1__1346->SetBinError(21,0.006087073);
   K_down1__1346->SetBinError(22,0.006087242);
   K_down1__1346->SetBinError(23,0.006087523);
   K_down1__1346->SetBinError(24,0.006088084);
   K_down1__1346->SetBinError(25,0.006088927);
   K_down1__1346->SetBinError(26,0.006092053);
   K_down1__1346->SetBinError(27,0.00610913);
   K_down1__1346->SetBinError(28,0.006195733);
   K_down1__1346->SetBinError(29,0.006730137);
   K_down1__1346->SetBinError(30,0.00894274);
   K_down1__1346->SetBinError(31,0.01603466);
   K_down1__1346->SetBinError(32,0.03389897);
   K_down1__1346->SetEntries(33);
   K_down1__1346->SetFillColor(38);
   K_down1__1346->SetLineColor(38);
   K_down1__1346->SetMarkerColor(38);
   K_down1__1346->SetMarkerStyle(21);
   K_down1__1346->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__1346->GetXaxis()->SetRange(1,400);
   K_down1__1346->GetXaxis()->SetLabelFont(42);
   K_down1__1346->GetXaxis()->SetLabelSize(0.035);
   K_down1__1346->GetXaxis()->SetTitleSize(0.035);
   K_down1__1346->GetXaxis()->SetTitleFont(42);
   K_down1__1346->GetYaxis()->SetTitle("Events / bin");
   K_down1__1346->GetYaxis()->SetLabelFont(42);
   K_down1__1346->GetYaxis()->SetLabelSize(0.035);
   K_down1__1346->GetYaxis()->SetTitleSize(0.035);
   K_down1__1346->GetYaxis()->SetTitleOffset(0);
   K_down1__1346->GetYaxis()->SetTitleFont(42);
   K_down1__1346->GetZaxis()->SetLabelFont(42);
   K_down1__1346->GetZaxis()->SetLabelSize(0.035);
   K_down1__1346->GetZaxis()->SetTitleSize(0.035);
   K_down1__1346->GetZaxis()->SetTitleFont(42);
   K_down1__1346->Draw("E0 same");
   Double_t xAxis1048[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__1347 = new TH1F("K_up1__1347","",32, xAxis1048);
   K_up1__1347->SetBinContent(0,0.9999999);
   K_up1__1347->SetBinContent(1,0.9999999);
   K_up1__1347->SetBinContent(2,0.9999999);
   K_up1__1347->SetBinContent(3,0.9999999);
   K_up1__1347->SetBinContent(4,0.9999999);
   K_up1__1347->SetBinContent(5,0.9999999);
   K_up1__1347->SetBinContent(6,0.9999999);
   K_up1__1347->SetBinContent(7,0.9999999);
   K_up1__1347->SetBinContent(8,0.9999999);
   K_up1__1347->SetBinContent(9,0.9999999);
   K_up1__1347->SetBinContent(10,0.9999999);
   K_up1__1347->SetBinContent(11,0.9999999);
   K_up1__1347->SetBinContent(12,0.9999999);
   K_up1__1347->SetBinContent(13,0.9999999);
   K_up1__1347->SetBinContent(14,0.9999999);
   K_up1__1347->SetBinContent(15,0.9999999);
   K_up1__1347->SetBinContent(16,0.9999999);
   K_up1__1347->SetBinContent(17,0.9999999);
   K_up1__1347->SetBinContent(18,0.9999999);
   K_up1__1347->SetBinContent(19,0.9999999);
   K_up1__1347->SetBinContent(20,0.9999999);
   K_up1__1347->SetBinContent(21,0.9999999);
   K_up1__1347->SetBinContent(22,0.9999999);
   K_up1__1347->SetBinContent(23,0.9999999);
   K_up1__1347->SetBinContent(24,0.9999999);
   K_up1__1347->SetBinContent(25,0.9999999);
   K_up1__1347->SetBinContent(26,1.000038);
   K_up1__1347->SetBinContent(27,1.00013);
   K_up1__1347->SetBinContent(28,1.001097);
   K_up1__1347->SetBinContent(29,1.00444);
   K_up1__1347->SetBinContent(30,1.012487);
   K_up1__1347->SetBinContent(31,1.017074);
   K_up1__1347->SetBinContent(32,1.014755);
   K_up1__1347->SetBinError(0,0.006086186);
   K_up1__1347->SetBinError(1,0.006086186);
   K_up1__1347->SetBinError(2,0.006086186);
   K_up1__1347->SetBinError(3,0.006086186);
   K_up1__1347->SetBinError(4,0.006086186);
   K_up1__1347->SetBinError(5,0.006086186);
   K_up1__1347->SetBinError(6,0.006086186);
   K_up1__1347->SetBinError(7,0.006086186);
   K_up1__1347->SetBinError(8,0.006086186);
   K_up1__1347->SetBinError(9,0.006086186);
   K_up1__1347->SetBinError(10,0.006086242);
   K_up1__1347->SetBinError(11,0.006086242);
   K_up1__1347->SetBinError(12,0.006086242);
   K_up1__1347->SetBinError(13,0.006086242);
   K_up1__1347->SetBinError(14,0.006086242);
   K_up1__1347->SetBinError(15,0.006086353);
   K_up1__1347->SetBinError(16,0.006086353);
   K_up1__1347->SetBinError(17,0.006086353);
   K_up1__1347->SetBinError(18,0.00608652);
   K_up1__1347->SetBinError(19,0.006086575);
   K_up1__1347->SetBinError(20,0.006086744);
   K_up1__1347->SetBinError(21,0.006087073);
   K_up1__1347->SetBinError(22,0.006087242);
   K_up1__1347->SetBinError(23,0.006087523);
   K_up1__1347->SetBinError(24,0.006088084);
   K_up1__1347->SetBinError(25,0.006088927);
   K_up1__1347->SetBinError(26,0.00609276);
   K_up1__1347->SetBinError(27,0.00611128);
   K_up1__1347->SetBinError(28,0.006212013);
   K_up1__1347->SetBinError(29,0.006801077);
   K_up1__1347->SetBinError(30,0.009213427);
   K_up1__1347->SetBinError(31,0.01671191);
   K_up1__1347->SetBinError(32,0.03535899);
   K_up1__1347->SetEntries(33);
   K_up1__1347->SetFillColor(46);
   K_up1__1347->SetLineColor(46);
   K_up1__1347->SetMarkerColor(46);
   K_up1__1347->SetMarkerStyle(21);
   K_up1__1347->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__1347->GetXaxis()->SetRange(1,400);
   K_up1__1347->GetXaxis()->SetLabelFont(42);
   K_up1__1347->GetXaxis()->SetLabelSize(0.035);
   K_up1__1347->GetXaxis()->SetTitleSize(0.035);
   K_up1__1347->GetXaxis()->SetTitleFont(42);
   K_up1__1347->GetYaxis()->SetTitle("Events / bin");
   K_up1__1347->GetYaxis()->SetLabelFont(42);
   K_up1__1347->GetYaxis()->SetLabelSize(0.035);
   K_up1__1347->GetYaxis()->SetTitleSize(0.035);
   K_up1__1347->GetYaxis()->SetTitleOffset(0);
   K_up1__1347->GetYaxis()->SetTitleFont(42);
   K_up1__1347->GetZaxis()->SetLabelFont(42);
   K_up1__1347->GetZaxis()->SetLabelSize(0.035);
   K_up1__1347->GetZaxis()->SetTitleSize(0.035);
   K_up1__1347->GetZaxis()->SetTitleFont(42);
   K_up1__1347->Draw("E0 same");
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
   
   TH1D *frameR2__1348 = new TH1D("frameR2__1348","",1,0,2000);
   frameR2__1348->SetMinimum(0.5);
   frameR2__1348->SetMaximum(1.5);
   frameR2__1348->SetStats(0);
   frameR2__1348->SetLineStyle(0);
   frameR2__1348->SetMarkerStyle(20);
   frameR2__1348->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1348->GetXaxis()->SetRange(1,1);
   frameR2__1348->GetXaxis()->SetLabelFont(43);
   frameR2__1348->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1348->GetXaxis()->SetLabelSize(16);
   frameR2__1348->GetXaxis()->SetTitleSize(24);
   frameR2__1348->GetXaxis()->SetTitleOffset(5);
   frameR2__1348->GetXaxis()->SetTitleFont(43);
   frameR2__1348->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1348->GetYaxis()->SetNdivisions(205);
   frameR2__1348->GetYaxis()->SetLabelFont(43);
   frameR2__1348->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1348->GetYaxis()->SetLabelSize(20);
   frameR2__1348->GetYaxis()->SetTitleSize(20);
   frameR2__1348->GetYaxis()->SetTitleOffset(2);
   frameR2__1348->GetYaxis()->SetTitleFont(43);
   frameR2__1348->GetZaxis()->SetLabelFont(42);
   frameR2__1348->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1348->GetZaxis()->SetLabelSize(0.05);
   frameR2__1348->GetZaxis()->SetTitleSize(0.06);
   frameR2__1348->GetZaxis()->SetTitleFont(42);
   frameR2__1348->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1049[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1349 = new TH1F("nominal__1349","",32, xAxis1049);
   nominal__1349->SetBinContent(9,1);
   nominal__1349->SetBinContent(14,1);
   nominal__1349->SetBinContent(17,1);
   nominal__1349->SetBinContent(18,1);
   nominal__1349->SetBinContent(19,1.496666);
   nominal__1349->SetBinContent(20,0.8567084);
   nominal__1349->SetBinContent(21,1);
   nominal__1349->SetBinContent(22,1);
   nominal__1349->SetBinContent(23,1);
   nominal__1349->SetBinContent(24,1);
   nominal__1349->SetBinContent(25,1.049916);
   nominal__1349->SetBinContent(26,1.016777);
   nominal__1349->SetBinContent(27,1.028419);
   nominal__1349->SetBinContent(28,1.01495);
   nominal__1349->SetBinContent(29,1.005534);
   nominal__1349->SetBinContent(30,0.9906808);
   nominal__1349->SetBinContent(31,0.9847981);
   nominal__1349->SetBinContent(32,0.9810473);
   nominal__1349->SetBinError(9,1.414214);
   nominal__1349->SetBinError(14,1.003067);
   nominal__1349->SetBinError(17,0.8212318);
   nominal__1349->SetBinError(18,1.414214);
   nominal__1349->SetBinError(19,1.366353);
   nominal__1349->SetBinError(20,0.4829159);
   nominal__1349->SetBinError(21,0.8168882);
   nominal__1349->SetBinError(22,0.6337714);
   nominal__1349->SetBinError(23,0.4483599);
   nominal__1349->SetBinError(24,0.3657496);
   nominal__1349->SetBinError(25,0.1897669);
   nominal__1349->SetBinError(26,0.08151045);
   nominal__1349->SetBinError(27,0.03664002);
   nominal__1349->SetBinError(28,0.01586149);
   nominal__1349->SetBinError(29,0.01021912);
   nominal__1349->SetBinError(30,0.01077383);
   nominal__1349->SetBinError(31,0.01819935);
   nominal__1349->SetBinError(32,0.02734495);
   nominal__1349->SetMinimum(5e-05);
   nominal__1349->SetMaximum(11556.87);
   nominal__1349->SetEntries(36.78206);
   nominal__1349->SetFillColor(38);
   nominal__1349->SetLineColor(38);
   nominal__1349->SetMarkerColor(38);
   nominal__1349->SetMarkerStyle(21);
   nominal__1349->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1349->GetXaxis()->SetRange(1,32);
   nominal__1349->GetXaxis()->SetLabelFont(42);
   nominal__1349->GetXaxis()->SetLabelSize(0.035);
   nominal__1349->GetXaxis()->SetTitleSize(0.035);
   nominal__1349->GetXaxis()->SetTitleFont(42);
   nominal__1349->GetYaxis()->SetTitle("Tracks");
   nominal__1349->GetYaxis()->SetLabelFont(42);
   nominal__1349->GetYaxis()->SetLabelSize(0.05);
   nominal__1349->GetYaxis()->SetTitleSize(0.07);
   nominal__1349->GetYaxis()->SetTitleOffset(0);
   nominal__1349->GetYaxis()->SetTitleFont(42);
   nominal__1349->GetZaxis()->SetLabelFont(42);
   nominal__1349->GetZaxis()->SetLabelSize(0.035);
   nominal__1349->GetZaxis()->SetTitleSize(0.035);
   nominal__1349->GetZaxis()->SetTitleFont(42);
   nominal__1349->Draw("E0 same");
   Double_t xAxis1050[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1350 = new TH1F("nominal__1350","",32, xAxis1050);
   nominal__1350->SetBinContent(9,1);
   nominal__1350->SetBinContent(14,1);
   nominal__1350->SetBinContent(17,1);
   nominal__1350->SetBinContent(18,1);
   nominal__1350->SetBinContent(19,1);
   nominal__1350->SetBinContent(20,1);
   nominal__1350->SetBinContent(21,1);
   nominal__1350->SetBinContent(22,1);
   nominal__1350->SetBinContent(23,1);
   nominal__1350->SetBinContent(24,1);
   nominal__1350->SetBinContent(25,0.9680643);
   nominal__1350->SetBinContent(26,0.9846157);
   nominal__1350->SetBinContent(27,0.9696408);
   nominal__1350->SetBinContent(28,0.9837263);
   nominal__1350->SetBinContent(29,0.9945287);
   nominal__1350->SetBinContent(30,1.010462);
   nominal__1350->SetBinContent(31,1.017738);
   nominal__1350->SetBinContent(32,1.014755);
   nominal__1350->SetBinError(9,1.414214);
   nominal__1350->SetBinError(14,1.003067);
   nominal__1350->SetBinError(17,0.8212318);
   nominal__1350->SetBinError(18,1.414214);
   nominal__1350->SetBinError(19,0.8165409);
   nominal__1350->SetBinError(20,0.5855034);
   nominal__1350->SetBinError(21,0.8168882);
   nominal__1350->SetBinError(22,0.6337714);
   nominal__1350->SetBinError(23,0.4483599);
   nominal__1350->SetBinError(24,0.3657495);
   nominal__1350->SetBinError(25,0.1714915);
   nominal__1350->SetBinError(26,0.07830637);
   nominal__1350->SetBinError(27,0.0340443);
   nominal__1350->SetBinError(28,0.01525438);
   nominal__1350->SetBinError(29,0.0100797);
   nominal__1350->SetBinError(30,0.01104314);
   nominal__1350->SetBinError(31,0.01896367);
   nominal__1350->SetBinError(32,0.0285537);
   nominal__1350->SetMinimum(5e-05);
   nominal__1350->SetMaximum(11556.87);
   nominal__1350->SetEntries(39.59192);
   nominal__1350->SetFillColor(46);
   nominal__1350->SetLineColor(46);
   nominal__1350->SetMarkerColor(46);
   nominal__1350->SetMarkerStyle(21);
   nominal__1350->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1350->GetXaxis()->SetRange(1,32);
   nominal__1350->GetXaxis()->SetLabelFont(42);
   nominal__1350->GetXaxis()->SetLabelSize(0.035);
   nominal__1350->GetXaxis()->SetTitleSize(0.035);
   nominal__1350->GetXaxis()->SetTitleFont(42);
   nominal__1350->GetYaxis()->SetTitle("Tracks");
   nominal__1350->GetYaxis()->SetLabelFont(42);
   nominal__1350->GetYaxis()->SetLabelSize(0.05);
   nominal__1350->GetYaxis()->SetTitleSize(0.07);
   nominal__1350->GetYaxis()->SetTitleOffset(0);
   nominal__1350->GetYaxis()->SetTitleFont(42);
   nominal__1350->GetZaxis()->SetLabelFont(42);
   nominal__1350->GetZaxis()->SetLabelSize(0.035);
   nominal__1350->GetZaxis()->SetTitleSize(0.035);
   nominal__1350->GetZaxis()->SetTitleFont(42);
   nominal__1350->Draw("E0 same");
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
