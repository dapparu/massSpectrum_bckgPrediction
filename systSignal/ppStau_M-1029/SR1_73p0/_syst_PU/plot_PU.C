void plot_PU()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:08 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.32537,2428.571,0.5426446);
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
   Double_t xAxis946[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1216 = new TH1F("nominal__1216","",32, xAxis946);
   nominal__1216->SetBinContent(2,1.563494e-05);
   nominal__1216->SetBinContent(4,1.486327e-05);
   nominal__1216->SetBinContent(6,1.603025e-05);
   nominal__1216->SetBinContent(9,4.928684e-05);
   nominal__1216->SetBinContent(11,1.636705e-05);
   nominal__1216->SetBinContent(12,3.200473e-05);
   nominal__1216->SetBinContent(14,4.609524e-05);
   nominal__1216->SetBinContent(15,1.663444e-05);
   nominal__1216->SetBinContent(16,1.623319e-05);
   nominal__1216->SetBinContent(17,1.603025e-05);
   nominal__1216->SetBinContent(18,8.742952e-05);
   nominal__1216->SetBinContent(20,4.660431e-05);
   nominal__1216->SetBinContent(21,5.074686e-05);
   nominal__1216->SetBinContent(22,0.0002131457);
   nominal__1216->SetBinContent(23,0.0004816079);
   nominal__1216->SetBinContent(24,0.001495932);
   nominal__1216->SetBinContent(25,0.004833851);
   nominal__1216->SetBinContent(26,0.02893999);
   nominal__1216->SetBinContent(27,0.1298854);
   nominal__1216->SetBinContent(28,0.3298409);
   nominal__1216->SetBinContent(29,0.261889);
   nominal__1216->SetBinContent(30,0.08305414);
   nominal__1216->SetBinContent(31,0.01510153);
   nominal__1216->SetBinContent(32,0.004513945);
   nominal__1216->SetBinError(2,1.563494e-05);
   nominal__1216->SetBinError(4,1.486327e-05);
   nominal__1216->SetBinError(6,1.603025e-05);
   nominal__1216->SetBinError(9,2.866041e-05);
   nominal__1216->SetBinError(11,1.636705e-05);
   nominal__1216->SetBinError(12,2.26308e-05);
   nominal__1216->SetBinError(14,2.671456e-05);
   nominal__1216->SetBinError(15,1.663444e-05);
   nominal__1216->SetBinError(16,1.623319e-05);
   nominal__1216->SetBinError(17,1.603025e-05);
   nominal__1216->SetBinError(18,3.943322e-05);
   nominal__1216->SetBinError(20,2.695298e-05);
   nominal__1216->SetBinError(21,2.936192e-05);
   nominal__1216->SetBinError(22,5.921806e-05);
   nominal__1216->SetBinError(23,8.803226e-05);
   nominal__1216->SetBinError(24,0.0001554508);
   nominal__1216->SetBinError(25,0.00028052);
   nominal__1216->SetBinError(26,0.0006857738);
   nominal__1216->SetBinError(27,0.001450656);
   nominal__1216->SetBinError(28,0.002312163);
   nominal__1216->SetBinError(29,0.002060617);
   nominal__1216->SetBinError(30,0.001162456);
   nominal__1216->SetBinError(31,0.0004959074);
   nominal__1216->SetBinError(32,0.0002099092);
   nominal__1216->SetBinError(33,0.0001709742);
   nominal__1216->SetMinimum(5e-05);
   nominal__1216->SetMaximum(3.298409);
   nominal__1216->SetEntries(53280);
   nominal__1216->SetFillColor(1);
   nominal__1216->SetMarkerStyle(20);
   nominal__1216->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1216->GetXaxis()->SetRange(1,32);
   nominal__1216->GetXaxis()->SetLabelFont(42);
   nominal__1216->GetXaxis()->SetLabelSize(0.035);
   nominal__1216->GetXaxis()->SetTitleSize(0.035);
   nominal__1216->GetXaxis()->SetTitleFont(42);
   nominal__1216->GetYaxis()->SetTitle("Tracks");
   nominal__1216->GetYaxis()->SetLabelFont(42);
   nominal__1216->GetYaxis()->SetLabelSize(0.05);
   nominal__1216->GetYaxis()->SetTitleSize(0.07);
   nominal__1216->GetYaxis()->SetTitleOffset(0);
   nominal__1216->GetYaxis()->SetTitleFont(42);
   nominal__1216->GetZaxis()->SetLabelFont(42);
   nominal__1216->GetZaxis()->SetLabelSize(0.035);
   nominal__1216->GetZaxis()->SetTitleSize(0.035);
   nominal__1216->GetZaxis()->SetTitleFont(42);
   nominal__1216->Draw("");
   Double_t xAxis947[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__1217 = new TH1F("PU_down__1217","",32, xAxis947);
   PU_down__1217->SetBinContent(2,1.34287e-05);
   PU_down__1217->SetBinContent(4,1.227001e-05);
   PU_down__1217->SetBinContent(6,1.433982e-05);
   PU_down__1217->SetBinContent(9,6.079708e-05);
   PU_down__1217->SetBinContent(11,1.595097e-05);
   PU_down__1217->SetBinContent(12,2.899158e-05);
   PU_down__1217->SetBinContent(14,4.195845e-05);
   PU_down__1217->SetBinContent(15,1.735631e-05);
   PU_down__1217->SetBinContent(16,1.517298e-05);
   PU_down__1217->SetBinContent(17,1.433982e-05);
   PU_down__1217->SetBinContent(18,8.775706e-05);
   PU_down__1217->SetBinContent(20,4.221857e-05);
   PU_down__1217->SetBinContent(21,5.838148e-05);
   PU_down__1217->SetBinContent(22,0.0002224387);
   PU_down__1217->SetBinContent(23,0.0004598979);
   PU_down__1217->SetBinContent(24,0.001481393);
   PU_down__1217->SetBinContent(25,0.004905493);
   PU_down__1217->SetBinContent(26,0.02895089);
   PU_down__1217->SetBinContent(27,0.1293565);
   PU_down__1217->SetBinContent(28,0.3288034);
   PU_down__1217->SetBinContent(29,0.2610845);
   PU_down__1217->SetBinContent(30,0.0838796);
   PU_down__1217->SetBinContent(31,0.01523133);
   PU_down__1217->SetBinContent(32,0.004521313);
   PU_down__1217->SetBinError(2,1.34287e-05);
   PU_down__1217->SetBinError(4,1.227001e-05);
   PU_down__1217->SetBinError(6,1.433982e-05);
   PU_down__1217->SetBinError(9,3.934021e-05);
   PU_down__1217->SetBinError(11,1.595097e-05);
   PU_down__1217->SetBinError(12,2.050133e-05);
   PU_down__1217->SetBinError(14,2.455977e-05);
   PU_down__1217->SetBinError(15,1.735631e-05);
   PU_down__1217->SetBinError(16,1.517298e-05);
   PU_down__1217->SetBinError(17,1.433982e-05);
   PU_down__1217->SetBinError(18,4.080923e-05);
   PU_down__1217->SetBinError(20,2.459322e-05);
   PU_down__1217->SetBinError(21,3.562654e-05);
   PU_down__1217->SetBinError(22,6.543674e-05);
   PU_down__1217->SetBinError(23,8.474138e-05);
   PU_down__1217->SetBinError(24,0.0001597314);
   PU_down__1217->SetBinError(25,0.0002938312);
   PU_down__1217->SetBinError(26,0.0007069478);
   PU_down__1217->SetBinError(27,0.001492538);
   PU_down__1217->SetBinError(28,0.002378223);
   PU_down__1217->SetBinError(29,0.002117611);
   PU_down__1217->SetBinError(30,0.001215101);
   PU_down__1217->SetBinError(31,0.0005164778);
   PU_down__1217->SetBinError(32,0.0002148141);
   PU_down__1217->SetBinError(33,0.0001773969);
   PU_down__1217->SetEntries(53280);
   PU_down__1217->SetFillColor(38);
   PU_down__1217->SetLineColor(38);
   PU_down__1217->SetMarkerColor(38);
   PU_down__1217->SetMarkerStyle(21);
   PU_down__1217->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__1217->GetXaxis()->SetRange(1,400);
   PU_down__1217->GetXaxis()->SetLabelFont(42);
   PU_down__1217->GetXaxis()->SetLabelSize(0.035);
   PU_down__1217->GetXaxis()->SetTitleSize(0.035);
   PU_down__1217->GetXaxis()->SetTitleFont(42);
   PU_down__1217->GetYaxis()->SetTitle("Events / bin");
   PU_down__1217->GetYaxis()->SetLabelFont(42);
   PU_down__1217->GetYaxis()->SetLabelSize(0.035);
   PU_down__1217->GetYaxis()->SetTitleSize(0.035);
   PU_down__1217->GetYaxis()->SetTitleOffset(0);
   PU_down__1217->GetYaxis()->SetTitleFont(42);
   PU_down__1217->GetZaxis()->SetLabelFont(42);
   PU_down__1217->GetZaxis()->SetLabelSize(0.035);
   PU_down__1217->GetZaxis()->SetTitleSize(0.035);
   PU_down__1217->GetZaxis()->SetTitleFont(42);
   PU_down__1217->Draw("same");
   Double_t xAxis948[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__1218 = new TH1F("PU_up__1218","",32, xAxis948);
   PU_up__1218->SetBinContent(2,1.792045e-05);
   PU_up__1218->SetBinContent(4,1.784777e-05);
   PU_up__1218->SetBinContent(6,1.770971e-05);
   PU_up__1218->SetBinContent(9,4.421497e-05);
   PU_up__1218->SetBinContent(11,1.621567e-05);
   PU_up__1218->SetBinContent(12,3.508212e-05);
   PU_up__1218->SetBinContent(14,5.037156e-05);
   PU_up__1218->SetBinContent(15,1.490026e-05);
   PU_up__1218->SetBinContent(16,1.722009e-05);
   PU_up__1218->SetBinContent(17,1.770971e-05);
   PU_up__1218->SetBinContent(18,8.939747e-05);
   PU_up__1218->SetBinContent(20,5.126142e-05);
   PU_up__1218->SetBinContent(21,4.516428e-05);
   PU_up__1218->SetBinContent(22,0.0002077501);
   PU_up__1218->SetBinContent(23,0.0004969265);
   PU_up__1218->SetBinContent(24,0.001522057);
   PU_up__1218->SetBinContent(25,0.004748026);
   PU_up__1218->SetBinContent(26,0.02888137);
   PU_up__1218->SetBinContent(27,0.1306304);
   PU_up__1218->SetBinContent(28,0.3309993);
   PU_up__1218->SetBinContent(29,0.2625431);
   PU_up__1218->SetBinContent(30,0.0823637);
   PU_up__1218->SetBinContent(31,0.01496002);
   PU_up__1218->SetBinContent(32,0.004482363);
   PU_up__1218->SetBinError(2,1.792045e-05);
   PU_up__1218->SetBinError(4,1.784777e-05);
   PU_up__1218->SetBinError(6,1.770971e-05);
   PU_up__1218->SetBinError(9,2.644902e-05);
   PU_up__1218->SetBinError(11,1.621567e-05);
   PU_up__1218->SetBinError(12,2.480795e-05);
   PU_up__1218->SetBinError(14,2.909057e-05);
   PU_up__1218->SetBinError(15,1.490026e-05);
   PU_up__1218->SetBinError(16,1.722009e-05);
   PU_up__1218->SetBinError(17,1.770971e-05);
   PU_up__1218->SetBinError(18,4.141363e-05);
   PU_up__1218->SetBinError(20,2.961679e-05);
   PU_up__1218->SetBinError(21,2.67378e-05);
   PU_up__1218->SetBinError(22,5.836411e-05);
   PU_up__1218->SetBinError(23,9.101763e-05);
   PU_up__1218->SetBinError(24,0.0001590062);
   PU_up__1218->SetBinError(25,0.0002775618);
   PU_up__1218->SetBinError(26,0.0006889034);
   PU_up__1218->SetBinError(27,0.001468293);
   PU_up__1218->SetBinError(28,0.00233494);
   PU_up__1218->SetBinError(29,0.002078755);
   PU_up__1218->SetBinError(30,0.001161179);
   PU_up__1218->SetBinError(31,0.0004947309);
   PU_up__1218->SetBinError(32,0.0002115017);
   PU_up__1218->SetBinError(33,0.0001686259);
   PU_up__1218->SetEntries(53280);
   PU_up__1218->SetFillColor(46);
   PU_up__1218->SetLineColor(46);
   PU_up__1218->SetMarkerColor(46);
   PU_up__1218->SetMarkerStyle(21);
   PU_up__1218->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__1218->GetXaxis()->SetRange(1,400);
   PU_up__1218->GetXaxis()->SetLabelFont(42);
   PU_up__1218->GetXaxis()->SetLabelSize(0.035);
   PU_up__1218->GetXaxis()->SetTitleSize(0.035);
   PU_up__1218->GetXaxis()->SetTitleFont(42);
   PU_up__1218->GetYaxis()->SetTitle("Events / bin");
   PU_up__1218->GetYaxis()->SetLabelFont(42);
   PU_up__1218->GetYaxis()->SetLabelSize(0.035);
   PU_up__1218->GetYaxis()->SetTitleSize(0.035);
   PU_up__1218->GetYaxis()->SetTitleOffset(0);
   PU_up__1218->GetYaxis()->SetTitleFont(42);
   PU_up__1218->GetZaxis()->SetLabelFont(42);
   PU_up__1218->GetZaxis()->SetLabelSize(0.035);
   PU_up__1218->GetZaxis()->SetTitleSize(0.035);
   PU_up__1218->GetZaxis()->SetTitleFont(42);
   PU_up__1218->Draw("same");
   TLine *line = new TLine(1730,0,1730,3.298409);
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
   
   TH1D *frameR2__1219 = new TH1D("frameR2__1219","",1,0,2000);
   frameR2__1219->SetMinimum(0.5);
   frameR2__1219->SetMaximum(1.5);
   frameR2__1219->SetStats(0);
   frameR2__1219->SetLineStyle(0);
   frameR2__1219->SetMarkerStyle(20);
   frameR2__1219->GetXaxis()->SetRange(1,1);
   frameR2__1219->GetXaxis()->SetLabelFont(43);
   frameR2__1219->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1219->GetXaxis()->SetLabelSize(16);
   frameR2__1219->GetXaxis()->SetTitleSize(24);
   frameR2__1219->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1219->GetXaxis()->SetTitleFont(43);
   frameR2__1219->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1219->GetYaxis()->SetNdivisions(205);
   frameR2__1219->GetYaxis()->SetLabelFont(43);
   frameR2__1219->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1219->GetYaxis()->SetLabelSize(20);
   frameR2__1219->GetYaxis()->SetTitleSize(20);
   frameR2__1219->GetYaxis()->SetTitleOffset(2);
   frameR2__1219->GetYaxis()->SetTitleFont(43);
   frameR2__1219->GetZaxis()->SetLabelFont(42);
   frameR2__1219->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1219->GetZaxis()->SetLabelSize(0.05);
   frameR2__1219->GetZaxis()->SetTitleSize(0.06);
   frameR2__1219->GetZaxis()->SetTitleFont(42);
   frameR2__1219->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis949[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_down__1220 = new TH1F("PU_down__1220","",32, xAxis949);
   PU_down__1220->SetBinContent(0,1.001575);
   PU_down__1220->SetBinContent(1,1.001575);
   PU_down__1220->SetBinContent(2,1.001575);
   PU_down__1220->SetBinContent(3,1.001573);
   PU_down__1220->SetBinContent(4,1.001573);
   PU_down__1220->SetBinContent(5,1.00157);
   PU_down__1220->SetBinContent(6,1.00157);
   PU_down__1220->SetBinContent(7,1.001568);
   PU_down__1220->SetBinContent(8,1.001568);
   PU_down__1220->SetBinContent(9,1.001568);
   PU_down__1220->SetBinContent(10,1.001581);
   PU_down__1220->SetBinContent(11,1.001581);
   PU_down__1220->SetBinContent(12,1.001581);
   PU_down__1220->SetBinContent(13,1.001577);
   PU_down__1220->SetBinContent(14,1.001577);
   PU_down__1220->SetBinContent(15,1.001573);
   PU_down__1220->SetBinContent(16,1.001574);
   PU_down__1220->SetBinContent(17,1.001572);
   PU_down__1220->SetBinContent(18,1.00157);
   PU_down__1220->SetBinContent(19,1.001571);
   PU_down__1220->SetBinContent(20,1.001571);
   PU_down__1220->SetBinContent(21,1.001566);
   PU_down__1220->SetBinContent(22,1.001575);
   PU_down__1220->SetBinContent(23,1.001586);
   PU_down__1220->SetBinContent(24,1.001562);
   PU_down__1220->SetBinContent(25,1.001547);
   PU_down__1220->SetBinContent(26,1.00164);
   PU_down__1220->SetBinContent(27,1.001711);
   PU_down__1220->SetBinContent(28,1.001268);
   PU_down__1220->SetBinContent(29,0.9995664);
   PU_down__1220->SetBinContent(30,0.990711);
   PU_down__1220->SetBinContent(31,0.9930555);
   PU_down__1220->SetBinContent(32,0.9983705);
   PU_down__1220->SetBinError(0,0.00624781);
   PU_down__1220->SetBinError(1,0.00624781);
   PU_down__1220->SetBinError(2,0.00624781);
   PU_down__1220->SetBinError(3,0.006247853);
   PU_down__1220->SetBinError(4,0.006247853);
   PU_down__1220->SetBinError(5,0.006247892);
   PU_down__1220->SetBinError(6,0.006247892);
   PU_down__1220->SetBinError(7,0.00624794);
   PU_down__1220->SetBinError(8,0.00624794);
   PU_down__1220->SetBinError(9,0.00624794);
   PU_down__1220->SetBinError(10,0.006248168);
   PU_down__1220->SetBinError(11,0.006248168);
   PU_down__1220->SetBinError(12,0.006248226);
   PU_down__1220->SetBinError(13,0.006248325);
   PU_down__1220->SetBinError(14,0.006248325);
   PU_down__1220->SetBinError(15,0.006248471);
   PU_down__1220->SetBinError(16,0.006248538);
   PU_down__1220->SetBinError(17,0.006248591);
   PU_down__1220->SetBinError(18,0.006248638);
   PU_down__1220->SetBinError(19,0.006248929);
   PU_down__1220->SetBinError(20,0.006248929);
   PU_down__1220->SetBinError(21,0.006249075);
   PU_down__1220->SetBinError(22,0.006249297);
   PU_down__1220->SetBinError(23,0.006250106);
   PU_down__1220->SetBinError(24,0.006251753);
   PU_down__1220->SetBinError(25,0.006257102);
   PU_down__1220->SetBinError(26,0.006275311);
   PU_down__1220->SetBinError(27,0.006384916);
   PU_down__1220->SetBinError(28,0.006953306);
   PU_down__1220->SetBinError(29,0.009583283);
   PU_down__1220->SetBinError(30,0.01794312);
   PU_down__1220->SetBinError(31,0.04109196);
   PU_down__1220->SetBinError(32,0.08584772);
   PU_down__1220->SetEntries(33);
   PU_down__1220->SetFillColor(38);
   PU_down__1220->SetLineColor(38);
   PU_down__1220->SetMarkerColor(38);
   PU_down__1220->SetMarkerStyle(21);
   PU_down__1220->GetXaxis()->SetTitle("Mass [GeV]");
   PU_down__1220->GetXaxis()->SetRange(1,400);
   PU_down__1220->GetXaxis()->SetLabelFont(42);
   PU_down__1220->GetXaxis()->SetLabelSize(0.035);
   PU_down__1220->GetXaxis()->SetTitleSize(0.035);
   PU_down__1220->GetXaxis()->SetTitleFont(42);
   PU_down__1220->GetYaxis()->SetTitle("Events / bin");
   PU_down__1220->GetYaxis()->SetLabelFont(42);
   PU_down__1220->GetYaxis()->SetLabelSize(0.035);
   PU_down__1220->GetYaxis()->SetTitleSize(0.035);
   PU_down__1220->GetYaxis()->SetTitleOffset(0);
   PU_down__1220->GetYaxis()->SetTitleFont(42);
   PU_down__1220->GetZaxis()->SetLabelFont(42);
   PU_down__1220->GetZaxis()->SetLabelSize(0.035);
   PU_down__1220->GetZaxis()->SetTitleSize(0.035);
   PU_down__1220->GetZaxis()->SetTitleFont(42);
   PU_down__1220->Draw("E0 same");
   Double_t xAxis950[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *PU_up__1221 = new TH1F("PU_up__1221","",32, xAxis950);
   PU_up__1221->SetBinContent(0,0.9981483);
   PU_up__1221->SetBinContent(1,0.9981483);
   PU_up__1221->SetBinContent(2,0.9981483);
   PU_up__1221->SetBinContent(3,0.9981509);
   PU_up__1221->SetBinContent(4,0.9981509);
   PU_up__1221->SetBinContent(5,0.9981543);
   PU_up__1221->SetBinContent(6,0.9981543);
   PU_up__1221->SetBinContent(7,0.9981562);
   PU_up__1221->SetBinContent(8,0.9981562);
   PU_up__1221->SetBinContent(9,0.9981562);
   PU_up__1221->SetBinContent(10,0.9981503);
   PU_up__1221->SetBinContent(11,0.9981503);
   PU_up__1221->SetBinContent(12,0.9981501);
   PU_up__1221->SetBinContent(13,0.9981536);
   PU_up__1221->SetBinContent(14,0.9981536);
   PU_up__1221->SetBinContent(15,0.9981584);
   PU_up__1221->SetBinContent(16,0.9981564);
   PU_up__1221->SetBinContent(17,0.9981575);
   PU_up__1221->SetBinContent(18,0.9981594);
   PU_up__1221->SetBinContent(19,0.9981615);
   PU_up__1221->SetBinContent(20,0.9981615);
   PU_up__1221->SetBinContent(21,0.9981668);
   PU_up__1221->SetBinContent(22,0.9981602);
   PU_up__1221->SetBinContent(23,0.9981535);
   PU_up__1221->SetBinContent(24,0.9981702);
   PU_up__1221->SetBinContent(25,0.9981974);
   PU_up__1221->SetBinContent(26,0.998087);
   PU_up__1221->SetBinContent(27,0.9979491);
   PU_up__1221->SetBinContent(28,0.9986353);
   PU_up__1221->SetBinContent(29,1.000575);
   PU_up__1221->SetBinContent(30,1.008482);
   PU_up__1221->SetBinContent(31,1.008903);
   PU_up__1221->SetBinContent(32,1.007046);
   PU_up__1221->SetBinError(0,0.006147384);
   PU_up__1221->SetBinError(1,0.006147384);
   PU_up__1221->SetBinError(2,0.006147384);
   PU_up__1221->SetBinError(3,0.006147458);
   PU_up__1221->SetBinError(4,0.006147458);
   PU_up__1221->SetBinError(5,0.006147537);
   PU_up__1221->SetBinError(6,0.006147537);
   PU_up__1221->SetBinError(7,0.006147606);
   PU_up__1221->SetBinError(8,0.006147606);
   PU_up__1221->SetBinError(9,0.006147606);
   PU_up__1221->SetBinError(10,0.006147737);
   PU_up__1221->SetBinError(11,0.006147737);
   PU_up__1221->SetBinError(12,0.006147794);
   PU_up__1221->SetBinError(13,0.006147932);
   PU_up__1221->SetBinError(14,0.006147932);
   PU_up__1221->SetBinError(15,0.006148136);
   PU_up__1221->SetBinError(16,0.006148181);
   PU_up__1221->SetBinError(17,0.006148246);
   PU_up__1221->SetBinError(18,0.006148316);
   PU_up__1221->SetBinError(19,0.006148603);
   PU_up__1221->SetBinError(20,0.006148603);
   PU_up__1221->SetBinError(21,0.00614881);
   PU_up__1221->SetBinError(22,0.00614894);
   PU_up__1221->SetBinError(23,0.006149646);
   PU_up__1221->SetBinError(24,0.006151492);
   PU_up__1221->SetBinError(25,0.006157047);
   PU_up__1221->SetBinError(26,0.006173697);
   PU_up__1221->SetBinError(27,0.006280032);
   PU_up__1221->SetBinError(28,0.006847511);
   PU_up__1221->SetBinError(29,0.009472333);
   PU_up__1221->SetBinError(30,0.01801897);
   PU_up__1221->SetBinError(31,0.04123813);
   PU_up__1221->SetBinError(32,0.08568083);
   PU_up__1221->SetEntries(33);
   PU_up__1221->SetFillColor(46);
   PU_up__1221->SetLineColor(46);
   PU_up__1221->SetMarkerColor(46);
   PU_up__1221->SetMarkerStyle(21);
   PU_up__1221->GetXaxis()->SetTitle("Mass [GeV]");
   PU_up__1221->GetXaxis()->SetRange(1,400);
   PU_up__1221->GetXaxis()->SetLabelFont(42);
   PU_up__1221->GetXaxis()->SetLabelSize(0.035);
   PU_up__1221->GetXaxis()->SetTitleSize(0.035);
   PU_up__1221->GetXaxis()->SetTitleFont(42);
   PU_up__1221->GetYaxis()->SetTitle("Events / bin");
   PU_up__1221->GetYaxis()->SetLabelFont(42);
   PU_up__1221->GetYaxis()->SetLabelSize(0.035);
   PU_up__1221->GetYaxis()->SetTitleSize(0.035);
   PU_up__1221->GetYaxis()->SetTitleOffset(0);
   PU_up__1221->GetYaxis()->SetTitleFont(42);
   PU_up__1221->GetZaxis()->SetLabelFont(42);
   PU_up__1221->GetZaxis()->SetLabelSize(0.035);
   PU_up__1221->GetZaxis()->SetTitleSize(0.035);
   PU_up__1221->GetZaxis()->SetTitleFont(42);
   PU_up__1221->Draw("E0 same");
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
   
   TH1D *frameR2__1222 = new TH1D("frameR2__1222","",1,0,2000);
   frameR2__1222->SetMinimum(0.5);
   frameR2__1222->SetMaximum(1.5);
   frameR2__1222->SetStats(0);
   frameR2__1222->SetLineStyle(0);
   frameR2__1222->SetMarkerStyle(20);
   frameR2__1222->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1222->GetXaxis()->SetRange(1,1);
   frameR2__1222->GetXaxis()->SetLabelFont(43);
   frameR2__1222->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1222->GetXaxis()->SetLabelSize(16);
   frameR2__1222->GetXaxis()->SetTitleSize(24);
   frameR2__1222->GetXaxis()->SetTitleOffset(5);
   frameR2__1222->GetXaxis()->SetTitleFont(43);
   frameR2__1222->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1222->GetYaxis()->SetNdivisions(205);
   frameR2__1222->GetYaxis()->SetLabelFont(43);
   frameR2__1222->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1222->GetYaxis()->SetLabelSize(20);
   frameR2__1222->GetYaxis()->SetTitleSize(20);
   frameR2__1222->GetYaxis()->SetTitleOffset(2);
   frameR2__1222->GetYaxis()->SetTitleFont(43);
   frameR2__1222->GetZaxis()->SetLabelFont(42);
   frameR2__1222->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1222->GetZaxis()->SetLabelSize(0.05);
   frameR2__1222->GetZaxis()->SetTitleSize(0.06);
   frameR2__1222->GetZaxis()->SetTitleFont(42);
   frameR2__1222->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis951[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1223 = new TH1F("nominal__1223","",32, xAxis951);
   nominal__1223->SetBinContent(2,1.164293);
   nominal__1223->SetBinContent(4,1.21135);
   nominal__1223->SetBinContent(6,1.117884);
   nominal__1223->SetBinContent(9,0.8106778);
   nominal__1223->SetBinContent(11,1.026085);
   nominal__1223->SetBinContent(12,1.103932);
   nominal__1223->SetBinContent(14,1.098593);
   nominal__1223->SetBinContent(15,0.9584086);
   nominal__1223->SetBinContent(16,1.069875);
   nominal__1223->SetBinContent(17,1.117884);
   nominal__1223->SetBinContent(18,0.9962677);
   nominal__1223->SetBinContent(20,1.103882);
   nominal__1223->SetBinContent(21,0.8692289);
   nominal__1223->SetBinContent(22,0.958222);
   nominal__1223->SetBinContent(23,1.047206);
   nominal__1223->SetBinContent(24,1.009815);
   nominal__1223->SetBinContent(25,0.9853957);
   nominal__1223->SetBinContent(26,0.9996234);
   nominal__1223->SetBinContent(27,1.004089);
   nominal__1223->SetBinContent(28,1.003155);
   nominal__1223->SetBinContent(29,1.003081);
   nominal__1223->SetBinContent(30,0.9901589);
   nominal__1223->SetBinContent(31,0.9914778);
   nominal__1223->SetBinContent(32,0.9983705);
   nominal__1223->SetBinError(2,1.646559);
   nominal__1223->SetBinError(4,1.713107);
   nominal__1223->SetBinError(6,1.580926);
   nominal__1223->SetBinError(9,0.7052662);
   nominal__1223->SetBinError(11,1.451104);
   nominal__1223->SetBinError(12,1.103965);
   nominal__1223->SetBinError(14,0.9049212);
   nominal__1223->SetBinError(15,1.355394);
   nominal__1223->SetBinError(16,1.513031);
   nominal__1223->SetBinError(17,1.580926);
   nominal__1223->SetBinError(18,0.6454056);
   nominal__1223->SetBinError(20,0.9061278);
   nominal__1223->SetBinError(21,0.7309601);
   nominal__1223->SetBinError(22,0.3877308);
   nominal__1223->SetBinError(23,0.2717975);
   nominal__1223->SetBinError(24,0.1512186);
   nominal__1223->SetBinError(25,0.08218211);
   nominal__1223->SetBinError(26,0.03401365);
   nominal__1223->SetBinError(27,0.016124);
   nominal__1223->SetBinError(28,0.01010427);
   nominal__1223->SetBinError(29,0.01133506);
   nominal__1223->SetBinError(30,0.019945);
   nominal__1223->SetBinError(31,0.04680115);
   nominal__1223->SetBinError(32,0.06637331);
   nominal__1223->SetMinimum(5e-05);
   nominal__1223->SetMaximum(3.298409);
   nominal__1223->SetEntries(28.30726);
   nominal__1223->SetFillColor(38);
   nominal__1223->SetLineColor(38);
   nominal__1223->SetMarkerColor(38);
   nominal__1223->SetMarkerStyle(21);
   nominal__1223->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1223->GetXaxis()->SetRange(1,32);
   nominal__1223->GetXaxis()->SetLabelFont(42);
   nominal__1223->GetXaxis()->SetLabelSize(0.035);
   nominal__1223->GetXaxis()->SetTitleSize(0.035);
   nominal__1223->GetXaxis()->SetTitleFont(42);
   nominal__1223->GetYaxis()->SetTitle("Tracks");
   nominal__1223->GetYaxis()->SetLabelFont(42);
   nominal__1223->GetYaxis()->SetLabelSize(0.05);
   nominal__1223->GetYaxis()->SetTitleSize(0.07);
   nominal__1223->GetYaxis()->SetTitleOffset(0);
   nominal__1223->GetYaxis()->SetTitleFont(42);
   nominal__1223->GetZaxis()->SetLabelFont(42);
   nominal__1223->GetZaxis()->SetLabelSize(0.035);
   nominal__1223->GetZaxis()->SetTitleSize(0.035);
   nominal__1223->GetZaxis()->SetTitleFont(42);
   nominal__1223->Draw("E0 same");
   Double_t xAxis952[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1224 = new TH1F("nominal__1224","",32, xAxis952);
   nominal__1224->SetBinContent(2,0.8724635);
   nominal__1224->SetBinContent(4,0.8327799);
   nominal__1224->SetBinContent(6,0.9051676);
   nominal__1224->SetBinContent(9,1.114709);
   nominal__1224->SetBinContent(11,1.009335);
   nominal__1224->SetBinContent(12,0.9122803);
   nominal__1224->SetBinContent(14,0.9151044);
   nominal__1224->SetBinContent(15,1.116385);
   nominal__1224->SetBinContent(16,0.9426885);
   nominal__1224->SetBinContent(17,0.9051676);
   nominal__1224->SetBinContent(18,0.9779865);
   nominal__1224->SetBinContent(20,0.9091498);
   nominal__1224->SetBinContent(21,1.123606);
   nominal__1224->SetBinContent(22,1.025972);
   nominal__1224->SetBinContent(23,0.9691734);
   nominal__1224->SetBinContent(24,0.9828354);
   nominal__1224->SetBinContent(25,1.018076);
   nominal__1224->SetBinContent(26,1.00203);
   nominal__1224->SetBinContent(27,0.9942966);
   nominal__1224->SetBinContent(28,0.9965003);
   nominal__1224->SetBinContent(29,0.9975086);
   nominal__1224->SetBinContent(30,1.008383);
   nominal__1224->SetBinContent(31,1.009459);
   nominal__1224->SetBinContent(32,1.007046);
   nominal__1224->SetBinError(2,1.23385);
   nominal__1224->SetBinError(4,1.177729);
   nominal__1224->SetBinError(6,1.2801);
   nominal__1224->SetBinError(9,0.9299495);
   nominal__1224->SetBinError(11,1.427415);
   nominal__1224->SetBinError(12,0.912302);
   nominal__1224->SetBinError(14,0.7487148);
   nominal__1224->SetBinError(15,1.578807);
   nominal__1224->SetBinError(16,1.333163);
   nominal__1224->SetBinError(17,1.2801);
   nominal__1224->SetBinError(18,0.6323195);
   nominal__1224->SetBinError(20,0.7432152);
   nominal__1224->SetBinError(21,0.9301201);
   nominal__1224->SetBinError(22,0.405373);
   nominal__1224->SetBinError(23,0.2507886);
   nominal__1224->SetBinError(24,0.1448208);
   nominal__1224->SetBinError(25,0.08386091);
   nominal__1224->SetBinError(26,0.03369083);
   nominal__1224->SetBinError(27,0.01575512);
   nominal__1224->SetBinError(28,0.009910101);
   nominal__1224->SetBinError(29,0.01113467);
   nominal__1224->SetBinError(30,0.02003251);
   nominal__1224->SetBinError(31,0.04704545);
   nominal__1224->SetBinError(32,0.06671574);
   nominal__1224->SetMinimum(5e-05);
   nominal__1224->SetMaximum(3.298409);
   nominal__1224->SetEntries(32.94182);
   nominal__1224->SetFillColor(46);
   nominal__1224->SetLineColor(46);
   nominal__1224->SetMarkerColor(46);
   nominal__1224->SetMarkerStyle(21);
   nominal__1224->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1224->GetXaxis()->SetRange(1,32);
   nominal__1224->GetXaxis()->SetLabelFont(42);
   nominal__1224->GetXaxis()->SetLabelSize(0.035);
   nominal__1224->GetXaxis()->SetTitleSize(0.035);
   nominal__1224->GetXaxis()->SetTitleFont(42);
   nominal__1224->GetYaxis()->SetTitle("Tracks");
   nominal__1224->GetYaxis()->SetLabelFont(42);
   nominal__1224->GetYaxis()->SetLabelSize(0.05);
   nominal__1224->GetYaxis()->SetTitleSize(0.07);
   nominal__1224->GetYaxis()->SetTitleOffset(0);
   nominal__1224->GetYaxis()->SetTitleFont(42);
   nominal__1224->GetZaxis()->SetLabelFont(42);
   nominal__1224->GetZaxis()->SetLabelSize(0.035);
   nominal__1224->GetZaxis()->SetTitleSize(0.035);
   nominal__1224->GetZaxis()->SetTitleFont(42);
   nominal__1224->Draw("E0 same");
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
