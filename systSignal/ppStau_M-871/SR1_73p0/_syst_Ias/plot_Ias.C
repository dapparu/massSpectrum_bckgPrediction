void plot_Ias()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:07 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.336977,2428.571,2.852504);
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
   Double_t xAxis897[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1153 = new TH1F("nominal__1153","",32, xAxis897);
   nominal__1153->SetBinContent(4,3.81416e-05);
   nominal__1153->SetBinContent(6,3.640635e-05);
   nominal__1153->SetBinContent(9,7.50132e-05);
   nominal__1153->SetBinContent(12,3.875062e-05);
   nominal__1153->SetBinContent(13,0.0001160943);
   nominal__1153->SetBinContent(14,0.0001246904);
   nominal__1153->SetBinContent(15,0.0001486052);
   nominal__1153->SetBinContent(16,7.575763e-05);
   nominal__1153->SetBinContent(17,7.503619e-05);
   nominal__1153->SetBinContent(18,0.0001129955);
   nominal__1153->SetBinContent(19,7.272686e-05);
   nominal__1153->SetBinContent(20,0.0003719724);
   nominal__1153->SetBinContent(21,0.0003015966);
   nominal__1153->SetBinContent(22,0.0009290634);
   nominal__1153->SetBinContent(23,0.004710929);
   nominal__1153->SetBinContent(24,0.01662294);
   nominal__1153->SetBinContent(25,0.061173);
   nominal__1153->SetBinContent(26,0.287656);
   nominal__1153->SetBinContent(27,0.6554763);
   nominal__1153->SetBinContent(28,0.4825398);
   nominal__1153->SetBinContent(29,0.125146);
   nominal__1153->SetBinContent(30,0.02389516);
   nominal__1153->SetBinContent(31,0.00418203);
   nominal__1153->SetBinContent(32,0.00216105);
   nominal__1153->SetBinError(4,3.81416e-05);
   nominal__1153->SetBinError(6,3.640634e-05);
   nominal__1153->SetBinError(9,5.304995e-05);
   nominal__1153->SetBinError(12,3.875062e-05);
   nominal__1153->SetBinError(13,6.703851e-05);
   nominal__1153->SetBinError(14,7.263418e-05);
   nominal__1153->SetBinError(15,7.433219e-05);
   nominal__1153->SetBinError(16,5.356886e-05);
   nominal__1153->SetBinError(17,5.305876e-05);
   nominal__1153->SetBinError(18,6.525134e-05);
   nominal__1153->SetBinError(19,5.142694e-05);
   nominal__1153->SetBinError(20,0.0001180085);
   nominal__1153->SetBinError(21,0.00010683);
   nominal__1153->SetBinError(22,0.0001860388);
   nominal__1153->SetBinError(23,0.0004186011);
   nominal__1153->SetBinError(24,0.0007905417);
   nominal__1153->SetBinError(25,0.001512552);
   nominal__1153->SetBinError(26,0.003282171);
   nominal__1153->SetBinError(27,0.004952964);
   nominal__1153->SetBinError(28,0.004252134);
   nominal__1153->SetBinError(29,0.002170874);
   nominal__1153->SetBinError(30,0.0009478426);
   nominal__1153->SetBinError(31,0.0003993255);
   nominal__1153->SetBinError(32,0.0002089647);
   nominal__1153->SetBinError(33,0.0001964304);
   nominal__1153->SetMinimum(5e-05);
   nominal__1153->SetMaximum(655.4763);
   nominal__1153->SetEntries(44645);
   nominal__1153->SetFillColor(1);
   nominal__1153->SetMarkerStyle(20);
   nominal__1153->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1153->GetXaxis()->SetRange(1,32);
   nominal__1153->GetXaxis()->SetLabelFont(42);
   nominal__1153->GetXaxis()->SetLabelSize(0.035);
   nominal__1153->GetXaxis()->SetTitleSize(0.035);
   nominal__1153->GetXaxis()->SetTitleFont(42);
   nominal__1153->GetYaxis()->SetTitle("Tracks");
   nominal__1153->GetYaxis()->SetLabelFont(42);
   nominal__1153->GetYaxis()->SetLabelSize(0.05);
   nominal__1153->GetYaxis()->SetTitleSize(0.07);
   nominal__1153->GetYaxis()->SetTitleOffset(0);
   nominal__1153->GetYaxis()->SetTitleFont(42);
   nominal__1153->GetZaxis()->SetLabelFont(42);
   nominal__1153->GetZaxis()->SetLabelSize(0.035);
   nominal__1153->GetZaxis()->SetTitleSize(0.035);
   nominal__1153->GetZaxis()->SetTitleFont(42);
   nominal__1153->Draw("");
   Double_t xAxis898[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Ias_down__1154 = new TH1F("Ias_down__1154","",32, xAxis898);
   Ias_down__1154->SetBinContent(4,3.81416e-05);
   Ias_down__1154->SetBinContent(6,3.640635e-05);
   Ias_down__1154->SetBinContent(9,7.50132e-05);
   Ias_down__1154->SetBinContent(12,3.875062e-05);
   Ias_down__1154->SetBinContent(13,0.0001160943);
   Ias_down__1154->SetBinContent(14,0.0001246904);
   Ias_down__1154->SetBinContent(15,0.0001486052);
   Ias_down__1154->SetBinContent(16,7.575763e-05);
   Ias_down__1154->SetBinContent(17,7.503619e-05);
   Ias_down__1154->SetBinContent(18,0.0001129955);
   Ias_down__1154->SetBinContent(19,7.272686e-05);
   Ias_down__1154->SetBinContent(20,0.0003343634);
   Ias_down__1154->SetBinContent(21,0.0003015966);
   Ias_down__1154->SetBinContent(22,0.0009290634);
   Ias_down__1154->SetBinContent(23,0.004491641);
   Ias_down__1154->SetBinContent(24,0.01631834);
   Ias_down__1154->SetBinContent(25,0.06076199);
   Ias_down__1154->SetBinContent(26,0.2868877);
   Ias_down__1154->SetBinContent(27,0.654951);
   Ias_down__1154->SetBinContent(28,0.4823872);
   Ias_down__1154->SetBinContent(29,0.1251028);
   Ias_down__1154->SetBinContent(30,0.02378404);
   Ias_down__1154->SetBinContent(31,0.00410767);
   Ias_down__1154->SetBinContent(32,0.002124444);
   Ias_down__1154->SetBinError(4,3.81416e-05);
   Ias_down__1154->SetBinError(6,3.640634e-05);
   Ias_down__1154->SetBinError(9,5.304995e-05);
   Ias_down__1154->SetBinError(12,3.875062e-05);
   Ias_down__1154->SetBinError(13,6.703851e-05);
   Ias_down__1154->SetBinError(14,7.263418e-05);
   Ias_down__1154->SetBinError(15,7.433219e-05);
   Ias_down__1154->SetBinError(16,5.356886e-05);
   Ias_down__1154->SetBinError(17,5.305876e-05);
   Ias_down__1154->SetBinError(18,6.525134e-05);
   Ias_down__1154->SetBinError(19,5.142694e-05);
   Ias_down__1154->SetBinError(20,0.0001118552);
   Ias_down__1154->SetBinError(21,0.00010683);
   Ias_down__1154->SetBinError(22,0.0001860388);
   Ias_down__1154->SetBinError(23,0.0004088945);
   Ias_down__1154->SetBinError(24,0.0007831662);
   Ias_down__1154->SetBinError(25,0.001507443);
   Ias_down__1154->SetBinError(26,0.003277873);
   Ias_down__1154->SetBinError(27,0.004950971);
   Ias_down__1154->SetBinError(28,0.004251446);
   Ias_down__1154->SetBinError(29,0.002170445);
   Ias_down__1154->SetBinError(30,0.0009456668);
   Ias_down__1154->SetBinError(31,0.0003958484);
   Ias_down__1154->SetBinError(32,0.0002057334);
   Ias_down__1154->SetBinError(33,0.0001964304);
   Ias_down__1154->SetEntries(44573);
   Ias_down__1154->SetFillColor(38);
   Ias_down__1154->SetLineColor(38);
   Ias_down__1154->SetMarkerColor(38);
   Ias_down__1154->SetMarkerStyle(21);
   Ias_down__1154->GetXaxis()->SetTitle("Mass [GeV]");
   Ias_down__1154->GetXaxis()->SetRange(1,400);
   Ias_down__1154->GetXaxis()->SetLabelFont(42);
   Ias_down__1154->GetXaxis()->SetLabelSize(0.035);
   Ias_down__1154->GetXaxis()->SetTitleSize(0.035);
   Ias_down__1154->GetXaxis()->SetTitleFont(42);
   Ias_down__1154->GetYaxis()->SetTitle("Events / bin");
   Ias_down__1154->GetYaxis()->SetLabelFont(42);
   Ias_down__1154->GetYaxis()->SetLabelSize(0.035);
   Ias_down__1154->GetYaxis()->SetTitleSize(0.035);
   Ias_down__1154->GetYaxis()->SetTitleOffset(0);
   Ias_down__1154->GetYaxis()->SetTitleFont(42);
   Ias_down__1154->GetZaxis()->SetLabelFont(42);
   Ias_down__1154->GetZaxis()->SetLabelSize(0.035);
   Ias_down__1154->GetZaxis()->SetTitleSize(0.035);
   Ias_down__1154->GetZaxis()->SetTitleFont(42);
   Ias_down__1154->Draw("same");
   Double_t xAxis899[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Ias_up__1155 = new TH1F("Ias_up__1155","",32, xAxis899);
   Ias_up__1155->SetBinContent(4,3.81416e-05);
   Ias_up__1155->SetBinContent(6,3.640635e-05);
   Ias_up__1155->SetBinContent(9,7.50132e-05);
   Ias_up__1155->SetBinContent(12,3.875062e-05);
   Ias_up__1155->SetBinContent(13,0.0001160943);
   Ias_up__1155->SetBinContent(14,0.0001246904);
   Ias_up__1155->SetBinContent(15,0.0001486052);
   Ias_up__1155->SetBinContent(16,7.575763e-05);
   Ias_up__1155->SetBinContent(17,7.503619e-05);
   Ias_up__1155->SetBinContent(18,0.0001129955);
   Ias_up__1155->SetBinContent(19,7.272686e-05);
   Ias_up__1155->SetBinContent(20,0.0003719724);
   Ias_up__1155->SetBinContent(21,0.0003015966);
   Ias_up__1155->SetBinContent(22,0.001008813);
   Ias_up__1155->SetBinContent(23,0.004748225);
   Ias_up__1155->SetBinContent(24,0.01695061);
   Ias_up__1155->SetBinContent(25,0.06150794);
   Ias_up__1155->SetBinContent(26,0.2883292);
   Ias_up__1155->SetBinContent(27,0.6556628);
   Ias_up__1155->SetBinContent(28,0.4826182);
   Ias_up__1155->SetBinContent(29,0.1251829);
   Ias_up__1155->SetBinContent(30,0.02389516);
   Ias_up__1155->SetBinContent(31,0.00421739);
   Ias_up__1155->SetBinContent(32,0.00216105);
   Ias_up__1155->SetBinError(4,3.81416e-05);
   Ias_up__1155->SetBinError(6,3.640634e-05);
   Ias_up__1155->SetBinError(9,5.304995e-05);
   Ias_up__1155->SetBinError(12,3.875062e-05);
   Ias_up__1155->SetBinError(13,6.703851e-05);
   Ias_up__1155->SetBinError(14,7.263418e-05);
   Ias_up__1155->SetBinError(15,7.433219e-05);
   Ias_up__1155->SetBinError(16,5.356886e-05);
   Ias_up__1155->SetBinError(17,5.305876e-05);
   Ias_up__1155->SetBinError(18,6.525134e-05);
   Ias_up__1155->SetBinError(19,5.142694e-05);
   Ias_up__1155->SetBinError(20,0.0001180085);
   Ias_up__1155->SetBinError(21,0.00010683);
   Ias_up__1155->SetBinError(22,0.0001944525);
   Ias_up__1155->SetBinError(23,0.0004202592);
   Ias_up__1155->SetBinError(24,0.0007980652);
   Ias_up__1155->SetBinError(25,0.001516699);
   Ias_up__1155->SetBinError(26,0.003286006);
   Ias_up__1155->SetBinError(27,0.004953666);
   Ias_up__1155->SetBinError(28,0.004252496);
   Ias_up__1155->SetBinError(29,0.002171187);
   Ias_up__1155->SetBinError(30,0.0009478426);
   Ias_up__1155->SetBinError(31,0.0004008879);
   Ias_up__1155->SetBinError(32,0.0002089647);
   Ias_up__1155->SetBinError(33,0.0001964304);
   Ias_up__1155->SetEntries(44693);
   Ias_up__1155->SetFillColor(46);
   Ias_up__1155->SetLineColor(46);
   Ias_up__1155->SetMarkerColor(46);
   Ias_up__1155->SetMarkerStyle(21);
   Ias_up__1155->GetXaxis()->SetTitle("Mass [GeV]");
   Ias_up__1155->GetXaxis()->SetRange(1,400);
   Ias_up__1155->GetXaxis()->SetLabelFont(42);
   Ias_up__1155->GetXaxis()->SetLabelSize(0.035);
   Ias_up__1155->GetXaxis()->SetTitleSize(0.035);
   Ias_up__1155->GetXaxis()->SetTitleFont(42);
   Ias_up__1155->GetYaxis()->SetTitle("Events / bin");
   Ias_up__1155->GetYaxis()->SetLabelFont(42);
   Ias_up__1155->GetYaxis()->SetLabelSize(0.035);
   Ias_up__1155->GetYaxis()->SetTitleSize(0.035);
   Ias_up__1155->GetYaxis()->SetTitleOffset(0);
   Ias_up__1155->GetYaxis()->SetTitleFont(42);
   Ias_up__1155->GetZaxis()->SetLabelFont(42);
   Ias_up__1155->GetZaxis()->SetLabelSize(0.035);
   Ias_up__1155->GetZaxis()->SetTitleSize(0.035);
   Ias_up__1155->GetZaxis()->SetTitleFont(42);
   Ias_up__1155->Draw("same");
   TLine *line = new TLine(1730,0,1730,655.4763);
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
   entry=leg->AddEntry("Ias_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Ias_up","Up","PE1");
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
   
   TH1D *frameR2__1156 = new TH1D("frameR2__1156","",1,0,2000);
   frameR2__1156->SetMinimum(0.5);
   frameR2__1156->SetMaximum(1.5);
   frameR2__1156->SetStats(0);
   frameR2__1156->SetLineStyle(0);
   frameR2__1156->SetMarkerStyle(20);
   frameR2__1156->GetXaxis()->SetRange(1,1);
   frameR2__1156->GetXaxis()->SetLabelFont(43);
   frameR2__1156->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1156->GetXaxis()->SetLabelSize(16);
   frameR2__1156->GetXaxis()->SetTitleSize(24);
   frameR2__1156->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1156->GetXaxis()->SetTitleFont(43);
   frameR2__1156->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1156->GetYaxis()->SetNdivisions(205);
   frameR2__1156->GetYaxis()->SetLabelFont(43);
   frameR2__1156->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1156->GetYaxis()->SetLabelSize(20);
   frameR2__1156->GetYaxis()->SetTitleSize(20);
   frameR2__1156->GetYaxis()->SetTitleOffset(2);
   frameR2__1156->GetYaxis()->SetTitleFont(43);
   frameR2__1156->GetZaxis()->SetLabelFont(42);
   frameR2__1156->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1156->GetZaxis()->SetLabelSize(0.05);
   frameR2__1156->GetZaxis()->SetTitleSize(0.06);
   frameR2__1156->GetZaxis()->SetTitleFont(42);
   frameR2__1156->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis900[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Ias_down__1157 = new TH1F("Ias_down__1157","",32, xAxis900);
   Ias_down__1157->SetBinContent(0,1.001613);
   Ias_down__1157->SetBinContent(1,1.001613);
   Ias_down__1157->SetBinContent(2,1.001613);
   Ias_down__1157->SetBinContent(3,1.001613);
   Ias_down__1157->SetBinContent(4,1.001613);
   Ias_down__1157->SetBinContent(5,1.001613);
   Ias_down__1157->SetBinContent(6,1.001613);
   Ias_down__1157->SetBinContent(7,1.001613);
   Ias_down__1157->SetBinContent(8,1.001613);
   Ias_down__1157->SetBinContent(9,1.001613);
   Ias_down__1157->SetBinContent(10,1.001614);
   Ias_down__1157->SetBinContent(11,1.001614);
   Ias_down__1157->SetBinContent(12,1.001614);
   Ias_down__1157->SetBinContent(13,1.001614);
   Ias_down__1157->SetBinContent(14,1.001614);
   Ias_down__1157->SetBinContent(15,1.001614);
   Ias_down__1157->SetBinContent(16,1.001614);
   Ias_down__1157->SetBinContent(17,1.001614);
   Ias_down__1157->SetBinContent(18,1.001614);
   Ias_down__1157->SetBinContent(19,1.001614);
   Ias_down__1157->SetBinContent(20,1.001614);
   Ias_down__1157->SetBinContent(21,1.001592);
   Ias_down__1157->SetBinContent(22,1.001592);
   Ias_down__1157->SetBinContent(23,1.001593);
   Ias_down__1157->SetBinContent(24,1.001465);
   Ias_down__1157->SetBinContent(25,1.001294);
   Ias_down__1157->SetBinContent(26,1.001083);
   Ias_down__1157->SetBinContent(27,1.00073);
   Ias_down__1157->SetBinContent(28,1.000655);
   Ias_down__1157->SetBinContent(29,1.00171);
   Ias_down__1157->SetBinContent(30,1.007399);
   Ias_down__1157->SetBinContent(31,1.017805);
   Ias_down__1157->SetBinContent(32,1.017231);
   Ias_down__1157->SetBinError(0,0.006719851);
   Ias_down__1157->SetBinError(1,0.006719851);
   Ias_down__1157->SetBinError(2,0.006719851);
   Ias_down__1157->SetBinError(3,0.006719851);
   Ias_down__1157->SetBinError(4,0.006719851);
   Ias_down__1157->SetBinError(5,0.006719927);
   Ias_down__1157->SetBinError(6,0.006719927);
   Ias_down__1157->SetBinError(7,0.006720003);
   Ias_down__1157->SetBinError(8,0.006720003);
   Ias_down__1157->SetBinError(9,0.006720003);
   Ias_down__1157->SetBinError(10,0.006720155);
   Ias_down__1157->SetBinError(11,0.006720155);
   Ias_down__1157->SetBinError(12,0.006720155);
   Ias_down__1157->SetBinError(13,0.00672023);
   Ias_down__1157->SetBinError(14,0.006720458);
   Ias_down__1157->SetBinError(15,0.006720678);
   Ias_down__1157->SetBinError(16,0.006720981);
   Ias_down__1157->SetBinError(17,0.006721133);
   Ias_down__1157->SetBinError(18,0.006721285);
   Ias_down__1157->SetBinError(19,0.006721512);
   Ias_down__1157->SetBinError(20,0.006721664);
   Ias_down__1157->SetBinError(21,0.006722228);
   Ias_down__1157->SetBinError(22,0.006722834);
   Ias_down__1157->SetBinError(23,0.006724729);
   Ias_down__1157->SetBinError(24,0.006733258);
   Ias_down__1157->SetBinError(25,0.006765691);
   Ias_down__1157->SetBinError(26,0.006893758);
   Ias_down__1157->SetBinError(27,0.007618908);
   Ias_down__1157->SetBinError(28,0.01085433);
   Ias_down__1157->SetBinError(29,0.02206561);
   Ias_down__1157->SetBinError(30,0.05040236);
   Ias_down__1157->SetBinError(31,0.1120754);
   Ias_down__1157->SetBinError(32,0.1917677);
   Ias_down__1157->SetEntries(33);
   Ias_down__1157->SetFillColor(38);
   Ias_down__1157->SetLineColor(38);
   Ias_down__1157->SetMarkerColor(38);
   Ias_down__1157->SetMarkerStyle(21);
   Ias_down__1157->GetXaxis()->SetTitle("Mass [GeV]");
   Ias_down__1157->GetXaxis()->SetRange(1,400);
   Ias_down__1157->GetXaxis()->SetLabelFont(42);
   Ias_down__1157->GetXaxis()->SetLabelSize(0.035);
   Ias_down__1157->GetXaxis()->SetTitleSize(0.035);
   Ias_down__1157->GetXaxis()->SetTitleFont(42);
   Ias_down__1157->GetYaxis()->SetTitle("Events / bin");
   Ias_down__1157->GetYaxis()->SetLabelFont(42);
   Ias_down__1157->GetYaxis()->SetLabelSize(0.035);
   Ias_down__1157->GetYaxis()->SetTitleSize(0.035);
   Ias_down__1157->GetYaxis()->SetTitleOffset(0);
   Ias_down__1157->GetYaxis()->SetTitleFont(42);
   Ias_down__1157->GetZaxis()->SetLabelFont(42);
   Ias_down__1157->GetZaxis()->SetLabelSize(0.035);
   Ias_down__1157->GetZaxis()->SetTitleSize(0.035);
   Ias_down__1157->GetZaxis()->SetTitleFont(42);
   Ias_down__1157->Draw("E0 same");
   Double_t xAxis901[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Ias_up__1158 = new TH1F("Ias_up__1158","",32, xAxis901);
   Ias_up__1158->SetBinContent(0,0.9989268);
   Ias_up__1158->SetBinContent(1,0.9989268);
   Ias_up__1158->SetBinContent(2,0.9989268);
   Ias_up__1158->SetBinContent(3,0.9989268);
   Ias_up__1158->SetBinContent(4,0.9989268);
   Ias_up__1158->SetBinContent(5,0.9989268);
   Ias_up__1158->SetBinContent(6,0.9989268);
   Ias_up__1158->SetBinContent(7,0.9989267);
   Ias_up__1158->SetBinContent(8,0.9989267);
   Ias_up__1158->SetBinContent(9,0.9989267);
   Ias_up__1158->SetBinContent(10,0.9989267);
   Ias_up__1158->SetBinContent(11,0.9989267);
   Ias_up__1158->SetBinContent(12,0.9989267);
   Ias_up__1158->SetBinContent(13,0.9989266);
   Ias_up__1158->SetBinContent(14,0.9989266);
   Ias_up__1158->SetBinContent(15,0.9989265);
   Ias_up__1158->SetBinContent(16,0.9989264);
   Ias_up__1158->SetBinContent(17,0.9989263);
   Ias_up__1158->SetBinContent(18,0.9989263);
   Ias_up__1158->SetBinContent(19,0.9989262);
   Ias_up__1158->SetBinContent(20,0.9989262);
   Ias_up__1158->SetBinContent(21,0.9989259);
   Ias_up__1158->SetBinContent(22,0.9989257);
   Ias_up__1158->SetBinContent(23,0.998973);
   Ias_up__1158->SetBinContent(24,0.9989925);
   Ias_up__1158->SetBinContent(25,0.9991815);
   Ias_up__1158->SetBinContent(26,0.9993614);
   Ias_up__1158->SetBinContent(27,0.9997394);
   Ias_up__1158->SetBinContent(28,0.9997638);
   Ias_up__1158->SetBinContent(29,0.9995352);
   Ias_up__1158->SetBinContent(30,0.998832);
   Ias_up__1158->SetBinContent(31,0.9944564);
   Ias_up__1158->SetBinContent(32,1);
   Ias_up__1158->SetBinError(0,0.006697318);
   Ias_up__1158->SetBinError(1,0.006697318);
   Ias_up__1158->SetBinError(2,0.006697318);
   Ias_up__1158->SetBinError(3,0.006697318);
   Ias_up__1158->SetBinError(4,0.006697318);
   Ias_up__1158->SetBinError(5,0.006697394);
   Ias_up__1158->SetBinError(6,0.006697394);
   Ias_up__1158->SetBinError(7,0.006697469);
   Ias_up__1158->SetBinError(8,0.006697469);
   Ias_up__1158->SetBinError(9,0.006697469);
   Ias_up__1158->SetBinError(10,0.006697619);
   Ias_up__1158->SetBinError(11,0.006697619);
   Ias_up__1158->SetBinError(12,0.006697619);
   Ias_up__1158->SetBinError(13,0.006697694);
   Ias_up__1158->SetBinError(14,0.006697919);
   Ias_up__1158->SetBinError(15,0.006698136);
   Ias_up__1158->SetBinError(16,0.006698437);
   Ias_up__1158->SetBinError(17,0.006698587);
   Ias_up__1158->SetBinError(18,0.006698737);
   Ias_up__1158->SetBinError(19,0.006698963);
   Ias_up__1158->SetBinError(20,0.006699113);
   Ias_up__1158->SetBinError(21,0.00669986);
   Ias_up__1158->SetBinError(22,0.006700459);
   Ias_up__1158->SetBinError(23,0.006702731);
   Ias_up__1158->SetBinError(24,0.006712474);
   Ias_up__1158->SetBinError(25,0.006747847);
   Ias_up__1158->SetBinError(26,0.006878929);
   Ias_up__1158->SetBinError(27,0.007609499);
   Ias_up__1158->SetBinError(28,0.01084228);
   Ias_up__1158->SetBinError(29,0.02200574);
   Ias_up__1158->SetBinError(30,0.04986445);
   Ias_up__1158->SetBinError(31,0.1088438);
   Ias_up__1158->SetBinError(32,0.1876814);
   Ias_up__1158->SetEntries(33);
   Ias_up__1158->SetFillColor(46);
   Ias_up__1158->SetLineColor(46);
   Ias_up__1158->SetMarkerColor(46);
   Ias_up__1158->SetMarkerStyle(21);
   Ias_up__1158->GetXaxis()->SetTitle("Mass [GeV]");
   Ias_up__1158->GetXaxis()->SetRange(1,400);
   Ias_up__1158->GetXaxis()->SetLabelFont(42);
   Ias_up__1158->GetXaxis()->SetLabelSize(0.035);
   Ias_up__1158->GetXaxis()->SetTitleSize(0.035);
   Ias_up__1158->GetXaxis()->SetTitleFont(42);
   Ias_up__1158->GetYaxis()->SetTitle("Events / bin");
   Ias_up__1158->GetYaxis()->SetLabelFont(42);
   Ias_up__1158->GetYaxis()->SetLabelSize(0.035);
   Ias_up__1158->GetYaxis()->SetTitleSize(0.035);
   Ias_up__1158->GetYaxis()->SetTitleOffset(0);
   Ias_up__1158->GetYaxis()->SetTitleFont(42);
   Ias_up__1158->GetZaxis()->SetLabelFont(42);
   Ias_up__1158->GetZaxis()->SetLabelSize(0.035);
   Ias_up__1158->GetZaxis()->SetTitleSize(0.035);
   Ias_up__1158->GetZaxis()->SetTitleFont(42);
   Ias_up__1158->Draw("E0 same");
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
   
   TH1D *frameR2__1159 = new TH1D("frameR2__1159","",1,0,2000);
   frameR2__1159->SetMinimum(0.5);
   frameR2__1159->SetMaximum(1.5);
   frameR2__1159->SetStats(0);
   frameR2__1159->SetLineStyle(0);
   frameR2__1159->SetMarkerStyle(20);
   frameR2__1159->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1159->GetXaxis()->SetRange(1,1);
   frameR2__1159->GetXaxis()->SetLabelFont(43);
   frameR2__1159->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1159->GetXaxis()->SetLabelSize(16);
   frameR2__1159->GetXaxis()->SetTitleSize(24);
   frameR2__1159->GetXaxis()->SetTitleOffset(5);
   frameR2__1159->GetXaxis()->SetTitleFont(43);
   frameR2__1159->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1159->GetYaxis()->SetNdivisions(205);
   frameR2__1159->GetYaxis()->SetLabelFont(43);
   frameR2__1159->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1159->GetYaxis()->SetLabelSize(20);
   frameR2__1159->GetYaxis()->SetTitleSize(20);
   frameR2__1159->GetYaxis()->SetTitleOffset(2);
   frameR2__1159->GetYaxis()->SetTitleFont(43);
   frameR2__1159->GetZaxis()->SetLabelFont(42);
   frameR2__1159->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1159->GetZaxis()->SetLabelSize(0.05);
   frameR2__1159->GetZaxis()->SetTitleSize(0.06);
   frameR2__1159->GetZaxis()->SetTitleFont(42);
   frameR2__1159->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis902[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1160 = new TH1F("nominal__1160","",32, xAxis902);
   nominal__1160->SetBinContent(4,1);
   nominal__1160->SetBinContent(6,1);
   nominal__1160->SetBinContent(9,1);
   nominal__1160->SetBinContent(12,1);
   nominal__1160->SetBinContent(13,1);
   nominal__1160->SetBinContent(14,1);
   nominal__1160->SetBinContent(15,1);
   nominal__1160->SetBinContent(16,1);
   nominal__1160->SetBinContent(17,1);
   nominal__1160->SetBinContent(18,1);
   nominal__1160->SetBinContent(19,1);
   nominal__1160->SetBinContent(20,1.112479);
   nominal__1160->SetBinContent(21,1);
   nominal__1160->SetBinContent(22,1);
   nominal__1160->SetBinContent(23,1.048822);
   nominal__1160->SetBinContent(24,1.018666);
   nominal__1160->SetBinContent(25,1.006764);
   nominal__1160->SetBinContent(26,1.002678);
   nominal__1160->SetBinContent(27,1.000802);
   nominal__1160->SetBinContent(28,1.000316);
   nominal__1160->SetBinContent(29,1.000345);
   nominal__1160->SetBinContent(30,1.004672);
   nominal__1160->SetBinContent(31,1.018103);
   nominal__1160->SetBinContent(32,1.017231);
   nominal__1160->SetBinError(4,1.414214);
   nominal__1160->SetBinError(6,1.414214);
   nominal__1160->SetBinError(9,1.000143);
   nominal__1160->SetBinError(12,1.414214);
   nominal__1160->SetBinError(13,0.816636);
   nominal__1160->SetBinError(14,0.823802);
   nominal__1160->SetBinError(15,0.7073883);
   nominal__1160->SetBinError(16,1.000002);
   nominal__1160->SetBinError(17,1.000003);
   nominal__1160->SetBinError(18,0.8166635);
   nominal__1160->SetBinError(19,1.000025);
   nominal__1160->SetBinError(20,0.5128996);
   nominal__1160->SetBinError(21,0.5009353);
   nominal__1160->SetBinError(22,0.283187);
   nominal__1160->SetBinError(23,0.1334228);
   nominal__1160->SetBinError(24,0.06882612);
   nominal__1160->SetBinError(25,0.03526337);
   nominal__1160->SetBinError(26,0.01619051);
   nominal__1160->SetBinError(27,0.0106969);
   nominal__1160->SetBinError(28,0.01246693);
   nominal__1160->SetBinError(29,0.02454225);
   nominal__1160->SetBinError(30,0.056426);
   nominal__1160->SetBinError(31,0.1381186);
   nominal__1160->SetBinError(32,0.1392094);
   nominal__1160->SetMinimum(5e-05);
   nominal__1160->SetMaximum(655.4763);
   nominal__1160->SetEntries(44.5678);
   nominal__1160->SetFillColor(38);
   nominal__1160->SetLineColor(38);
   nominal__1160->SetMarkerColor(38);
   nominal__1160->SetMarkerStyle(21);
   nominal__1160->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1160->GetXaxis()->SetRange(1,32);
   nominal__1160->GetXaxis()->SetLabelFont(42);
   nominal__1160->GetXaxis()->SetLabelSize(0.035);
   nominal__1160->GetXaxis()->SetTitleSize(0.035);
   nominal__1160->GetXaxis()->SetTitleFont(42);
   nominal__1160->GetYaxis()->SetTitle("Tracks");
   nominal__1160->GetYaxis()->SetLabelFont(42);
   nominal__1160->GetYaxis()->SetLabelSize(0.05);
   nominal__1160->GetYaxis()->SetTitleSize(0.07);
   nominal__1160->GetYaxis()->SetTitleOffset(0);
   nominal__1160->GetYaxis()->SetTitleFont(42);
   nominal__1160->GetZaxis()->SetLabelFont(42);
   nominal__1160->GetZaxis()->SetLabelSize(0.035);
   nominal__1160->GetZaxis()->SetTitleSize(0.035);
   nominal__1160->GetZaxis()->SetTitleFont(42);
   nominal__1160->Draw("E0 same");
   Double_t xAxis903[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1161 = new TH1F("nominal__1161","",32, xAxis903);
   nominal__1161->SetBinContent(4,1);
   nominal__1161->SetBinContent(6,1);
   nominal__1161->SetBinContent(9,1);
   nominal__1161->SetBinContent(12,1);
   nominal__1161->SetBinContent(13,1);
   nominal__1161->SetBinContent(14,1);
   nominal__1161->SetBinContent(15,1);
   nominal__1161->SetBinContent(16,1);
   nominal__1161->SetBinContent(17,1);
   nominal__1161->SetBinContent(18,1);
   nominal__1161->SetBinContent(19,1);
   nominal__1161->SetBinContent(20,1);
   nominal__1161->SetBinContent(21,1);
   nominal__1161->SetBinContent(22,0.9209474);
   nominal__1161->SetBinContent(23,0.9921454);
   nominal__1161->SetBinContent(24,0.9806687);
   nominal__1161->SetBinContent(25,0.9945545);
   nominal__1161->SetBinContent(26,0.997665);
   nominal__1161->SetBinContent(27,0.9997156);
   nominal__1161->SetBinContent(28,0.9998375);
   nominal__1161->SetBinContent(29,0.9997053);
   nominal__1161->SetBinContent(30,1);
   nominal__1161->SetBinContent(31,0.9916158);
   nominal__1161->SetBinContent(32,1);
   nominal__1161->SetBinError(4,1.414214);
   nominal__1161->SetBinError(6,1.414214);
   nominal__1161->SetBinError(9,1.000143);
   nominal__1161->SetBinError(12,1.414214);
   nominal__1161->SetBinError(13,0.816636);
   nominal__1161->SetBinError(14,0.823802);
   nominal__1161->SetBinError(15,0.7073883);
   nominal__1161->SetBinError(16,1.000002);
   nominal__1161->SetBinError(17,1.000003);
   nominal__1161->SetBinError(18,0.8166635);
   nominal__1161->SetBinError(19,1.000025);
   nominal__1161->SetBinError(20,0.4486604);
   nominal__1161->SetBinError(21,0.5009353);
   nominal__1161->SetBinError(22,0.2559695);
   nominal__1161->SetBinError(23,0.1244319);
   nominal__1161->SetBinError(24,0.0656271);
   nominal__1161->SetBinError(25,0.03472992);
   nominal__1161->SetBinError(26,0.01608917);
   nominal__1161->SetBinError(27,0.0106824);
   nominal__1161->SetBinError(28,0.01245952);
   nominal__1161->SetBinError(29,0.02452291);
   nominal__1161->SetBinError(30,0.05609721);
   nominal__1161->SetBinError(31,0.1336042);
   nominal__1161->SetBinError(32,0.1367486);
   nominal__1161->SetMinimum(5e-05);
   nominal__1161->SetMaximum(655.4763);
   nominal__1161->SetEntries(43.54332);
   nominal__1161->SetFillColor(46);
   nominal__1161->SetLineColor(46);
   nominal__1161->SetMarkerColor(46);
   nominal__1161->SetMarkerStyle(21);
   nominal__1161->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1161->GetXaxis()->SetRange(1,32);
   nominal__1161->GetXaxis()->SetLabelFont(42);
   nominal__1161->GetXaxis()->SetLabelSize(0.035);
   nominal__1161->GetXaxis()->SetTitleSize(0.035);
   nominal__1161->GetXaxis()->SetTitleFont(42);
   nominal__1161->GetYaxis()->SetTitle("Tracks");
   nominal__1161->GetYaxis()->SetLabelFont(42);
   nominal__1161->GetYaxis()->SetLabelSize(0.05);
   nominal__1161->GetYaxis()->SetTitleSize(0.07);
   nominal__1161->GetYaxis()->SetTitleOffset(0);
   nominal__1161->GetYaxis()->SetTitleFont(42);
   nominal__1161->GetZaxis()->SetLabelFont(42);
   nominal__1161->GetZaxis()->SetLabelSize(0.035);
   nominal__1161->GetZaxis()->SetTitleSize(0.035);
   nominal__1161->GetZaxis()->SetTitleFont(42);
   nominal__1161->Draw("E0 same");
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
