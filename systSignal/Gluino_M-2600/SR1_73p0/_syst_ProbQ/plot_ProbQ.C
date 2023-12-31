void plot_ProbQ()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:02 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.329041,2428.571,1.273083);
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
   Double_t xAxis449[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__577 = new TH1F("nominal__577","",32, xAxis449);
   nominal__577->SetBinContent(20,4.716178e-05);
   nominal__577->SetBinContent(21,4.750374e-05);
   nominal__577->SetBinContent(26,4.856034e-05);
   nominal__577->SetBinContent(27,5.207413e-05);
   nominal__577->SetBinContent(28,0.0003452141);
   nominal__577->SetBinContent(29,0.0005223806);
   nominal__577->SetBinContent(30,0.003358804);
   nominal__577->SetBinContent(31,0.01355104);
   nominal__577->SetBinContent(32,0.1758217);
   nominal__577->SetBinError(20,4.716178e-05);
   nominal__577->SetBinError(21,4.750374e-05);
   nominal__577->SetBinError(26,4.856034e-05);
   nominal__577->SetBinError(27,5.207413e-05);
   nominal__577->SetBinError(28,0.0001306065);
   nominal__577->SetBinError(29,0.0001582745);
   nominal__577->SetBinError(30,0.0004023621);
   nominal__577->SetBinError(31,0.0008105387);
   nominal__577->SetBinError(32,0.001129738);
   nominal__577->SetBinError(33,0.002682408);
   nominal__577->SetMinimum(5e-05);
   nominal__577->SetMaximum(17.58217);
   nominal__577->SetEntries(4041);
   nominal__577->SetFillColor(1);
   nominal__577->SetMarkerStyle(20);
   nominal__577->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__577->GetXaxis()->SetRange(1,32);
   nominal__577->GetXaxis()->SetLabelFont(42);
   nominal__577->GetXaxis()->SetLabelSize(0.035);
   nominal__577->GetXaxis()->SetTitleSize(0.035);
   nominal__577->GetXaxis()->SetTitleFont(42);
   nominal__577->GetYaxis()->SetTitle("Tracks");
   nominal__577->GetYaxis()->SetLabelFont(42);
   nominal__577->GetYaxis()->SetLabelSize(0.05);
   nominal__577->GetYaxis()->SetTitleSize(0.07);
   nominal__577->GetYaxis()->SetTitleOffset(0);
   nominal__577->GetYaxis()->SetTitleFont(42);
   nominal__577->GetZaxis()->SetLabelFont(42);
   nominal__577->GetZaxis()->SetLabelSize(0.035);
   nominal__577->GetZaxis()->SetTitleSize(0.035);
   nominal__577->GetZaxis()->SetTitleFont(42);
   nominal__577->Draw("");
   Double_t xAxis450[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *ProbQ_down__578 = new TH1F("ProbQ_down__578","",32, xAxis450);
   ProbQ_down__578->SetBinContent(20,4.716178e-05);
   ProbQ_down__578->SetBinContent(21,4.750374e-05);
   ProbQ_down__578->SetBinContent(26,4.856034e-05);
   ProbQ_down__578->SetBinContent(27,5.207413e-05);
   ProbQ_down__578->SetBinContent(28,0.0003452141);
   ProbQ_down__578->SetBinContent(29,0.0005223806);
   ProbQ_down__578->SetBinContent(30,0.003358804);
   ProbQ_down__578->SetBinContent(31,0.01355104);
   ProbQ_down__578->SetBinContent(32,0.1758217);
   ProbQ_down__578->SetBinError(20,4.716178e-05);
   ProbQ_down__578->SetBinError(21,4.750374e-05);
   ProbQ_down__578->SetBinError(26,4.856034e-05);
   ProbQ_down__578->SetBinError(27,5.207413e-05);
   ProbQ_down__578->SetBinError(28,0.0001306065);
   ProbQ_down__578->SetBinError(29,0.0001582745);
   ProbQ_down__578->SetBinError(30,0.0004023621);
   ProbQ_down__578->SetBinError(31,0.0008105387);
   ProbQ_down__578->SetBinError(32,0.001129738);
   ProbQ_down__578->SetBinError(33,0.002682408);
   ProbQ_down__578->SetEntries(4041);
   ProbQ_down__578->SetFillColor(38);
   ProbQ_down__578->SetLineColor(38);
   ProbQ_down__578->SetMarkerColor(38);
   ProbQ_down__578->SetMarkerStyle(21);
   ProbQ_down__578->GetXaxis()->SetTitle("Mass [GeV]");
   ProbQ_down__578->GetXaxis()->SetRange(1,400);
   ProbQ_down__578->GetXaxis()->SetLabelFont(42);
   ProbQ_down__578->GetXaxis()->SetLabelSize(0.035);
   ProbQ_down__578->GetXaxis()->SetTitleSize(0.035);
   ProbQ_down__578->GetXaxis()->SetTitleFont(42);
   ProbQ_down__578->GetYaxis()->SetTitle("Events / bin");
   ProbQ_down__578->GetYaxis()->SetLabelFont(42);
   ProbQ_down__578->GetYaxis()->SetLabelSize(0.035);
   ProbQ_down__578->GetYaxis()->SetTitleSize(0.035);
   ProbQ_down__578->GetYaxis()->SetTitleOffset(0);
   ProbQ_down__578->GetYaxis()->SetTitleFont(42);
   ProbQ_down__578->GetZaxis()->SetLabelFont(42);
   ProbQ_down__578->GetZaxis()->SetLabelSize(0.035);
   ProbQ_down__578->GetZaxis()->SetTitleSize(0.035);
   ProbQ_down__578->GetZaxis()->SetTitleFont(42);
   ProbQ_down__578->Draw("same");
   Double_t xAxis451[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *ProbQ_up__579 = new TH1F("ProbQ_up__579","",32, xAxis451);
   ProbQ_up__579->SetBinContent(20,4.716178e-05);
   ProbQ_up__579->SetBinContent(21,4.750374e-05);
   ProbQ_up__579->SetBinContent(26,4.856034e-05);
   ProbQ_up__579->SetBinContent(27,5.207413e-05);
   ProbQ_up__579->SetBinContent(28,0.0003452141);
   ProbQ_up__579->SetBinContent(29,0.0005223806);
   ProbQ_up__579->SetBinContent(30,0.003358804);
   ProbQ_up__579->SetBinContent(31,0.01355104);
   ProbQ_up__579->SetBinContent(32,0.1758217);
   ProbQ_up__579->SetBinError(20,4.716178e-05);
   ProbQ_up__579->SetBinError(21,4.750374e-05);
   ProbQ_up__579->SetBinError(26,4.856034e-05);
   ProbQ_up__579->SetBinError(27,5.207413e-05);
   ProbQ_up__579->SetBinError(28,0.0001306065);
   ProbQ_up__579->SetBinError(29,0.0001582745);
   ProbQ_up__579->SetBinError(30,0.0004023621);
   ProbQ_up__579->SetBinError(31,0.0008105387);
   ProbQ_up__579->SetBinError(32,0.001129738);
   ProbQ_up__579->SetBinError(33,0.002682408);
   ProbQ_up__579->SetEntries(4041);
   ProbQ_up__579->SetFillColor(46);
   ProbQ_up__579->SetLineColor(46);
   ProbQ_up__579->SetMarkerColor(46);
   ProbQ_up__579->SetMarkerStyle(21);
   ProbQ_up__579->GetXaxis()->SetTitle("Mass [GeV]");
   ProbQ_up__579->GetXaxis()->SetRange(1,400);
   ProbQ_up__579->GetXaxis()->SetLabelFont(42);
   ProbQ_up__579->GetXaxis()->SetLabelSize(0.035);
   ProbQ_up__579->GetXaxis()->SetTitleSize(0.035);
   ProbQ_up__579->GetXaxis()->SetTitleFont(42);
   ProbQ_up__579->GetYaxis()->SetTitle("Events / bin");
   ProbQ_up__579->GetYaxis()->SetLabelFont(42);
   ProbQ_up__579->GetYaxis()->SetLabelSize(0.035);
   ProbQ_up__579->GetYaxis()->SetTitleSize(0.035);
   ProbQ_up__579->GetYaxis()->SetTitleOffset(0);
   ProbQ_up__579->GetYaxis()->SetTitleFont(42);
   ProbQ_up__579->GetZaxis()->SetLabelFont(42);
   ProbQ_up__579->GetZaxis()->SetLabelSize(0.035);
   ProbQ_up__579->GetZaxis()->SetTitleSize(0.035);
   ProbQ_up__579->GetZaxis()->SetTitleFont(42);
   ProbQ_up__579->Draw("same");
   TLine *line = new TLine(1730,0,1730,17.58217);
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
   entry=leg->AddEntry("ProbQ_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("ProbQ_up","Up","PE1");
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
   
   TH1D *frameR2__580 = new TH1D("frameR2__580","",1,0,2000);
   frameR2__580->SetMinimum(0.5);
   frameR2__580->SetMaximum(1.5);
   frameR2__580->SetStats(0);
   frameR2__580->SetLineStyle(0);
   frameR2__580->SetMarkerStyle(20);
   frameR2__580->GetXaxis()->SetRange(1,1);
   frameR2__580->GetXaxis()->SetLabelFont(43);
   frameR2__580->GetXaxis()->SetLabelOffset(0.007);
   frameR2__580->GetXaxis()->SetLabelSize(16);
   frameR2__580->GetXaxis()->SetTitleSize(24);
   frameR2__580->GetXaxis()->SetTitleOffset(3.75);
   frameR2__580->GetXaxis()->SetTitleFont(43);
   frameR2__580->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__580->GetYaxis()->SetNdivisions(205);
   frameR2__580->GetYaxis()->SetLabelFont(43);
   frameR2__580->GetYaxis()->SetLabelOffset(0.007);
   frameR2__580->GetYaxis()->SetLabelSize(20);
   frameR2__580->GetYaxis()->SetTitleSize(20);
   frameR2__580->GetYaxis()->SetTitleOffset(2);
   frameR2__580->GetYaxis()->SetTitleFont(43);
   frameR2__580->GetZaxis()->SetLabelFont(42);
   frameR2__580->GetZaxis()->SetLabelOffset(0.007);
   frameR2__580->GetZaxis()->SetLabelSize(0.05);
   frameR2__580->GetZaxis()->SetTitleSize(0.06);
   frameR2__580->GetZaxis()->SetTitleFont(42);
   frameR2__580->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis452[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *ProbQ_down__581 = new TH1F("ProbQ_down__581","",32, xAxis452);
   ProbQ_down__581->SetBinContent(0,1);
   ProbQ_down__581->SetBinContent(1,1);
   ProbQ_down__581->SetBinContent(2,1);
   ProbQ_down__581->SetBinContent(3,1);
   ProbQ_down__581->SetBinContent(4,1);
   ProbQ_down__581->SetBinContent(5,1);
   ProbQ_down__581->SetBinContent(6,1);
   ProbQ_down__581->SetBinContent(7,1);
   ProbQ_down__581->SetBinContent(8,1);
   ProbQ_down__581->SetBinContent(9,1);
   ProbQ_down__581->SetBinContent(10,1);
   ProbQ_down__581->SetBinContent(11,1);
   ProbQ_down__581->SetBinContent(12,1);
   ProbQ_down__581->SetBinContent(13,1);
   ProbQ_down__581->SetBinContent(14,1);
   ProbQ_down__581->SetBinContent(15,1);
   ProbQ_down__581->SetBinContent(16,1);
   ProbQ_down__581->SetBinContent(17,1);
   ProbQ_down__581->SetBinContent(18,1);
   ProbQ_down__581->SetBinContent(19,1);
   ProbQ_down__581->SetBinContent(20,1);
   ProbQ_down__581->SetBinContent(21,1);
   ProbQ_down__581->SetBinContent(22,1);
   ProbQ_down__581->SetBinContent(23,1);
   ProbQ_down__581->SetBinContent(24,1);
   ProbQ_down__581->SetBinContent(25,1);
   ProbQ_down__581->SetBinContent(26,1);
   ProbQ_down__581->SetBinContent(27,1);
   ProbQ_down__581->SetBinContent(28,1);
   ProbQ_down__581->SetBinContent(29,1);
   ProbQ_down__581->SetBinContent(30,1);
   ProbQ_down__581->SetBinContent(31,1);
   ProbQ_down__581->SetBinContent(32,1);
   ProbQ_down__581->SetBinError(0,0.02230473);
   ProbQ_down__581->SetBinError(1,0.02230473);
   ProbQ_down__581->SetBinError(2,0.02230473);
   ProbQ_down__581->SetBinError(3,0.02230473);
   ProbQ_down__581->SetBinError(4,0.02230473);
   ProbQ_down__581->SetBinError(5,0.02230473);
   ProbQ_down__581->SetBinError(6,0.02230473);
   ProbQ_down__581->SetBinError(7,0.02230473);
   ProbQ_down__581->SetBinError(8,0.02230473);
   ProbQ_down__581->SetBinError(9,0.02230473);
   ProbQ_down__581->SetBinError(10,0.02230473);
   ProbQ_down__581->SetBinError(11,0.02230473);
   ProbQ_down__581->SetBinError(12,0.02230473);
   ProbQ_down__581->SetBinError(13,0.02230473);
   ProbQ_down__581->SetBinError(14,0.02230473);
   ProbQ_down__581->SetBinError(15,0.02230473);
   ProbQ_down__581->SetBinError(16,0.02230473);
   ProbQ_down__581->SetBinError(17,0.02230473);
   ProbQ_down__581->SetBinError(18,0.02230473);
   ProbQ_down__581->SetBinError(19,0.02230473);
   ProbQ_down__581->SetBinError(20,0.02230473);
   ProbQ_down__581->SetBinError(21,0.02230751);
   ProbQ_down__581->SetBinError(22,0.02231028);
   ProbQ_down__581->SetBinError(23,0.02231028);
   ProbQ_down__581->SetBinError(24,0.02231028);
   ProbQ_down__581->SetBinError(25,0.02231028);
   ProbQ_down__581->SetBinError(26,0.02231028);
   ProbQ_down__581->SetBinError(27,0.02231306);
   ProbQ_down__581->SetBinError(28,0.02231582);
   ProbQ_down__581->SetBinError(29,0.02233524);
   ProbQ_down__581->SetBinError(30,0.02236564);
   ProbQ_down__581->SetBinError(31,0.02256313);
   ProbQ_down__581->SetBinError(32,0.02341131);
   ProbQ_down__581->SetEntries(33);
   ProbQ_down__581->SetFillColor(38);
   ProbQ_down__581->SetLineColor(38);
   ProbQ_down__581->SetMarkerColor(38);
   ProbQ_down__581->SetMarkerStyle(21);
   ProbQ_down__581->GetXaxis()->SetTitle("Mass [GeV]");
   ProbQ_down__581->GetXaxis()->SetRange(1,400);
   ProbQ_down__581->GetXaxis()->SetLabelFont(42);
   ProbQ_down__581->GetXaxis()->SetLabelSize(0.035);
   ProbQ_down__581->GetXaxis()->SetTitleSize(0.035);
   ProbQ_down__581->GetXaxis()->SetTitleFont(42);
   ProbQ_down__581->GetYaxis()->SetTitle("Events / bin");
   ProbQ_down__581->GetYaxis()->SetLabelFont(42);
   ProbQ_down__581->GetYaxis()->SetLabelSize(0.035);
   ProbQ_down__581->GetYaxis()->SetTitleSize(0.035);
   ProbQ_down__581->GetYaxis()->SetTitleOffset(0);
   ProbQ_down__581->GetYaxis()->SetTitleFont(42);
   ProbQ_down__581->GetZaxis()->SetLabelFont(42);
   ProbQ_down__581->GetZaxis()->SetLabelSize(0.035);
   ProbQ_down__581->GetZaxis()->SetTitleSize(0.035);
   ProbQ_down__581->GetZaxis()->SetTitleFont(42);
   ProbQ_down__581->Draw("E0 same");
   Double_t xAxis453[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *ProbQ_up__582 = new TH1F("ProbQ_up__582","",32, xAxis453);
   ProbQ_up__582->SetBinContent(0,1);
   ProbQ_up__582->SetBinContent(1,1);
   ProbQ_up__582->SetBinContent(2,1);
   ProbQ_up__582->SetBinContent(3,1);
   ProbQ_up__582->SetBinContent(4,1);
   ProbQ_up__582->SetBinContent(5,1);
   ProbQ_up__582->SetBinContent(6,1);
   ProbQ_up__582->SetBinContent(7,1);
   ProbQ_up__582->SetBinContent(8,1);
   ProbQ_up__582->SetBinContent(9,1);
   ProbQ_up__582->SetBinContent(10,1);
   ProbQ_up__582->SetBinContent(11,1);
   ProbQ_up__582->SetBinContent(12,1);
   ProbQ_up__582->SetBinContent(13,1);
   ProbQ_up__582->SetBinContent(14,1);
   ProbQ_up__582->SetBinContent(15,1);
   ProbQ_up__582->SetBinContent(16,1);
   ProbQ_up__582->SetBinContent(17,1);
   ProbQ_up__582->SetBinContent(18,1);
   ProbQ_up__582->SetBinContent(19,1);
   ProbQ_up__582->SetBinContent(20,1);
   ProbQ_up__582->SetBinContent(21,1);
   ProbQ_up__582->SetBinContent(22,1);
   ProbQ_up__582->SetBinContent(23,1);
   ProbQ_up__582->SetBinContent(24,1);
   ProbQ_up__582->SetBinContent(25,1);
   ProbQ_up__582->SetBinContent(26,1);
   ProbQ_up__582->SetBinContent(27,1);
   ProbQ_up__582->SetBinContent(28,1);
   ProbQ_up__582->SetBinContent(29,1);
   ProbQ_up__582->SetBinContent(30,1);
   ProbQ_up__582->SetBinContent(31,1);
   ProbQ_up__582->SetBinContent(32,1);
   ProbQ_up__582->SetBinError(0,0.02230473);
   ProbQ_up__582->SetBinError(1,0.02230473);
   ProbQ_up__582->SetBinError(2,0.02230473);
   ProbQ_up__582->SetBinError(3,0.02230473);
   ProbQ_up__582->SetBinError(4,0.02230473);
   ProbQ_up__582->SetBinError(5,0.02230473);
   ProbQ_up__582->SetBinError(6,0.02230473);
   ProbQ_up__582->SetBinError(7,0.02230473);
   ProbQ_up__582->SetBinError(8,0.02230473);
   ProbQ_up__582->SetBinError(9,0.02230473);
   ProbQ_up__582->SetBinError(10,0.02230473);
   ProbQ_up__582->SetBinError(11,0.02230473);
   ProbQ_up__582->SetBinError(12,0.02230473);
   ProbQ_up__582->SetBinError(13,0.02230473);
   ProbQ_up__582->SetBinError(14,0.02230473);
   ProbQ_up__582->SetBinError(15,0.02230473);
   ProbQ_up__582->SetBinError(16,0.02230473);
   ProbQ_up__582->SetBinError(17,0.02230473);
   ProbQ_up__582->SetBinError(18,0.02230473);
   ProbQ_up__582->SetBinError(19,0.02230473);
   ProbQ_up__582->SetBinError(20,0.02230473);
   ProbQ_up__582->SetBinError(21,0.02230751);
   ProbQ_up__582->SetBinError(22,0.02231028);
   ProbQ_up__582->SetBinError(23,0.02231028);
   ProbQ_up__582->SetBinError(24,0.02231028);
   ProbQ_up__582->SetBinError(25,0.02231028);
   ProbQ_up__582->SetBinError(26,0.02231028);
   ProbQ_up__582->SetBinError(27,0.02231306);
   ProbQ_up__582->SetBinError(28,0.02231582);
   ProbQ_up__582->SetBinError(29,0.02233524);
   ProbQ_up__582->SetBinError(30,0.02236564);
   ProbQ_up__582->SetBinError(31,0.02256313);
   ProbQ_up__582->SetBinError(32,0.02341131);
   ProbQ_up__582->SetEntries(33);
   ProbQ_up__582->SetFillColor(46);
   ProbQ_up__582->SetLineColor(46);
   ProbQ_up__582->SetMarkerColor(46);
   ProbQ_up__582->SetMarkerStyle(21);
   ProbQ_up__582->GetXaxis()->SetTitle("Mass [GeV]");
   ProbQ_up__582->GetXaxis()->SetRange(1,400);
   ProbQ_up__582->GetXaxis()->SetLabelFont(42);
   ProbQ_up__582->GetXaxis()->SetLabelSize(0.035);
   ProbQ_up__582->GetXaxis()->SetTitleSize(0.035);
   ProbQ_up__582->GetXaxis()->SetTitleFont(42);
   ProbQ_up__582->GetYaxis()->SetTitle("Events / bin");
   ProbQ_up__582->GetYaxis()->SetLabelFont(42);
   ProbQ_up__582->GetYaxis()->SetLabelSize(0.035);
   ProbQ_up__582->GetYaxis()->SetTitleSize(0.035);
   ProbQ_up__582->GetYaxis()->SetTitleOffset(0);
   ProbQ_up__582->GetYaxis()->SetTitleFont(42);
   ProbQ_up__582->GetZaxis()->SetLabelFont(42);
   ProbQ_up__582->GetZaxis()->SetLabelSize(0.035);
   ProbQ_up__582->GetZaxis()->SetTitleSize(0.035);
   ProbQ_up__582->GetZaxis()->SetTitleFont(42);
   ProbQ_up__582->Draw("E0 same");
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
   
   TH1D *frameR2__583 = new TH1D("frameR2__583","",1,0,2000);
   frameR2__583->SetMinimum(0.5);
   frameR2__583->SetMaximum(1.5);
   frameR2__583->SetStats(0);
   frameR2__583->SetLineStyle(0);
   frameR2__583->SetMarkerStyle(20);
   frameR2__583->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__583->GetXaxis()->SetRange(1,1);
   frameR2__583->GetXaxis()->SetLabelFont(43);
   frameR2__583->GetXaxis()->SetLabelOffset(0.007);
   frameR2__583->GetXaxis()->SetLabelSize(16);
   frameR2__583->GetXaxis()->SetTitleSize(24);
   frameR2__583->GetXaxis()->SetTitleOffset(5);
   frameR2__583->GetXaxis()->SetTitleFont(43);
   frameR2__583->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__583->GetYaxis()->SetNdivisions(205);
   frameR2__583->GetYaxis()->SetLabelFont(43);
   frameR2__583->GetYaxis()->SetLabelOffset(0.007);
   frameR2__583->GetYaxis()->SetLabelSize(20);
   frameR2__583->GetYaxis()->SetTitleSize(20);
   frameR2__583->GetYaxis()->SetTitleOffset(2);
   frameR2__583->GetYaxis()->SetTitleFont(43);
   frameR2__583->GetZaxis()->SetLabelFont(42);
   frameR2__583->GetZaxis()->SetLabelOffset(0.007);
   frameR2__583->GetZaxis()->SetLabelSize(0.05);
   frameR2__583->GetZaxis()->SetTitleSize(0.06);
   frameR2__583->GetZaxis()->SetTitleFont(42);
   frameR2__583->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis454[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__584 = new TH1F("nominal__584","",32, xAxis454);
   nominal__584->SetBinContent(20,1);
   nominal__584->SetBinContent(21,1);
   nominal__584->SetBinContent(26,1);
   nominal__584->SetBinContent(27,1);
   nominal__584->SetBinContent(28,1);
   nominal__584->SetBinContent(29,1);
   nominal__584->SetBinContent(30,1);
   nominal__584->SetBinContent(31,1);
   nominal__584->SetBinContent(32,1);
   nominal__584->SetBinError(20,1.414214);
   nominal__584->SetBinError(21,1.414214);
   nominal__584->SetBinError(26,1.414214);
   nominal__584->SetBinError(27,1.414214);
   nominal__584->SetBinError(28,0.5350461);
   nominal__584->SetBinError(29,0.4284881);
   nominal__584->SetBinError(30,0.1694132);
   nominal__584->SetBinError(31,0.08458943);
   nominal__584->SetBinError(32,0.009086997);
   nominal__584->SetMinimum(5e-05);
   nominal__584->SetMaximum(17.58217);
   nominal__584->SetEntries(9.522897);
   nominal__584->SetFillColor(38);
   nominal__584->SetLineColor(38);
   nominal__584->SetMarkerColor(38);
   nominal__584->SetMarkerStyle(21);
   nominal__584->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__584->GetXaxis()->SetRange(1,32);
   nominal__584->GetXaxis()->SetLabelFont(42);
   nominal__584->GetXaxis()->SetLabelSize(0.035);
   nominal__584->GetXaxis()->SetTitleSize(0.035);
   nominal__584->GetXaxis()->SetTitleFont(42);
   nominal__584->GetYaxis()->SetTitle("Tracks");
   nominal__584->GetYaxis()->SetLabelFont(42);
   nominal__584->GetYaxis()->SetLabelSize(0.05);
   nominal__584->GetYaxis()->SetTitleSize(0.07);
   nominal__584->GetYaxis()->SetTitleOffset(0);
   nominal__584->GetYaxis()->SetTitleFont(42);
   nominal__584->GetZaxis()->SetLabelFont(42);
   nominal__584->GetZaxis()->SetLabelSize(0.035);
   nominal__584->GetZaxis()->SetTitleSize(0.035);
   nominal__584->GetZaxis()->SetTitleFont(42);
   nominal__584->Draw("E0 same");
   Double_t xAxis455[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__585 = new TH1F("nominal__585","",32, xAxis455);
   nominal__585->SetBinContent(20,1);
   nominal__585->SetBinContent(21,1);
   nominal__585->SetBinContent(26,1);
   nominal__585->SetBinContent(27,1);
   nominal__585->SetBinContent(28,1);
   nominal__585->SetBinContent(29,1);
   nominal__585->SetBinContent(30,1);
   nominal__585->SetBinContent(31,1);
   nominal__585->SetBinContent(32,1);
   nominal__585->SetBinError(20,1.414214);
   nominal__585->SetBinError(21,1.414214);
   nominal__585->SetBinError(26,1.414214);
   nominal__585->SetBinError(27,1.414214);
   nominal__585->SetBinError(28,0.5350461);
   nominal__585->SetBinError(29,0.4284881);
   nominal__585->SetBinError(30,0.1694132);
   nominal__585->SetBinError(31,0.08458943);
   nominal__585->SetBinError(32,0.009086997);
   nominal__585->SetMinimum(5e-05);
   nominal__585->SetMaximum(17.58217);
   nominal__585->SetEntries(9.522897);
   nominal__585->SetFillColor(46);
   nominal__585->SetLineColor(46);
   nominal__585->SetMarkerColor(46);
   nominal__585->SetMarkerStyle(21);
   nominal__585->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__585->GetXaxis()->SetRange(1,32);
   nominal__585->GetXaxis()->SetLabelFont(42);
   nominal__585->GetXaxis()->SetLabelSize(0.035);
   nominal__585->GetXaxis()->SetTitleSize(0.035);
   nominal__585->GetXaxis()->SetTitleFont(42);
   nominal__585->GetYaxis()->SetTitle("Tracks");
   nominal__585->GetYaxis()->SetLabelFont(42);
   nominal__585->GetYaxis()->SetLabelSize(0.05);
   nominal__585->GetYaxis()->SetTitleSize(0.07);
   nominal__585->GetYaxis()->SetTitleOffset(0);
   nominal__585->GetYaxis()->SetTitleFont(42);
   nominal__585->GetZaxis()->SetLabelFont(42);
   nominal__585->GetZaxis()->SetLabelSize(0.035);
   nominal__585->GetZaxis()->SetTitleSize(0.035);
   nominal__585->GetZaxis()->SetTitleFont(42);
   nominal__585->Draw("E0 same");
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
