void plot_C2()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:10 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.358204,2428.571,7.076637);
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
   Double_t xAxis1191[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1531 = new TH1F("nominal__1531","",32, xAxis1191);
   nominal__1531->SetBinContent(5,9.581986e-07);
   nominal__1531->SetBinContent(9,1.028137e-06);
   nominal__1531->SetBinContent(12,9.955889e-07);
   nominal__1531->SetBinContent(15,1.905316e-06);
   nominal__1531->SetBinContent(16,9.408592e-07);
   nominal__1531->SetBinContent(17,9.408592e-07);
   nominal__1531->SetBinContent(19,9.648228e-07);
   nominal__1531->SetBinContent(23,2.764097e-06);
   nominal__1531->SetBinContent(24,1.89891e-06);
   nominal__1531->SetBinContent(25,4.893618e-06);
   nominal__1531->SetBinContent(26,2.108757e-05);
   nominal__1531->SetBinContent(27,5.90942e-05);
   nominal__1531->SetBinContent(28,0.0004133201);
   nominal__1531->SetBinContent(29,0.002132682);
   nominal__1531->SetBinContent(30,0.007451629);
   nominal__1531->SetBinContent(31,0.01045835);
   nominal__1531->SetBinContent(32,0.008744266);
   nominal__1531->SetBinError(5,9.581986e-07);
   nominal__1531->SetBinError(9,1.028137e-06);
   nominal__1531->SetBinError(12,9.95589e-07);
   nominal__1531->SetBinError(15,1.347372e-06);
   nominal__1531->SetBinError(16,9.408592e-07);
   nominal__1531->SetBinError(17,9.408592e-07);
   nominal__1531->SetBinError(19,9.648228e-07);
   nominal__1531->SetBinError(23,1.598504e-06);
   nominal__1531->SetBinError(24,1.342734e-06);
   nominal__1531->SetBinError(25,2.189295e-06);
   nominal__1531->SetBinError(26,4.501391e-06);
   nominal__1531->SetBinError(27,7.583449e-06);
   nominal__1531->SetBinError(28,1.999415e-05);
   nominal__1531->SetBinError(29,4.53378e-05);
   nominal__1531->SetBinError(30,8.467707e-05);
   nominal__1531->SetBinError(31,0.0001003342);
   nominal__1531->SetBinError(32,6.921264e-05);
   nominal__1531->SetBinError(33,6.023362e-05);
   nominal__1531->SetMinimum(5e-05);
   nominal__1531->SetMaximum(1.045835e+07);
   nominal__1531->SetEntries(30555);
   nominal__1531->SetFillColor(1);
   nominal__1531->SetMarkerStyle(20);
   nominal__1531->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1531->GetXaxis()->SetRange(1,32);
   nominal__1531->GetXaxis()->SetLabelFont(42);
   nominal__1531->GetXaxis()->SetLabelSize(0.035);
   nominal__1531->GetXaxis()->SetTitleSize(0.035);
   nominal__1531->GetXaxis()->SetTitleFont(42);
   nominal__1531->GetYaxis()->SetTitle("Tracks");
   nominal__1531->GetYaxis()->SetLabelFont(42);
   nominal__1531->GetYaxis()->SetLabelSize(0.05);
   nominal__1531->GetYaxis()->SetTitleSize(0.07);
   nominal__1531->GetYaxis()->SetTitleOffset(0);
   nominal__1531->GetYaxis()->SetTitleFont(42);
   nominal__1531->GetZaxis()->SetLabelFont(42);
   nominal__1531->GetZaxis()->SetLabelSize(0.035);
   nominal__1531->GetZaxis()->SetTitleSize(0.035);
   nominal__1531->GetZaxis()->SetTitleFont(42);
   nominal__1531->Draw("");
   Double_t xAxis1192[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_down2__1532 = new TH1F("C_down2__1532","",32, xAxis1192);
   C_down2__1532->SetBinContent(6,9.581986e-07);
   C_down2__1532->SetBinContent(9,1.028137e-06);
   C_down2__1532->SetBinContent(14,9.955889e-07);
   C_down2__1532->SetBinContent(15,9.404936e-07);
   C_down2__1532->SetBinContent(16,1.905682e-06);
   C_down2__1532->SetBinContent(17,9.408592e-07);
   C_down2__1532->SetBinContent(19,9.648228e-07);
   C_down2__1532->SetBinContent(23,2.764097e-06);
   C_down2__1532->SetBinContent(24,1.89891e-06);
   C_down2__1532->SetBinContent(25,2.937386e-06);
   C_down2__1532->SetBinContent(26,2.033346e-05);
   C_down2__1532->SetBinContent(27,5.41404e-05);
   C_down2__1532->SetBinContent(28,0.0003864389);
   C_down2__1532->SetBinContent(29,0.002042578);
   C_down2__1532->SetBinContent(30,0.007332467);
   C_down2__1532->SetBinContent(31,0.01052138);
   C_down2__1532->SetBinContent(32,0.008925045);
   C_down2__1532->SetBinError(6,9.581986e-07);
   C_down2__1532->SetBinError(9,1.028137e-06);
   C_down2__1532->SetBinError(14,9.95589e-07);
   C_down2__1532->SetBinError(15,9.404936e-07);
   C_down2__1532->SetBinError(16,1.347627e-06);
   C_down2__1532->SetBinError(17,9.408592e-07);
   C_down2__1532->SetBinError(19,9.648228e-07);
   C_down2__1532->SetBinError(23,1.598504e-06);
   C_down2__1532->SetBinError(24,1.342734e-06);
   C_down2__1532->SetBinError(25,1.696507e-06);
   C_down2__1532->SetBinError(26,4.438891e-06);
   C_down2__1532->SetBinError(27,7.254992e-06);
   C_down2__1532->SetBinError(28,1.933677e-05);
   C_down2__1532->SetBinError(29,4.437056e-05);
   C_down2__1532->SetBinError(30,8.40006e-05);
   C_down2__1532->SetBinError(31,0.0001006355);
   C_down2__1532->SetBinError(32,6.983334e-05);
   C_down2__1532->SetBinError(33,6.095328e-05);
   C_down2__1532->SetEntries(30555);
   C_down2__1532->SetFillColor(38);
   C_down2__1532->SetLineColor(38);
   C_down2__1532->SetMarkerColor(38);
   C_down2__1532->SetMarkerStyle(21);
   C_down2__1532->GetXaxis()->SetTitle("Mass [GeV]");
   C_down2__1532->GetXaxis()->SetRange(1,400);
   C_down2__1532->GetXaxis()->SetLabelFont(42);
   C_down2__1532->GetXaxis()->SetLabelSize(0.035);
   C_down2__1532->GetXaxis()->SetTitleSize(0.035);
   C_down2__1532->GetXaxis()->SetTitleFont(42);
   C_down2__1532->GetYaxis()->SetTitle("Events / bin");
   C_down2__1532->GetYaxis()->SetLabelFont(42);
   C_down2__1532->GetYaxis()->SetLabelSize(0.035);
   C_down2__1532->GetYaxis()->SetTitleSize(0.035);
   C_down2__1532->GetYaxis()->SetTitleOffset(0);
   C_down2__1532->GetYaxis()->SetTitleFont(42);
   C_down2__1532->GetZaxis()->SetLabelFont(42);
   C_down2__1532->GetZaxis()->SetLabelSize(0.035);
   C_down2__1532->GetZaxis()->SetTitleSize(0.035);
   C_down2__1532->GetZaxis()->SetTitleFont(42);
   C_down2__1532->Draw("same");
   Double_t xAxis1193[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_up2__1533 = new TH1F("C_up2__1533","",32, xAxis1193);
   C_up2__1533->SetBinContent(4,9.581986e-07);
   C_up2__1533->SetBinContent(8,1.028137e-06);
   C_up2__1533->SetBinContent(10,9.955889e-07);
   C_up2__1533->SetBinContent(14,9.404936e-07);
   C_up2__1533->SetBinContent(15,1.905682e-06);
   C_up2__1533->SetBinContent(17,9.408592e-07);
   C_up2__1533->SetBinContent(19,9.648228e-07);
   C_up2__1533->SetBinContent(23,2.764097e-06);
   C_up2__1533->SetBinContent(24,1.89891e-06);
   C_up2__1533->SetBinContent(25,4.893618e-06);
   C_up2__1533->SetBinContent(26,2.203488e-05);
   C_up2__1533->SetBinContent(27,6.763284e-05);
   C_up2__1533->SetBinContent(28,0.0004372059);
   C_up2__1533->SetBinContent(29,0.00220814);
   C_up2__1533->SetBinContent(30,0.0075742);
   C_up2__1533->SetBinContent(31,0.01040039);
   C_up2__1533->SetBinContent(32,0.008570824);
   C_up2__1533->SetBinError(4,9.581986e-07);
   C_up2__1533->SetBinError(8,1.028137e-06);
   C_up2__1533->SetBinError(10,9.95589e-07);
   C_up2__1533->SetBinError(14,9.404936e-07);
   C_up2__1533->SetBinError(15,1.347627e-06);
   C_up2__1533->SetBinError(17,9.408592e-07);
   C_up2__1533->SetBinError(19,9.648228e-07);
   C_up2__1533->SetBinError(23,1.598504e-06);
   C_up2__1533->SetBinError(24,1.342734e-06);
   C_up2__1533->SetBinError(25,2.189295e-06);
   C_up2__1533->SetBinError(26,4.599992e-06);
   C_up2__1533->SetBinError(27,8.101699e-06);
   C_up2__1533->SetBinError(28,2.055978e-05);
   C_up2__1533->SetBinError(29,4.613876e-05);
   C_up2__1533->SetBinError(30,8.53686e-05);
   C_up2__1533->SetBinError(31,0.0001000588);
   C_up2__1533->SetBinError(32,6.865055e-05);
   C_up2__1533->SetBinError(33,5.948051e-05);
   C_up2__1533->SetEntries(30555);
   C_up2__1533->SetFillColor(46);
   C_up2__1533->SetLineColor(46);
   C_up2__1533->SetMarkerColor(46);
   C_up2__1533->SetMarkerStyle(21);
   C_up2__1533->GetXaxis()->SetTitle("Mass [GeV]");
   C_up2__1533->GetXaxis()->SetRange(1,400);
   C_up2__1533->GetXaxis()->SetLabelFont(42);
   C_up2__1533->GetXaxis()->SetLabelSize(0.035);
   C_up2__1533->GetXaxis()->SetTitleSize(0.035);
   C_up2__1533->GetXaxis()->SetTitleFont(42);
   C_up2__1533->GetYaxis()->SetTitle("Events / bin");
   C_up2__1533->GetYaxis()->SetLabelFont(42);
   C_up2__1533->GetYaxis()->SetLabelSize(0.035);
   C_up2__1533->GetYaxis()->SetTitleSize(0.035);
   C_up2__1533->GetYaxis()->SetTitleOffset(0);
   C_up2__1533->GetYaxis()->SetTitleFont(42);
   C_up2__1533->GetZaxis()->SetLabelFont(42);
   C_up2__1533->GetZaxis()->SetLabelSize(0.035);
   C_up2__1533->GetZaxis()->SetTitleSize(0.035);
   C_up2__1533->GetZaxis()->SetTitleFont(42);
   C_up2__1533->Draw("same");
   TLine *line = new TLine(1730,0,1730,1.045835e+07);
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
   
   TH1D *frameR2__1534 = new TH1D("frameR2__1534","",1,0,2000);
   frameR2__1534->SetMinimum(0.5);
   frameR2__1534->SetMaximum(1.5);
   frameR2__1534->SetStats(0);
   frameR2__1534->SetLineStyle(0);
   frameR2__1534->SetMarkerStyle(20);
   frameR2__1534->GetXaxis()->SetRange(1,1);
   frameR2__1534->GetXaxis()->SetLabelFont(43);
   frameR2__1534->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1534->GetXaxis()->SetLabelSize(16);
   frameR2__1534->GetXaxis()->SetTitleSize(24);
   frameR2__1534->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1534->GetXaxis()->SetTitleFont(43);
   frameR2__1534->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1534->GetYaxis()->SetNdivisions(205);
   frameR2__1534->GetYaxis()->SetLabelFont(43);
   frameR2__1534->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1534->GetYaxis()->SetLabelSize(20);
   frameR2__1534->GetYaxis()->SetTitleSize(20);
   frameR2__1534->GetYaxis()->SetTitleOffset(2);
   frameR2__1534->GetYaxis()->SetTitleFont(43);
   frameR2__1534->GetZaxis()->SetLabelFont(42);
   frameR2__1534->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1534->GetZaxis()->SetLabelSize(0.05);
   frameR2__1534->GetZaxis()->SetTitleSize(0.06);
   frameR2__1534->GetZaxis()->SetTitleFont(42);
   frameR2__1534->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1194[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_down2__1535 = new TH1F("C_down2__1535","",32, xAxis1194);
   C_down2__1535->SetBinContent(0,1);
   C_down2__1535->SetBinContent(1,1);
   C_down2__1535->SetBinContent(2,1);
   C_down2__1535->SetBinContent(3,1);
   C_down2__1535->SetBinContent(4,1);
   C_down2__1535->SetBinContent(5,1);
   C_down2__1535->SetBinContent(6,0.9999673);
   C_down2__1535->SetBinContent(7,1);
   C_down2__1535->SetBinContent(8,1);
   C_down2__1535->SetBinContent(9,1);
   C_down2__1535->SetBinContent(10,1);
   C_down2__1535->SetBinContent(11,1);
   C_down2__1535->SetBinContent(12,1);
   C_down2__1535->SetBinContent(13,0.999966);
   C_down2__1535->SetBinContent(14,0.999966);
   C_down2__1535->SetBinContent(15,1);
   C_down2__1535->SetBinContent(16,0.9999671);
   C_down2__1535->SetBinContent(17,1);
   C_down2__1535->SetBinContent(18,1);
   C_down2__1535->SetBinContent(19,1);
   C_down2__1535->SetBinContent(20,1);
   C_down2__1535->SetBinContent(21,1);
   C_down2__1535->SetBinContent(22,1);
   C_down2__1535->SetBinContent(23,1);
   C_down2__1535->SetBinContent(24,1);
   C_down2__1535->SetBinContent(25,1);
   C_down2__1535->SetBinContent(26,0.9999332);
   C_down2__1535->SetBinContent(27,0.9999074);
   C_down2__1535->SetBinContent(28,0.9997376);
   C_down2__1535->SetBinContent(29,0.9988014);
   C_down2__1535->SetBinContent(30,0.9953453);
   C_down2__1535->SetBinContent(31,0.9874625);
   C_down2__1535->SetBinContent(32,0.9797447);
   C_down2__1535->SetBinError(0,0.008106748);
   C_down2__1535->SetBinError(1,0.008106748);
   C_down2__1535->SetBinError(2,0.008106748);
   C_down2__1535->SetBinError(3,0.008106748);
   C_down2__1535->SetBinError(4,0.008106748);
   C_down2__1535->SetBinError(5,0.008106748);
   C_down2__1535->SetBinError(6,0.00810655);
   C_down2__1535->SetBinError(7,0.008106881);
   C_down2__1535->SetBinError(8,0.008106881);
   C_down2__1535->SetBinError(9,0.008106881);
   C_down2__1535->SetBinError(10,0.008107014);
   C_down2__1535->SetBinError(11,0.008107014);
   C_down2__1535->SetBinError(12,0.008107014);
   C_down2__1535->SetBinError(13,0.008106805);
   C_down2__1535->SetBinError(14,0.008106805);
   C_down2__1535->SetBinError(15,0.008107147);
   C_down2__1535->SetBinError(16,0.00810708);
   C_down2__1535->SetBinError(17,0.008107547);
   C_down2__1535->SetBinError(18,0.00810768);
   C_down2__1535->SetBinError(19,0.00810768);
   C_down2__1535->SetBinError(20,0.008107813);
   C_down2__1535->SetBinError(21,0.008107813);
   C_down2__1535->SetBinError(22,0.008107813);
   C_down2__1535->SetBinError(23,0.008107813);
   C_down2__1535->SetBinError(24,0.008108211);
   C_down2__1535->SetBinError(25,0.008108477);
   C_down2__1535->SetBinError(26,0.008108468);
   C_down2__1535->SetBinError(27,0.008111121);
   C_down2__1535->SetBinError(28,0.00811752);
   C_down2__1535->SetBinError(29,0.008165743);
   C_down2__1535->SetBinError(30,0.008449047);
   C_down2__1535->SetBinError(31,0.009856455);
   C_down2__1535->SetBinError(32,0.01446453);
   C_down2__1535->SetEntries(33);
   C_down2__1535->SetFillColor(38);
   C_down2__1535->SetLineColor(38);
   C_down2__1535->SetMarkerColor(38);
   C_down2__1535->SetMarkerStyle(21);
   C_down2__1535->GetXaxis()->SetTitle("Mass [GeV]");
   C_down2__1535->GetXaxis()->SetRange(1,400);
   C_down2__1535->GetXaxis()->SetLabelFont(42);
   C_down2__1535->GetXaxis()->SetLabelSize(0.035);
   C_down2__1535->GetXaxis()->SetTitleSize(0.035);
   C_down2__1535->GetXaxis()->SetTitleFont(42);
   C_down2__1535->GetYaxis()->SetTitle("Events / bin");
   C_down2__1535->GetYaxis()->SetLabelFont(42);
   C_down2__1535->GetYaxis()->SetLabelSize(0.035);
   C_down2__1535->GetYaxis()->SetTitleSize(0.035);
   C_down2__1535->GetYaxis()->SetTitleOffset(0);
   C_down2__1535->GetYaxis()->SetTitleFont(42);
   C_down2__1535->GetZaxis()->SetLabelFont(42);
   C_down2__1535->GetZaxis()->SetLabelSize(0.035);
   C_down2__1535->GetZaxis()->SetTitleSize(0.035);
   C_down2__1535->GetZaxis()->SetTitleFont(42);
   C_down2__1535->Draw("E0 same");
   Double_t xAxis1195[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *C_up2__1536 = new TH1F("C_up2__1536","",32, xAxis1195);
   C_up2__1536->SetBinContent(0,1);
   C_up2__1536->SetBinContent(1,1);
   C_up2__1536->SetBinContent(2,1);
   C_up2__1536->SetBinContent(3,1);
   C_up2__1536->SetBinContent(4,1);
   C_up2__1536->SetBinContent(5,1.000033);
   C_up2__1536->SetBinContent(6,1);
   C_up2__1536->SetBinContent(7,1);
   C_up2__1536->SetBinContent(8,1);
   C_up2__1536->SetBinContent(9,1.000035);
   C_up2__1536->SetBinContent(10,1);
   C_up2__1536->SetBinContent(11,1.000034);
   C_up2__1536->SetBinContent(12,1.000034);
   C_up2__1536->SetBinContent(13,1);
   C_up2__1536->SetBinContent(14,1);
   C_up2__1536->SetBinContent(15,1.000032);
   C_up2__1536->SetBinContent(16,1.000032);
   C_up2__1536->SetBinContent(17,1);
   C_up2__1536->SetBinContent(18,1);
   C_up2__1536->SetBinContent(19,1);
   C_up2__1536->SetBinContent(20,1);
   C_up2__1536->SetBinContent(21,1);
   C_up2__1536->SetBinContent(22,1);
   C_up2__1536->SetBinContent(23,1);
   C_up2__1536->SetBinContent(24,1);
   C_up2__1536->SetBinContent(25,1);
   C_up2__1536->SetBinContent(26,1);
   C_up2__1536->SetBinContent(27,1.000032);
   C_up2__1536->SetBinContent(28,1.000325);
   C_up2__1536->SetBinContent(29,1.001161);
   C_up2__1536->SetBinContent(30,1.0041);
   C_up2__1536->SetBinContent(31,1.012197);
   C_up2__1536->SetBinContent(32,1.020236);
   C_up2__1536->SetBinError(0,0.008106748);
   C_up2__1536->SetBinError(1,0.008106748);
   C_up2__1536->SetBinError(2,0.008106748);
   C_up2__1536->SetBinError(3,0.008106748);
   C_up2__1536->SetBinError(4,0.008106748);
   C_up2__1536->SetBinError(5,0.008107079);
   C_up2__1536->SetBinError(6,0.008106881);
   C_up2__1536->SetBinError(7,0.008106881);
   C_up2__1536->SetBinError(8,0.008106881);
   C_up2__1536->SetBinError(9,0.008107232);
   C_up2__1536->SetBinError(10,0.008107014);
   C_up2__1536->SetBinError(11,0.008107356);
   C_up2__1536->SetBinError(12,0.008107356);
   C_up2__1536->SetBinError(13,0.008107147);
   C_up2__1536->SetBinError(14,0.008107147);
   C_up2__1536->SetBinError(15,0.008107473);
   C_up2__1536->SetBinError(16,0.00810774);
   C_up2__1536->SetBinError(17,0.008107546);
   C_up2__1536->SetBinError(18,0.008107679);
   C_up2__1536->SetBinError(19,0.008107679);
   C_up2__1536->SetBinError(20,0.008107812);
   C_up2__1536->SetBinError(21,0.008107812);
   C_up2__1536->SetBinError(22,0.008107812);
   C_up2__1536->SetBinError(23,0.008107812);
   C_up2__1536->SetBinError(24,0.00810821);
   C_up2__1536->SetBinError(25,0.008108477);
   C_up2__1536->SetBinError(26,0.008109143);
   C_up2__1536->SetBinError(27,0.008112399);
   C_up2__1536->SetBinError(28,0.008123491);
   C_up2__1536->SetBinError(29,0.008189875);
   C_up2__1536->SetBinError(30,0.008541997);
   C_up2__1536->SetBinError(31,0.0101661);
   C_up2__1536->SetBinError(32,0.01521551);
   C_up2__1536->SetEntries(33);
   C_up2__1536->SetFillColor(46);
   C_up2__1536->SetLineColor(46);
   C_up2__1536->SetMarkerColor(46);
   C_up2__1536->SetMarkerStyle(21);
   C_up2__1536->GetXaxis()->SetTitle("Mass [GeV]");
   C_up2__1536->GetXaxis()->SetRange(1,400);
   C_up2__1536->GetXaxis()->SetLabelFont(42);
   C_up2__1536->GetXaxis()->SetLabelSize(0.035);
   C_up2__1536->GetXaxis()->SetTitleSize(0.035);
   C_up2__1536->GetXaxis()->SetTitleFont(42);
   C_up2__1536->GetYaxis()->SetTitle("Events / bin");
   C_up2__1536->GetYaxis()->SetLabelFont(42);
   C_up2__1536->GetYaxis()->SetLabelSize(0.035);
   C_up2__1536->GetYaxis()->SetTitleSize(0.035);
   C_up2__1536->GetYaxis()->SetTitleOffset(0);
   C_up2__1536->GetYaxis()->SetTitleFont(42);
   C_up2__1536->GetZaxis()->SetLabelFont(42);
   C_up2__1536->GetZaxis()->SetLabelSize(0.035);
   C_up2__1536->GetZaxis()->SetTitleSize(0.035);
   C_up2__1536->GetZaxis()->SetTitleFont(42);
   C_up2__1536->Draw("E0 same");
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
   
   TH1D *frameR2__1537 = new TH1D("frameR2__1537","",1,0,2000);
   frameR2__1537->SetMinimum(0.5);
   frameR2__1537->SetMaximum(1.5);
   frameR2__1537->SetStats(0);
   frameR2__1537->SetLineStyle(0);
   frameR2__1537->SetMarkerStyle(20);
   frameR2__1537->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1537->GetXaxis()->SetRange(1,1);
   frameR2__1537->GetXaxis()->SetLabelFont(43);
   frameR2__1537->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1537->GetXaxis()->SetLabelSize(16);
   frameR2__1537->GetXaxis()->SetTitleSize(24);
   frameR2__1537->GetXaxis()->SetTitleOffset(5);
   frameR2__1537->GetXaxis()->SetTitleFont(43);
   frameR2__1537->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1537->GetYaxis()->SetNdivisions(205);
   frameR2__1537->GetYaxis()->SetLabelFont(43);
   frameR2__1537->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1537->GetYaxis()->SetLabelSize(20);
   frameR2__1537->GetYaxis()->SetTitleSize(20);
   frameR2__1537->GetYaxis()->SetTitleOffset(2);
   frameR2__1537->GetYaxis()->SetTitleFont(43);
   frameR2__1537->GetZaxis()->SetLabelFont(42);
   frameR2__1537->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1537->GetZaxis()->SetLabelSize(0.05);
   frameR2__1537->GetZaxis()->SetTitleSize(0.06);
   frameR2__1537->GetZaxis()->SetTitleFont(42);
   frameR2__1537->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1196[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1538 = new TH1F("nominal__1538","",32, xAxis1196);
   nominal__1538->SetBinContent(9,1);
   nominal__1538->SetBinContent(15,2.025869);
   nominal__1538->SetBinContent(16,0.4937126);
   nominal__1538->SetBinContent(17,1);
   nominal__1538->SetBinContent(19,1);
   nominal__1538->SetBinContent(23,1);
   nominal__1538->SetBinContent(24,1);
   nominal__1538->SetBinContent(25,1.665977);
   nominal__1538->SetBinContent(26,1.037087);
   nominal__1538->SetBinContent(27,1.091499);
   nominal__1538->SetBinContent(28,1.069561);
   nominal__1538->SetBinContent(29,1.044113);
   nominal__1538->SetBinContent(30,1.016251);
   nominal__1538->SetBinContent(31,0.9940093);
   nominal__1538->SetBinContent(32,0.9797447);
   nominal__1538->SetBinError(9,1.414214);
   nominal__1538->SetBinError(15,2.48124);
   nominal__1538->SetBinError(16,0.6046879);
   nominal__1538->SetBinError(17,1.414214);
   nominal__1538->SetBinError(19,1.414214);
   nominal__1538->SetBinError(23,0.8178533);
   nominal__1538->SetBinError(24,1.000002);
   nominal__1538->SetBinError(25,1.217097);
   nominal__1538->SetBinError(26,0.316648);
   nominal__1538->SetBinError(27,0.2025165);
   nominal__1538->SetBinError(28,0.07443971);
   nominal__1538->SetBinError(29,0.03173499);
   nominal__1538->SetBinError(30,0.01639822);
   nominal__1538->SetBinError(31,0.01346599);
   nominal__1538->SetBinError(32,0.01090435);
   nominal__1538->SetMinimum(5e-05);
   nominal__1538->SetMaximum(1.045835e+07);
   nominal__1538->SetEntries(17.03737);
   nominal__1538->SetFillColor(38);
   nominal__1538->SetLineColor(38);
   nominal__1538->SetMarkerColor(38);
   nominal__1538->SetMarkerStyle(21);
   nominal__1538->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1538->GetXaxis()->SetRange(1,32);
   nominal__1538->GetXaxis()->SetLabelFont(42);
   nominal__1538->GetXaxis()->SetLabelSize(0.035);
   nominal__1538->GetXaxis()->SetTitleSize(0.035);
   nominal__1538->GetXaxis()->SetTitleFont(42);
   nominal__1538->GetYaxis()->SetTitle("Tracks");
   nominal__1538->GetYaxis()->SetLabelFont(42);
   nominal__1538->GetYaxis()->SetLabelSize(0.05);
   nominal__1538->GetYaxis()->SetTitleSize(0.07);
   nominal__1538->GetYaxis()->SetTitleOffset(0);
   nominal__1538->GetYaxis()->SetTitleFont(42);
   nominal__1538->GetZaxis()->SetLabelFont(42);
   nominal__1538->GetZaxis()->SetLabelSize(0.035);
   nominal__1538->GetZaxis()->SetTitleSize(0.035);
   nominal__1538->GetZaxis()->SetTitleFont(42);
   nominal__1538->Draw("E0 same");
   Double_t xAxis1197[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1539 = new TH1F("nominal__1539","",32, xAxis1197);
   nominal__1539->SetBinContent(15,0.9998081);
   nominal__1539->SetBinContent(17,1);
   nominal__1539->SetBinContent(19,1);
   nominal__1539->SetBinContent(23,1);
   nominal__1539->SetBinContent(24,1);
   nominal__1539->SetBinContent(25,1);
   nominal__1539->SetBinContent(26,0.9570085);
   nominal__1539->SetBinContent(27,0.8737501);
   nominal__1539->SetBinContent(28,0.9453671);
   nominal__1539->SetBinContent(29,0.9658273);
   nominal__1539->SetBinContent(30,0.9838173);
   nominal__1539->SetBinContent(31,1.005573);
   nominal__1539->SetBinContent(32,1.020236);
   nominal__1539->SetBinError(15,0.9998884);
   nominal__1539->SetBinError(17,1.414214);
   nominal__1539->SetBinError(19,1.414214);
   nominal__1539->SetBinError(23,0.8178533);
   nominal__1539->SetBinError(24,1.000002);
   nominal__1539->SetBinError(25,0.6326876);
   nominal__1539->SetBinError(26,0.2857379);
   nominal__1539->SetBinError(27,0.1533864);
   nominal__1539->SetBinError(28,0.06377888);
   nominal__1539->SetBinError(29,0.02878947);
   nominal__1539->SetBinError(30,0.01574616);
   nominal__1539->SetBinError(31,0.01366234);
   nominal__1539->SetBinError(32,0.01148875);
   nominal__1539->SetMinimum(5e-05);
   nominal__1539->SetMaximum(1.045835e+07);
   nominal__1539->SetEntries(22.64722);
   nominal__1539->SetFillColor(46);
   nominal__1539->SetLineColor(46);
   nominal__1539->SetMarkerColor(46);
   nominal__1539->SetMarkerStyle(21);
   nominal__1539->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1539->GetXaxis()->SetRange(1,32);
   nominal__1539->GetXaxis()->SetLabelFont(42);
   nominal__1539->GetXaxis()->SetLabelSize(0.035);
   nominal__1539->GetXaxis()->SetTitleSize(0.035);
   nominal__1539->GetXaxis()->SetTitleFont(42);
   nominal__1539->GetYaxis()->SetTitle("Tracks");
   nominal__1539->GetYaxis()->SetLabelFont(42);
   nominal__1539->GetYaxis()->SetLabelSize(0.05);
   nominal__1539->GetYaxis()->SetTitleSize(0.07);
   nominal__1539->GetYaxis()->SetTitleOffset(0);
   nominal__1539->GetYaxis()->SetTitleFont(42);
   nominal__1539->GetZaxis()->SetLabelFont(42);
   nominal__1539->GetZaxis()->SetLabelSize(0.035);
   nominal__1539->GetZaxis()->SetTitleSize(0.035);
   nominal__1539->GetZaxis()->SetTitleFont(42);
   nominal__1539->Draw("E0 same");
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
