void plot_Pt()
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
   t1->Range(-428.5714,-4.332952,2428.571,2.051385);
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
   Double_t xAxis1156[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1486 = new TH1F("nominal__1486","",32, xAxis1156);
   nominal__1486->SetBinContent(5,9.581986e-07);
   nominal__1486->SetBinContent(9,1.028137e-06);
   nominal__1486->SetBinContent(12,9.955889e-07);
   nominal__1486->SetBinContent(15,1.905316e-06);
   nominal__1486->SetBinContent(16,9.408592e-07);
   nominal__1486->SetBinContent(17,9.408592e-07);
   nominal__1486->SetBinContent(19,9.648228e-07);
   nominal__1486->SetBinContent(23,2.764097e-06);
   nominal__1486->SetBinContent(24,1.89891e-06);
   nominal__1486->SetBinContent(25,4.893618e-06);
   nominal__1486->SetBinContent(26,2.108757e-05);
   nominal__1486->SetBinContent(27,5.90942e-05);
   nominal__1486->SetBinContent(28,0.0004133201);
   nominal__1486->SetBinContent(29,0.002132682);
   nominal__1486->SetBinContent(30,0.007451629);
   nominal__1486->SetBinContent(31,0.01045835);
   nominal__1486->SetBinContent(32,0.008744266);
   nominal__1486->SetBinError(5,9.581986e-07);
   nominal__1486->SetBinError(9,1.028137e-06);
   nominal__1486->SetBinError(12,9.95589e-07);
   nominal__1486->SetBinError(15,1.347372e-06);
   nominal__1486->SetBinError(16,9.408592e-07);
   nominal__1486->SetBinError(17,9.408592e-07);
   nominal__1486->SetBinError(19,9.648228e-07);
   nominal__1486->SetBinError(23,1.598504e-06);
   nominal__1486->SetBinError(24,1.342734e-06);
   nominal__1486->SetBinError(25,2.189295e-06);
   nominal__1486->SetBinError(26,4.501391e-06);
   nominal__1486->SetBinError(27,7.583449e-06);
   nominal__1486->SetBinError(28,1.999415e-05);
   nominal__1486->SetBinError(29,4.53378e-05);
   nominal__1486->SetBinError(30,8.467707e-05);
   nominal__1486->SetBinError(31,0.0001003342);
   nominal__1486->SetBinError(32,6.921264e-05);
   nominal__1486->SetBinError(33,6.023362e-05);
   nominal__1486->SetMinimum(5e-05);
   nominal__1486->SetMaximum(104.5835);
   nominal__1486->SetEntries(30555);
   nominal__1486->SetFillColor(1);
   nominal__1486->SetMarkerStyle(20);
   nominal__1486->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1486->GetXaxis()->SetRange(1,32);
   nominal__1486->GetXaxis()->SetLabelFont(42);
   nominal__1486->GetXaxis()->SetLabelSize(0.035);
   nominal__1486->GetXaxis()->SetTitleSize(0.035);
   nominal__1486->GetXaxis()->SetTitleFont(42);
   nominal__1486->GetYaxis()->SetTitle("Tracks");
   nominal__1486->GetYaxis()->SetLabelFont(42);
   nominal__1486->GetYaxis()->SetLabelSize(0.05);
   nominal__1486->GetYaxis()->SetTitleSize(0.07);
   nominal__1486->GetYaxis()->SetTitleOffset(0);
   nominal__1486->GetYaxis()->SetTitleFont(42);
   nominal__1486->GetZaxis()->SetLabelFont(42);
   nominal__1486->GetZaxis()->SetLabelSize(0.035);
   nominal__1486->GetZaxis()->SetTitleSize(0.035);
   nominal__1486->GetZaxis()->SetTitleFont(42);
   nominal__1486->Draw("");
   Double_t xAxis1157[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_down__1487 = new TH1F("Pt_down__1487","",32, xAxis1157);
   Pt_down__1487->SetBinContent(5,9.581986e-07);
   Pt_down__1487->SetBinContent(9,1.028137e-06);
   Pt_down__1487->SetBinContent(12,9.955889e-07);
   Pt_down__1487->SetBinContent(15,1.905316e-06);
   Pt_down__1487->SetBinContent(16,9.408592e-07);
   Pt_down__1487->SetBinContent(17,9.408592e-07);
   Pt_down__1487->SetBinContent(19,9.648228e-07);
   Pt_down__1487->SetBinContent(23,2.764097e-06);
   Pt_down__1487->SetBinContent(24,1.89891e-06);
   Pt_down__1487->SetBinContent(25,4.893618e-06);
   Pt_down__1487->SetBinContent(26,2.108757e-05);
   Pt_down__1487->SetBinContent(27,5.90942e-05);
   Pt_down__1487->SetBinContent(28,0.0004133201);
   Pt_down__1487->SetBinContent(29,0.002132682);
   Pt_down__1487->SetBinContent(30,0.007451629);
   Pt_down__1487->SetBinContent(31,0.01045835);
   Pt_down__1487->SetBinContent(32,0.008744266);
   Pt_down__1487->SetBinError(5,9.581986e-07);
   Pt_down__1487->SetBinError(9,1.028137e-06);
   Pt_down__1487->SetBinError(12,9.95589e-07);
   Pt_down__1487->SetBinError(15,1.347372e-06);
   Pt_down__1487->SetBinError(16,9.408592e-07);
   Pt_down__1487->SetBinError(17,9.408592e-07);
   Pt_down__1487->SetBinError(19,9.648228e-07);
   Pt_down__1487->SetBinError(23,1.598504e-06);
   Pt_down__1487->SetBinError(24,1.342734e-06);
   Pt_down__1487->SetBinError(25,2.189295e-06);
   Pt_down__1487->SetBinError(26,4.501391e-06);
   Pt_down__1487->SetBinError(27,7.583449e-06);
   Pt_down__1487->SetBinError(28,1.999415e-05);
   Pt_down__1487->SetBinError(29,4.53378e-05);
   Pt_down__1487->SetBinError(30,8.467707e-05);
   Pt_down__1487->SetBinError(31,0.0001003342);
   Pt_down__1487->SetBinError(32,6.921264e-05);
   Pt_down__1487->SetBinError(33,6.023362e-05);
   Pt_down__1487->SetEntries(30555);
   Pt_down__1487->SetFillColor(38);
   Pt_down__1487->SetLineColor(38);
   Pt_down__1487->SetMarkerColor(38);
   Pt_down__1487->SetMarkerStyle(21);
   Pt_down__1487->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_down__1487->GetXaxis()->SetRange(1,400);
   Pt_down__1487->GetXaxis()->SetLabelFont(42);
   Pt_down__1487->GetXaxis()->SetLabelSize(0.035);
   Pt_down__1487->GetXaxis()->SetTitleSize(0.035);
   Pt_down__1487->GetXaxis()->SetTitleFont(42);
   Pt_down__1487->GetYaxis()->SetTitle("Events / bin");
   Pt_down__1487->GetYaxis()->SetLabelFont(42);
   Pt_down__1487->GetYaxis()->SetLabelSize(0.035);
   Pt_down__1487->GetYaxis()->SetTitleSize(0.035);
   Pt_down__1487->GetYaxis()->SetTitleOffset(0);
   Pt_down__1487->GetYaxis()->SetTitleFont(42);
   Pt_down__1487->GetZaxis()->SetLabelFont(42);
   Pt_down__1487->GetZaxis()->SetLabelSize(0.035);
   Pt_down__1487->GetZaxis()->SetTitleSize(0.035);
   Pt_down__1487->GetZaxis()->SetTitleFont(42);
   Pt_down__1487->Draw("same");
   Double_t xAxis1158[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_up__1488 = new TH1F("Pt_up__1488","",32, xAxis1158);
   Pt_up__1488->SetBinContent(5,9.581986e-07);
   Pt_up__1488->SetBinContent(9,1.028137e-06);
   Pt_up__1488->SetBinContent(12,9.955889e-07);
   Pt_up__1488->SetBinContent(15,1.905316e-06);
   Pt_up__1488->SetBinContent(16,9.408592e-07);
   Pt_up__1488->SetBinContent(17,9.408592e-07);
   Pt_up__1488->SetBinContent(19,9.648228e-07);
   Pt_up__1488->SetBinContent(23,2.764097e-06);
   Pt_up__1488->SetBinContent(24,1.89891e-06);
   Pt_up__1488->SetBinContent(25,4.893618e-06);
   Pt_up__1488->SetBinContent(26,2.108757e-05);
   Pt_up__1488->SetBinContent(27,5.90942e-05);
   Pt_up__1488->SetBinContent(28,0.0004133201);
   Pt_up__1488->SetBinContent(29,0.002132682);
   Pt_up__1488->SetBinContent(30,0.007451629);
   Pt_up__1488->SetBinContent(31,0.01045835);
   Pt_up__1488->SetBinContent(32,0.008744266);
   Pt_up__1488->SetBinError(5,9.581986e-07);
   Pt_up__1488->SetBinError(9,1.028137e-06);
   Pt_up__1488->SetBinError(12,9.95589e-07);
   Pt_up__1488->SetBinError(15,1.347372e-06);
   Pt_up__1488->SetBinError(16,9.408592e-07);
   Pt_up__1488->SetBinError(17,9.408592e-07);
   Pt_up__1488->SetBinError(19,9.648228e-07);
   Pt_up__1488->SetBinError(23,1.598504e-06);
   Pt_up__1488->SetBinError(24,1.342734e-06);
   Pt_up__1488->SetBinError(25,2.189295e-06);
   Pt_up__1488->SetBinError(26,4.501391e-06);
   Pt_up__1488->SetBinError(27,7.583449e-06);
   Pt_up__1488->SetBinError(28,1.999415e-05);
   Pt_up__1488->SetBinError(29,4.53378e-05);
   Pt_up__1488->SetBinError(30,8.467707e-05);
   Pt_up__1488->SetBinError(31,0.0001003342);
   Pt_up__1488->SetBinError(32,6.921264e-05);
   Pt_up__1488->SetBinError(33,6.023362e-05);
   Pt_up__1488->SetEntries(30555);
   Pt_up__1488->SetFillColor(46);
   Pt_up__1488->SetLineColor(46);
   Pt_up__1488->SetMarkerColor(46);
   Pt_up__1488->SetMarkerStyle(21);
   Pt_up__1488->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_up__1488->GetXaxis()->SetRange(1,400);
   Pt_up__1488->GetXaxis()->SetLabelFont(42);
   Pt_up__1488->GetXaxis()->SetLabelSize(0.035);
   Pt_up__1488->GetXaxis()->SetTitleSize(0.035);
   Pt_up__1488->GetXaxis()->SetTitleFont(42);
   Pt_up__1488->GetYaxis()->SetTitle("Events / bin");
   Pt_up__1488->GetYaxis()->SetLabelFont(42);
   Pt_up__1488->GetYaxis()->SetLabelSize(0.035);
   Pt_up__1488->GetYaxis()->SetTitleSize(0.035);
   Pt_up__1488->GetYaxis()->SetTitleOffset(0);
   Pt_up__1488->GetYaxis()->SetTitleFont(42);
   Pt_up__1488->GetZaxis()->SetLabelFont(42);
   Pt_up__1488->GetZaxis()->SetLabelSize(0.035);
   Pt_up__1488->GetZaxis()->SetTitleSize(0.035);
   Pt_up__1488->GetZaxis()->SetTitleFont(42);
   Pt_up__1488->Draw("same");
   TLine *line = new TLine(1730,0,1730,104.5835);
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
   entry=leg->AddEntry("Pt_down","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Pt_up","Up","PE1");
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
   
   TH1D *frameR2__1489 = new TH1D("frameR2__1489","",1,0,2000);
   frameR2__1489->SetMinimum(0.5);
   frameR2__1489->SetMaximum(1.5);
   frameR2__1489->SetStats(0);
   frameR2__1489->SetLineStyle(0);
   frameR2__1489->SetMarkerStyle(20);
   frameR2__1489->GetXaxis()->SetRange(1,1);
   frameR2__1489->GetXaxis()->SetLabelFont(43);
   frameR2__1489->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1489->GetXaxis()->SetLabelSize(16);
   frameR2__1489->GetXaxis()->SetTitleSize(24);
   frameR2__1489->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1489->GetXaxis()->SetTitleFont(43);
   frameR2__1489->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1489->GetYaxis()->SetNdivisions(205);
   frameR2__1489->GetYaxis()->SetLabelFont(43);
   frameR2__1489->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1489->GetYaxis()->SetLabelSize(20);
   frameR2__1489->GetYaxis()->SetTitleSize(20);
   frameR2__1489->GetYaxis()->SetTitleOffset(2);
   frameR2__1489->GetYaxis()->SetTitleFont(43);
   frameR2__1489->GetZaxis()->SetLabelFont(42);
   frameR2__1489->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1489->GetZaxis()->SetLabelSize(0.05);
   frameR2__1489->GetZaxis()->SetTitleSize(0.06);
   frameR2__1489->GetZaxis()->SetTitleFont(42);
   frameR2__1489->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1159[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_down__1490 = new TH1F("Pt_down__1490","",32, xAxis1159);
   Pt_down__1490->SetBinContent(0,1);
   Pt_down__1490->SetBinContent(1,1);
   Pt_down__1490->SetBinContent(2,1);
   Pt_down__1490->SetBinContent(3,1);
   Pt_down__1490->SetBinContent(4,1);
   Pt_down__1490->SetBinContent(5,1);
   Pt_down__1490->SetBinContent(6,1);
   Pt_down__1490->SetBinContent(7,1);
   Pt_down__1490->SetBinContent(8,1);
   Pt_down__1490->SetBinContent(9,1);
   Pt_down__1490->SetBinContent(10,1);
   Pt_down__1490->SetBinContent(11,1);
   Pt_down__1490->SetBinContent(12,1);
   Pt_down__1490->SetBinContent(13,1);
   Pt_down__1490->SetBinContent(14,1);
   Pt_down__1490->SetBinContent(15,1);
   Pt_down__1490->SetBinContent(16,1);
   Pt_down__1490->SetBinContent(17,1);
   Pt_down__1490->SetBinContent(18,1);
   Pt_down__1490->SetBinContent(19,1);
   Pt_down__1490->SetBinContent(20,1);
   Pt_down__1490->SetBinContent(21,1);
   Pt_down__1490->SetBinContent(22,1);
   Pt_down__1490->SetBinContent(23,1);
   Pt_down__1490->SetBinContent(24,1);
   Pt_down__1490->SetBinContent(25,1);
   Pt_down__1490->SetBinContent(26,1);
   Pt_down__1490->SetBinContent(27,1);
   Pt_down__1490->SetBinContent(28,1);
   Pt_down__1490->SetBinContent(29,1);
   Pt_down__1490->SetBinContent(30,1);
   Pt_down__1490->SetBinContent(31,1);
   Pt_down__1490->SetBinContent(32,1);
   Pt_down__1490->SetBinError(0,0.008106748);
   Pt_down__1490->SetBinError(1,0.008106748);
   Pt_down__1490->SetBinError(2,0.008106748);
   Pt_down__1490->SetBinError(3,0.008106748);
   Pt_down__1490->SetBinError(4,0.008106748);
   Pt_down__1490->SetBinError(5,0.008106748);
   Pt_down__1490->SetBinError(6,0.008106881);
   Pt_down__1490->SetBinError(7,0.008106881);
   Pt_down__1490->SetBinError(8,0.008106881);
   Pt_down__1490->SetBinError(9,0.008106881);
   Pt_down__1490->SetBinError(10,0.008107014);
   Pt_down__1490->SetBinError(11,0.008107014);
   Pt_down__1490->SetBinError(12,0.008107014);
   Pt_down__1490->SetBinError(13,0.008107147);
   Pt_down__1490->SetBinError(14,0.008107147);
   Pt_down__1490->SetBinError(15,0.008107147);
   Pt_down__1490->SetBinError(16,0.008107413);
   Pt_down__1490->SetBinError(17,0.008107546);
   Pt_down__1490->SetBinError(18,0.00810768);
   Pt_down__1490->SetBinError(19,0.00810768);
   Pt_down__1490->SetBinError(20,0.008107813);
   Pt_down__1490->SetBinError(21,0.008107813);
   Pt_down__1490->SetBinError(22,0.008107813);
   Pt_down__1490->SetBinError(23,0.008107813);
   Pt_down__1490->SetBinError(24,0.008108211);
   Pt_down__1490->SetBinError(25,0.008108477);
   Pt_down__1490->SetBinError(26,0.008109143);
   Pt_down__1490->SetBinError(27,0.00811207);
   Pt_down__1490->SetBinError(28,0.008120185);
   Pt_down__1490->SetBinError(29,0.008177997);
   Pt_down__1490->SetBinError(30,0.008498454);
   Pt_down__1490->SetBinError(31,0.01001312);
   Pt_down__1490->SetBinError(32,0.01483913);
   Pt_down__1490->SetEntries(33);
   Pt_down__1490->SetFillColor(38);
   Pt_down__1490->SetLineColor(38);
   Pt_down__1490->SetMarkerColor(38);
   Pt_down__1490->SetMarkerStyle(21);
   Pt_down__1490->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_down__1490->GetXaxis()->SetRange(1,400);
   Pt_down__1490->GetXaxis()->SetLabelFont(42);
   Pt_down__1490->GetXaxis()->SetLabelSize(0.035);
   Pt_down__1490->GetXaxis()->SetTitleSize(0.035);
   Pt_down__1490->GetXaxis()->SetTitleFont(42);
   Pt_down__1490->GetYaxis()->SetTitle("Events / bin");
   Pt_down__1490->GetYaxis()->SetLabelFont(42);
   Pt_down__1490->GetYaxis()->SetLabelSize(0.035);
   Pt_down__1490->GetYaxis()->SetTitleSize(0.035);
   Pt_down__1490->GetYaxis()->SetTitleOffset(0);
   Pt_down__1490->GetYaxis()->SetTitleFont(42);
   Pt_down__1490->GetZaxis()->SetLabelFont(42);
   Pt_down__1490->GetZaxis()->SetLabelSize(0.035);
   Pt_down__1490->GetZaxis()->SetTitleSize(0.035);
   Pt_down__1490->GetZaxis()->SetTitleFont(42);
   Pt_down__1490->Draw("E0 same");
   Double_t xAxis1160[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Pt_up__1491 = new TH1F("Pt_up__1491","",32, xAxis1160);
   Pt_up__1491->SetBinContent(0,1);
   Pt_up__1491->SetBinContent(1,1);
   Pt_up__1491->SetBinContent(2,1);
   Pt_up__1491->SetBinContent(3,1);
   Pt_up__1491->SetBinContent(4,1);
   Pt_up__1491->SetBinContent(5,1);
   Pt_up__1491->SetBinContent(6,1);
   Pt_up__1491->SetBinContent(7,1);
   Pt_up__1491->SetBinContent(8,1);
   Pt_up__1491->SetBinContent(9,1);
   Pt_up__1491->SetBinContent(10,1);
   Pt_up__1491->SetBinContent(11,1);
   Pt_up__1491->SetBinContent(12,1);
   Pt_up__1491->SetBinContent(13,1);
   Pt_up__1491->SetBinContent(14,1);
   Pt_up__1491->SetBinContent(15,1);
   Pt_up__1491->SetBinContent(16,1);
   Pt_up__1491->SetBinContent(17,1);
   Pt_up__1491->SetBinContent(18,1);
   Pt_up__1491->SetBinContent(19,1);
   Pt_up__1491->SetBinContent(20,1);
   Pt_up__1491->SetBinContent(21,1);
   Pt_up__1491->SetBinContent(22,1);
   Pt_up__1491->SetBinContent(23,1);
   Pt_up__1491->SetBinContent(24,1);
   Pt_up__1491->SetBinContent(25,1);
   Pt_up__1491->SetBinContent(26,1);
   Pt_up__1491->SetBinContent(27,1);
   Pt_up__1491->SetBinContent(28,1);
   Pt_up__1491->SetBinContent(29,1);
   Pt_up__1491->SetBinContent(30,1);
   Pt_up__1491->SetBinContent(31,1);
   Pt_up__1491->SetBinContent(32,1);
   Pt_up__1491->SetBinError(0,0.008106748);
   Pt_up__1491->SetBinError(1,0.008106748);
   Pt_up__1491->SetBinError(2,0.008106748);
   Pt_up__1491->SetBinError(3,0.008106748);
   Pt_up__1491->SetBinError(4,0.008106748);
   Pt_up__1491->SetBinError(5,0.008106748);
   Pt_up__1491->SetBinError(6,0.008106881);
   Pt_up__1491->SetBinError(7,0.008106881);
   Pt_up__1491->SetBinError(8,0.008106881);
   Pt_up__1491->SetBinError(9,0.008106881);
   Pt_up__1491->SetBinError(10,0.008107014);
   Pt_up__1491->SetBinError(11,0.008107014);
   Pt_up__1491->SetBinError(12,0.008107014);
   Pt_up__1491->SetBinError(13,0.008107147);
   Pt_up__1491->SetBinError(14,0.008107147);
   Pt_up__1491->SetBinError(15,0.008107147);
   Pt_up__1491->SetBinError(16,0.008107413);
   Pt_up__1491->SetBinError(17,0.008107546);
   Pt_up__1491->SetBinError(18,0.00810768);
   Pt_up__1491->SetBinError(19,0.00810768);
   Pt_up__1491->SetBinError(20,0.008107813);
   Pt_up__1491->SetBinError(21,0.008107813);
   Pt_up__1491->SetBinError(22,0.008107813);
   Pt_up__1491->SetBinError(23,0.008107813);
   Pt_up__1491->SetBinError(24,0.008108211);
   Pt_up__1491->SetBinError(25,0.008108477);
   Pt_up__1491->SetBinError(26,0.008109143);
   Pt_up__1491->SetBinError(27,0.00811207);
   Pt_up__1491->SetBinError(28,0.008120185);
   Pt_up__1491->SetBinError(29,0.008177997);
   Pt_up__1491->SetBinError(30,0.008498454);
   Pt_up__1491->SetBinError(31,0.01001312);
   Pt_up__1491->SetBinError(32,0.01483913);
   Pt_up__1491->SetEntries(33);
   Pt_up__1491->SetFillColor(46);
   Pt_up__1491->SetLineColor(46);
   Pt_up__1491->SetMarkerColor(46);
   Pt_up__1491->SetMarkerStyle(21);
   Pt_up__1491->GetXaxis()->SetTitle("Mass [GeV]");
   Pt_up__1491->GetXaxis()->SetRange(1,400);
   Pt_up__1491->GetXaxis()->SetLabelFont(42);
   Pt_up__1491->GetXaxis()->SetLabelSize(0.035);
   Pt_up__1491->GetXaxis()->SetTitleSize(0.035);
   Pt_up__1491->GetXaxis()->SetTitleFont(42);
   Pt_up__1491->GetYaxis()->SetTitle("Events / bin");
   Pt_up__1491->GetYaxis()->SetLabelFont(42);
   Pt_up__1491->GetYaxis()->SetLabelSize(0.035);
   Pt_up__1491->GetYaxis()->SetTitleSize(0.035);
   Pt_up__1491->GetYaxis()->SetTitleOffset(0);
   Pt_up__1491->GetYaxis()->SetTitleFont(42);
   Pt_up__1491->GetZaxis()->SetLabelFont(42);
   Pt_up__1491->GetZaxis()->SetLabelSize(0.035);
   Pt_up__1491->GetZaxis()->SetTitleSize(0.035);
   Pt_up__1491->GetZaxis()->SetTitleFont(42);
   Pt_up__1491->Draw("E0 same");
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
   
   TH1D *frameR2__1492 = new TH1D("frameR2__1492","",1,0,2000);
   frameR2__1492->SetMinimum(0.5);
   frameR2__1492->SetMaximum(1.5);
   frameR2__1492->SetStats(0);
   frameR2__1492->SetLineStyle(0);
   frameR2__1492->SetMarkerStyle(20);
   frameR2__1492->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1492->GetXaxis()->SetRange(1,1);
   frameR2__1492->GetXaxis()->SetLabelFont(43);
   frameR2__1492->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1492->GetXaxis()->SetLabelSize(16);
   frameR2__1492->GetXaxis()->SetTitleSize(24);
   frameR2__1492->GetXaxis()->SetTitleOffset(5);
   frameR2__1492->GetXaxis()->SetTitleFont(43);
   frameR2__1492->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1492->GetYaxis()->SetNdivisions(205);
   frameR2__1492->GetYaxis()->SetLabelFont(43);
   frameR2__1492->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1492->GetYaxis()->SetLabelSize(20);
   frameR2__1492->GetYaxis()->SetTitleSize(20);
   frameR2__1492->GetYaxis()->SetTitleOffset(2);
   frameR2__1492->GetYaxis()->SetTitleFont(43);
   frameR2__1492->GetZaxis()->SetLabelFont(42);
   frameR2__1492->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1492->GetZaxis()->SetLabelSize(0.05);
   frameR2__1492->GetZaxis()->SetTitleSize(0.06);
   frameR2__1492->GetZaxis()->SetTitleFont(42);
   frameR2__1492->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1161[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1493 = new TH1F("nominal__1493","",32, xAxis1161);
   nominal__1493->SetBinContent(5,1);
   nominal__1493->SetBinContent(9,1);
   nominal__1493->SetBinContent(12,1);
   nominal__1493->SetBinContent(15,1);
   nominal__1493->SetBinContent(16,1);
   nominal__1493->SetBinContent(17,1);
   nominal__1493->SetBinContent(19,1);
   nominal__1493->SetBinContent(23,1);
   nominal__1493->SetBinContent(24,1);
   nominal__1493->SetBinContent(25,1);
   nominal__1493->SetBinContent(26,1);
   nominal__1493->SetBinContent(27,1);
   nominal__1493->SetBinContent(28,1);
   nominal__1493->SetBinContent(29,1);
   nominal__1493->SetBinContent(30,1);
   nominal__1493->SetBinContent(31,1);
   nominal__1493->SetBinContent(32,1);
   nominal__1493->SetBinError(5,1.414214);
   nominal__1493->SetBinError(9,1.414214);
   nominal__1493->SetBinError(12,1.414214);
   nominal__1493->SetBinError(15,1.000082);
   nominal__1493->SetBinError(16,1.414214);
   nominal__1493->SetBinError(17,1.414214);
   nominal__1493->SetBinError(19,1.414214);
   nominal__1493->SetBinError(23,0.8178533);
   nominal__1493->SetBinError(24,1.000002);
   nominal__1493->SetBinError(25,0.6326876);
   nominal__1493->SetBinError(26,0.3018806);
   nominal__1493->SetBinError(27,0.1814834);
   nominal__1493->SetBinError(28,0.06841188);
   nominal__1493->SetBinError(29,0.03006418);
   nominal__1493->SetBinError(30,0.01607051);
   nominal__1493->SetBinError(31,0.01356753);
   nominal__1493->SetBinError(32,0.01119379);
   nominal__1493->SetMinimum(5e-05);
   nominal__1493->SetMaximum(104.5835);
   nominal__1493->SetEntries(19.0137);
   nominal__1493->SetFillColor(38);
   nominal__1493->SetLineColor(38);
   nominal__1493->SetMarkerColor(38);
   nominal__1493->SetMarkerStyle(21);
   nominal__1493->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1493->GetXaxis()->SetRange(1,32);
   nominal__1493->GetXaxis()->SetLabelFont(42);
   nominal__1493->GetXaxis()->SetLabelSize(0.035);
   nominal__1493->GetXaxis()->SetTitleSize(0.035);
   nominal__1493->GetXaxis()->SetTitleFont(42);
   nominal__1493->GetYaxis()->SetTitle("Tracks");
   nominal__1493->GetYaxis()->SetLabelFont(42);
   nominal__1493->GetYaxis()->SetLabelSize(0.05);
   nominal__1493->GetYaxis()->SetTitleSize(0.07);
   nominal__1493->GetYaxis()->SetTitleOffset(0);
   nominal__1493->GetYaxis()->SetTitleFont(42);
   nominal__1493->GetZaxis()->SetLabelFont(42);
   nominal__1493->GetZaxis()->SetLabelSize(0.035);
   nominal__1493->GetZaxis()->SetTitleSize(0.035);
   nominal__1493->GetZaxis()->SetTitleFont(42);
   nominal__1493->Draw("E0 same");
   Double_t xAxis1162[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1494 = new TH1F("nominal__1494","",32, xAxis1162);
   nominal__1494->SetBinContent(5,1);
   nominal__1494->SetBinContent(9,1);
   nominal__1494->SetBinContent(12,1);
   nominal__1494->SetBinContent(15,1);
   nominal__1494->SetBinContent(16,1);
   nominal__1494->SetBinContent(17,1);
   nominal__1494->SetBinContent(19,1);
   nominal__1494->SetBinContent(23,1);
   nominal__1494->SetBinContent(24,1);
   nominal__1494->SetBinContent(25,1);
   nominal__1494->SetBinContent(26,1);
   nominal__1494->SetBinContent(27,1);
   nominal__1494->SetBinContent(28,1);
   nominal__1494->SetBinContent(29,1);
   nominal__1494->SetBinContent(30,1);
   nominal__1494->SetBinContent(31,1);
   nominal__1494->SetBinContent(32,1);
   nominal__1494->SetBinError(5,1.414214);
   nominal__1494->SetBinError(9,1.414214);
   nominal__1494->SetBinError(12,1.414214);
   nominal__1494->SetBinError(15,1.000082);
   nominal__1494->SetBinError(16,1.414214);
   nominal__1494->SetBinError(17,1.414214);
   nominal__1494->SetBinError(19,1.414214);
   nominal__1494->SetBinError(23,0.8178533);
   nominal__1494->SetBinError(24,1.000002);
   nominal__1494->SetBinError(25,0.6326876);
   nominal__1494->SetBinError(26,0.3018806);
   nominal__1494->SetBinError(27,0.1814834);
   nominal__1494->SetBinError(28,0.06841188);
   nominal__1494->SetBinError(29,0.03006418);
   nominal__1494->SetBinError(30,0.01607051);
   nominal__1494->SetBinError(31,0.01356753);
   nominal__1494->SetBinError(32,0.01119379);
   nominal__1494->SetMinimum(5e-05);
   nominal__1494->SetMaximum(104.5835);
   nominal__1494->SetEntries(19.0137);
   nominal__1494->SetFillColor(46);
   nominal__1494->SetLineColor(46);
   nominal__1494->SetMarkerColor(46);
   nominal__1494->SetMarkerStyle(21);
   nominal__1494->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1494->GetXaxis()->SetRange(1,32);
   nominal__1494->GetXaxis()->SetLabelFont(42);
   nominal__1494->GetXaxis()->SetLabelSize(0.035);
   nominal__1494->GetXaxis()->SetTitleSize(0.035);
   nominal__1494->GetXaxis()->SetTitleFont(42);
   nominal__1494->GetYaxis()->SetTitle("Tracks");
   nominal__1494->GetYaxis()->SetLabelFont(42);
   nominal__1494->GetYaxis()->SetLabelSize(0.05);
   nominal__1494->GetYaxis()->SetTitleSize(0.07);
   nominal__1494->GetYaxis()->SetTitleOffset(0);
   nominal__1494->GetYaxis()->SetTitleFont(42);
   nominal__1494->GetZaxis()->SetLabelFont(42);
   nominal__1494->GetZaxis()->SetLabelSize(0.035);
   nominal__1494->GetZaxis()->SetTitleSize(0.035);
   nominal__1494->GetZaxis()->SetTitleFont(42);
   nominal__1494->Draw("E0 same");
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
