void plot_K1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:00 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.359683,2428.571,7.370992);
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
   Double_t xAxis162[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__208 = new TH1F("nominal__208","",32, xAxis162);
   nominal__208->SetBinContent(17,0.02156111);
   nominal__208->SetBinContent(22,0.04150828);
   nominal__208->SetBinContent(24,0.02087854);
   nominal__208->SetBinContent(25,0.05551993);
   nominal__208->SetBinContent(26,0.04206394);
   nominal__208->SetBinContent(27,0.1826455);
   nominal__208->SetBinContent(28,1.447503);
   nominal__208->SetBinContent(29,4.932075);
   nominal__208->SetBinContent(30,15.40064);
   nominal__208->SetBinContent(31,20.52763);
   nominal__208->SetBinContent(32,17.68969);
   nominal__208->SetBinError(17,0.02156111);
   nominal__208->SetBinError(22,0.02938644);
   nominal__208->SetBinError(24,0.02087854);
   nominal__208->SetBinError(25,0.03209129);
   nominal__208->SetBinError(26,0.02980124);
   nominal__208->SetBinError(27,0.06093327);
   nominal__208->SetBinError(28,0.1708915);
   nominal__208->SetBinError(29,0.3163469);
   nominal__208->SetBinError(30,0.557835);
   nominal__208->SetBinError(31,0.64459);
   nominal__208->SetBinError(32,0.4431387);
   nominal__208->SetBinError(33,0.4002011);
   nominal__208->SetMinimum(5e-05);
   nominal__208->SetMaximum(2.052763e+07);
   nominal__208->SetEntries(3002);
   nominal__208->SetFillColor(1);
   nominal__208->SetMarkerStyle(20);
   nominal__208->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__208->GetXaxis()->SetRange(1,32);
   nominal__208->GetXaxis()->SetLabelFont(42);
   nominal__208->GetXaxis()->SetLabelSize(0.035);
   nominal__208->GetXaxis()->SetTitleSize(0.035);
   nominal__208->GetXaxis()->SetTitleFont(42);
   nominal__208->GetYaxis()->SetTitle("Tracks");
   nominal__208->GetYaxis()->SetLabelFont(42);
   nominal__208->GetYaxis()->SetLabelSize(0.05);
   nominal__208->GetYaxis()->SetTitleSize(0.07);
   nominal__208->GetYaxis()->SetTitleOffset(0);
   nominal__208->GetYaxis()->SetTitleFont(42);
   nominal__208->GetZaxis()->SetLabelFont(42);
   nominal__208->GetZaxis()->SetLabelSize(0.035);
   nominal__208->GetZaxis()->SetTitleSize(0.035);
   nominal__208->GetZaxis()->SetTitleFont(42);
   nominal__208->Draw("");
   Double_t xAxis163[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__209 = new TH1F("K_down1__209","",32, xAxis163);
   K_down1__209->SetBinContent(17,0.02156111);
   K_down1__209->SetBinContent(22,0.04150828);
   K_down1__209->SetBinContent(24,0.02087854);
   K_down1__209->SetBinContent(25,0.05551993);
   K_down1__209->SetBinContent(26,0.04206394);
   K_down1__209->SetBinContent(27,0.1826455);
   K_down1__209->SetBinContent(28,1.428452);
   K_down1__209->SetBinContent(29,4.910481);
   K_down1__209->SetBinContent(30,15.18106);
   K_down1__209->SetBinContent(31,20.50463);
   K_down1__209->SetBinContent(32,17.97292);
   K_down1__209->SetBinError(17,0.02156111);
   K_down1__209->SetBinError(22,0.02938644);
   K_down1__209->SetBinError(24,0.02087854);
   K_down1__209->SetBinError(25,0.03209129);
   K_down1__209->SetBinError(26,0.02980124);
   K_down1__209->SetBinError(27,0.06093327);
   K_down1__209->SetBinError(28,0.1698262);
   K_down1__209->SetBinError(29,0.315614);
   K_down1__209->SetBinError(30,0.5538787);
   K_down1__209->SetBinError(31,0.6441938);
   K_down1__209->SetBinError(32,0.4469296);
   K_down1__209->SetBinError(33,0.4031398);
   K_down1__209->SetEntries(3002);
   K_down1__209->SetFillColor(38);
   K_down1__209->SetLineColor(38);
   K_down1__209->SetMarkerColor(38);
   K_down1__209->SetMarkerStyle(21);
   K_down1__209->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__209->GetXaxis()->SetRange(1,400);
   K_down1__209->GetXaxis()->SetLabelFont(42);
   K_down1__209->GetXaxis()->SetLabelSize(0.035);
   K_down1__209->GetXaxis()->SetTitleSize(0.035);
   K_down1__209->GetXaxis()->SetTitleFont(42);
   K_down1__209->GetYaxis()->SetTitle("Events / bin");
   K_down1__209->GetYaxis()->SetLabelFont(42);
   K_down1__209->GetYaxis()->SetLabelSize(0.035);
   K_down1__209->GetYaxis()->SetTitleSize(0.035);
   K_down1__209->GetYaxis()->SetTitleOffset(0);
   K_down1__209->GetYaxis()->SetTitleFont(42);
   K_down1__209->GetZaxis()->SetLabelFont(42);
   K_down1__209->GetZaxis()->SetLabelSize(0.035);
   K_down1__209->GetZaxis()->SetTitleSize(0.035);
   K_down1__209->GetZaxis()->SetTitleFont(42);
   K_down1__209->Draw("same");
   Double_t xAxis164[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__210 = new TH1F("K_up1__210","",32, xAxis164);
   K_up1__210->SetBinContent(17,0.02156111);
   K_up1__210->SetBinContent(22,0.04150828);
   K_up1__210->SetBinContent(24,0.04063199);
   K_up1__210->SetBinContent(25,0.03576648);
   K_up1__210->SetBinContent(26,0.04206394);
   K_up1__210->SetBinContent(27,0.1826455);
   K_up1__210->SetBinContent(28,1.488366);
   K_up1__210->SetBinContent(29,5.047035);
   K_up1__210->SetBinContent(30,15.45234);
   K_up1__210->SetBinContent(31,20.5787);
   K_up1__210->SetBinContent(32,17.4311);
   K_up1__210->SetBinError(17,0.02156111);
   K_up1__210->SetBinError(22,0.02938644);
   K_up1__210->SetBinError(24,0.02874216);
   K_up1__210->SetBinError(25,0.02529135);
   K_up1__210->SetBinError(26,0.02980124);
   K_up1__210->SetBinError(27,0.06093327);
   K_up1__210->SetBinError(28,0.1733175);
   K_up1__210->SetBinError(29,0.3198274);
   K_up1__210->SetBinError(30,0.5589785);
   K_up1__210->SetBinError(31,0.6452566);
   K_up1__210->SetBinError(32,0.4397972);
   K_up1__210->SetBinError(33,0.3973979);
   K_up1__210->SetEntries(3002);
   K_up1__210->SetFillColor(46);
   K_up1__210->SetLineColor(46);
   K_up1__210->SetMarkerColor(46);
   K_up1__210->SetMarkerStyle(21);
   K_up1__210->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__210->GetXaxis()->SetRange(1,400);
   K_up1__210->GetXaxis()->SetLabelFont(42);
   K_up1__210->GetXaxis()->SetLabelSize(0.035);
   K_up1__210->GetXaxis()->SetTitleSize(0.035);
   K_up1__210->GetXaxis()->SetTitleFont(42);
   K_up1__210->GetYaxis()->SetTitle("Events / bin");
   K_up1__210->GetYaxis()->SetLabelFont(42);
   K_up1__210->GetYaxis()->SetLabelSize(0.035);
   K_up1__210->GetYaxis()->SetTitleSize(0.035);
   K_up1__210->GetYaxis()->SetTitleOffset(0);
   K_up1__210->GetYaxis()->SetTitleFont(42);
   K_up1__210->GetZaxis()->SetLabelFont(42);
   K_up1__210->GetZaxis()->SetLabelSize(0.035);
   K_up1__210->GetZaxis()->SetTitleSize(0.035);
   K_up1__210->GetZaxis()->SetTitleFont(42);
   K_up1__210->Draw("same");
   TLine *line = new TLine(1730,0,1730,2.052763e+07);
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
   
   TH1D *frameR2__211 = new TH1D("frameR2__211","",1,0,2000);
   frameR2__211->SetMinimum(0.5);
   frameR2__211->SetMaximum(1.5);
   frameR2__211->SetStats(0);
   frameR2__211->SetLineStyle(0);
   frameR2__211->SetMarkerStyle(20);
   frameR2__211->GetXaxis()->SetRange(1,1);
   frameR2__211->GetXaxis()->SetLabelFont(43);
   frameR2__211->GetXaxis()->SetLabelOffset(0.007);
   frameR2__211->GetXaxis()->SetLabelSize(16);
   frameR2__211->GetXaxis()->SetTitleSize(24);
   frameR2__211->GetXaxis()->SetTitleOffset(3.75);
   frameR2__211->GetXaxis()->SetTitleFont(43);
   frameR2__211->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__211->GetYaxis()->SetNdivisions(205);
   frameR2__211->GetYaxis()->SetLabelFont(43);
   frameR2__211->GetYaxis()->SetLabelOffset(0.007);
   frameR2__211->GetYaxis()->SetLabelSize(20);
   frameR2__211->GetYaxis()->SetTitleSize(20);
   frameR2__211->GetYaxis()->SetTitleOffset(2);
   frameR2__211->GetYaxis()->SetTitleFont(43);
   frameR2__211->GetZaxis()->SetLabelFont(42);
   frameR2__211->GetZaxis()->SetLabelOffset(0.007);
   frameR2__211->GetZaxis()->SetLabelSize(0.05);
   frameR2__211->GetZaxis()->SetTitleSize(0.06);
   frameR2__211->GetZaxis()->SetTitleFont(42);
   frameR2__211->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis165[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down1__212 = new TH1F("K_down1__212","",32, xAxis165);
   K_down1__212->SetBinContent(0,1);
   K_down1__212->SetBinContent(1,1);
   K_down1__212->SetBinContent(2,1);
   K_down1__212->SetBinContent(3,1);
   K_down1__212->SetBinContent(4,1);
   K_down1__212->SetBinContent(5,1);
   K_down1__212->SetBinContent(6,1);
   K_down1__212->SetBinContent(7,1);
   K_down1__212->SetBinContent(8,1);
   K_down1__212->SetBinContent(9,1);
   K_down1__212->SetBinContent(10,1);
   K_down1__212->SetBinContent(11,1);
   K_down1__212->SetBinContent(12,1);
   K_down1__212->SetBinContent(13,1);
   K_down1__212->SetBinContent(14,1);
   K_down1__212->SetBinContent(15,1);
   K_down1__212->SetBinContent(16,1);
   K_down1__212->SetBinContent(17,1);
   K_down1__212->SetBinContent(18,1);
   K_down1__212->SetBinContent(19,1);
   K_down1__212->SetBinContent(20,1);
   K_down1__212->SetBinContent(21,1);
   K_down1__212->SetBinContent(22,1);
   K_down1__212->SetBinContent(23,1);
   K_down1__212->SetBinContent(24,1);
   K_down1__212->SetBinContent(25,1);
   K_down1__212->SetBinContent(26,1);
   K_down1__212->SetBinContent(27,1);
   K_down1__212->SetBinContent(28,1);
   K_down1__212->SetBinContent(29,0.9996747);
   K_down1__212->SetBinContent(30,0.9992425);
   K_down1__212->SetBinContent(31,0.993237);
   K_down1__212->SetBinContent(32,0.9842416);
   K_down1__212->SetBinError(0,0.02587692);
   K_down1__212->SetBinError(1,0.02587692);
   K_down1__212->SetBinError(2,0.02587692);
   K_down1__212->SetBinError(3,0.02587692);
   K_down1__212->SetBinError(4,0.02587692);
   K_down1__212->SetBinError(5,0.02587692);
   K_down1__212->SetBinError(6,0.02587692);
   K_down1__212->SetBinError(7,0.02587692);
   K_down1__212->SetBinError(8,0.02587692);
   K_down1__212->SetBinError(9,0.02587692);
   K_down1__212->SetBinError(10,0.02587692);
   K_down1__212->SetBinError(11,0.02587692);
   K_down1__212->SetBinError(12,0.02587692);
   K_down1__212->SetBinError(13,0.02587692);
   K_down1__212->SetBinError(14,0.02587692);
   K_down1__212->SetBinError(15,0.02587692);
   K_down1__212->SetBinError(16,0.02587692);
   K_down1__212->SetBinError(17,0.02587692);
   K_down1__212->SetBinError(18,0.02588123);
   K_down1__212->SetBinError(19,0.02588123);
   K_down1__212->SetBinError(20,0.02588123);
   K_down1__212->SetBinError(21,0.02588123);
   K_down1__212->SetBinError(22,0.02588123);
   K_down1__212->SetBinError(23,0.02588988);
   K_down1__212->SetBinError(24,0.02588988);
   K_down1__212->SetBinError(25,0.02589421);
   K_down1__212->SetBinError(26,0.02590713);
   K_down1__212->SetBinError(27,0.02591577);
   K_down1__212->SetBinError(28,0.02595496);
   K_down1__212->SetBinError(29,0.02626357);
   K_down1__212->SetBinError(30,0.02742521);
   K_down1__212->SetBinError(31,0.03223921);
   K_down1__212->SetBinError(32,0.04679898);
   K_down1__212->SetEntries(33);
   K_down1__212->SetFillColor(38);
   K_down1__212->SetLineColor(38);
   K_down1__212->SetMarkerColor(38);
   K_down1__212->SetMarkerStyle(21);
   K_down1__212->GetXaxis()->SetTitle("Mass [GeV]");
   K_down1__212->GetXaxis()->SetRange(1,400);
   K_down1__212->GetXaxis()->SetLabelFont(42);
   K_down1__212->GetXaxis()->SetLabelSize(0.035);
   K_down1__212->GetXaxis()->SetTitleSize(0.035);
   K_down1__212->GetXaxis()->SetTitleFont(42);
   K_down1__212->GetYaxis()->SetTitle("Events / bin");
   K_down1__212->GetYaxis()->SetLabelFont(42);
   K_down1__212->GetYaxis()->SetLabelSize(0.035);
   K_down1__212->GetYaxis()->SetTitleSize(0.035);
   K_down1__212->GetYaxis()->SetTitleOffset(0);
   K_down1__212->GetYaxis()->SetTitleFont(42);
   K_down1__212->GetZaxis()->SetLabelFont(42);
   K_down1__212->GetZaxis()->SetLabelSize(0.035);
   K_down1__212->GetZaxis()->SetTitleSize(0.035);
   K_down1__212->GetZaxis()->SetTitleFont(42);
   K_down1__212->Draw("E0 same");
   Double_t xAxis166[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up1__213 = new TH1F("K_up1__213","",32, xAxis166);
   K_up1__213->SetBinContent(0,1);
   K_up1__213->SetBinContent(1,1);
   K_up1__213->SetBinContent(2,1);
   K_up1__213->SetBinContent(3,1);
   K_up1__213->SetBinContent(4,1);
   K_up1__213->SetBinContent(5,1);
   K_up1__213->SetBinContent(6,1);
   K_up1__213->SetBinContent(7,1);
   K_up1__213->SetBinContent(8,1);
   K_up1__213->SetBinContent(9,1);
   K_up1__213->SetBinContent(10,1);
   K_up1__213->SetBinContent(11,1);
   K_up1__213->SetBinContent(12,1);
   K_up1__213->SetBinContent(13,1);
   K_up1__213->SetBinContent(14,1);
   K_up1__213->SetBinContent(15,1);
   K_up1__213->SetBinContent(16,1);
   K_up1__213->SetBinContent(17,1);
   K_up1__213->SetBinContent(18,1);
   K_up1__213->SetBinContent(19,1);
   K_up1__213->SetBinContent(20,1);
   K_up1__213->SetBinContent(21,1);
   K_up1__213->SetBinContent(22,1);
   K_up1__213->SetBinContent(23,1);
   K_up1__213->SetBinContent(24,1);
   K_up1__213->SetBinContent(25,1.000328);
   K_up1__213->SetBinContent(26,1);
   K_up1__213->SetBinContent(27,1);
   K_up1__213->SetBinContent(28,1);
   K_up1__213->SetBinContent(29,1.000698);
   K_up1__213->SetBinContent(30,1.002915);
   K_up1__213->SetBinContent(31,1.00546);
   K_up1__213->SetBinContent(32,1.014835);
   K_up1__213->SetBinError(0,0.02587692);
   K_up1__213->SetBinError(1,0.02587692);
   K_up1__213->SetBinError(2,0.02587692);
   K_up1__213->SetBinError(3,0.02587692);
   K_up1__213->SetBinError(4,0.02587692);
   K_up1__213->SetBinError(5,0.02587692);
   K_up1__213->SetBinError(6,0.02587692);
   K_up1__213->SetBinError(7,0.02587692);
   K_up1__213->SetBinError(8,0.02587692);
   K_up1__213->SetBinError(9,0.02587692);
   K_up1__213->SetBinError(10,0.02587692);
   K_up1__213->SetBinError(11,0.02587692);
   K_up1__213->SetBinError(12,0.02587692);
   K_up1__213->SetBinError(13,0.02587692);
   K_up1__213->SetBinError(14,0.02587692);
   K_up1__213->SetBinError(15,0.02587692);
   K_up1__213->SetBinError(16,0.02587692);
   K_up1__213->SetBinError(17,0.02587692);
   K_up1__213->SetBinError(18,0.02588124);
   K_up1__213->SetBinError(19,0.02588124);
   K_up1__213->SetBinError(20,0.02588124);
   K_up1__213->SetBinError(21,0.02588124);
   K_up1__213->SetBinError(22,0.02588124);
   K_up1__213->SetBinError(23,0.02588988);
   K_up1__213->SetBinError(24,0.02588988);
   K_up1__213->SetBinError(25,0.02590487);
   K_up1__213->SetBinError(26,0.02590713);
   K_up1__213->SetBinError(27,0.02591578);
   K_up1__213->SetBinError(28,0.02595496);
   K_up1__213->SetBinError(29,0.02629727);
   K_up1__213->SetBinError(30,0.02755186);
   K_up1__213->SetBinError(31,0.03273476);
   K_up1__213->SetBinError(32,0.04862423);
   K_up1__213->SetEntries(33);
   K_up1__213->SetFillColor(46);
   K_up1__213->SetLineColor(46);
   K_up1__213->SetMarkerColor(46);
   K_up1__213->SetMarkerStyle(21);
   K_up1__213->GetXaxis()->SetTitle("Mass [GeV]");
   K_up1__213->GetXaxis()->SetRange(1,400);
   K_up1__213->GetXaxis()->SetLabelFont(42);
   K_up1__213->GetXaxis()->SetLabelSize(0.035);
   K_up1__213->GetXaxis()->SetTitleSize(0.035);
   K_up1__213->GetXaxis()->SetTitleFont(42);
   K_up1__213->GetYaxis()->SetTitle("Events / bin");
   K_up1__213->GetYaxis()->SetLabelFont(42);
   K_up1__213->GetYaxis()->SetLabelSize(0.035);
   K_up1__213->GetYaxis()->SetTitleSize(0.035);
   K_up1__213->GetYaxis()->SetTitleOffset(0);
   K_up1__213->GetYaxis()->SetTitleFont(42);
   K_up1__213->GetZaxis()->SetLabelFont(42);
   K_up1__213->GetZaxis()->SetLabelSize(0.035);
   K_up1__213->GetZaxis()->SetTitleSize(0.035);
   K_up1__213->GetZaxis()->SetTitleFont(42);
   K_up1__213->Draw("E0 same");
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
   
   TH1D *frameR2__214 = new TH1D("frameR2__214","",1,0,2000);
   frameR2__214->SetMinimum(0.5);
   frameR2__214->SetMaximum(1.5);
   frameR2__214->SetStats(0);
   frameR2__214->SetLineStyle(0);
   frameR2__214->SetMarkerStyle(20);
   frameR2__214->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__214->GetXaxis()->SetRange(1,1);
   frameR2__214->GetXaxis()->SetLabelFont(43);
   frameR2__214->GetXaxis()->SetLabelOffset(0.007);
   frameR2__214->GetXaxis()->SetLabelSize(16);
   frameR2__214->GetXaxis()->SetTitleSize(24);
   frameR2__214->GetXaxis()->SetTitleOffset(5);
   frameR2__214->GetXaxis()->SetTitleFont(43);
   frameR2__214->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__214->GetYaxis()->SetNdivisions(205);
   frameR2__214->GetYaxis()->SetLabelFont(43);
   frameR2__214->GetYaxis()->SetLabelOffset(0.007);
   frameR2__214->GetYaxis()->SetLabelSize(20);
   frameR2__214->GetYaxis()->SetTitleSize(20);
   frameR2__214->GetYaxis()->SetTitleOffset(2);
   frameR2__214->GetYaxis()->SetTitleFont(43);
   frameR2__214->GetZaxis()->SetLabelFont(42);
   frameR2__214->GetZaxis()->SetLabelOffset(0.007);
   frameR2__214->GetZaxis()->SetLabelSize(0.05);
   frameR2__214->GetZaxis()->SetTitleSize(0.06);
   frameR2__214->GetZaxis()->SetTitleFont(42);
   frameR2__214->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis167[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__215 = new TH1F("nominal__215","",32, xAxis167);
   nominal__215->SetBinContent(17,1);
   nominal__215->SetBinContent(22,1);
   nominal__215->SetBinContent(24,1);
   nominal__215->SetBinContent(25,1);
   nominal__215->SetBinContent(26,1);
   nominal__215->SetBinContent(27,1);
   nominal__215->SetBinContent(28,1.013337);
   nominal__215->SetBinContent(29,1.004397);
   nominal__215->SetBinContent(30,1.014464);
   nominal__215->SetBinContent(31,1.001122);
   nominal__215->SetBinContent(32,0.9842416);
   nominal__215->SetBinError(17,1.414214);
   nominal__215->SetBinError(22,1.001215);
   nominal__215->SetBinError(24,1.414214);
   nominal__215->SetBinError(25,0.8174352);
   nominal__215->SetBinError(26,1.001935);
   nominal__215->SetBinError(27,0.4718029);
   nominal__215->SetBinError(28,0.1697831);
   nominal__215->SetBinError(29,0.09120196);
   nominal__215->SetBinError(30,0.05215514);
   nominal__215->SetBinError(31,0.04446893);
   nominal__215->SetBinError(32,0.03474102);
   nominal__215->SetMinimum(5e-05);
   nominal__215->SetMaximum(2.052763e+07);
   nominal__215->SetEntries(17.49049);
   nominal__215->SetFillColor(38);
   nominal__215->SetLineColor(38);
   nominal__215->SetMarkerColor(38);
   nominal__215->SetMarkerStyle(21);
   nominal__215->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__215->GetXaxis()->SetRange(1,32);
   nominal__215->GetXaxis()->SetLabelFont(42);
   nominal__215->GetXaxis()->SetLabelSize(0.035);
   nominal__215->GetXaxis()->SetTitleSize(0.035);
   nominal__215->GetXaxis()->SetTitleFont(42);
   nominal__215->GetYaxis()->SetTitle("Tracks");
   nominal__215->GetYaxis()->SetLabelFont(42);
   nominal__215->GetYaxis()->SetLabelSize(0.05);
   nominal__215->GetYaxis()->SetTitleSize(0.07);
   nominal__215->GetYaxis()->SetTitleOffset(0);
   nominal__215->GetYaxis()->SetTitleFont(42);
   nominal__215->GetZaxis()->SetLabelFont(42);
   nominal__215->GetZaxis()->SetLabelSize(0.035);
   nominal__215->GetZaxis()->SetTitleSize(0.035);
   nominal__215->GetZaxis()->SetTitleFont(42);
   nominal__215->Draw("E0 same");
   Double_t xAxis168[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__216 = new TH1F("nominal__216","",32, xAxis168);
   nominal__216->SetBinContent(17,1);
   nominal__216->SetBinContent(22,1);
   nominal__216->SetBinContent(24,0.5138448);
   nominal__216->SetBinContent(25,1.552289);
   nominal__216->SetBinContent(26,1);
   nominal__216->SetBinContent(27,1);
   nominal__216->SetBinContent(28,0.972545);
   nominal__216->SetBinContent(29,0.9772221);
   nominal__216->SetBinContent(30,0.9966547);
   nominal__216->SetBinContent(31,0.9975181);
   nominal__216->SetBinContent(32,1.014835);
   nominal__216->SetBinError(17,1.414214);
   nominal__216->SetBinError(22,1.001215);
   nominal__216->SetBinError(24,0.6294092);
   nominal__216->SetBinError(25,1.417713);
   nominal__216->SetBinError(26,1.001935);
   nominal__216->SetBinError(27,0.4718029);
   nominal__216->SetBinError(28,0.1612731);
   nominal__216->SetBinError(29,0.08811115);
   nominal__216->SetBinError(30,0.0510204);
   nominal__216->SetBinError(31,0.04426553);
   nominal__216->SetBinError(32,0.03608192);
   nominal__216->SetMinimum(5e-05);
   nominal__216->SetMaximum(2.052763e+07);
   nominal__216->SetEntries(18.21061);
   nominal__216->SetFillColor(46);
   nominal__216->SetLineColor(46);
   nominal__216->SetMarkerColor(46);
   nominal__216->SetMarkerStyle(21);
   nominal__216->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__216->GetXaxis()->SetRange(1,32);
   nominal__216->GetXaxis()->SetLabelFont(42);
   nominal__216->GetXaxis()->SetLabelSize(0.035);
   nominal__216->GetXaxis()->SetTitleSize(0.035);
   nominal__216->GetXaxis()->SetTitleFont(42);
   nominal__216->GetYaxis()->SetTitle("Tracks");
   nominal__216->GetYaxis()->SetLabelFont(42);
   nominal__216->GetYaxis()->SetLabelSize(0.05);
   nominal__216->GetYaxis()->SetTitleSize(0.07);
   nominal__216->GetYaxis()->SetTitleOffset(0);
   nominal__216->GetYaxis()->SetTitleFont(42);
   nominal__216->GetZaxis()->SetLabelFont(42);
   nominal__216->GetZaxis()->SetLabelSize(0.035);
   nominal__216->GetZaxis()->SetTitleSize(0.035);
   nominal__216->GetZaxis()->SetTitleFont(42);
   nominal__216->Draw("E0 same");
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
