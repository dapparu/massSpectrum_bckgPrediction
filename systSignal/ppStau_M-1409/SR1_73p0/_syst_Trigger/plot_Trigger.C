void plot_Trigger()
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
   t1->Range(-428.5714,-4.338199,2428.571,3.095511);
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
   Double_t xAxis1100[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1414 = new TH1F("nominal__1414","",32, xAxis1100);
   nominal__1414->SetBinContent(4,1.133859e-06);
   nominal__1414->SetBinContent(12,1.172288e-06);
   nominal__1414->SetBinContent(14,1.076732e-06);
   nominal__1414->SetBinContent(16,1.187564e-06);
   nominal__1414->SetBinContent(17,2.410127e-06);
   nominal__1414->SetBinContent(18,5.643037e-07);
   nominal__1414->SetBinContent(19,5.450868e-07);
   nominal__1414->SetBinContent(20,5.651711e-07);
   nominal__1414->SetBinContent(21,5.778977e-07);
   nominal__1414->SetBinContent(22,1.810394e-06);
   nominal__1414->SetBinContent(23,2.834062e-06);
   nominal__1414->SetBinContent(24,2.286669e-06);
   nominal__1414->SetBinContent(25,1.175712e-05);
   nominal__1414->SetBinContent(26,4.076618e-05);
   nominal__1414->SetBinContent(27,0.0002393376);
   nominal__1414->SetBinContent(28,0.00138461);
   nominal__1414->SetBinContent(29,0.005642689);
   nominal__1414->SetBinContent(30,0.0114378);
   nominal__1414->SetBinContent(31,0.008162949);
   nominal__1414->SetBinContent(32,0.003875451);
   nominal__1414->SetBinError(4,8.017631e-07);
   nominal__1414->SetBinError(12,8.290847e-07);
   nominal__1414->SetBinError(14,7.625711e-07);
   nominal__1414->SetBinError(16,8.403703e-07);
   nominal__1414->SetBinError(17,1.206807e-06);
   nominal__1414->SetBinError(18,5.643037e-07);
   nominal__1414->SetBinError(19,5.450868e-07);
   nominal__1414->SetBinError(20,5.651712e-07);
   nominal__1414->SetBinError(21,5.778977e-07);
   nominal__1414->SetBinError(22,1.046944e-06);
   nominal__1414->SetBinError(23,1.267576e-06);
   nominal__1414->SetBinError(24,1.143746e-06);
   nominal__1414->SetBinError(25,2.631857e-06);
   nominal__1414->SetBinError(26,4.849418e-06);
   nominal__1414->SetBinError(27,1.173014e-05);
   nominal__1414->SetBinError(28,2.826415e-05);
   nominal__1414->SetBinError(29,5.706391e-05);
   nominal__1414->SetBinError(30,8.126949e-05);
   nominal__1414->SetBinError(31,6.867795e-05);
   nominal__1414->SetBinError(32,3.765616e-05);
   nominal__1414->SetBinError(33,2.867408e-05);
   nominal__1414->SetMinimum(5e-05);
   nominal__1414->SetMaximum(1143.78);
   nominal__1414->SetEntries(53566);
   nominal__1414->SetFillColor(1);
   nominal__1414->SetMarkerStyle(20);
   nominal__1414->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1414->GetXaxis()->SetRange(1,32);
   nominal__1414->GetXaxis()->SetLabelFont(42);
   nominal__1414->GetXaxis()->SetLabelSize(0.035);
   nominal__1414->GetXaxis()->SetTitleSize(0.035);
   nominal__1414->GetXaxis()->SetTitleFont(42);
   nominal__1414->GetYaxis()->SetTitle("Tracks");
   nominal__1414->GetYaxis()->SetLabelFont(42);
   nominal__1414->GetYaxis()->SetLabelSize(0.05);
   nominal__1414->GetYaxis()->SetTitleSize(0.07);
   nominal__1414->GetYaxis()->SetTitleOffset(0);
   nominal__1414->GetYaxis()->SetTitleFont(42);
   nominal__1414->GetZaxis()->SetLabelFont(42);
   nominal__1414->GetZaxis()->SetLabelSize(0.035);
   nominal__1414->GetZaxis()->SetTitleSize(0.035);
   nominal__1414->GetZaxis()->SetTitleFont(42);
   nominal__1414->Draw("");
   Double_t xAxis1101[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__1415 = new TH1F("Trigger_down__1415","",32, xAxis1101);
   Trigger_down__1415->SetBinContent(4,8.448831e-07);
   Trigger_down__1415->SetBinContent(12,1.017647e-06);
   Trigger_down__1415->SetBinContent(14,9.928707e-07);
   Trigger_down__1415->SetBinContent(16,9.924699e-07);
   Trigger_down__1415->SetBinContent(17,2.062383e-06);
   Trigger_down__1415->SetBinContent(18,5.360885e-07);
   Trigger_down__1415->SetBinContent(19,5.069308e-07);
   Trigger_down__1415->SetBinContent(20,5.369126e-07);
   Trigger_down__1415->SetBinContent(21,5.201079e-07);
   Trigger_down__1415->SetBinContent(22,1.617161e-06);
   Trigger_down__1415->SetBinContent(23,2.736883e-06);
   Trigger_down__1415->SetBinContent(24,2.182959e-06);
   Trigger_down__1415->SetBinContent(25,1.065884e-05);
   Trigger_down__1415->SetBinContent(26,3.871518e-05);
   Trigger_down__1415->SetBinContent(27,0.0002249033);
   Trigger_down__1415->SetBinContent(28,0.001265271);
   Trigger_down__1415->SetBinContent(29,0.004963682);
   Trigger_down__1415->SetBinContent(30,0.009791534);
   Trigger_down__1415->SetBinContent(31,0.006968718);
   Trigger_down__1415->SetBinContent(32,0.003379913);
   Trigger_down__1415->SetBinError(4,5.986255e-07);
   Trigger_down__1415->SetBinError(12,7.228947e-07);
   Trigger_down__1415->SetBinError(14,7.05228e-07);
   Trigger_down__1415->SetBinError(16,7.047785e-07);
   Trigger_down__1415->SetBinError(17,1.040177e-06);
   Trigger_down__1415->SetBinError(18,5.360885e-07);
   Trigger_down__1415->SetBinError(19,5.069308e-07);
   Trigger_down__1415->SetBinError(20,5.369126e-07);
   Trigger_down__1415->SetBinError(21,5.201079e-07);
   Trigger_down__1415->SetBinError(22,9.392593e-07);
   Trigger_down__1415->SetBinError(23,1.229265e-06);
   Trigger_down__1415->SetBinError(24,1.100633e-06);
   Trigger_down__1415->SetBinError(25,2.404522e-06);
   Trigger_down__1415->SetBinError(26,4.625778e-06);
   Trigger_down__1415->SetBinError(27,1.108953e-05);
   Trigger_down__1415->SetBinError(28,2.600124e-05);
   Trigger_down__1415->SetBinError(29,5.054822e-05);
   Trigger_down__1415->SetBinError(30,7.005622e-05);
   Trigger_down__1415->SetBinError(31,5.900572e-05);
   Trigger_down__1415->SetBinError(32,3.28036e-05);
   Trigger_down__1415->SetBinError(33,2.546178e-05);
   Trigger_down__1415->SetEntries(53566);
   Trigger_down__1415->SetFillColor(38);
   Trigger_down__1415->SetLineColor(38);
   Trigger_down__1415->SetMarkerColor(38);
   Trigger_down__1415->SetMarkerStyle(21);
   Trigger_down__1415->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__1415->GetXaxis()->SetRange(1,400);
   Trigger_down__1415->GetXaxis()->SetLabelFont(42);
   Trigger_down__1415->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__1415->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__1415->GetXaxis()->SetTitleFont(42);
   Trigger_down__1415->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__1415->GetYaxis()->SetLabelFont(42);
   Trigger_down__1415->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__1415->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__1415->GetYaxis()->SetTitleOffset(0);
   Trigger_down__1415->GetYaxis()->SetTitleFont(42);
   Trigger_down__1415->GetZaxis()->SetLabelFont(42);
   Trigger_down__1415->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__1415->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__1415->GetZaxis()->SetTitleFont(42);
   Trigger_down__1415->Draw("same");
   Double_t xAxis1102[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__1416 = new TH1F("Trigger_up__1416","",32, xAxis1102);
   Trigger_up__1416->SetBinContent(4,1.207542e-06);
   Trigger_up__1416->SetBinContent(12,1.143093e-06);
   Trigger_up__1416->SetBinContent(14,1.127307e-06);
   Trigger_up__1416->SetBinContent(16,1.157182e-06);
   Trigger_up__1416->SetBinContent(17,2.461391e-06);
   Trigger_up__1416->SetBinContent(18,5.078733e-07);
   Trigger_up__1416->SetBinContent(19,5.39636e-07);
   Trigger_up__1416->SetBinContent(20,5.08654e-07);
   Trigger_up__1416->SetBinContent(21,5.778977e-07);
   Trigger_up__1416->SetBinContent(22,1.878844e-06);
   Trigger_up__1416->SetBinContent(23,2.839254e-06);
   Trigger_up__1416->SetBinContent(24,2.263766e-06);
   Trigger_up__1416->SetBinContent(25,1.179543e-05);
   Trigger_up__1416->SetBinContent(26,4.049716e-05);
   Trigger_up__1416->SetBinContent(27,0.0002388607);
   Trigger_up__1416->SetBinContent(28,0.001378276);
   Trigger_up__1416->SetBinContent(29,0.005652313);
   Trigger_up__1416->SetBinContent(30,0.01146764);
   Trigger_up__1416->SetBinContent(31,0.008153494);
   Trigger_up__1416->SetBinContent(32,0.003850132);
   Trigger_up__1416->SetBinError(4,8.538623e-07);
   Trigger_up__1416->SetBinError(12,8.085268e-07);
   Trigger_up__1416->SetBinError(14,8.003233e-07);
   Trigger_up__1416->SetBinError(16,8.184778e-07);
   Trigger_up__1416->SetBinError(17,1.232454e-06);
   Trigger_up__1416->SetBinError(18,5.078733e-07);
   Trigger_up__1416->SetBinError(19,5.39636e-07);
   Trigger_up__1416->SetBinError(20,5.08654e-07);
   Trigger_up__1416->SetBinError(21,5.778977e-07);
   Trigger_up__1416->SetBinError(22,1.087933e-06);
   Trigger_up__1416->SetBinError(23,1.2709e-06);
   Trigger_up__1416->SetBinError(24,1.132424e-06);
   Trigger_up__1416->SetBinError(25,2.641828e-06);
   Trigger_up__1416->SetBinError(26,4.820964e-06);
   Trigger_up__1416->SetBinError(27,1.171444e-05);
   Trigger_up__1416->SetBinError(28,2.81523e-05);
   Trigger_up__1416->SetBinError(29,5.721988e-05);
   Trigger_up__1416->SetBinError(30,8.156727e-05);
   Trigger_up__1416->SetBinError(31,6.866319e-05);
   Trigger_up__1416->SetBinError(32,3.747025e-05);
   Trigger_up__1416->SetBinError(33,2.845966e-05);
   Trigger_up__1416->SetEntries(53566);
   Trigger_up__1416->SetFillColor(46);
   Trigger_up__1416->SetLineColor(46);
   Trigger_up__1416->SetMarkerColor(46);
   Trigger_up__1416->SetMarkerStyle(21);
   Trigger_up__1416->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__1416->GetXaxis()->SetRange(1,400);
   Trigger_up__1416->GetXaxis()->SetLabelFont(42);
   Trigger_up__1416->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__1416->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__1416->GetXaxis()->SetTitleFont(42);
   Trigger_up__1416->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__1416->GetYaxis()->SetLabelFont(42);
   Trigger_up__1416->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__1416->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__1416->GetYaxis()->SetTitleOffset(0);
   Trigger_up__1416->GetYaxis()->SetTitleFont(42);
   Trigger_up__1416->GetZaxis()->SetLabelFont(42);
   Trigger_up__1416->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__1416->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__1416->GetZaxis()->SetTitleFont(42);
   Trigger_up__1416->Draw("same");
   TLine *line = new TLine(1730,0,1730,1143.78);
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
   
   TH1D *frameR2__1417 = new TH1D("frameR2__1417","",1,0,2000);
   frameR2__1417->SetMinimum(0.5);
   frameR2__1417->SetMaximum(1.5);
   frameR2__1417->SetStats(0);
   frameR2__1417->SetLineStyle(0);
   frameR2__1417->SetMarkerStyle(20);
   frameR2__1417->GetXaxis()->SetRange(1,1);
   frameR2__1417->GetXaxis()->SetLabelFont(43);
   frameR2__1417->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1417->GetXaxis()->SetLabelSize(16);
   frameR2__1417->GetXaxis()->SetTitleSize(24);
   frameR2__1417->GetXaxis()->SetTitleOffset(3.75);
   frameR2__1417->GetXaxis()->SetTitleFont(43);
   frameR2__1417->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__1417->GetYaxis()->SetNdivisions(205);
   frameR2__1417->GetYaxis()->SetLabelFont(43);
   frameR2__1417->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1417->GetYaxis()->SetLabelSize(20);
   frameR2__1417->GetYaxis()->SetTitleSize(20);
   frameR2__1417->GetYaxis()->SetTitleOffset(2);
   frameR2__1417->GetYaxis()->SetTitleFont(43);
   frameR2__1417->GetZaxis()->SetLabelFont(42);
   frameR2__1417->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1417->GetZaxis()->SetLabelSize(0.05);
   frameR2__1417->GetZaxis()->SetTitleSize(0.06);
   frameR2__1417->GetZaxis()->SetTitleFont(42);
   frameR2__1417->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1103[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_down__1418 = new TH1F("Trigger_down__1418","",32, xAxis1103);
   Trigger_down__1418->SetBinContent(0,1.15581);
   Trigger_down__1418->SetBinContent(1,1.15581);
   Trigger_down__1418->SetBinContent(2,1.15581);
   Trigger_down__1418->SetBinContent(3,1.15581);
   Trigger_down__1418->SetBinContent(4,1.15581);
   Trigger_down__1418->SetBinContent(5,1.155805);
   Trigger_down__1418->SetBinContent(6,1.155805);
   Trigger_down__1418->SetBinContent(7,1.155805);
   Trigger_down__1418->SetBinContent(8,1.155805);
   Trigger_down__1418->SetBinContent(9,1.155805);
   Trigger_down__1418->SetBinContent(10,1.155805);
   Trigger_down__1418->SetBinContent(11,1.155805);
   Trigger_down__1418->SetBinContent(12,1.155805);
   Trigger_down__1418->SetBinContent(13,1.155805);
   Trigger_down__1418->SetBinContent(14,1.155805);
   Trigger_down__1418->SetBinContent(15,1.155807);
   Trigger_down__1418->SetBinContent(16,1.155807);
   Trigger_down__1418->SetBinContent(17,1.155806);
   Trigger_down__1418->SetBinContent(18,1.155805);
   Trigger_down__1418->SetBinContent(19,1.155807);
   Trigger_down__1418->SetBinContent(20,1.155808);
   Trigger_down__1418->SetBinContent(21,1.15581);
   Trigger_down__1418->SetBinContent(22,1.155811);
   Trigger_down__1418->SetBinContent(23,1.155814);
   Trigger_down__1418->SetBinContent(24,1.155826);
   Trigger_down__1418->SetBinContent(25,1.155835);
   Trigger_down__1418->SetBinContent(26,1.155856);
   Trigger_down__1418->SetBinContent(27,1.156006);
   Trigger_down__1418->SetBinContent(28,1.156789);
   Trigger_down__1418->SetBinContent(29,1.159938);
   Trigger_down__1418->SetBinContent(30,1.165641);
   Trigger_down__1418->SetBinContent(31,1.163284);
   Trigger_down__1418->SetBinContent(32,1.146613);
   Trigger_down__1418->SetBinError(0,0.007100892);
   Trigger_down__1418->SetBinError(1,0.007100892);
   Trigger_down__1418->SetBinError(2,0.007100892);
   Trigger_down__1418->SetBinError(3,0.007100892);
   Trigger_down__1418->SetBinError(4,0.007100892);
   Trigger_down__1418->SetBinError(5,0.007100988);
   Trigger_down__1418->SetBinError(6,0.007100988);
   Trigger_down__1418->SetBinError(7,0.007100988);
   Trigger_down__1418->SetBinError(8,0.007100988);
   Trigger_down__1418->SetBinError(9,0.007100988);
   Trigger_down__1418->SetBinError(10,0.007100988);
   Trigger_down__1418->SetBinError(11,0.007100988);
   Trigger_down__1418->SetBinError(12,0.007100988);
   Trigger_down__1418->SetBinError(13,0.007101122);
   Trigger_down__1418->SetBinError(14,0.007101122);
   Trigger_down__1418->SetBinError(15,0.007101272);
   Trigger_down__1418->SetBinError(16,0.007101272);
   Trigger_down__1418->SetBinError(17,0.007101395);
   Trigger_down__1418->SetBinError(18,0.007101654);
   Trigger_down__1418->SetBinError(19,0.007101734);
   Trigger_down__1418->SetBinError(20,0.00710181);
   Trigger_down__1418->SetBinError(21,0.00710189);
   Trigger_down__1418->SetBinError(22,0.007101962);
   Trigger_down__1418->SetBinError(23,0.007102175);
   Trigger_down__1418->SetBinError(24,0.007102583);
   Trigger_down__1418->SetBinError(25,0.007102902);
   Trigger_down__1418->SetBinError(26,0.007104357);
   Trigger_down__1418->SetBinError(27,0.007109983);
   Trigger_down__1418->SetBinError(28,0.007142716);
   Trigger_down__1418->SetBinError(29,0.007330366);
   Trigger_down__1418->SetBinError(30,0.008204397);
   Trigger_down__1418->SetBinError(31,0.0114343);
   Trigger_down__1418->SetBinError(32,0.01986328);
   Trigger_down__1418->SetEntries(33);
   Trigger_down__1418->SetFillColor(38);
   Trigger_down__1418->SetLineColor(38);
   Trigger_down__1418->SetMarkerColor(38);
   Trigger_down__1418->SetMarkerStyle(21);
   Trigger_down__1418->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_down__1418->GetXaxis()->SetRange(1,400);
   Trigger_down__1418->GetXaxis()->SetLabelFont(42);
   Trigger_down__1418->GetXaxis()->SetLabelSize(0.035);
   Trigger_down__1418->GetXaxis()->SetTitleSize(0.035);
   Trigger_down__1418->GetXaxis()->SetTitleFont(42);
   Trigger_down__1418->GetYaxis()->SetTitle("Events / bin");
   Trigger_down__1418->GetYaxis()->SetLabelFont(42);
   Trigger_down__1418->GetYaxis()->SetLabelSize(0.035);
   Trigger_down__1418->GetYaxis()->SetTitleSize(0.035);
   Trigger_down__1418->GetYaxis()->SetTitleOffset(0);
   Trigger_down__1418->GetYaxis()->SetTitleFont(42);
   Trigger_down__1418->GetZaxis()->SetLabelFont(42);
   Trigger_down__1418->GetZaxis()->SetLabelSize(0.035);
   Trigger_down__1418->GetZaxis()->SetTitleSize(0.035);
   Trigger_down__1418->GetZaxis()->SetTitleFont(42);
   Trigger_down__1418->Draw("E0 same");
   Double_t xAxis1104[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *Trigger_up__1419 = new TH1F("Trigger_up__1419","",32, xAxis1104);
   Trigger_up__1419->SetBinContent(0,1.000075);
   Trigger_up__1419->SetBinContent(1,1.000075);
   Trigger_up__1419->SetBinContent(2,1.000075);
   Trigger_up__1419->SetBinContent(3,1.000075);
   Trigger_up__1419->SetBinContent(4,1.000075);
   Trigger_up__1419->SetBinContent(5,1.000077);
   Trigger_up__1419->SetBinContent(6,1.000077);
   Trigger_up__1419->SetBinContent(7,1.000077);
   Trigger_up__1419->SetBinContent(8,1.000077);
   Trigger_up__1419->SetBinContent(9,1.000077);
   Trigger_up__1419->SetBinContent(10,1.000077);
   Trigger_up__1419->SetBinContent(11,1.000077);
   Trigger_up__1419->SetBinContent(12,1.000077);
   Trigger_up__1419->SetBinContent(13,1.000076);
   Trigger_up__1419->SetBinContent(14,1.000076);
   Trigger_up__1419->SetBinContent(15,1.000078);
   Trigger_up__1419->SetBinContent(16,1.000078);
   Trigger_up__1419->SetBinContent(17,1.000077);
   Trigger_up__1419->SetBinContent(18,1.000079);
   Trigger_up__1419->SetBinContent(19,1.000077);
   Trigger_up__1419->SetBinContent(20,1.000077);
   Trigger_up__1419->SetBinContent(21,1.000075);
   Trigger_up__1419->SetBinContent(22,1.000075);
   Trigger_up__1419->SetBinContent(23,1.000077);
   Trigger_up__1419->SetBinContent(24,1.000077);
   Trigger_up__1419->SetBinContent(25,1.000076);
   Trigger_up__1419->SetBinContent(26,1.000078);
   Trigger_up__1419->SetBinContent(27,1.000069);
   Trigger_up__1419->SetBinContent(28,1.000054);
   Trigger_up__1419->SetBinContent(29,0.9998391);
   Trigger_up__1419->SetBinContent(30,1.00021);
   Trigger_up__1419->SetBinContent(31,1.002897);
   Trigger_up__1419->SetBinContent(32,1.006576);
   Trigger_up__1419->SetBinError(0,0.006126035);
   Trigger_up__1419->SetBinError(1,0.006126035);
   Trigger_up__1419->SetBinError(2,0.006126035);
   Trigger_up__1419->SetBinError(3,0.006126035);
   Trigger_up__1419->SetBinError(4,0.006126035);
   Trigger_up__1419->SetBinError(5,0.006126164);
   Trigger_up__1419->SetBinError(6,0.006126164);
   Trigger_up__1419->SetBinError(7,0.006126164);
   Trigger_up__1419->SetBinError(8,0.006126164);
   Trigger_up__1419->SetBinError(9,0.006126164);
   Trigger_up__1419->SetBinError(10,0.006126164);
   Trigger_up__1419->SetBinError(11,0.006126164);
   Trigger_up__1419->SetBinError(12,0.006126164);
   Trigger_up__1419->SetBinError(13,0.006126273);
   Trigger_up__1419->SetBinError(14,0.006126273);
   Trigger_up__1419->SetBinError(15,0.006126397);
   Trigger_up__1419->SetBinError(16,0.006126397);
   Trigger_up__1419->SetBinError(17,0.006126506);
   Trigger_up__1419->SetBinError(18,0.006126745);
   Trigger_up__1419->SetBinError(19,0.006126791);
   Trigger_up__1419->SetBinError(20,0.006126847);
   Trigger_up__1419->SetBinError(21,0.006126893);
   Trigger_up__1419->SetBinError(22,0.00612695);
   Trigger_up__1419->SetBinError(23,0.006127135);
   Trigger_up__1419->SetBinError(24,0.006127423);
   Trigger_up__1419->SetBinError(25,0.006127648);
   Trigger_up__1419->SetBinError(26,0.006128803);
   Trigger_up__1419->SetBinError(27,0.006132817);
   Trigger_up__1419->SetBinError(28,0.006156853);
   Trigger_up__1419->SetBinError(29,0.006300375);
   Trigger_up__1419->SetBinError(30,0.007020068);
   Trigger_up__1419->SetBinError(31,0.009830947);
   Trigger_up__1419->SetBinError(32,0.01739109);
   Trigger_up__1419->SetEntries(33);
   Trigger_up__1419->SetFillColor(46);
   Trigger_up__1419->SetLineColor(46);
   Trigger_up__1419->SetMarkerColor(46);
   Trigger_up__1419->SetMarkerStyle(21);
   Trigger_up__1419->GetXaxis()->SetTitle("Mass [GeV]");
   Trigger_up__1419->GetXaxis()->SetRange(1,400);
   Trigger_up__1419->GetXaxis()->SetLabelFont(42);
   Trigger_up__1419->GetXaxis()->SetLabelSize(0.035);
   Trigger_up__1419->GetXaxis()->SetTitleSize(0.035);
   Trigger_up__1419->GetXaxis()->SetTitleFont(42);
   Trigger_up__1419->GetYaxis()->SetTitle("Events / bin");
   Trigger_up__1419->GetYaxis()->SetLabelFont(42);
   Trigger_up__1419->GetYaxis()->SetLabelSize(0.035);
   Trigger_up__1419->GetYaxis()->SetTitleSize(0.035);
   Trigger_up__1419->GetYaxis()->SetTitleOffset(0);
   Trigger_up__1419->GetYaxis()->SetTitleFont(42);
   Trigger_up__1419->GetZaxis()->SetLabelFont(42);
   Trigger_up__1419->GetZaxis()->SetLabelSize(0.035);
   Trigger_up__1419->GetZaxis()->SetTitleSize(0.035);
   Trigger_up__1419->GetZaxis()->SetTitleFont(42);
   Trigger_up__1419->Draw("E0 same");
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
   
   TH1D *frameR2__1420 = new TH1D("frameR2__1420","",1,0,2000);
   frameR2__1420->SetMinimum(0.5);
   frameR2__1420->SetMaximum(1.5);
   frameR2__1420->SetStats(0);
   frameR2__1420->SetLineStyle(0);
   frameR2__1420->SetMarkerStyle(20);
   frameR2__1420->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__1420->GetXaxis()->SetRange(1,1);
   frameR2__1420->GetXaxis()->SetLabelFont(43);
   frameR2__1420->GetXaxis()->SetLabelOffset(0.007);
   frameR2__1420->GetXaxis()->SetLabelSize(16);
   frameR2__1420->GetXaxis()->SetTitleSize(24);
   frameR2__1420->GetXaxis()->SetTitleOffset(5);
   frameR2__1420->GetXaxis()->SetTitleFont(43);
   frameR2__1420->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__1420->GetYaxis()->SetNdivisions(205);
   frameR2__1420->GetYaxis()->SetLabelFont(43);
   frameR2__1420->GetYaxis()->SetLabelOffset(0.007);
   frameR2__1420->GetYaxis()->SetLabelSize(20);
   frameR2__1420->GetYaxis()->SetTitleSize(20);
   frameR2__1420->GetYaxis()->SetTitleOffset(2);
   frameR2__1420->GetYaxis()->SetTitleFont(43);
   frameR2__1420->GetZaxis()->SetLabelFont(42);
   frameR2__1420->GetZaxis()->SetLabelOffset(0.007);
   frameR2__1420->GetZaxis()->SetLabelSize(0.05);
   frameR2__1420->GetZaxis()->SetTitleSize(0.06);
   frameR2__1420->GetZaxis()->SetTitleFont(42);
   frameR2__1420->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis1105[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1421 = new TH1F("nominal__1421","",32, xAxis1105);
   nominal__1421->SetBinContent(4,1.342031);
   nominal__1421->SetBinContent(12,1.151959);
   nominal__1421->SetBinContent(14,1.084463);
   nominal__1421->SetBinContent(16,1.196574);
   nominal__1421->SetBinContent(17,1.168613);
   nominal__1421->SetBinContent(18,1.052632);
   nominal__1421->SetBinContent(19,1.075269);
   nominal__1421->SetBinContent(20,1.052632);
   nominal__1421->SetBinContent(21,1.111111);
   nominal__1421->SetBinContent(22,1.119489);
   nominal__1421->SetBinContent(23,1.035507);
   nominal__1421->SetBinContent(24,1.047509);
   nominal__1421->SetBinContent(25,1.103039);
   nominal__1421->SetBinContent(26,1.052977);
   nominal__1421->SetBinContent(27,1.06418);
   nominal__1421->SetBinContent(28,1.094318);
   nominal__1421->SetBinContent(29,1.136795);
   nominal__1421->SetBinContent(30,1.168132);
   nominal__1421->SetBinContent(31,1.17137);
   nominal__1421->SetBinContent(32,1.146613);
   nominal__1421->SetBinError(4,1.343386);
   nominal__1421->SetBinError(12,1.154716);
   nominal__1421->SetBinError(14,1.087766);
   nominal__1421->SetBinError(16,1.199583);
   nominal__1421->SetBinError(17,0.830537);
   nominal__1421->SetBinError(18,1.488646);
   nominal__1421->SetBinError(19,1.52066);
   nominal__1421->SetBinError(20,1.488646);
   nominal__1421->SetBinError(21,1.571349);
   nominal__1421->SetBinError(22,0.9175469);
   nominal__1421->SetBinError(23,0.6563676);
   nominal__1421->SetBinError(24,0.7439458);
   nominal__1421->SetBinError(25,0.3505518);
   nominal__1421->SetBinError(26,0.1775344);
   nominal__1421->SetBinError(27,0.07398416);
   nominal__1421->SetBinError(28,0.03169735);
   nominal__1421->SetBinError(29,0.01631515);
   nominal__1421->SetBinError(30,0.01177884);
   nominal__1421->SetBinError(31,0.013982);
   nominal__1421->SetBinError(32,0.01574696);
   nominal__1421->SetMinimum(5e-05);
   nominal__1421->SetMaximum(1143.78);
   nominal__1421->SetEntries(28.36279);
   nominal__1421->SetFillColor(38);
   nominal__1421->SetLineColor(38);
   nominal__1421->SetMarkerColor(38);
   nominal__1421->SetMarkerStyle(21);
   nominal__1421->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1421->GetXaxis()->SetRange(1,32);
   nominal__1421->GetXaxis()->SetLabelFont(42);
   nominal__1421->GetXaxis()->SetLabelSize(0.035);
   nominal__1421->GetXaxis()->SetTitleSize(0.035);
   nominal__1421->GetXaxis()->SetTitleFont(42);
   nominal__1421->GetYaxis()->SetTitle("Tracks");
   nominal__1421->GetYaxis()->SetLabelFont(42);
   nominal__1421->GetYaxis()->SetLabelSize(0.05);
   nominal__1421->GetYaxis()->SetTitleSize(0.07);
   nominal__1421->GetYaxis()->SetTitleOffset(0);
   nominal__1421->GetYaxis()->SetTitleFont(42);
   nominal__1421->GetZaxis()->SetLabelFont(42);
   nominal__1421->GetZaxis()->SetLabelSize(0.035);
   nominal__1421->GetZaxis()->SetTitleSize(0.035);
   nominal__1421->GetZaxis()->SetTitleFont(42);
   nominal__1421->Draw("E0 same");
   Double_t xAxis1106[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__1422 = new TH1F("nominal__1422","",32, xAxis1106);
   nominal__1422->SetBinContent(4,0.9389809);
   nominal__1422->SetBinContent(12,1.02554);
   nominal__1422->SetBinContent(14,0.955136);
   nominal__1422->SetBinContent(16,1.026255);
   nominal__1422->SetBinContent(17,0.9791727);
   nominal__1422->SetBinContent(18,1.111111);
   nominal__1422->SetBinContent(19,1.010101);
   nominal__1422->SetBinContent(20,1.111111);
   nominal__1422->SetBinContent(21,1);
   nominal__1422->SetBinContent(22,0.9635679);
   nominal__1422->SetBinContent(23,0.9981716);
   nominal__1422->SetBinContent(24,1.010117);
   nominal__1422->SetBinContent(25,0.9967516);
   nominal__1422->SetBinContent(26,1.006643);
   nominal__1422->SetBinContent(27,1.001997);
   nominal__1422->SetBinContent(28,1.004595);
   nominal__1422->SetBinContent(29,0.9982975);
   nominal__1422->SetBinContent(30,0.9973982);
   nominal__1422->SetBinContent(31,1.00116);
   nominal__1422->SetBinContent(32,1.006576);
   nominal__1422->SetBinError(4,0.9389837);
   nominal__1422->SetBinError(12,1.025785);
   nominal__1422->SetBinError(14,0.9578087);
   nominal__1422->SetBinError(16,1.026786);
   nominal__1422->SetBinError(17,0.6933753);
   nominal__1422->SetBinError(18,1.571348);
   nominal__1422->SetBinError(19,1.428498);
   nominal__1422->SetBinError(20,1.571349);
   nominal__1422->SetBinError(21,1.414214);
   nominal__1422->SetBinError(22,0.7885482);
   nominal__1422->SetBinError(23,0.6316205);
   nominal__1422->SetBinError(24,0.7145599);
   nominal__1422->SetBinError(25,0.3156298);
   nominal__1422->SetBinError(26,0.1694104);
   nominal__1422->SetBinError(27,0.06947298);
   nominal__1422->SetBinError(28,0.0290101);
   nominal__1422->SetBinError(29,0.01428477);
   nominal__1422->SetBinError(30,0.0100276);
   nominal__1422->SetBinError(31,0.01191773);
   nominal__1422->SetBinError(32,0.01384281);
   nominal__1422->SetMinimum(5e-05);
   nominal__1422->SetMaximum(1143.78);
   nominal__1422->SetEntries(26.9923);
   nominal__1422->SetFillColor(46);
   nominal__1422->SetLineColor(46);
   nominal__1422->SetMarkerColor(46);
   nominal__1422->SetMarkerStyle(21);
   nominal__1422->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__1422->GetXaxis()->SetRange(1,32);
   nominal__1422->GetXaxis()->SetLabelFont(42);
   nominal__1422->GetXaxis()->SetLabelSize(0.035);
   nominal__1422->GetXaxis()->SetTitleSize(0.035);
   nominal__1422->GetXaxis()->SetTitleFont(42);
   nominal__1422->GetYaxis()->SetTitle("Tracks");
   nominal__1422->GetYaxis()->SetLabelFont(42);
   nominal__1422->GetYaxis()->SetLabelSize(0.05);
   nominal__1422->GetYaxis()->SetTitleSize(0.07);
   nominal__1422->GetYaxis()->SetTitleOffset(0);
   nominal__1422->GetYaxis()->SetTitleFont(42);
   nominal__1422->GetZaxis()->SetLabelFont(42);
   nominal__1422->GetZaxis()->SetLabelSize(0.035);
   nominal__1422->GetZaxis()->SetTitleSize(0.035);
   nominal__1422->GetZaxis()->SetTitleFont(42);
   nominal__1422->Draw("E0 same");
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
