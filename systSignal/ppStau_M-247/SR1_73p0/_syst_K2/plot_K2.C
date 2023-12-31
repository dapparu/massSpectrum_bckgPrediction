void plot_K2()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr 10 23:18:03 2023) by ROOT version 6.14/09
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
   t1->Range(-428.5714,-4.373865,2428.571,10.19308);
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
   Double_t xAxis554[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__712 = new TH1F("nominal__712","",32, xAxis554);
   nominal__712->SetBinContent(2,0.01441054);
   nominal__712->SetBinContent(3,0.01441054);
   nominal__712->SetBinContent(4,0.01391778);
   nominal__712->SetBinContent(5,0.05432382);
   nominal__712->SetBinContent(6,0.08749244);
   nominal__712->SetBinContent(7,0.2495225);
   nominal__712->SetBinContent(8,1.359236);
   nominal__712->SetBinContent(9,7.188493);
   nominal__712->SetBinContent(10,26.18366);
   nominal__712->SetBinContent(11,67.32104);
   nominal__712->SetBinContent(12,130.176);
   nominal__712->SetBinContent(13,131.8988);
   nominal__712->SetBinContent(14,61.8816);
   nominal__712->SetBinContent(15,19.00125);
   nominal__712->SetBinContent(16,6.612836);
   nominal__712->SetBinContent(17,2.733715);
   nominal__712->SetBinContent(18,1.220673);
   nominal__712->SetBinContent(19,0.8732852);
   nominal__712->SetBinContent(20,0.5433682);
   nominal__712->SetBinContent(21,0.3258789);
   nominal__712->SetBinContent(22,0.1552988);
   nominal__712->SetBinContent(23,0.1287887);
   nominal__712->SetBinContent(24,0.1276015);
   nominal__712->SetBinContent(25,0.05518857);
   nominal__712->SetBinContent(26,0.05718816);
   nominal__712->SetBinContent(28,0.0142559);
   nominal__712->SetBinError(2,0.01441054);
   nominal__712->SetBinError(3,0.01441054);
   nominal__712->SetBinError(4,0.01391778);
   nominal__712->SetBinError(5,0.02720134);
   nominal__712->SetBinError(6,0.03574617);
   nominal__712->SetBinError(7,0.0589242);
   nominal__712->SetBinError(8,0.138277);
   nominal__712->SetBinError(9,0.3183696);
   nominal__712->SetBinError(10,0.6083438);
   nominal__712->SetBinError(11,0.9773891);
   nominal__712->SetBinError(12,1.358925);
   nominal__712->SetBinError(13,1.369211);
   nominal__712->SetBinError(14,0.9385086);
   nominal__712->SetBinError(15,0.5204361);
   nominal__712->SetBinError(16,0.3075156);
   nominal__712->SetBinError(17,0.1981458);
   nominal__712->SetBinError(18,0.1311492);
   nominal__712->SetBinError(19,0.1111456);
   nominal__712->SetBinError(20,0.0882359);
   nominal__712->SetBinError(21,0.06819582);
   nominal__712->SetBinError(22,0.04700373);
   nominal__712->SetBinError(23,0.04295133);
   nominal__712->SetBinError(24,0.04257366);
   nominal__712->SetBinError(25,0.02766176);
   nominal__712->SetBinError(26,0.0286151);
   nominal__712->SetBinError(28,0.0142559);
   nominal__712->SetMinimum(5e-05);
   nominal__712->SetMaximum(1.318988e+10);
   nominal__712->SetEntries(32401);
   nominal__712->SetFillColor(1);
   nominal__712->SetMarkerStyle(20);
   nominal__712->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__712->GetXaxis()->SetRange(1,32);
   nominal__712->GetXaxis()->SetLabelFont(42);
   nominal__712->GetXaxis()->SetLabelSize(0.035);
   nominal__712->GetXaxis()->SetTitleSize(0.035);
   nominal__712->GetXaxis()->SetTitleFont(42);
   nominal__712->GetYaxis()->SetTitle("Tracks");
   nominal__712->GetYaxis()->SetLabelFont(42);
   nominal__712->GetYaxis()->SetLabelSize(0.05);
   nominal__712->GetYaxis()->SetTitleSize(0.07);
   nominal__712->GetYaxis()->SetTitleOffset(0);
   nominal__712->GetYaxis()->SetTitleFont(42);
   nominal__712->GetZaxis()->SetLabelFont(42);
   nominal__712->GetZaxis()->SetLabelSize(0.035);
   nominal__712->GetZaxis()->SetTitleSize(0.035);
   nominal__712->GetZaxis()->SetTitleFont(42);
   nominal__712->Draw("");
   Double_t xAxis555[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down2__713 = new TH1F("K_down2__713","",32, xAxis555);
   K_down2__713->SetBinContent(2,0.01441054);
   K_down2__713->SetBinContent(3,0.01441054);
   K_down2__713->SetBinContent(4,0.01391778);
   K_down2__713->SetBinContent(5,0.05432382);
   K_down2__713->SetBinContent(6,0.07341424);
   K_down2__713->SetBinContent(7,0.1924725);
   K_down2__713->SetBinContent(8,0.9834347);
   K_down2__713->SetBinContent(9,5.31412);
   K_down2__713->SetBinContent(10,20.33568);
   K_down2__713->SetBinContent(11,53.3503);
   K_down2__713->SetBinContent(12,113.1008);
   K_down2__713->SetBinContent(13,138.9485);
   K_down2__713->SetBinContent(14,81.00983);
   K_down2__713->SetBinContent(15,27.84955);
   K_down2__713->SetBinContent(16,8.99421);
   K_down2__713->SetBinContent(17,3.569766);
   K_down2__713->SetBinContent(18,1.815161);
   K_down2__713->SetBinContent(19,0.8401939);
   K_down2__713->SetBinContent(20,0.7633467);
   K_down2__713->SetBinContent(21,0.3992448);
   K_down2__713->SetBinContent(22,0.2442695);
   K_down2__713->SetBinContent(23,0.1145205);
   K_down2__713->SetBinContent(24,0.1555641);
   K_down2__713->SetBinContent(25,0.06939666);
   K_down2__713->SetBinContent(26,0.04298007);
   K_down2__713->SetBinContent(27,0.01420809);
   K_down2__713->SetBinContent(28,0.0142559);
   K_down2__713->SetBinError(2,0.01441054);
   K_down2__713->SetBinError(3,0.01441054);
   K_down2__713->SetBinError(4,0.01391778);
   K_down2__713->SetBinError(5,0.02720134);
   K_down2__713->SetBinError(6,0.03285716);
   K_down2__713->SetBinError(7,0.05155934);
   K_down2__713->SetBinError(8,0.1177443);
   K_down2__713->SetBinError(9,0.2738805);
   K_down2__713->SetBinError(10,0.5355635);
   K_down2__713->SetBinError(11,0.8696686);
   K_down2__713->SetBinError(12,1.266625);
   K_down2__713->SetBinError(13,1.405062);
   K_down2__713->SetBinError(14,1.073934);
   K_down2__713->SetBinError(15,0.6296799);
   K_down2__713->SetBinError(16,0.3581567);
   K_down2__713->SetBinError(17,0.226214);
   K_down2__713->SetBinError(18,0.1607064);
   K_down2__713->SetBinError(19,0.1087035);
   K_down2__713->SetBinError(20,0.1040304);
   K_down2__713->SetBinError(21,0.07562298);
   K_down2__713->SetBinError(22,0.05949146);
   K_down2__713->SetBinError(23,0.04051758);
   K_down2__713->SetBinError(24,0.04694274);
   K_down2__713->SetBinError(25,0.03109731);
   K_down2__713->SetBinError(26,0.02483856);
   K_down2__713->SetBinError(27,0.01420809);
   K_down2__713->SetBinError(28,0.0142559);
   K_down2__713->SetEntries(32401);
   K_down2__713->SetFillColor(38);
   K_down2__713->SetLineColor(38);
   K_down2__713->SetMarkerColor(38);
   K_down2__713->SetMarkerStyle(21);
   K_down2__713->GetXaxis()->SetTitle("Mass [GeV]");
   K_down2__713->GetXaxis()->SetRange(1,400);
   K_down2__713->GetXaxis()->SetLabelFont(42);
   K_down2__713->GetXaxis()->SetLabelSize(0.035);
   K_down2__713->GetXaxis()->SetTitleSize(0.035);
   K_down2__713->GetXaxis()->SetTitleFont(42);
   K_down2__713->GetYaxis()->SetTitle("Events / bin");
   K_down2__713->GetYaxis()->SetLabelFont(42);
   K_down2__713->GetYaxis()->SetLabelSize(0.035);
   K_down2__713->GetYaxis()->SetTitleSize(0.035);
   K_down2__713->GetYaxis()->SetTitleOffset(0);
   K_down2__713->GetYaxis()->SetTitleFont(42);
   K_down2__713->GetZaxis()->SetLabelFont(42);
   K_down2__713->GetZaxis()->SetLabelSize(0.035);
   K_down2__713->GetZaxis()->SetTitleSize(0.035);
   K_down2__713->GetZaxis()->SetTitleFont(42);
   K_down2__713->Draw("same");
   Double_t xAxis556[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up2__714 = new TH1F("K_up2__714","",32, xAxis556);
   K_up2__714->SetBinContent(2,0.01441054);
   K_up2__714->SetBinContent(3,0.01441054);
   K_up2__714->SetBinContent(4,0.01391778);
   K_up2__714->SetBinContent(5,0.06880312);
   K_up2__714->SetBinContent(6,0.08570232);
   K_up2__714->SetBinContent(7,0.3087298);
   K_up2__714->SetBinContent(8,1.834902);
   K_up2__714->SetBinContent(9,9.63912);
   K_up2__714->SetBinContent(10,32.80276);
   K_up2__714->SetBinContent(11,82.92493);
   K_up2__714->SetBinContent(12,144.3911);
   K_up2__714->SetBinContent(13,116.4386);
   K_up2__714->SetBinContent(14,45.83369);
   K_up2__714->SetBinContent(15,13.77269);
   K_up2__714->SetBinContent(16,5.086434);
   K_up2__714->SetBinContent(17,2.110933);
   K_up2__714->SetBinContent(18,1.021967);
   K_up2__714->SetBinContent(19,0.7317405);
   K_up2__714->SetBinContent(20,0.4026227);
   K_up2__714->SetBinContent(21,0.2999752);
   K_up2__714->SetBinContent(22,0.1542709);
   K_up2__714->SetBinContent(23,0.1575206);
   K_up2__714->SetBinContent(24,0.06859891);
   K_up2__714->SetBinContent(25,0.04304202);
   K_up2__714->SetBinContent(26,0.05718816);
   K_up2__714->SetBinContent(28,0.0142559);
   K_up2__714->SetBinError(2,0.01441054);
   K_up2__714->SetBinError(3,0.01441054);
   K_up2__714->SetBinError(4,0.01391778);
   K_up2__714->SetBinError(5,0.03081498);
   K_up2__714->SetBinError(6,0.03505929);
   K_up2__714->SetBinError(7,0.06597384);
   K_up2__714->SetBinError(8,0.1606242);
   K_up2__714->SetBinError(9,0.368854);
   K_up2__714->SetBinError(10,0.6810687);
   K_up2__714->SetBinError(11,1.084679);
   K_up2__714->SetBinError(12,1.431518);
   K_up2__714->SetBinError(13,1.286844);
   K_up2__714->SetBinError(14,0.8079512);
   K_up2__714->SetBinError(15,0.4430945);
   K_up2__714->SetBinError(16,0.2696954);
   K_up2__714->SetBinError(17,0.1737815);
   K_up2__714->SetBinError(18,0.1199019);
   K_up2__714->SetBinError(19,0.1016305);
   K_up2__714->SetBinError(20,0.07619337);
   K_up2__714->SetBinError(21,0.06575443);
   K_up2__714->SetBinError(22,0.04656907);
   K_up2__714->SetBinError(23,0.04753844);
   K_up2__714->SetBinError(24,0.03073952);
   K_up2__714->SetBinError(25,0.02485224);
   K_up2__714->SetBinError(26,0.0286151);
   K_up2__714->SetBinError(28,0.0142559);
   K_up2__714->SetEntries(32401);
   K_up2__714->SetFillColor(46);
   K_up2__714->SetLineColor(46);
   K_up2__714->SetMarkerColor(46);
   K_up2__714->SetMarkerStyle(21);
   K_up2__714->GetXaxis()->SetTitle("Mass [GeV]");
   K_up2__714->GetXaxis()->SetRange(1,400);
   K_up2__714->GetXaxis()->SetLabelFont(42);
   K_up2__714->GetXaxis()->SetLabelSize(0.035);
   K_up2__714->GetXaxis()->SetTitleSize(0.035);
   K_up2__714->GetXaxis()->SetTitleFont(42);
   K_up2__714->GetYaxis()->SetTitle("Events / bin");
   K_up2__714->GetYaxis()->SetLabelFont(42);
   K_up2__714->GetYaxis()->SetLabelSize(0.035);
   K_up2__714->GetYaxis()->SetTitleSize(0.035);
   K_up2__714->GetYaxis()->SetTitleOffset(0);
   K_up2__714->GetYaxis()->SetTitleFont(42);
   K_up2__714->GetZaxis()->SetLabelFont(42);
   K_up2__714->GetZaxis()->SetLabelSize(0.035);
   K_up2__714->GetZaxis()->SetTitleSize(0.035);
   K_up2__714->GetZaxis()->SetTitleFont(42);
   K_up2__714->Draw("same");
   TLine *line = new TLine(1730,0,1730,1.318988e+10);
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
   entry=leg->AddEntry("K_down2","Down","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("K_up2","Up","PE1");
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
   
   TH1D *frameR2__715 = new TH1D("frameR2__715","",1,0,2000);
   frameR2__715->SetMinimum(0.5);
   frameR2__715->SetMaximum(1.5);
   frameR2__715->SetStats(0);
   frameR2__715->SetLineStyle(0);
   frameR2__715->SetMarkerStyle(20);
   frameR2__715->GetXaxis()->SetRange(1,1);
   frameR2__715->GetXaxis()->SetLabelFont(43);
   frameR2__715->GetXaxis()->SetLabelOffset(0.007);
   frameR2__715->GetXaxis()->SetLabelSize(16);
   frameR2__715->GetXaxis()->SetTitleSize(24);
   frameR2__715->GetXaxis()->SetTitleOffset(3.75);
   frameR2__715->GetXaxis()->SetTitleFont(43);
   frameR2__715->GetYaxis()->SetTitle("Ratio #int_{m}^{#infty}");
   frameR2__715->GetYaxis()->SetNdivisions(205);
   frameR2__715->GetYaxis()->SetLabelFont(43);
   frameR2__715->GetYaxis()->SetLabelOffset(0.007);
   frameR2__715->GetYaxis()->SetLabelSize(20);
   frameR2__715->GetYaxis()->SetTitleSize(20);
   frameR2__715->GetYaxis()->SetTitleOffset(2);
   frameR2__715->GetYaxis()->SetTitleFont(43);
   frameR2__715->GetZaxis()->SetLabelFont(42);
   frameR2__715->GetZaxis()->SetLabelOffset(0.007);
   frameR2__715->GetZaxis()->SetLabelSize(0.05);
   frameR2__715->GetZaxis()->SetTitleSize(0.06);
   frameR2__715->GetZaxis()->SetTitleFont(42);
   frameR2__715->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis557[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_down2__716 = new TH1F("K_down2__716","",32, xAxis557);
   K_down2__716->SetBinContent(0,1);
   K_down2__716->SetBinContent(1,1);
   K_down2__716->SetBinContent(2,1);
   K_down2__716->SetBinContent(3,1);
   K_down2__716->SetBinContent(4,1);
   K_down2__716->SetBinContent(5,1);
   K_down2__716->SetBinContent(6,1);
   K_down2__716->SetBinContent(7,0.9999693);
   K_down2__716->SetBinContent(8,0.9998447);
   K_down2__716->SetBinContent(9,0.9990219);
   K_down2__716->SetBinContent(10,0.9948602);
   K_down2__716->SetBinContent(11,0.9810588);
   K_down2__716->SetBinContent(12,0.9414201);
   K_down2__716->SetBinContent(13,0.8519312);
   K_down2__716->SetBinContent(14,0.7445078);
   K_down2__716->SetBinContent(15,0.7095498);
   K_down2__716->SetBinContent(16,0.7541227);
   K_down2__716->SetBinContent(17,0.7752472);
   K_down2__716->SetBinContent(18,0.7827891);
   K_down2__716->SetBinContent(19,0.8581154);
   K_down2__716->SetBinContent(20,0.7743312);
   K_down2__716->SetBinContent(21,0.8195826);
   K_down2__716->SetBinContent(22,0.8216206);
   K_down2__716->SetBinContent(23,0.9320984);
   K_down2__716->SetBinContent(24,0.857726);
   K_down2__716->SetBinContent(25,0.8991194);
   K_down2__716->SetBinContent(26,1);
   K_down2__716->SetBinContent(27,0.5008398);
   K_down2__716->SetBinContent(28,1);
   K_down2__716->SetBinError(0,0.007872498);
   K_down2__716->SetBinError(1,0.007872498);
   K_down2__716->SetBinError(2,0.007872498);
   K_down2__716->SetBinError(3,0.007872619);
   K_down2__716->SetBinError(4,0.007872741);
   K_down2__716->SetBinError(5,0.007872863);
   K_down2__716->SetBinError(6,0.007873349);
   K_down2__716->SetBinError(7,0.007873777);
   K_down2__716->SetBinError(8,0.00787474);
   K_down2__716->SetBinError(9,0.007878431);
   K_down2__716->SetBinError(10,0.007900255);
   K_down2__716->SetBinError(11,0.008001528);
   K_down2__716->SetBinError(12,0.008290118);
   K_down2__716->SetBinError(13,0.009205448);
   K_down2__716->SetBinError(14,0.01212159);
   K_down2__716->SetBinError(15,0.01962884);
   K_down2__716->SetBinError(16,0.03328961);
   K_down2__716->SetBinError(17,0.04939945);
   K_down2__716->SetBinError(18,0.06653923);
   K_down2__716->SetBinError(19,0.09233048);
   K_down2__716->SetBinError(20,0.1038377);
   K_down2__716->SetBinError(21,0.1420735);
   K_down2__716->SetBinError(22,0.1804717);
   K_down2__716->SetBinError(23,0.2495234);
   K_down2__716->SetBinError(24,0.2758308);
   K_down2__716->SetBinError(25,0.4137208);
   K_down2__716->SetBinError(26,0.6328283);
   K_down2__716->SetBinError(27,0.6134013);
   K_down2__716->SetBinError(28,1.414214);
   K_down2__716->SetEntries(33);
   K_down2__716->SetFillColor(38);
   K_down2__716->SetLineColor(38);
   K_down2__716->SetMarkerColor(38);
   K_down2__716->SetMarkerStyle(21);
   K_down2__716->GetXaxis()->SetTitle("Mass [GeV]");
   K_down2__716->GetXaxis()->SetRange(1,400);
   K_down2__716->GetXaxis()->SetLabelFont(42);
   K_down2__716->GetXaxis()->SetLabelSize(0.035);
   K_down2__716->GetXaxis()->SetTitleSize(0.035);
   K_down2__716->GetXaxis()->SetTitleFont(42);
   K_down2__716->GetYaxis()->SetTitle("Events / bin");
   K_down2__716->GetYaxis()->SetLabelFont(42);
   K_down2__716->GetYaxis()->SetLabelSize(0.035);
   K_down2__716->GetYaxis()->SetTitleSize(0.035);
   K_down2__716->GetYaxis()->SetTitleOffset(0);
   K_down2__716->GetYaxis()->SetTitleFont(42);
   K_down2__716->GetZaxis()->SetLabelFont(42);
   K_down2__716->GetZaxis()->SetLabelSize(0.035);
   K_down2__716->GetZaxis()->SetTitleSize(0.035);
   K_down2__716->GetZaxis()->SetTitleFont(42);
   K_down2__716->Draw("E0 same");
   Double_t xAxis558[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *K_up2__717 = new TH1F("K_up2__717","",32, xAxis558);
   K_up2__717->SetBinContent(0,1);
   K_up2__717->SetBinContent(1,1);
   K_up2__717->SetBinContent(2,1);
   K_up2__717->SetBinContent(3,1);
   K_up2__717->SetBinContent(4,1);
   K_up2__717->SetBinContent(5,1);
   K_up2__717->SetBinContent(6,1.000032);
   K_up2__717->SetBinContent(7,1.000028);
   K_up2__717->SetBinContent(8,1.000157);
   K_up2__717->SetBinContent(9,1.001201);
   K_up2__717->SetBinContent(10,1.006718);
   K_up2__717->SetBinContent(11,1.023258);
   K_up2__717->SetBinContent(12,1.076293);
   K_up2__717->SetBinContent(13,1.211802);
   K_up2__717->SetBinContent(14,1.343718);
   K_up2__717->SetBinContent(15,1.331425);
   K_up2__717->SetBinContent(16,1.266001);
   K_up2__717->SetBinContent(17,1.231747);
   K_up2__717->SetBinContent(18,1.186483);
   K_up2__717->SetBinContent(19,1.182271);
   K_up2__717->SetBinContent(20,1.175448);
   K_up2__717->SetBinContent(21,1.087247);
   K_up2__717->SetBinContent(22,1.08779);
   K_up2__717->SetBinContent(23,1.124535);
   K_up2__717->SetBinContent(24,1.388613);
   K_up2__717->SetBinContent(25,1.106096);
   K_up2__717->SetBinContent(26,1);
   K_up2__717->SetBinContent(27,1);
   K_up2__717->SetBinContent(28,1);
   K_up2__717->SetBinError(0,0.007872498);
   K_up2__717->SetBinError(1,0.007872498);
   K_up2__717->SetBinError(2,0.007872498);
   K_up2__717->SetBinError(3,0.00787262);
   K_up2__717->SetBinError(4,0.007872742);
   K_up2__717->SetBinError(5,0.007872864);
   K_up2__717->SetBinError(6,0.007873659);
   K_up2__717->SetBinError(7,0.007874358);
   K_up2__717->SetBinError(8,0.007877809);
   K_up2__717->SetBinError(9,0.007899958);
   K_up2__717->SetBinError(10,0.008018305);
   K_up2__717->SetBinError(11,0.008434356);
   K_up2__717->SetBinError(12,0.009801949);
   K_up2__717->SetBinError(13,0.01431193);
   K_up2__717->SetBinError(14,0.02536154);
   K_up2__717->SetBinError(15,0.04301976);
   K_up2__717->SetBinError(16,0.06351758);
   K_up2__717->SetBinError(17,0.08795027);
   K_up2__717->SetBinError(18,0.1116419);
   K_up2__717->SetBinError(19,0.1379632);
   K_up2__717->SetBinError(20,0.1747222);
   K_up2__717->SetBinError(21,0.2017115);
   K_up2__717->SetBinError(22,0.2552219);
   K_up2__717->SetBinError(23,0.3158132);
   K_up2__717->SetBinError(24,0.5060337);
   K_up2__717->SetBinError(25,0.5379676);
   K_up2__717->SetBinError(26,0.6328284);
   K_up2__717->SetBinError(27,1.414214);
   K_up2__717->SetBinError(28,1.414214);
   K_up2__717->SetEntries(33);
   K_up2__717->SetFillColor(46);
   K_up2__717->SetLineColor(46);
   K_up2__717->SetMarkerColor(46);
   K_up2__717->SetMarkerStyle(21);
   K_up2__717->GetXaxis()->SetTitle("Mass [GeV]");
   K_up2__717->GetXaxis()->SetRange(1,400);
   K_up2__717->GetXaxis()->SetLabelFont(42);
   K_up2__717->GetXaxis()->SetLabelSize(0.035);
   K_up2__717->GetXaxis()->SetTitleSize(0.035);
   K_up2__717->GetXaxis()->SetTitleFont(42);
   K_up2__717->GetYaxis()->SetTitle("Events / bin");
   K_up2__717->GetYaxis()->SetLabelFont(42);
   K_up2__717->GetYaxis()->SetLabelSize(0.035);
   K_up2__717->GetYaxis()->SetTitleSize(0.035);
   K_up2__717->GetYaxis()->SetTitleOffset(0);
   K_up2__717->GetYaxis()->SetTitleFont(42);
   K_up2__717->GetZaxis()->SetLabelFont(42);
   K_up2__717->GetZaxis()->SetLabelSize(0.035);
   K_up2__717->GetZaxis()->SetTitleSize(0.035);
   K_up2__717->GetZaxis()->SetTitleFont(42);
   K_up2__717->Draw("E0 same");
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
   
   TH1D *frameR2__718 = new TH1D("frameR2__718","",1,0,2000);
   frameR2__718->SetMinimum(0.5);
   frameR2__718->SetMaximum(1.5);
   frameR2__718->SetStats(0);
   frameR2__718->SetLineStyle(0);
   frameR2__718->SetMarkerStyle(20);
   frameR2__718->GetXaxis()->SetTitle("Mass (GeV)");
   frameR2__718->GetXaxis()->SetRange(1,1);
   frameR2__718->GetXaxis()->SetLabelFont(43);
   frameR2__718->GetXaxis()->SetLabelOffset(0.007);
   frameR2__718->GetXaxis()->SetLabelSize(16);
   frameR2__718->GetXaxis()->SetTitleSize(24);
   frameR2__718->GetXaxis()->SetTitleOffset(5);
   frameR2__718->GetXaxis()->SetTitleFont(43);
   frameR2__718->GetYaxis()->SetTitle("#frac{nominal}{var}");
   frameR2__718->GetYaxis()->SetNdivisions(205);
   frameR2__718->GetYaxis()->SetLabelFont(43);
   frameR2__718->GetYaxis()->SetLabelOffset(0.007);
   frameR2__718->GetYaxis()->SetLabelSize(20);
   frameR2__718->GetYaxis()->SetTitleSize(20);
   frameR2__718->GetYaxis()->SetTitleOffset(2);
   frameR2__718->GetYaxis()->SetTitleFont(43);
   frameR2__718->GetZaxis()->SetLabelFont(42);
   frameR2__718->GetZaxis()->SetLabelOffset(0.007);
   frameR2__718->GetZaxis()->SetLabelSize(0.05);
   frameR2__718->GetZaxis()->SetTitleSize(0.06);
   frameR2__718->GetZaxis()->SetTitleFont(42);
   frameR2__718->Draw("AXIS");
   line = new TLine(0,1,2000,1);
   line->SetLineStyle(3);
   line->Draw();
   line = new TLine(0,1.2,2000,1.2);
   line->SetLineStyle(4);
   line->Draw();
   line = new TLine(0,0.8,2000,0.8);
   line->SetLineStyle(4);
   line->Draw();
   Double_t xAxis559[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__719 = new TH1F("nominal__719","",32, xAxis559);
   nominal__719->SetBinContent(2,1);
   nominal__719->SetBinContent(3,1);
   nominal__719->SetBinContent(4,1);
   nominal__719->SetBinContent(5,1);
   nominal__719->SetBinContent(6,1.191764);
   nominal__719->SetBinContent(7,1.296406);
   nominal__719->SetBinContent(8,1.382132);
   nominal__719->SetBinContent(9,1.352715);
   nominal__719->SetBinContent(10,1.287572);
   nominal__719->SetBinContent(11,1.261868);
   nominal__719->SetBinContent(12,1.150974);
   nominal__719->SetBinContent(13,0.9492639);
   nominal__719->SetBinContent(14,0.7638777);
   nominal__719->SetBinContent(15,0.6822821);
   nominal__719->SetBinContent(16,0.7352325);
   nominal__719->SetBinContent(17,0.7657967);
   nominal__719->SetBinContent(18,0.6724871);
   nominal__719->SetBinContent(19,1.039385);
   nominal__719->SetBinContent(20,0.7118236);
   nominal__719->SetBinContent(21,0.8162382);
   nominal__719->SetBinContent(22,0.6357682);
   nominal__719->SetBinContent(23,1.124591);
   nominal__719->SetBinContent(24,0.8202502);
   nominal__719->SetBinContent(25,0.7952626);
   nominal__719->SetBinContent(26,1.330574);
   nominal__719->SetBinContent(28,1);
   nominal__719->SetBinError(2,1.414214);
   nominal__719->SetBinError(3,1.414214);
   nominal__719->SetBinError(4,1.414214);
   nominal__719->SetBinError(5,0.7081332);
   nominal__719->SetBinError(6,0.7222052);
   nominal__719->SetBinError(7,0.4629549);
   nominal__719->SetBinError(8,0.2171486);
   nominal__719->SetBinError(9,0.09192185);
   nominal__719->SetBinError(10,0.04521925);
   nominal__719->SetBinError(11,0.02754539);
   nominal__719->SetBinError(12,0.01762137);
   nominal__719->SetBinError(13,0.01375663);
   nominal__719->SetBinError(14,0.01538711);
   nominal__719->SetBinError(15,0.02423209);
   nominal__719->SetBinError(16,0.04501286);
   nominal__719->SetBinError(17,0.073729);
   nominal__719->SetBinError(18,0.093623);
   nominal__719->SetBinError(19,0.1886344);
   nominal__719->SetBinError(20,0.1509037);
   nominal__719->SetBinError(21,0.2303917);
   nominal__719->SetBinError(22,0.2469882);
   nominal__719->SetBinError(23,0.5467868);
   nominal__719->SetBinError(24,0.3690008);
   nominal__719->SetBinError(25,0.5346781);
   nominal__719->SetBinError(26,1.017125);
   nominal__719->SetBinError(28,1.414214);
   nominal__719->SetMinimum(5e-05);
   nominal__719->SetMaximum(1.318988e+10);
   nominal__719->SetEntries(59.05895);
   nominal__719->SetFillColor(38);
   nominal__719->SetLineColor(38);
   nominal__719->SetMarkerColor(38);
   nominal__719->SetMarkerStyle(21);
   nominal__719->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__719->GetXaxis()->SetRange(1,32);
   nominal__719->GetXaxis()->SetLabelFont(42);
   nominal__719->GetXaxis()->SetLabelSize(0.035);
   nominal__719->GetXaxis()->SetTitleSize(0.035);
   nominal__719->GetXaxis()->SetTitleFont(42);
   nominal__719->GetYaxis()->SetTitle("Tracks");
   nominal__719->GetYaxis()->SetLabelFont(42);
   nominal__719->GetYaxis()->SetLabelSize(0.05);
   nominal__719->GetYaxis()->SetTitleSize(0.07);
   nominal__719->GetYaxis()->SetTitleOffset(0);
   nominal__719->GetYaxis()->SetTitleFont(42);
   nominal__719->GetZaxis()->SetLabelFont(42);
   nominal__719->GetZaxis()->SetLabelSize(0.035);
   nominal__719->GetZaxis()->SetTitleSize(0.035);
   nominal__719->GetZaxis()->SetTitleFont(42);
   nominal__719->Draw("E0 same");
   Double_t xAxis560[33] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 405, 435, 475, 525, 585, 660, 755, 875, 1025, 1210, 1440, 1730, 2000}; 
   
   TH1F *nominal__720 = new TH1F("nominal__720","",32, xAxis560);
   nominal__720->SetBinContent(2,1);
   nominal__720->SetBinContent(3,1);
   nominal__720->SetBinContent(4,1);
   nominal__720->SetBinContent(5,0.7895547);
   nominal__720->SetBinContent(6,1.020888);
   nominal__720->SetBinContent(7,0.8082228);
   nominal__720->SetBinContent(8,0.7407677);
   nominal__720->SetBinContent(9,0.7457623);
   nominal__720->SetBinContent(10,0.798215);
   nominal__720->SetBinContent(11,0.8118311);
   nominal__720->SetBinContent(12,0.901552);
   nominal__720->SetBinContent(13,1.132776);
   nominal__720->SetBinContent(14,1.350133);
   nominal__720->SetBinContent(15,1.379632);
   nominal__720->SetBinContent(16,1.300093);
   nominal__720->SetBinContent(17,1.295027);
   nominal__720->SetBinContent(18,1.194435);
   nominal__720->SetBinContent(19,1.193436);
   nominal__720->SetBinContent(20,1.349572);
   nominal__720->SetBinContent(21,1.086353);
   nominal__720->SetBinContent(22,1.006663);
   nominal__720->SetBinContent(23,0.8175989);
   nominal__720->SetBinContent(24,1.86011);
   nominal__720->SetBinContent(25,1.282202);
   nominal__720->SetBinContent(26,1);
   nominal__720->SetBinContent(28,1);
   nominal__720->SetBinError(2,1.414214);
   nominal__720->SetBinError(3,1.414214);
   nominal__720->SetBinError(4,1.414214);
   nominal__720->SetBinError(5,0.530423);
   nominal__720->SetBinError(6,0.590239);
   nominal__720->SetBinError(7,0.2574048);
   nominal__720->SetBinError(8,0.09941818);
   nominal__720->SetBinError(9,0.04364979);
   nominal__720->SetBinError(10,0.02487165);
   nominal__720->SetBinError(11,0.0158645);
   nominal__720->SetBinError(12,0.01297942);
   nominal__720->SetBinError(13,0.01717567);
   nominal__720->SetBinError(14,0.03139622);
   nominal__720->SetBinError(15,0.0582921);
   nominal__720->SetBinError(16,0.09169016);
   nominal__720->SetBinError(17,0.1420462);
   nominal__720->SetBinError(18,0.1900182);
   nominal__720->SetBinError(19,0.224824);
   nominal__720->SetBinError(20,0.3365343);
   nominal__720->SetBinError(21,0.3292227);
   nominal__720->SetBinError(22,0.4303173);
   nominal__720->SetBinError(23,0.3677397);
   nominal__720->SetBinError(24,1.039196);
   nominal__720->SetBinError(25,0.9803684);
   nominal__720->SetBinError(26,0.7076266);
   nominal__720->SetBinError(28,1.414214);
   nominal__720->SetMinimum(5e-05);
   nominal__720->SetMaximum(1.318988e+10);
   nominal__720->SetEntries(65.17893);
   nominal__720->SetFillColor(46);
   nominal__720->SetLineColor(46);
   nominal__720->SetMarkerColor(46);
   nominal__720->SetMarkerStyle(21);
   nominal__720->GetXaxis()->SetTitle("Mass (GeV)");
   nominal__720->GetXaxis()->SetRange(1,32);
   nominal__720->GetXaxis()->SetLabelFont(42);
   nominal__720->GetXaxis()->SetLabelSize(0.035);
   nominal__720->GetXaxis()->SetTitleSize(0.035);
   nominal__720->GetXaxis()->SetTitleFont(42);
   nominal__720->GetYaxis()->SetTitle("Tracks");
   nominal__720->GetYaxis()->SetLabelFont(42);
   nominal__720->GetYaxis()->SetLabelSize(0.05);
   nominal__720->GetYaxis()->SetTitleSize(0.07);
   nominal__720->GetYaxis()->SetTitleOffset(0);
   nominal__720->GetYaxis()->SetTitleFont(42);
   nominal__720->GetZaxis()->SetLabelFont(42);
   nominal__720->GetZaxis()->SetLabelSize(0.035);
   nominal__720->GetZaxis()->SetTitleSize(0.035);
   nominal__720->GetZaxis()->SetTitleFont(42);
   nominal__720->Draw("E0 same");
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
