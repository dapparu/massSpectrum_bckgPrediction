// Usage:
// root -l -q -b step2_backgroundPrediction.C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "SUSYBSMAnalysis/Analyzer/interface/Regions.h"

using namespace std;



void step2_backgroundPrediction(){
    ifstream infile;
    infile.open("configFile_readHisto_toLaunch.txt");
    std::string line;
    std::string filename;
    std::string st_sample;
    std::string dirname;
    int nPE, cutIndex;
    int rebineta,rebinih,rebinp,rebinmass;
    bool rebin;
    bool corrTemplateIh, corrTemplateP;
    int fitIh, fitP;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> st_sample >> dirname >> nPE >> cutIndex >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass >> corrTemplateIh >> corrTemplateP >> fitIh >> fitP;
    }

    std::string commandDir = "mkdir -p "+dirname;
    system(commandDir.c_str());

    std::string outfilename_;
    if(!rebin)outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_analysed";
    else outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass);

    if(corrTemplateIh) outfilename_ += "_corrTemplateIh";
    if(corrTemplateP) outfilename_ += "_corrTemplateP";
    if(fitIh==0) outfilename_ += "_fitIhDown";
    if(fitIh==2) outfilename_ += "_fitIhUp";
    if(fitP==0) outfilename_ += "_fitPDown";
    if(fitP==2) outfilename_ += "_fitPUp";

    outfilename_ += "_nPE"+to_string(nPE);
    //outfilename_ += "_11april_wFit_ih0sigmapt0iso0";
    //outfilename_ += "_14april_ih0sigmapt2iso1";
    outfilename_ += "_UnB_v3_Data_v2";

    std::cout << outfilename_ << std::endl;

    bool bool_rebin=rebin;
    
    TFile* ifile = new TFile((filename+".root").c_str());

    // histograms used for the mass prediction
    //------------

/*
    std::string dir = "analyzer/BaseName/";
    //std::string dir = "HSCParticleAnalyzer/BaseName/";
    TH2F* eta_cutIndex_regA = (TH2F*)ifile->Get((dir+"Pred_EtaB").c_str())->Clone(); 
    TH2F* eta_cutIndex_regB = (TH2F*)ifile->Get((dir+"Pred_EtaS").c_str())->Clone(); 
    TH3F* ih_eta_cutIndex_regB = (TH3F*)ifile->Get((dir+"Pred_EtaI").c_str())->Clone(); 
    TH3F* eta_p_cutIndex_regC = (TH3F*)ifile->Get((dir+"Pred_EtaP").c_str())->Clone(); 
    TH1F* H_A = (TH1F*)ifile->Get((dir+"H_A").c_str())->Clone();
    TH1F* H_B = (TH1F*)ifile->Get((dir+"H_B").c_str())->Clone();
    TH1F* H_C = (TH1F*)ifile->Get((dir+"H_C").c_str())->Clone();
    TH2F* mass_cutIndex = (TH2F*)ifile->Get((dir+"Mass").c_str())->Clone();*/

    //------------

    Region ra_ias50;
    Region rc_ias50;

    Region rb_50ias60;
    Region rb_60ias70;
    Region rb_70ias80;
    Region rb_80ias90;
    Region rb_50ias90;
    Region rb_50ias99;
    Region rb_50ias999;
    Region rb_90ias100;
    Region rb_99ias100;
    Region rb_999ias100;

    Region rd_50ias60;
    Region rd_60ias70;
    Region rd_70ias80;
    Region rd_80ias90;
    Region rd_50ias90;
    Region rd_50ias99;
    Region rd_50ias999;
    Region rd_90ias100;
    Region rd_99ias100;
    Region rd_999ias100;
    
    Region rbc_50ias60;
    Region rbc_60ias70;
    Region rbc_70ias80;
    Region rbc_80ias90;
    Region rbc_50ias90;
    Region rbc_50ias99;
    Region rbc_50ias999;
    Region rbc_90ias100;
    Region rbc_99ias100;
    Region rbc_999ias100;

    /*Region ra_1o4;
    Region ra_2o4;
    Region ra_3o4;
    Region ra_4o4;
    Region rb_1o2;
    Region rb_2o2;
    Region rc_1o2;
    Region rc_2o2;
    Region rd_1o1;*/
    Region rb_true;
    Region rc_true;

    Region rb_true_50ias60;
    Region rc_true_50ias60;

    Region rb_true_60ias70;
    Region rc_true_60ias70;
    
    Region rb_true_70ias80;
    Region rc_true_70ias80;

    Region rb_true_80ias90;
    Region rc_true_80ias90;

    Region rb_true_50ias90;
    Region rc_true_50ias90;

    Region rb_true_90ias100;
    Region rc_true_90ias100;

    /*Region rbc_1o2_b;
    Region rbc_2o2_b;
    Region rbc_1o2_c;
    Region rbc_2o2_c;*/
    Region rbc_true_b;
    Region rbc_true_c;
    
   
    // loading histograms used to validate the background estimate method in data --> base on Ias slices 
    // ------------------------------------------------------------------------------------------------------
    
    std::cout << "loading..." << std::endl;
/*    
    loadHistograms(ra_ias50,ifile,"regionA_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionA_50"); 
    loadHistograms(rc_ias50,ifile,"regionC_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionC_50"); 
    
    loadHistograms(rb_50ias60,ifile,"regionB_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_50"); 
    loadHistograms(rb_60ias70,ifile,"regionB_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_60"); 
    loadHistograms(rb_70ias80,ifile,"regionB_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_70"); 
    loadHistograms(rb_80ias90,ifile,"regionB_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_80"); 
    loadHistograms(rb_50ias90,ifile,"regionB_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_50_90"); 
    loadHistograms(rb_50ias99,ifile,"regionB_50_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_50_99"); 
    loadHistograms(rb_50ias999,ifile,"regionB_50_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_50_999"); 
    loadHistograms(rb_90ias100,ifile,"regionB_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_90_100"); 
    loadHistograms(rb_99ias100,ifile,"regionB_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_99_100"); 
    loadHistograms(rb_999ias100,ifile,"regionB_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionB_999_100"); 

    loadHistograms(rd_50ias60,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_50"); 
    loadHistograms(rd_60ias70,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_60"); 
    loadHistograms(rd_70ias80,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_70"); 
    loadHistograms(rd_80ias90,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_80"); 
    loadHistograms(rd_50ias90,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_50_90"); 
    loadHistograms(rd_50ias99,ifile,"regionD_50_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_50_99"); 
    loadHistograms(rd_50ias999,ifile,"regionD_50_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_50_999"); 
    loadHistograms(rd_90ias100,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_90_100");
    loadHistograms(rd_99ias100,ifile,"regionD_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_99_100");
    loadHistograms(rd_999ias100,ifile,"regionD_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionD_999_100");
     
    loadHistograms(rbc_50ias60,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_50"); 
    loadHistograms(rbc_60ias70,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_60"); 
    loadHistograms(rbc_70ias80,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_70"); 
    loadHistograms(rbc_80ias90,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_80"); 
    loadHistograms(rbc_50ias90,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_50_90"); 
    loadHistograms(rbc_50ias99,ifile,"regionD_50_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_50_99"); 
    loadHistograms(rbc_50ias999,ifile,"regionD_50_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_50_999"); 
    loadHistograms(rbc_90ias100,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_90_100");
    loadHistograms(rbc_99ias100,ifile,"regionD_99",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_99_100");
    loadHistograms(rbc_999ias100,ifile,"regionD_999",bool_rebin,rebineta,rebinp,rebinih,rebinmass,"regionBC_999_100");
*/
/*    loadHistograms(ra_ias50,ifile,"regionA_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regionC_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regionB_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
*/
/*
    loadHistograms(ra_ias50,ifile,"regionA_ias50_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regionC_ias50_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias90,ifile,"regionB_50ias90_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regionD_50ias90_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_50ias90,ifile,"regionD_50ias90_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_90ias100,ifile,"regionB_90ias100_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regionD_90ias100_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_90ias100,ifile,"regionD_90ias100_testIhPt",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

*/

//comment test Raph pT>200 GeV

    loadHistograms(ra_ias50,ifile,"regionA_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regionC_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regionB_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regionB_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regionB_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regionB_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regionB_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90ias100,ifile,"regionB_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_99ias100,ifile,"regionB_99ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_999ias100,ifile,"regionB_999ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rd_99ias100,ifile,"regionD_99ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rd_999ias100,ifile,"regionD_999ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_99ias100,ifile,"regionD_99ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_999ias100,ifile,"regionD_999ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);


    loadHistograms(rb_true,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 


/*
    loadHistograms(ra_1o4,ifile,"regA_1o4",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(ra_2o4,ifile,"regA_2o4",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(ra_3o4,ifile,"regA_3o4",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(ra_4o4,ifile,"regA_4o4",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_1o2,ifile,"regB_1o2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_2o2,ifile,"regB_2o2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_1o2,ifile,"regC_1o2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_2o2,ifile,"regC_2o2",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_1o1,ifile,"regD_1o1",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rb_true,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_true,ifile,"regionD_50_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
*/



/*    loadHistograms(ra_ias50,ifile,"regA_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regC_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regB_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regB_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regB_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regB_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regB_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90ias100,ifile,"regB_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90ias100,ifile,"regD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);*/
    
    // ------------------------------------------------------------------------------------------------------
    
    std::cout << "Regions loaded" << std::endl;

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving... " << std::endl;


    // estimate the background in different Ias slices, each containing 10% of the statistic 
    // ------------------------------------------------------------------------------------------------------
 

    //bckgEstimate(st_sample, dirname, rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    //bckgEstimate(st_sample, dirname, rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    //bckgEstimate(st_sample, dirname, rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    //bckgEstimate(st_sample, dirname, rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, "80ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);

    bckgEstimate(st_sample, dirname, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    
    bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_99ias100, rc_ias50, rbc_99ias100, ra_ias50, rd_99ias100, "99ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_999ias100, rc_ias50, rbc_999ias100, ra_ias50, rd_999ias100, "999ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);

/*
    vector<float> vX;
    vector<float> vXerr;

    vX.push_back(1);
    vX.push_back(2);
    vX.push_back(3);
    vX.push_back(4);
    vX.push_back(5);
    vX.push_back(7);

    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);

    bckgEstimate(st_sample, dirname, rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    bckgEstimate(st_sample, dirname, rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    bckgEstimate(st_sample, dirname, rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    bckgEstimate(st_sample, dirname, rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, "80ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    
    bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    
    bckgEstimate(st_sample, dirname, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);

    //bckgEstimate(st_sample, dirname, rb_true, rc_true, rbc_50ias90, ra_ias50, rd_50ias90, "True_BC", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    for(int i=0;i<v_par0.size();i++){
        cout << vX[i] << " " << vXerr[i] << " " << v_par0[i] << " " << v_err0[i] << " " << v_par1[i] << " " << v_err1[i] << endl;
    }
    TGraphErrors grPar0(vX.size(),vX.data(),v_par0.data(),vXerr.data(),v_err0.data());
    TGraphErrors grPar1(vX.size(),vX.data(),v_par1.data(),vXerr.data(),v_err1.data());

    grPar0.SetTitle("Par0");
    grPar1.SetTitle("Par1");

    char const *range[7] = {"D_{50-60}","D_{60-70}","D_{70-80}","D_{80-90}","D_{90-100}","","D_{50-90}"};

    for(int i=1;i<=7;i++){
        grPar0.GetXaxis()->SetBinLabel(i,range[i-1]);
        grPar1.GetXaxis()->SetBinLabel(i,range[i-1]);
    }

    TCanvas c;
    grPar0.Draw("AP*");
    c.SaveAs((dirname+"grPar0.pdf").c_str());
    c.SaveAs((dirname+"grPar0.root").c_str());
    c.SaveAs((dirname+"grPar0.C").c_str());
    grPar1.Draw("AP*");
    c.SaveAs((dirname+"grPar1.pdf").c_str());
    c.SaveAs((dirname+"grPar1.root").c_str());
    c.SaveAs((dirname+"grPar1.C").c_str());

    grPar0.Write();
    grPar1.Write();
*/
    /*
    vector<float> vX;
    vector<float> vXerr;

    vX.push_back(1);
    vX.push_back(2);
    vX.push_back(3);
    vX.push_back(4);
    vX.push_back(5);
    vX.push_back(7);

    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    vXerr.push_back(0);
    
    
    bckgEstimate(st_sample, dirname, rb_true_60ias70, rc_true_60ias70, rbc_60ias70, ra_ias50, rd_60ias70, "True_60ias70", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_true_50ias60, rc_true_50ias60, rbc_50ias60, ra_ias50, rd_50ias60, "True_50ias60", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_true_70ias80, rc_true_70ias80, rbc_70ias80, ra_ias50, rd_70ias80, "True_70ias80", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_true_80ias90, rc_true_80ias90, rbc_80ias90, ra_ias50, rd_80ias90, "True_80ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    bckgEstimate(st_sample, dirname, rb_true_90ias100, rc_true_90ias100, rbc_90ias100, ra_ias50, rd_90ias100, "True_90ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    bckgEstimate(st_sample, dirname, rb_true_50ias90, rc_true_50ias90, rbc_50ias90, ra_ias50, rd_50ias90, "True_50ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, false);
    

    for(int i=0;i<v_par0.size();i++){
        cout << vX[i] << " " << vXerr[i] << " " << v_par0[i] << " " << v_err0[i] << " " << v_par1[i] << " " << v_err1[i] << endl;
    }
    TGraphErrors grPar0(vX.size(),vX.data(),v_par0.data(),vXerr.data(),v_err0.data());
    TGraphErrors grPar1(vX.size(),vX.data(),v_par1.data(),vXerr.data(),v_err1.data());

    grPar0.SetTitle("Par0");
    grPar1.SetTitle("Par1");

    char const *range[7] = {"D_{50-60}","D_{60-70}","D_{70-80}","D_{80-90}","D_{90-100}","","D_{50-90}"};

    for(int i=1;i<=7;i++){
        grPar0.GetXaxis()->SetBinLabel(i,range[i-1]);
        grPar1.GetXaxis()->SetBinLabel(i,range[i-1]);
    }

    TCanvas c;
    grPar0.Draw("AP*");
    c.SaveAs((dirname+"grPar0.pdf").c_str());
    c.SaveAs((dirname+"grPar0.root").c_str());
    c.SaveAs((dirname+"grPar0.C").c_str());
    grPar1.Draw("AP*");
    c.SaveAs((dirname+"grPar1.pdf").c_str());
    c.SaveAs((dirname+"grPar1.root").c_str());
    c.SaveAs((dirname+"grPar1.C").c_str());

    grPar0.Write();
    grPar1.Write();
    
    
    
    */

    //bckgEstimate(st_sample, dirname, rd_90ias100, rd_90ias100, rbc_90ias100, ra_ias50, rd_90ias100, "True_90ias100", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, true);
    
    
    //bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100_nominal", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    //bckgEstimate(st_sample, dirname, rb_99ias100, rc_ias50, rbc_99ias100, ra_ias50, rd_99ias100, "99ias100_nominal", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    //bckgEstimate(st_sample, dirname, rb_999ias100, rc_ias50, rbc_999ias100, ra_ias50, rd_999ias100, "999ias100_nominal", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    
    //bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100_down", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    //bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100_up", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    

    //bckgEstimate(st_sample, dirname, rb_50ias99, rc_ias50, rbc_50ias99, ra_ias50, rd_50ias99, "50ias99", nPE);
    //bckgEstimate(st_sample, dirname, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", nPE);
    //bckgEstimate(st_sample, dirname, rb_99ias100, rc_ias50, rbc_99ias100, ra_ias50, rd_99ias100, "99ias100", nPE);
    //bckgEstimate(st_sample, dirname, rb_50ias999, rc_ias50, rbc_50ias999, ra_ias50, rd_50ias999, "50ias999", nPE);
    //bckgEstimate(st_sample, dirname, rb_999ias100, rc_ias50, rbc_999ias100, ra_ias50, rd_999ias100, "999ias100", nPE);


    //bckgEstimate(st_sample, dirname, rb_1o2, rc_1o2, rbc_50ias90, ra_1o4, rd_50ias90, "1o4", nPE);
    //bckgEstimate(st_sample, dirname, rb_2o2, rc_1o2, rbc_50ias90, ra_2o4, rd_50ias90, "2o4", nPE);
    //bckgEstimate(st_sample, dirname, rb_1o2, rc_2o2, rbc_50ias90, ra_3o4, rd_50ias90, "3o4", nPE);
    //bckgEstimate(st_sample, dirname, rb_2o2, rc_2o2, rbc_50ias90, ra_4o4, rd_50ias90, "4o4", nPE);

    //bckgEstimate(st_sample, dirname, rb_true, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "True_B", nPE);
    //bckgEstimate(st_sample, dirname, rb_50ias90, rc_true, rbc_50ias90, ra_ias50, rd_50ias90, "True_C", nPE);
    
    
    //bckgEstimate(st_sample, dirname, rb_true, rc_true, rbc_50ias90, ra_ias50, rd_50ias90, "True_BC", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP);
    
    //bckgEstimate(st_sample, dirname, rd_50ias90, rd_50ias90, rd_50ias90, rd_50ias90, rd_50ias90, "RegionD", nPE);
    
    // ------------------------------------------------------------------------------------------------------
   
    // bkg estimate for a selected cut index 
    // cutIndex = 3 --> pT > 60 GeV & Ias > 0.05
   

    //bckgEstimate_fromHistos(st_sample, dirname, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, nPE);

    delete ofile;
    //delete mass_cutIndex;

    return;
}
