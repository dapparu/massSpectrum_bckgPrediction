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
    outfilename_ += "_UnB_v3_Data_v2";

    std::cout << outfilename_ << std::endl;

    bool bool_rebin=rebin;
    
    TFile* ifile = new TFile((filename+".root").c_str());

    // histograms used for the mass prediction
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

    std::cout << "background estimation... " << std::endl;


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

    delete ofile;

    return;
}
