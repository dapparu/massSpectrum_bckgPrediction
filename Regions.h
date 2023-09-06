#ifndef SUSYBSMAnalysis_Analyzer_Regions_h
#define SUSYBSMAnalysis_Analyzer_Regions_h

#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TDirectory.h"
#include <TRatioPlot.h>
#include <THStack.h>

using namespace std::placeholders;



//data 2017
//float K_data2017 = 2.30;
//float C_data2017 = 3.17;
float K_data2017 = 2.54;
float C_data2017 = 3.14;
//data 2018
//float K_data2018 = 2.27;
//float C_data2018 = 3.16;
float K_data2018 = 2.55;
float C_data2018 = 3.14;
//MC 2017
//float K_mc2017 = 2.26;
//float C_mc2017 = 3.22;
float K_mc2017 = 2.48;
float C_mc2017 = 3.19;
//MC 2018
//float K_mc2018 = 2.27;
//float C_mc2018 = 3.22;
float K_mc2018 = 2.49;
float C_mc2018 = 3.19;

//Systematic error due to the background estimate method
float systErr_ = 0.; //set to 0 for systematic studies

void underflowAndOverflow(TH1F* h){
    h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(0,0);
    h->SetBinContent(h->GetNbinsX()+1,0);
    return h;
}

// Scale the 1D-histogram given to the unit 
void scale(TH1F* h){
    h->Scale(1./h->Integral(0,h->GetNbinsX()+1));
}

// normalizes an histogram to a given norm 
void massNormalisation(TH1F* h, const float& normalisation){
    for(int k=0;k<h->GetNbinsX()+1;k++){
        h->SetBinContent(k,h->GetBinContent(k)*normalisation);
    }
}

void corrIh(TH2F* ih_eta){
    TF1 f_correlationPtIh("f_correlationPtIh","pol1",3,8);
    //f_correlationPtIh.SetParameter(0,1.215);
    //f_correlationPtIh.SetParameter(1,-6.313e-2);
    f_correlationPtIh.SetParameter(0,1.2);
    f_correlationPtIh.SetParameter(1,-5.3e-2);
    for(int i=0; i<ih_eta->GetNbinsX(); i++){
        for(int j=0; j<ih_eta->GetNbinsY(); j++){
            ih_eta->SetBinContent(i,j,ih_eta->GetBinContent(i,j)/f_correlationPtIh.Eval(ih_eta->GetYaxis()->GetBinCenter(j)));
        }
    }
}

void corrP(TH2F* eta_p){
    TF1 f_correlationPIas("f_correlationPIas","pol1",0,200);
    f_correlationPIas.SetParameter(0,9.8e-1);
    f_correlationPIas.SetParameter(1,2.4e-4);
    for(int i=0; i<eta_p->GetNbinsX(); i++){
        for(int j=0; j<eta_p->GetNbinsY(); j++){
            eta_p->SetBinContent(i,j,eta_p->GetBinContent(i,j)/f_correlationPIas.Eval(eta_p->GetXaxis()->GetBinCenter(i)));
        }
    }
}

void blindMass(TH1F* h_m,float mass_value=300){
    for(int i=0; i<h_m->GetNbinsX()+2; i++){
        if(h_m->GetBinLowEdge(i)>=mass_value) {
            h_m->SetBinContent(i,0);
            h_m->SetBinError(i,0);
        }
    }
}

// class using to definite signal and control regions. 
class Region{
    public:
        Region();
        Region(TFileDirectory &dir,std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins);
        ~Region();
        void setSuffix(std::string suffix);
        void initHisto(TFileDirectory &dir,int etabins,int ihbins,int pbins,int massbins);
        void fill(float& eta, float&p, float& pt, float& pterr, float& ih, float& ias, float& m, float& tof, float& w);
        void fillPredMass(const std::string&,TF1&,TF1&,const int&,const int&,float weight_);
        void write();

        float K_;
        float C_;


        int np;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int nias;
        float iaslow;
        float iasup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        std::vector<double> VectOfBins_P_;
        std::string suffix_;
        TH3F* ih_p_eta;
        TH2F* eta_p;
        TH2F* ih_eta;
        TH2F* ih_p;
        TH2F* ias_p;
        TH2F* ias_pt;
        TH1F* mass;
        TH1F* pred_mass;
        TH1F* pred_mass_correction;
        TH1F* pred_mass_fitIh;
        TH1F* pred_mass_fitP;
        TH1F* pred_mass_fitIh_fitP;
        TH1F* pred_mass_noFit;
        TH2F* eta_p_rebinned;
        TH2F* pt_pterroverpt;
        TH1F* hTOF;
        TH2F* ih_p_cross1D;
        TH2F* ih_p_cross1D_fit;
        TH2F* ih_p_cross1D_corr;
};

Region::Region(){}

Region::Region(TFileDirectory &dir, std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins){
    suffix_ = suffix;
    initHisto(dir,etabins,ihbins,pbins,massbins);
} 

Region::~Region(){}

void Region::setSuffix(std::string suffix){
    suffix_ = suffix;
}

// Function which intializes the histograms with given binnings 
void Region::initHisto(TFileDirectory &dir,int etabins,int ihbins,int pbins,int massbins){
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH3::SetDefaultSumw2(kTRUE);
    np = pbins;
    plow = 0;
    pup = 10000;
    npt = pbins;
    ptlow = 0;
    ptup = 10000; 
    nih = ihbins;
    ihlow = 0;
    ihup = 20;
    nias = ihbins;
    iaslow = 0;
    iasup = 1;
    neta = etabins;
    etalow = -3;
    etaup = 3;
    nmass = massbins;
    masslow = 0;
    massup = 4000;
    std::string suffix = suffix_;
    ih_p_eta = dir.make<TH3F>(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup); 
    eta_p = dir.make<TH2F>(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup); 
    ih_eta = dir.make<TH2F>(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup); 
    ih_p = dir.make<TH2F>(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D = dir.make<TH2F>(("ih_p_cross1D"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D_fit = dir.make<TH2F>(("ih_p_cross1D_fit"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D_corr = dir.make<TH2F>(("ih_p_cross1D_corr"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ias_p = dir.make<TH2F>(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup); 
    ias_pt = dir.make<TH2F>(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    mass = dir.make<TH1F>(("mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); 
    pred_mass = dir.make<TH1F>(("pred_mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); 
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pred_mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pt_pterroverpt = dir.make<TH2F>(("pt_pterroverpt"+suffix).c_str(),";p_{T} [GeV];#frac{#sigma_{pT}}{p_{T}}",npt,ptlow,ptup,100,0,1); 
    hTOF    = dir.make<TH1F>(("hTOF_"+suffix).c_str(),";TOF",200,-10,10); 
}

// Function which fills histograms
void Region::fill(float& eta, float& p, float& pt, float& pterr, float& ih, float& ias, float& m, float& tof, float& w){
   ih_p_eta->Fill(eta,p,ih,w);
   eta_p->Fill(p,eta,w);
   ih_eta->Fill(eta,ih,w);
   ih_p->Fill(p,ih,w);
   ias_p->Fill(p,ias,w);
   ias_pt->Fill(pt,ias,w);
   mass->Fill(m,w);
   pt_pterroverpt->Fill(pt,pterr/pt,w);
   hTOF->Fill(tof,w);
}

// in order to compute properly the uncertainties we use the methods SetBinContent SetBinError instead of Fill
// as several couples of bins in (p,ih) can provide the same mass estimate we need to properly sum the entries and errors
// for a couple of bins in (p,ih) where the bin content were (N_p,N_ih) the associated quantities should be 
// content: (N_p * N_ih) / N_total, where N_total represents the total number of events in the region (integral of p, ih & mass distributions)
// error: content * sqrt( 1 / N_p + 1 / N_ih ) where we assume Poisson uncertainties in both distributions (independent distributions) and we neglect the uncertainty on N_total
// While combining the input for several couples leading to the same mass: 
// contents are added 
// errors: the sqrt of the squared uncertainties are added
void Region::fillPredMass(const std::string& st_sample,TF1& f_p,TF1& f_ih,const int& fit_ih_err=1,const int& fit_p_err=1,float weight_=-1) {

    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    float K=2.27, C=3.16;
    //float K=2.27, C=3.22; //MC

    if(st_sample=="data2017"){K=K_data2017;C=C_data2017;}
    else if(st_sample=="data2018"){K=K_data2018;C=C_data2018;}
    else if(st_sample=="mc2017"){K=K_mc2017;C=C_mc2017;}
    else if(st_sample=="mc2018"){K=K_mc2018;C=C_mc2018;}

    bool useFitIh=true;
    bool useFitP=true;
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);

    for(int i=1;i<eta->GetNbinsX()+1;i++)
    {
        useFitIh=true;
        useFitP=true;
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i,"e");
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i,"e");

        scale(p); //only scale one of the two distributions ih or p --> keep the information of the normalisation 
        //scale(ih);
        if(ih->GetEntries()<1) continue;
        if(p->GetEntries()<1) continue;
        float a=3., b=8, c=0, d=30; //range fits
        TFitResultPtr ptr1 = 0;
        if(useFitIh) ptr1 = ih->Fit(&f_ih,"QRS","",a,b);

        if(ptr1->Status()!=0) useFitIh=false;
        
        TF1* const f_ih2=&f_ih;
        float intFih=f_ih2->Integral(a,b);
        float intIh=ih->Integral(ih->FindBin(a),ih->FindBin(b));

        float SFih=intIh/intFih;
        if(SFih<0) useFitIh=false;

        TF1* f_ih3=f_ih2;
        const double* fit_ih_params = ptr1->GetParams();
        double* fit_ih_cov = ptr1->GetCovarianceMatrix().GetMatrixArray();
        auto covMatrix_ih = ptr1->GetCovarianceMatrix();

        float SFp=10;
        TF1* f_p3;
        int incrFit=0;
        
        TFitResultPtr ptr2 = 0;
        int statusFit=1;
        
        while(statusFit!=0 && incrFit<5){
        incrFit++;

        if(useFitP) ptr2 = p->Fit(&f_p,"QRS","",c,d);
        TF1* f_p2=&f_p;
        
        float intFp=1;
        if(useFitP){ 
            ROOT::Math::IntegratorOneDim intOneDim_p(*f_p2,ROOT::Math::IntegrationOneDim::kGAUSS);
            intFp=intOneDim_p.Integral(c,d);
        }
        
        float intP=p->Integral(p->FindBin(c),p->FindBin(d));

        SFp=intP/intFp;
        f_p3=f_p2;
        
        statusFit=ptr2->Status();
        
        }
        if(statusFit!=0) useFitP=false;

        ROOT::Math::IntegratorOneDim intOneDimFP3(*f_p3,ROOT::Math::IntegrationOneDim::kGAUSS);
        
        useFitIh=true;
        useFitP=true;

        for(int j=1;j<p->GetNbinsX()+2;j++)
        {
            for(int k=1;k<ih->GetNbinsX()+2;k++)
            {
                float mom = p->GetBinLowEdge(j);
                float dedx = ih->GetBinLowEdge(k);
                double c_p = p->GetBinContent(j);
                double c_ih = ih->GetBinContent(k);
                float pLowEdge = p->GetBinLowEdge(j);
                float pUpEdge = p->GetBinLowEdge(j+1);
                float dedxLowEdge = ih->GetBinLowEdge(k);
                float dedxUpEdge = ih->GetBinLowEdge(k+1);

                double weight = 0;

                float invMom = 0;
                float mass = -1;
                int bin_mass = 0;

                float dedx_sampling = (dedxUpEdge-dedxLowEdge)/5.;
                float mom_sampling = (pUpEdge-pLowEdge)/5.;

                ih_p_cross1D->SetBinContent(j,k,ih_p_cross1D->GetBinContent(j,k)+(p->GetBinContent(j)*ih->GetBinContent(k)));

                // use Ih fit
                if(c_ih<100 && dedx>3.5 && useFitIh) {
                    for(double divdedx=dedxLowEdge; divdedx<dedxUpEdge; divdedx+=dedx_sampling){
                        
                        
                        c_ih = f_ih3->Integral(divdedx,divdedx+dedx_sampling);
                        c_ih *= SFih;
                        if(c_ih==0) continue;
                        
                        if (fit_ih_err==0) c_ih -= (SFih*f_ih3->IntegralError(divdedx,divdedx+dedx_sampling,fit_ih_params,fit_ih_cov,5e-2));
                        if (fit_ih_err==2) c_ih += (SFih*f_ih3->IntegralError(divdedx,divdedx+dedx_sampling,fit_ih_params,fit_ih_cov,5e-2));
                       
                        if(mom<20 && mom>0 && useFitP){
                            for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                                
                                c_p=intOneDimFP3.Integral(divmom,divmom+mom_sampling);
                                c_p *= SFp;
                                if(c_p==0) continue;

                                if (fit_p_err==0) c_p -= (SFp*intOneDimFP3.Error());
                                if (fit_p_err==2) c_p += (SFp*intOneDimFP3.Error());
                                
                                weight = c_ih * c_p;
                                
                                dedx = divdedx+dedx_sampling/2.;
                                invMom = 10000./(divmom+mom_sampling/2.);
                                mass = GetMass(invMom,dedx,K,C);
                                bin_mass = pred_mass->FindBin(mass);
                                pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                                pred_mass_fitIh_fitP->SetBinContent(bin_mass,pred_mass_fitIh_fitP->GetBinContent(bin_mass)+weight);
                                ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                            }
                        }
                        else{
                            c_p = p->GetBinContent(j);
                            weight = c_ih * c_p;
                            dedx = divdedx+dedx_sampling/2.;
                            invMom = 10000./p->GetBinCenter(j);
                            mass = GetMass(invMom,dedx,K,C);
                            bin_mass = pred_mass->FindBin(mass);
                            pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_fitIh->SetBinContent(bin_mass,pred_mass_fitIh->GetBinContent(bin_mass)+weight);
                            ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                        }
                    }
                }
                else{
                    // use 1/p fit
                    if(mom<20 && mom>0 && useFitP){
                        for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                            c_p=intOneDimFP3.Integral(divmom,divmom+mom_sampling);
                            c_p *= SFp;
                            if(c_p==0) continue;
                            if (fit_p_err==0) c_p -= (SFp*intOneDimFP3.Error());
                            if (fit_p_err==2) c_p += (SFp*intOneDimFP3.Error());
                            
                            weight = c_ih * c_p;
                            dedx = ih->GetBinCenter(k);
                            invMom = 10000./(divmom+mom_sampling/2.);
                            mass = GetMass(invMom,dedx,K,C);
                            bin_mass = pred_mass->FindBin(mass);
                            pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_fitP->SetBinContent(bin_mass,pred_mass_fitP->GetBinContent(bin_mass)+weight);
                            ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                        }
                    }
                    else{
                        c_p = p->GetBinContent(j);
                        c_ih = ih->GetBinContent(k);
                        weight = c_ih * c_p;
                        
                        dedx = ih->GetBinCenter(k);
                        invMom = 10000./p->GetBinCenter(j);
                        mass = GetMass(invMom,dedx,K,C);
                        bin_mass = pred_mass->FindBin(mass);
                        pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                        pred_mass_noFit->SetBinContent(bin_mass,pred_mass_noFit->GetBinContent(bin_mass)+weight);
                        ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                    }
                }
            }
        }
        delete p;
        delete ih;
    }
    delete eta;
}

void Region::write(){
    ih_p_eta->Write();
    eta_p->Write();
    ih_eta->Write();
    ih_p->Write();
    ih_p_cross1D->Write();
    ih_p_cross1D_fit->Write();
    ias_p->Write();
    ias_pt->Write();
    mass->Write();
    pred_mass->Write();
    pred_mass_fitIh->Write();
    pred_mass_fitP->Write();
    pred_mass_fitIh_fitP->Write();
    pred_mass_noFit->Write();
    pt_pterroverpt->Write();
    hTOF->Write();
}

void loadHistograms(Region& r, TFile* f, const std::string& regionName, bool bool_rebin=true, int rebineta=1, int rebinp=1, int rebinih=1, int rebinmass=1, std::string labelTest=""){
    std::string dir = "analyzer/BaseName/";
    dir="";
    //dir = "HSCParticleAnalyzer/SigmaPt3_iso2_IhCut1_PtCut1/";
    //dir = "HSCParticleAnalyzer/ih0_sigmapt2_iso1/";
    dir = "HSCParticleAnalyzer/BaseName/";
    cout<<"loading region "<<regionName<<endl;
    //r.ih_p_eta                          = (TH3F*)f->Get((dir+"ih_p_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p_eta->Rebin3D(rebineta,rebinp,rebinih);
    r.ih_p_eta                          = NULL;
    r.eta_p                             = (TH2F*)f->Get((dir+"eta_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.eta_p->Rebin2D(rebinp,rebineta);
    //r.p_npv                             = (TH2F*)f->Get((dir+"p_npv_"+regionName).c_str())->Clone(); if(bool_rebin) r.p_npv->Rebin2D(1,rebinp);
    r.ih_eta                            = (TH2F*)f->Get((dir+"ih_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_eta->Rebin2D(rebineta,rebinih);
    r.ih_p                              = (TH2F*)f->Get((dir+"ih_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p->Rebin2D(rebinp,rebinih);
    r.ih_p_cross1D                      = (TH2F*)r.ih_p->Clone(); r.ih_p_cross1D->Reset(); r.ih_p_cross1D->SetName(("cross1D_"+regionName).c_str());
    r.ih_p_cross1D_fit                  = (TH2F*)r.ih_p->Clone(); r.ih_p_cross1D_fit->Reset(); r.ih_p_cross1D_fit->SetName(("cross1D_fit_"+regionName).c_str());
    r.ih_p_cross1D_corr                 = (TH2F*)r.ih_p->Clone(); r.ih_p_cross1D_corr->Reset(); r.ih_p_cross1D_corr->SetName(("cross1D_corr_"+regionName).c_str());
    r.ias_p                             = (TH2F*)f->Get((dir+"ias_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_p->Rebin2D(rebinp,rebinih);
    r.ias_pt                            = (TH2F*)f->Get((dir+"ias_pt_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_pt->Rebin2D(rebinp,rebinih);
    r.mass                              = (TH1F*)f->Get((dir+"mass_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    //r.mass                              = (TH1F*)f->Get((dir+"massFromTree_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    //r.pred_mass                         = (TH1F*)f->Get((dir+"pred_mass_"+regionName).c_str())->Clone(); r.pred_mass->Reset(); if(bool_rebin) r.pred_mass->Rebin(rebinmass);
    r.pred_mass                         = (TH1F*) r.mass->Clone(); r.pred_mass->SetName(("pred_mass_"+regionName).c_str()); r.pred_mass->Reset();
    r.pred_mass_fitIh                   = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitIh->SetName(("pred_mass_fitIh_"+regionName).c_str());
    r.pred_mass_fitP                    = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitP->SetName(("pred_mass_fitP_"+regionName).c_str());
    r.pred_mass_fitIh_fitP              = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitIh_fitP->SetName(("pred_mass_fitIh_fitP_"+regionName).c_str());
    r.pred_mass_noFit                   = (TH1F*) r.pred_mass->Clone(); r.pred_mass_noFit->SetName(("pred_mass_noFit_"+regionName).c_str());
    //r.pred_mass                         = (TH1F*)f->Get((dir+"massFrom1DTemplatesEtaBinning_"+regionName).c_str())->Clone(); r.pred_mass->Reset(); if(bool_rebin) r.pred_mass->Rebin(rebinmass);
    //r.pred_mass_correction              = (TH1F*)f->Get((dir+"massFrom1DTemplatesEtaBinning_"+regionName).c_str())->Clone(); r.pred_mass_correction->Reset(); if(bool_rebin) r.pred_mass_correction->Rebin(rebinmass);
}

// Return randomly select histo 
TH1F* poissonHisto(const TH1F& h,TRandom3* RNG){
    TH1F* hres = (TH1F*) h.Clone();
    for(int i=0;i<h.GetNbinsX()+1;i++){
        hres->SetBinContent(i,RNG->Poisson(h.GetBinContent(i)));
    }
    return hres;
}
TH2F* poissonHisto(const TH2F& h,TRandom3* RNG){
    TH2F* hres = (TH2F*) h.Clone();
    for(int i=0;i<h.GetNbinsX()+1;i++){
        for(int j=0;j<h.GetNbinsY()+1;j++){
            hres->SetBinContent(i,j,RNG->Poisson(h.GetBinContent(i,j)));
        }
    }
    return hres;
}

// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
void etaReweighingP(TH2F* eta_p_1, const TH1F* eta2_)
{
    TH1F* eta1 = (TH1F*) eta_p_1->ProjectionY(); 
    TH1F* eta2 = (TH1F*) eta2_->Clone();
    eta1->Scale(1./eta1->Integral(0,eta1->GetNbinsX()+1));
    eta2->Scale(1./eta2->Integral(0,eta2->GetNbinsX()+1));
    eta2->Divide(eta1);
    for(int i=0;i<eta_p_1->GetNbinsX()+2;i++)
    {
        for(int j=0;j<eta_p_1->GetNbinsY()+2;j++)
        {
            float val_ij = eta_p_1->GetBinContent(i,j);
            float err_ij = eta_p_1->GetBinError(i,j);
            
            eta_p_1->SetBinContent(i,j,val_ij*eta2->GetBinContent(j));
            eta_p_1->SetBinError(i,j,err_ij*eta2->GetBinContent(j));
        }
    }
}

// Function doing the nof. primary vertices reweighing between two 2D-histograms, as done with eta 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
void npvReweighingP(TH2F* npv_p_1, TH1F* npv2)
{
    TH1F* npv1 = (TH1F*) npv_p_1->ProjectionX(); 
    npv1->Scale(1./npv1->Integral());
    npv2->Scale(1./npv2->Integral());
    npv2->Divide(npv1);
    for(int i=0;i<npv_p_1->GetNbinsX()+1;i++)
    {
        for(int j=0;j<npv_p_1->GetNbinsY()+1;j++)
        {
            float val_ij = npv_p_1->GetBinContent(i,j);
            float err_ij = npv_p_1->GetBinError(i,j);
            npv_p_1->SetBinContent(i,j,val_ij*npv2->GetBinContent(j));
            npv_p_1->SetBinError(i,j,err_ij*npv2->GetBinContent(j));
        }
    }
}



// add the overflow bin to the last one
void overflowLastBin(TH1F* h){
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(h->GetNbinsX()+1,0);
}

// rebinning histogram according to an array of bins
TH1F* rebinHisto(TH1F* h){
    //overflowLastBin(h);
    double xbins[33]={0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,405.,435.,475.,525.,585.,660.,755.,875.,1025.,1210.,1440.,1730.,2000.};
    std::vector<double> xbins_v;
    for(double i=0.0;i<=1000.0;i+=50) xbins_v.push_back(i);
    std::string newname = h->GetName(); 
    newname += "_rebinned";
    TH1F* hres = (TH1F*) h->Rebin(32,newname.c_str(),xbins);
    //overflowLastBin(hres);
    return hres;
}

// Function returning the ratio of right integer (from x to infty) for two 1D-histograms
// This function is used in the Hscp data-driven background estimate to test the mass shape prediction
// The argument to use this type of ratio is that we're in case of cut & count experiment 
TH1F* ratioIntegral(TH1F* h1, TH1F* h2){    
    float SystError = systErr_;
    TH1F* res = (TH1F*) h1->Clone(); res->Reset();
    for(int i=1;i<h1->GetNbinsX()+1;i++)
    {   
        double Perr=0, Derr=0;
        double P=h1->IntegralAndError(i,h1->GetNbinsX()+1,Perr); if(P<=0) continue;
        double D=h2->IntegralAndError(i,h2->GetNbinsX()+1,Derr);
        Perr = sqrt(Perr*Perr + pow(P*SystError,2));
        res->SetBinContent(i,D/P);
        res->SetBinError(i,sqrt(pow(Derr*P,2)+pow(Perr*D,2))/pow(P,2));
    }
    return res;
}

TH1F* pull(TH1F* h1, TH1F* h2){
    float SystError = systErr_;
    TH1F* res = (TH1F*) h2->Clone(); //res->Reset();
    res->Divide(h1);

    return res;
}

void saveHistoRatio(TH1F* h1,TH1F* h2,std::string st1,std::string st2,std::string st3,bool rebin=false){
    h1->SetName(st1.c_str());
    h2->SetName(st2.c_str());
    if(rebin){
        h1=rebinHisto(h1);
        h2=rebinHisto(h2);
    }
    h1->Write();
    h2->Write();
    TH1F* R = (TH1F*) ratioIntegral(h2,h1)->Clone();
    if(rebin) st3+="_rebinned";
    R->SetName(st3.c_str());
    R->Write();
}

TH1F meanHistoPE(std::vector<TH1F> vPE){
    float SystError = systErr_;
    TH1F h=TH1F(vPE[0]);
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    for(int i=0;i<h.GetNbinsX()+1;i++){
        float mean=0, err=0;
        for(unsigned int pe=0;pe<vPE.size();pe++){
            mean += vPE[pe].GetBinContent(i);
        }
        mean /= vPE.size();
        for(unsigned int pe=0;pe<vPE.size();pe++){
            err += pow(mean - vPE[pe].GetBinContent(i),2);
        }
        float fact=1;
        if(vPE.size()>1) {err = sqrt(err/(vPE.size()-1));fact=vPE.size()/(vPE.size()-1);}
        else err = sqrt(err);
        h.SetBinContent(i,mean);
        h.SetBinError(i,err);
    }
    return h;
}

static std::vector<float> v_par0;
static std::vector<float> v_err0;
static std::vector<float> v_par1;
static std::vector<float> v_err1;

// Function returning a canvas divide in two windows
// The first one contains the two 1D-histograms, given as arguments, superposed.
// There also is a legend associated to this window where the names are defined as arguments.
// The second window contains the ratio of these 1D-histograms or the ratio of right integers of them. 
// We define which kind of ratio we want with tha 'ratioSimple' boolean.
// The 'name' given corresponds to the name of the canvas 
TCanvas* plotting(TH1F* h1, TH1F* h2, bool ratioSimple=true, std::string dirname="", std::string name="", std::string leg1="", std::string leg2="", bool rebin=false,bool logy=true, bool normalised=false){
    if(normalised) h1->Scale(1./h1->Integral());
    if(normalised) h2->Scale(1./h2->Integral());
    if(rebin) h1=rebinHisto(h1);
    if(rebin) h2=rebinHisto(h2);
    std::string canvName;
    if(rebin) {
      canvName = "plotting_"+name+"_rebin";
    } else {
      canvName = "plotting_"+name;
    }
    TCanvas* c1 = new TCanvas(canvName.c_str(),"", 800,800);
    c1->Divide(1,3);
    gStyle->SetOptStat(0);
    c1->cd(1);
    TPad* p1 = (TPad*)(c1->cd(1));
    if(logy) p1->SetLogy();
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h1,leg1.c_str(),"lep");
    leg->AddEntry(h2,leg2.c_str(),"lep");
    h1->SetStats(0);
    h1->Draw();
    if(logy) h1->GetYaxis()->SetRangeUser(1e-8,h1->GetMaximum()*1.5);
    else h1->GetYaxis()->SetRangeUser(h1->GetMinimum(),h1->GetMaximum());
    h1->SetLineColor(1);
    h2->SetLineColor(2);
    h2->SetStats(0);
    h2->Draw("esame");
    leg->Draw("same");
    h1->SetName((name+"_obs").c_str());
    h1->Write();
    h2->SetName((name+"_pred").c_str());
    h2->Write();
    cout<<"chi2: "<<h1->Chi2Test(h2)<<endl;
    c1->cd(2);
    TH1F* tmp = (TH1F*) h1->Clone(); tmp->Reset();
    if(ratioSimple){
        tmp = (TH1F*)h1->Clone();
        tmp->Divide(h2);
        tmp->GetYaxis()->SetTitle("#frac{N_{obs}}{N_{pred}}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }else{
        tmp=ratioIntegral(h2,h1);
        tmp->GetYaxis()->SetTitle("#int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }
    tmp->GetYaxis()->SetRangeUser(0,2);
    tmp->Draw();
    tmp->SetLineColor(1);
    tmp->SetMarkerColor(1);
    tmp->SetName((name+"_ratioInt").c_str());
    tmp->Write();
    c1->cd(3);
    TH1F* tmp2 = (TH1F*) pull(h2,h1)->Clone();
    tmp2->Draw("E1");
    tmp2->GetYaxis()->SetTitle("obs / pred");
    tmp2->GetYaxis()->SetRangeUser(0,2);
    tmp2->SetLineColor(1);
    tmp2->SetMarkerColor(1);
    TFitResultPtr ptr = tmp2->Fit("pol1","RS","same",25,300);
    TF1* fLin = tmp2->GetFunction("pol1");
    v_par0.push_back(fLin->GetParameter(0));
    v_err0.push_back(fLin->GetParError(0));
    v_par1.push_back(fLin->GetParameter(1));
    v_err1.push_back(fLin->GetParError(1));
    tmp2->SetName((name+"_ratio").c_str());
    tmp2->Write();
    c1->SaveAs((dirname+canvName+".pdf").c_str());
    c1->SaveAs((dirname+canvName+".root").c_str());
    return c1;
}

void bckgEstimate(const std::string& st_sample, const std::string& dirname, const Region& B, const Region& C, const Region& BC, const Region& A, const Region& D, const std::string& st, const int& nPE=200, const bool& corrTemplateIh=false, const bool& corrTemplateP=false, const int& fitIh=1, const int& fitP=1, bool blind=false, const int& rebinMass=1){

    std::vector<TH1F> vPE_;
    std::vector<TH1F> vPE_corr;
    std::vector<TH2F> vPE_cross1D;
    std::vector<TH2F> vPE_cross1D_corr;  
    ROOT::EnableImplicitMT(true);

    Region bc = BC;
    Region d = D;


    //underflowAndOverflow(d.mass);

    Region a = A;
    Region b = B;
    Region c = C;
    //scale(d.mass);
    TH2F a_ih_eta_base(*a.ih_eta);
    TH2F b_ih_eta_base(*b.ih_eta);
    TH2F c_ih_eta_base(*c.ih_eta);
    TH2F b_eta_p_base(*b.eta_p);
    TH2F c_eta_p_base(*c.eta_p);
    //TH2F b_p_npv_base(*b.p_npv);
    //TH2F c_p_npv_base(*c.p_npv);
    TH1F d_mass_base(*d.mass);


    auto workItem = [] (UInt_t workerID,const std::string& st_sample, const std::string& dirname, const Region& B, const Region& C, const Region& BC, const Region& A, const Region& D, const std::string& st, const int& nPE=200, const bool& corrTemplateIh=false, const bool& corrTemplateP=false, const int& fitIh=1, const int& fitP=1, bool blind=false, const int& rebinMass=1){

        cout<<"workerId: "<<workerID<<endl;

        TF1 f_p("f_p","[0]*([1]+erf((log(x)-[2])/[3]))",0,90);
        f_p.SetParameter(0,1);
        f_p.FixParameter(1,1.0);
        f_p.SetParameter(2,5.10116e+00);
        f_p.SetParameter(3,1.10152e+00);
        TF1 f_ih("f_ih","gaus",3.25,8);
        f_ih.SetParameter(0,7.05442e-03);
        f_ih.SetParameter(1,3.34e+00);
        f_ih.SetParameter(2,1.59e-01);

        Region bc = BC;
        Region d = D;

        //underflowAndOverflow(d.mass);

        Region a = A;
        Region b = B;
        Region c = C;
        //scale(d.mass);
        TH2F a_ih_eta_base(*a.ih_eta);
        TH2F b_ih_eta_base(*b.ih_eta);
        TH2F c_ih_eta_base(*c.ih_eta);
        TH2F b_eta_p_base(*b.eta_p);
        TH2F c_eta_p_base(*c.eta_p);
        TH1F d_mass_base(*d.mass);

        if(corrTemplateIh) corrIh(&b_ih_eta_base);
        if(corrTemplateP) corrP(&c_eta_p_base);

        float normalisationABC=B.ih_eta->GetEntries()*C.eta_p->GetEntries()/A.ih_eta->GetEntries();

        if(st_sample=="data2017"){bc.K_=K_data2017;bc.C_=C_data2017;}
        else if(st_sample=="data2018"){bc.K_=K_data2018;bc.C_=C_data2018;}
        else if(st_sample=="mc2017"){bc.K_=K_mc2017;bc.C_=C_mc2017;}
        else if(st_sample=="mc2018"){bc.K_=K_mc2018;bc.C_=C_mc2018;}

        TH1F* b_eta_base = (TH1F*)b_eta_p_base.ProjectionY();

        TRandom3* RNG = new TRandom3(workerID);
        bc.pred_mass->Reset();
        TH1F* d_mass = &d_mass_base;
        TH2F* a_ih_eta = poissonHisto(a_ih_eta_base,RNG);
        TH2F* b_ih_eta = poissonHisto(b_ih_eta_base,RNG);
        TH2F* b_eta_p = poissonHisto(b_eta_p_base,RNG);
        TH2F* c_eta_p = poissonHisto(c_eta_p_base,RNG);
        
        TH1F* b_ih = (TH1F*)b_ih_eta->ProjectionY();
        TH1F* b_eta = poissonHisto(*b_eta_base,RNG);

        etaReweighingP(c_eta_p,b_eta);

        bc.eta_p = c_eta_p;    bc.ih_eta = b_ih_eta;
        float normA = a_ih_eta->Integral(0,a_ih_eta->GetNbinsX()+1);
        float normB = b_ih_eta->Integral(0,b_ih_eta->GetNbinsX()+1);
        float normC = c_eta_p->Integral(0,c_eta_p->GetNbinsX()+1);

        bc.fillPredMass(st_sample,f_p,f_ih,fitIh,fitP);
        normalisationABC=normB*normC/normA;
        bc.pred_mass->Scale(normalisationABC/bc.pred_mass->Integral());
        float sfMass = d_mass->Integral()/bc.pred_mass->Integral();

        delete a_ih_eta;
        delete b_ih_eta;
        delete b_eta_p;
        delete c_eta_p;
        delete b_eta;
        return *bc.pred_mass;
    };

    
    // Create the pool of workers
    ROOT::TProcessExecutor workers(25);
    // Fill the pool with work
    auto workItemToRun = std::bind (workItem, _1, st_sample,dirname,B,C,BC,A,D,st,nPE,corrTemplateIh,corrTemplateP,fitIh,fitP,blind,rebinMass);
    auto vPE = workers.Map(workItemToRun, ROOT::TSeqI(nPE));

    TH1F h_temp = meanHistoPE(vPE);
    if(nPE>1) bc.pred_mass = &h_temp;

    TH1F* htemp = &h_temp;

    if(blind) {
        blindMass(d.mass,300);
    }
     
    saveHistoRatio(d.mass,bc.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str());
    //saveHistoRatio(d.mass,bc.pred_mass_correction,("mass_obs_"+st).c_str(),("mass_predBC_corr_"+st).c_str(),("mass_predBCR_corr_"+st).c_str());
    //saveHistoRatio(d.mass,bc.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),true);
    //saveHistoRatio(d.mass,bc.pred_mass_correction,("mass_obs_"+st).c_str(),("mass_predBC_corr_"+st).c_str(),("mass_predBCR_corr_"+st).c_str(),true);
    
    overflowLastBin(d.mass);
    overflowLastBin(bc.pred_mass);

    A.ih_eta->Write();
    A.eta_p->Write();
    B.ih_eta->Write();
    D.ih_eta->Write();
    C.eta_p->Write();
    D.eta_p->Write();

    b_ih_eta_base.Write();
    b_ih_eta_base.ProjectionY("_ih_py")->Write();
    c_eta_p_base.Write();
    c_eta_p_base.ProjectionX("_p_px")->Write();
    
    D.ih_eta->ProjectionY()->Write();
    B.ih_eta->ProjectionY()->Write();
    D.ih_eta->ProjectionX()->Write();
    B.ih_eta->ProjectionX()->Write();
    D.eta_p->ProjectionX()->Write();
    C.eta_p->ProjectionX()->Write();   
    C.eta_p->ProjectionY()->Write();   


    plotting(d.mass,bc.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",true,true,false)->Write();
    //plotting(d.mass,bc.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",false,true,true)->Write();
    //plotting(d.mass,bc.pred_mass_correction,false,dirname,("mass1D_regionBC_corr"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction")->Write();
    //plotting(d.mass,bc.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",true)->Write();
    //plotting(d.mass,bc.pred_mass_correction,false,dirname,("mass1D_regionBC_corr_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",true)->Write();
    //plotting((TH1F*)D.ih_eta->ProjectionY(),(TH1F*)BC.ih_eta->ProjectionY(),false,dirname,("ih_template_"+st).c_str(),"Region D","Template from B",true)->Write();
    //plotting((TH1F*)D.eta_p->ProjectionX(),(TH1F*)BC.eta_p->ProjectionX(),false,dirname,("p_template_"+st).c_str(),"Region D","Template from C",true)->Write();
    
}

#endif
