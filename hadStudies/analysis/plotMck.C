#include "TH1.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPrincipal.h"
#include "TDecompChol.h"
#include "TEfficiency.h"
enum Mode { kKO, kAbs };
struct FitEntry {
    int target;
    int category; // 0 = Sum, 1 = Difference
    int quantity; // 0 = Gamma, 1 = Std, 2 = Mean

    double hA_s, hA_i;
    double hN_s, hN_i;
    double incl_s, incl_i;
    double g4_s, g4_i;
};



// =======================================================
// Unified histogram builder
// =======================================================
TH1D* makeHistogram(TTree* t, std::string label, Mode mode, bool useThreshold,int A1, bool compound)
{
    int min=0; int max=40;
    if (A1<40){ min=0; max=A1+2;}
    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),
                          Form("hist_%s",label.c_str()),max,min,max);

    if (!t || t->GetEntries() < 1) return hist;

    int pdgf[100], nf, probe_fsi;
    double Ef[100], mh[100];
    int np; int nn;
    t->SetBranchAddress("pdgh",&pdgf);
    t->SetBranchAddress("nh",&nf);
    t->SetBranchAddress("Eh",&Ef);
    t->SetBranchAddress("mh",&mh);
    t->SetBranchAddress("np",&np);
    t->SetBranchAddress("nn",&nn);
    t->SetBranchAddress("probe_fsi",&probe_fsi);
    
    int required_fsi = (mode == kAbs) ? 5 : 6;
    double thresh = useThreshold ? 0.005 : 0.0;
    int numScatters=0;
    for (Long64_t j=0; j<t->GetEntries(); j++){
        t->GetEntry(j);
        if (probe_fsi!=5 && probe_fsi!=6) continue;
       int mult=0;
       int piMult=0;
    int numP=0;
       for(int i=0; i<nf; i++){
            int apdg = abs(pdgf[i]);
            double KE = Ef[i] - mh[i];
            if (apdg==211|| apdg==111){
                piMult++;
            }
            
            if (KE <= thresh) continue;

            if (apdg==2212 || apdg==2112){
                mult++;
                if (apdg==2212) numP++;
            }
           
            if (apdg > 1000000000){
                int Z = (apdg/10000) % 1000;
                int A = (apdg/10) % 1000;
                int N = A - Z;

                   if (Z > 2) continue;
             if (KE<=thresh*A) continue;

             if (compound==true)    mult += Z + N;
            }
            
        }
        if (required_fsi==6 && (mult<3 || piMult!=0)) continue;
        if (required_fsi==5 && (piMult!=0 || mult<2)) continue;
        if (mult==2) numScatters++;

        if (useThreshold==0 && compound==false) mult=np+nn;
        if (np==2 && mult==2) hist->Fill(mult);
    }
    hist->Scale(1.f/numScatters);
    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotMck(int probe, double ke, int target,
                  bool useThreshold = false, bool compound=false)
{
    std::string probeStr;
    Mode mode = (probe==2212 || probe==-2212 || probe==2112) ? kKO : kAbs;
    bool isPion = (probe==211 || probe==-211 || probe==111);
    // ---------------------------------------------------
    // Probe string (for filenames)
    // ---------------------------------------------------
    if (mode == kAbs){
        if (probe==211) probeStr="piPlus";
        if (probe==-211) probeStr="piMinus";
        if (probe==111) probeStr="pi0";
    } else {
        if (probe==2212) probeStr="protonPlus";
        if (probe==-2212) probeStr="protonMinus";
        if (probe==2112) probeStr="neutron";
    }

    std::string keStr = std::to_string(ke).substr(0,5);
    std::string keLow=std::to_string(ke-0.025).substr(0, 4);
    std::string keHigh=std::to_string(ke+0.025).substr(0, 4);
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.45);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);

    // ---------------------------------------------------
    // Load files
    // ---------------------------------------------------
    TFile *f1 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_hA2018_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f2 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_hN2018_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f3 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_HINCL_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f4 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_HG4BertCasc_100k.ginuke.root",probeStr.c_str(),target,keStr.c_str()));

    int Z = (target/10000) % 1000;
    int A = (target/10) % 1000;
    double min=-20; double max=20;
    if (A<40){ min=0; max=12;}
                              

    TTree *t1 = (TTree*)f1->Get("ginuke");
    TTree *t2 = (TTree*)f2->Get("ginuke");
    TTree *t3 = (TTree*)f3->Get("ginuke");
    TTree *t4 = (TTree*)f4->Get("ginuke");

    // ---------------------------------------------------
    // Histograms
    // ---------------------------------------------------
    TH1D* h1 = makeHistogram(t1,"20i",mode,useThreshold,A,compound);
    TH1D* h2 = makeHistogram(t2,"20j",mode,useThreshold,A,compound);
    TH1D* h3 = makeHistogram(t3,"20k",mode,useThreshold,A,compound);
    TH1D* h4 = makeHistogram(t4,"20l",mode,useThreshold,A,compound);


    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);   h2->SetLineStyle(2);
    h3->SetLineColor(kBlue);  h3->SetLineStyle(4);
    h4->SetLineColor(kGreen); h4->SetLineStyle(3);

    // ---------------------------------------------------
    // Labels and titles
    // ---------------------------------------------------
    std::string title;
    if (mode==kAbs)
        h1->SetTitle("Pion Absorption");
    else
        h1->SetTitle("Nucleon Knockout");


    if (useThreshold)
        h1->GetXaxis()->SetTitle("N_{p}+N_{n} (T_{p,n}>5 MeV)");
    else
        h1->GetXaxis()->SetTitle("N_{p}+N_{n}");
    
    h1->GetYaxis()->SetTitle("Fraction of Absorption Int.");
    if (mode!=kAbs)    h1->GetYaxis()->SetTitle("Fraction of Knockout Int.");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    if (A==12) h1->GetYaxis()->SetRangeUser(0,1);
    else h1->GetYaxis()->SetRangeUser(0,0.6);
    if (h4->GetMaximum()>0.5) h1->GetYaxis()->SetRangeUser(0,1);
    h1->GetYaxis()->SetTitleOffset(1.3);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    h4->Draw("HIST SAME");

    // ---------------------------------------------------
    // Legend header (correct particle + target)
    // ---------------------------------------------------
    std::string particleLabel;

    if (mode == kAbs){
        if (probe==211) particleLabel = "#pi^{+}";
        if (probe==-211) particleLabel = "#pi^{-}";
        if (probe==111) particleLabel = "#pi^{0}";
    } else {
        if (probe==2212) particleLabel = "p^{+}";
        if (probe==-2212) particleLabel = "p^{-}";
        if (probe==2112) particleLabel = "n";
    }

    std::string targetStr = "Carbon";
    if (target==1000822080) targetStr="Lead";
    if (target==1000260560) targetStr="Iron";
    if (target==1000080160) targetStr="Oxygen";
    if (target==1000180400) targetStr="Argon";


    TLegend *leg = new TLegend(0.55, 0.65, 0.45, 0.88);
leg->SetHeader(Form("%s on %s at %s GeV<T<%s GeV",
                    particleLabel.c_str(),
                    targetStr.c_str(),
                    keLow.c_str(),keHigh.c_str()));
leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.035);     // smaller, fits entries


leg->AddEntry(h1,"hA2018","l");
leg->AddEntry(h2,"hN2018","l");
leg->AddEntry(h3,"INCL++","l");
leg->AddEntry(h4,"Geant4","l");

leg->Draw();

    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.025);
    c1->Update();  // commit before drawing

    // ---------------------------------------------------
    // Output naming (NO "merged")
    // ---------------------------------------------------
    std::string channelStr="nucleonSumThreshPionAbsMck";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "nucleonSumThreshPionAbsMck" :
            "nucleonSumPionAbsMck";
    if (mode==kKO)
        channelStr = useThreshold ?
            "nucleonSumThreshKO" :
            "nucleonSumKO";

    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%d.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound));

    // ---------------------------------------------------
    // Fit
    // ---------------------------------------------------
    double edge1 = -1.0;
    double edge2 = -1.0;
    double edge3 = -1.0;
    double edge4 = -1.0;
    
    TH1* hists[] = {h1, h2, h3, h4};
    double* edges[] = {&edge1, &edge2, &edge3, &edge4};
    
    

    TF1* g1 = (mode==kAbs) ? new TF1("g1","gaus",2,max) : new TF1("g1","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g2 = (mode==kAbs) ? new TF1("g2","gaus",2,max) : new TF1("g2","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g3 = (mode==kAbs) ? new TF1("g3","gaus",2,max) : new TF1("g3","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g4 = (mode==kAbs) ? new TF1("g4","gaus",2,max) : new TF1("g4","[0]*exp(-(x-3)*[1])",3,max);

    if (mode==kKO){
    /*g1->SetParLimits(0, 0,1E4);
    g2->SetParLimits(0, 0,1E4);
    g3->SetParLimits(0, 0,1E4);
    g4->SetParLimits(0, 0,1E4);
    g1->SetParLimits(1, 0,10);
    g2->SetParLimits(1, 0,10);
    g3->SetParLimits(1, 0,10);
    g4->SetParLimits(1, 0,10);*/
    }
    else{

    g1->SetParLimits(0, 0,1E4);
    g2->SetParLimits(0, 0,1E4);
    g3->SetParLimits(0, 0,1E4);
    g4->SetParLimits(0, 0,1E4);
    g1->SetParLimits(1, 2,20);
    g2->SetParLimits(1, 2,20);
    g3->SetParLimits(1, 2,20);
    g4->SetParLimits(1, 2,20);
    g1->SetParLimits(2, 0,10);
    g2->SetParLimits(2, 0,10);
    g3->SetParLimits(2, 0,10);
    g4->SetParLimits(2, 0,10);



    }
    h1->Fit(g1,"LRQ");
    h2->Fit(g2,"LRQ");
    h3->Fit(g3,"LRQ");
    h4->Fit(g4,"LRQ");



}