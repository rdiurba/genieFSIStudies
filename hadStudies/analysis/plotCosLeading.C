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

// =======================================================
// Unified histogram builder
// =======================================================
TH1D* makeHistogram(TTree* t, std::string label, Mode mode, bool useThreshold,int A1, bool compound)
{

    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),
                          Form("hist_%s",label.c_str()),50,0.5,1);

    if (!t || t->GetEntries() < 1) return hist;

    int pdgf[100], nf, probe_fsi;
    double Ef[100], mh[100];
    int np; int nn; double cth[100];
    t->SetBranchAddress("pdgh",&pdgf);
    t->SetBranchAddress("nh",&nf);
    t->SetBranchAddress("Eh",&Ef);
    t->SetBranchAddress("mh",&mh);
    t->SetBranchAddress("np",&np);
    t->SetBranchAddress("nn",&nn);
    t->SetBranchAddress("cth",&cth);
    t->SetBranchAddress("probe_fsi",&probe_fsi);

    int required_fsi = (mode == kAbs) ? 5 : 6;
    double thresh = useThreshold ? 0.02 : 0.0;
    
   double maxKE=0.02; double maxCth=0;
    for (Long64_t j=0; j<t->GetEntries(); j++){
        t->GetEntry(j);
        if (probe_fsi!=5 && probe_fsi!=6) continue;
       int mult=0;
       int piMult=0;
       for(int i=0; i<nf; i++){
            int apdg = abs(pdgf[i]);
            double KE = Ef[i] - mh[i];

            if (apdg==211|| apdg==111){
                piMult++;
            }
            
            if (KE <= thresh) continue;

            if (apdg==2212 || apdg==2112){
                mult++;
            }
           
            if (apdg > 1000000000){
                int Z = (apdg/10000) % 1000;
                int A = (apdg/10) % 1000;
                int N = A - Z;

                   if (Z > 2) continue;
             if (KE<=thresh*A) continue;

             if (compound==true)    mult += Z + N;
            }
            if (KE>maxKE && apdg==2212){ maxKE=KE;
                    maxCth=cth[i];

                                       }
        }
        if (required_fsi==6 && (mult<3 || piMult!=0)) continue;
        if (required_fsi==5 && probe_fsi!=5 && piMult!=0 && mult<2) continue;

        if (useThreshold==0 && compound==false) mult=np+nn;
        if (maxKE>0.02) hist->Fill(maxCth);
    }

    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotCosLeading(int probe, double ke, int target,
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

    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.35);
    gStyle->SetOptFit(111);
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

    h1->Scale(1.0/h1->GetEntries());
    h2->Scale(1.0/h2->GetEntries());
    h3->Scale(1.0/h3->GetEntries());
    h4->Scale(1.0/h4->GetEntries());

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
        h1->GetXaxis()->SetTitle("cos(#theta_{p})");
    else
        h1->GetXaxis()->SetTitle("cos(#theta_{p})");
    
    h1->GetYaxis()->SetTitle("Fraction of Absorption Int.");
    if (mode!=kAbs)    h1->GetYaxis()->SetTitle("Fraction of Knockout Int.");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetYaxis()->SetRangeUser(0,1);

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


TLegend *leg = new TLegend(0.2, 0.6, 0.45, 0.85);
leg->SetHeader(Form("%s on %s with T=%s GeV",
                    particleLabel.c_str(),
                    targetStr.c_str(),
                    keStr.c_str()));
leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.03);     // smaller, fits entries


leg->AddEntry(h1,"hA2018","l");
leg->AddEntry(h2,"hN2018","l");
leg->AddEntry(h3,"INCL++","l");
leg->AddEntry(h4,"Geant4","l");

leg->Draw();


    // ---------------------------------------------------
    // Output naming (NO "merged")
    // ---------------------------------------------------
    std::string channelStr="cosLeadingThreshPionAbs";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "cosLeadingThreshPionAbs" :
            "cosLeadingPionAbs";
    if (mode==kKO)
        channelStr = useThreshold ?
            "cosLeadingThreshKO" :
            "cosLeadingSumKO";

    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%d.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound));


}