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

double GetVisEReweight(
    TH2D* hist_nom,
    TH2D* hist_alt,
    double KEini,
    double Ebias)
{
    if (!hist_nom || !hist_alt) return 1.;

    int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_Ebias = hist_nom->GetYaxis()->FindBin(Ebias);

    double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_Ebias);
    double weight_alt = hist_alt->GetBinContent(idx_KEini, idx_Ebias);
    std::cout<<weight_nom<<","<<weight_alt<<","<<KEini<<","<<Ebias<<std::endl;
    if (weight_nom == 0.) return 1.;


    return weight_alt / weight_nom;
}

// =======================================================
// Unified histogram builder
// =======================================================
TH1D* makeHistogram(TTree* t, std::string label, Mode mode, bool useThreshold,int A1, int target, bool compound, bool RW, std::string pstr="pionp",std::string reweightModel="INCL++")
{
    TH1D* hist = new TH1D(Form("hist_%s", label.c_str()),
                          Form("hist_%s", label.c_str()),
                          12, -0.2,1);
    TFile* f2 = TFile::Open("FSI_KOAbs_Evis_reweight_template.root", "READ");
        TH2D* hist_nom = (TH2D*)f2->Get(Form("hA2018_%s_KEini_vs_Ebias", pstr.c_str()));
    TH2D* hist_alt = (TH2D*)f2->Get(Form("%s_%s_KEini_vs_Ebias",     reweightModel.c_str(), pstr.c_str()));
    if (!t || t->GetEntries() < 1) return hist;
    int required_fsi = (mode == kAbs) ? 5 : 6;
    // Branch variables
    int pdgh[100], nh; 
    double Eh[100], mh[100];
    int probe_fsi;
    int npi0, npip, npim;
    double KEini, Eini;
    
    t->SetBranchAddress("pdgh",      pdgh);
    t->SetBranchAddress("npi0",   &npi0);
    t->SetBranchAddress("npip",   &npip);
    t->SetBranchAddress("npim",   &npim);
    t->SetBranchAddress("ke",       &KEini);
    t->SetBranchAddress("e",        &Eini);
    t->SetBranchAddress("nh",       &nh);
    t->SetBranchAddress("Eh",        Eh);
    t->SetBranchAddress("mh",        mh);
    t->SetBranchAddress("probe_fsi",&probe_fsi);

    double thresh = useThreshold ? 0.005 : 0.0;

    for (Long64_t j = 0; j < t->GetEntries(); j++) {
        t->GetEntry(j);
      if (Eh[0]==Eini) continue;

        int    mult   = 0;
        int    piMult = 0;
       for(int i=0; i<nh; i++){
            int apdg = abs(pdgh[i]);
            double KE = Eh[i] - mh[i];
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
            
        }
        if (required_fsi==6 && (mult<3 || piMult!=0)) continue;
        if (required_fsi==5 && (piMult!=0 || mult<2)) continue;
        double visE = 0.0;
        for (int i = 0; i < nh; ++i) {
            int    apdg = abs(pdgh[i]);
            double KE   = Eh[i] - mh[i];

            if (apdg == 211 || apdg == 111) {
                visE += Eh[i];          // pi±, pi0: total energy
            }
            else if (apdg == 22) {
                visE += Eh[i];          // gamma: total energy
            }
            else if (apdg == 2212) {
                if (KE > thresh) visE += KE;   // proton: KE only, threshold gated
            }
            else if (apdg == 2112) {
                ;                       // neutron: invisible, skip
            }
            else if (apdg > 1000000000) {
                int Z = (apdg / 10000) % 1000;
                int A = (apdg / 10)    % 1000;
                if (Z > 2)          continue;
                if (KE <= thresh * A) continue;
                if (compound) visE += KE;
            }
        }
        // -------------------------------------------------------------

        // Compute energy bias
        double Ehad = visE;                       // fixed: removed re-declaration of visE, named clearly
        double Ebias = 0.0;
        if (mode==kKO)
            Ebias = (KEini - Ehad) / KEini;
        else 
            Ebias = (Eini - Ehad) / Eini;
    double rw=1;
    if (RW) {


    double visERW=GetVisEReweight(hist_nom, hist_alt, KEini, Ebias);

     rw =visERW;

    }
        hist->Fill(Ebias,rw);
    }
    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotVisibleERW(int probe=211, double ke=0.575, int target=1000180400,
                   bool useThreshold = true,std::string reweightModel="INCL++")
{
    bool compound=true;
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
    std::string particle=probeStr;
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
    TH1D* h1 = makeHistogram(t1,"20i",mode,useThreshold,A,target,compound,1,particle,reweightModel);
    TH1D* h2 = makeHistogram(t2,"20j",mode,useThreshold,A,target,compound,0);
    TH1D* h3 = makeHistogram(t3,"20k",mode,useThreshold,A,target,compound,0);
    TH1D* h4 = makeHistogram(t4,"20l",mode,useThreshold,A,target,compound,0);

    h1->Scale(1.0/h1->Integral());
    h2->Scale(1.0/h2->Integral());
    h3->Scale(1.0/h3->Integral());
    h4->Scale(1.0/h4->Integral());

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


    h1->GetXaxis()->SetTitle("(E^{ini}-E^{vis})/E^{ini}");
    if (probe>211) h1->GetXaxis()->SetTitle("(T^{ini}-E^{vis})/T^{ini}");

    h1->GetYaxis()->SetTitle("Fraction of Absorption Int.");
    if (mode!=kAbs)    h1->GetYaxis()->SetTitle("Fraction of Knockout Int.");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetYaxis()->SetRangeUser(0,0.6);
    if (probe<212) h1->GetYaxis()->SetRangeUser(0,0.4);
    if (ke<0.350) h1->GetYaxis()->SetRangeUser(0,1);
    if (A>16) h1->GetYaxis()->SetRangeUser(0,0.6);
    h1->GetYaxis()->SetTitleOffset(1.3);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    h4->Draw("HIST SAME");

    // ---------------------------------------------------
    // Legend header (correct particle + target)
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

    std::string targetStr = "C";
    if (target==1000822080) targetStr="Pb";
    if (target==1000260560) targetStr="Fe";
    if (target==1000080160) targetStr="O";
    if (target==1000180400) targetStr="Ar";


    TLegend *leg = new TLegend(0.55, 0.65, 0.45, 0.88);
leg->SetHeader(Form("%s+%s (%s#font[122]{-}%s GeV)",
                    particleLabel.c_str(),
                    targetStr.c_str(),
                    keLow.c_str(), keHigh.c_str()));
leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.05);     // smaller, fits entries

leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.035);     // smaller, fits entries


leg->AddEntry(h1,Form("hA2018 Reweighted to %s",reweightModel.c_str()),"l");
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
    std::string channelStr="visESumThreshPionAbs";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "visESumThreshPionAbs" :
            "visESumPionAbs";
    if (mode==kKO)
        channelStr = useThreshold ?
            "visESumThreshKO" :
            "visESumKO";

    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%dRW%s.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound,reweightModel.c_str()));


}