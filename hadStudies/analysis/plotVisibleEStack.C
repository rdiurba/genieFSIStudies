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
TH1D* makeHistogram(TTree* t, std::string label, int mode, int probe,
                    bool useThreshold, bool compound)
{
    TH1D* hist = new TH1D(Form("hist_%s", label.c_str()),
                          Form("hist_%s", label.c_str()),
                          12, -0.2,1);

    if (!t || t->GetEntries() < 1) return hist;

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

        int    mult   = 0;
        int    piMult = 0;
        double visE   = 0.0;  // fixed: was uninitialised, accumulating across events

        for (int i = 0; i < nh; i++) {           // fixed: nf/pdgf/Ef -> nh/pdgh/Eh
            int    apdg = abs(pdgh[i]);
            double KE   = Eh[i] - mh[i];

            if (apdg == 211) {
                //if (KE <= 0.5 * thresh) continue;
                visE += Eh[i];                    // fixed: Ef -> Eh

            }
            if (apdg == 22) {
                //if (KE <= 0.5 * thresh) continue;
                visE += Eh[i];                    // fixed: Ef -> Eh

            }
            
            if (apdg == 111) {
                visE += Eh[i];                    // fixed: Ef -> Eh
            }
            if (KE <= thresh) continue;

            if (apdg == 2212 || apdg == 2112) {
                if (apdg == 2212) visE += KE;
            }

            if (apdg > 1000000000) {
                int Z = (apdg / 10000) % 1000;
                int A = (apdg / 10)    % 1000;
                int N = A - Z;

                if (Z > 2)          continue;
                if (KE <= thresh*A) continue;


                   if (compound==1){ 
                                        visE += KE;}
            }
        }


        // =========================
        // Recompute probe_fsi
        // =========================
        int fsi = -1;
        
        if (Eini==Eh[0]) {
            fsi = 1;
        }
        else if ((abs(probe) == 211 || abs(probe) == 111) &&
                 (piMult + npi0 == 0)) {
            fsi = 5; // Absorption
        }
        else if ((abs(probe) == 2212 || abs(probe) == 2112) &&
                 (piMult + npi0 == 0) && mult > 2) {  // fixed: missing opening {
            fsi = 6; // Knock-out
        }
        else if ((abs(probe)==211 && piMult + npi0 >1) || (abs(probe)==2212 && piMult+npi0>0)) {
            fsi = 7;
        }

        else if (mult + npi0 + piMult >= 2) {            
            fsi = 2;
            //std::cout<<""<<std::endl;
            for (int i = 0; i < nh; ++i) {
               // std::cout<<Eh[i]-mh[i]<<","<<pdgh[i]<<std::endl;
                if ((abs(probe) == 211 || abs(probe) == 111) && probe == pdgh[i]) {
                    fsi = 4; 
                }
                if ((abs(probe) == 2212 || abs(probe) == 2112)) {
                    fsi = 4; 
                }
            }
        }
        else if (abs(probe) == 211 || abs(probe) == 111) {  // fixed: bare else if{ -> else if (condition)
            for (int j = 0; j < nh; ++j) {
                if ((probe ==  211 && pdgh[j] == -211) ||
                    (probe == -211 && pdgh[j] ==  211)) {
                    fsi = 8;
                }
            }
        }
        else {                                    // fixed: mismatched braces, restructured
            if (((mult== 1 && probe>2100) || ((piMult==1 || npi0==1) && probe<212)) && ( fsi != 1))
                fsi = 3;                          // fixed: fsi==3 -> fsi=3
        }

        // --- enforce match ---
        if (compound==0 && useThreshold==0) fsi=probe_fsi;
        if (fsi != mode) continue;
        // Compute energy bias
        double Ehad = visE;                       // fixed: removed re-declaration of visE, named clearly
        double Ebias = 0.0;
        if (probe == 2212 || probe == 2112)
            Ebias = (KEini - Ehad) / KEini;
        else if (probe == 211 || probe == 111 || probe == -211)
            Ebias = (Eini - Ehad) / Eini;

        hist->Fill(Ebias);
    }
    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotVisibleEStack(int probe, double ke, int target,
                   bool useThreshold = false, bool compound = false, std::string model="hA2018")
{
    Mode mode = (probe == 2212 || probe == -2212 || probe == 2112) ? kKO : kAbs;

    // ---------------------------------------------------
    // Probe string (for filenames)
    // ---------------------------------------------------
    std::string probeStr;
    if (mode == kAbs) {
        if (probe ==  211) probeStr = "piPlus";
        if (probe == -211) probeStr = "piMinus";
        if (probe ==  111) probeStr = "pi0";
    } else {
        if (probe ==  2212) probeStr = "protonPlus";
        if (probe == -2212) probeStr = "protonMinus";
        if (probe ==  2112) probeStr = "neutron";
    }

    std::string keStr = std::to_string(ke).substr(0, 5);
    std::string keLow=std::to_string(ke-0.025).substr(0, 4);
    std::string keHigh=std::to_string(ke+0.025).substr(0, 4);
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.45);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
    c1->SetLeftMargin(0.18);
    // ---------------------------------------------------
    // Load file
    // ---------------------------------------------------
    TFile *f1 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_%s_100k.ginuke.root",
                               probeStr.c_str(), target, keStr.c_str(),model.c_str()));

    TTree *t1 = (TTree*)f1->Get("ginuke");
// ---------------------------------------------------
    // Build per-FSI-mode histograms for the stack
    //   fsi=2  quasi-elastic scatter
    //   fsi=3  elastic
    //   fsi=4  charge exchange
    //   fsi=5  absorption (pions)  /  fsi=6 knock-out (nucleons)
    //   fsi=7  pion production
    //   fsi=8  double charge exchange
    // ---------------------------------------------------

    const int fsiModes[]    = {  4,                 2,
                                 (mode==kAbs ? 5 : 6),
                                 7,                8  };
    const char* fsiLabels[] = { "Inelastic",
                                "Charge Exchange",
                                (mode==kAbs ? "Absorption" : "Knock-out"),
                                "Pion Production",
                                "Double Charge Exchange" };
    const int fsiColors[]   = { kRed+1, kAzure+1, kGreen+2,
                                 kOrange+1, kMagenta+1};
  
    const int nModes = 5;

    // Total histogram (for normalisation)
    TH1D* hTotal = new TH1D("hTotal", Form("%s;(E^{ini}-E^{vis}})/E^{ini};Fraction of Interactions",
             model.c_str()), 12, -0.2,1);
    if (probe>211) hTotal->GetXaxis()->SetTitle("(T^{ini}-E^{vis})/T^{ini}");
    hTotal->GetYaxis()->SetRangeUser(0,0.5);
    THStack* stack = new THStack("stack",
        Form("%s;(E^{#mathrm{vis}}-E^{initial})/E^{initial};Fraction of Interactions", model.c_str()));
    hTotal->SetLineWidth(0);
    hTotal->GetYaxis()->SetRangeUser(0,0.5);
    TH1D* hFSI[nModes];
    for (int m = 0; m < nModes; m++) {
        hFSI[m] = makeHistogram(t1,
                                Form("fsi%d", fsiModes[m]),
                                fsiModes[m],   // mode = FSI channel
                                probe,
                                useThreshold,
                                compound);
        hFSI[m]->SetFillColor(fsiColors[m]);
        hFSI[m]->SetLineColor(fsiColors[m]);
        hFSI[m]->SetLineWidth(0);
        hFSI[m]->SetFillStyle(1001);
        hTotal->Add(hFSI[m]);
    }

    // Normalise each slice to total then add to stack
    double totalEntries = hTotal->Integral();
    for (int m = 0; m < nModes; m++) {
        if (totalEntries > 0)
            hFSI[m]->Scale(1.0 / totalEntries);
        if (hFSI[m]->Integral() > 0)  // only add non-empty histograms
            stack->Add(hFSI[m]);
    }
    hTotal->Scale(1.0/totalEntries);
    hTotal->GetYaxis()->SetRangeUser(0,hTotal->GetMaximum()*1.2);
    if (hTotal->GetMaximum()*1.2<0.25)   hTotal->GetYaxis()->SetRangeUser(0,0.25);
    hTotal->SetTitle("GENIE hA2018");
    if (model!="hA2018") hTotal->SetTitle("INCL++");
    hTotal->GetYaxis()->SetTitleOffset(1.3);
    hTotal->Draw("HIST");
    stack->Draw("SAME HIST ][");
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.025);
    c1->Update();  // commit before drawing
     hTotal->GetXaxis()->CenterTitle();
     hTotal->GetYaxis()->CenterTitle();

    // ---------------------------------------------------
    // Legend
    // ---------------------------------------------------
    TLegend *leg = new TLegend(0.20, 0.65, 0.45, 0.88);
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

    
leg->SetHeader(Form("%s on %s at %s GeV<T<%s GeV",
                    particleLabel.c_str(),
                    targetStr.c_str(),
                    keLow.c_str(),keHigh.c_str()));

leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.035);     // smaller, fits entries
    for (int m = 0; m < nModes; m++){
      if (abs(probe)==211) leg->AddEntry(hFSI[m], fsiLabels[m], "f");
    
      if ((abs(probe)==2212 || abs(probe)==2112)){ if (m==1 || m==4) continue; leg->AddEntry(hFSI[m], fsiLabels[m], "f");}
    }
    leg->Draw();

    // ---------------------------------------------------
    // Output
    // ---------------------------------------------------
    std::string channelStr;
    if (mode == kAbs)
        channelStr = useThreshold ? "stackThreshPionAbs" : "stackPionAbs";
    else
        channelStr = useThreshold ? "stackThreshKO"      : "stackKO";

    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%d_%s.pdf",
                   channelStr.c_str(),
                   probeStr.c_str(),
                   keStr.c_str(),
                   target, useThreshold, compound,model.c_str()));
}