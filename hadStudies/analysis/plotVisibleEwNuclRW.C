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

struct LinearParams {
    double slope, intercept;
};

struct GammaParams {
    double slope, intercept, floor, exp;
};

struct ModelParamsLinear {
    LinearParams hA, hN, INCL, G4;
};

struct ModelParamsGamma {
    GammaParams hA, hN, INCL, G4;
};

struct PionTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsLinear SumMean;
    ModelParamsLinear SumStd;
};



std::map<int, PionTargetParams>    pionFitParams;


struct TargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsGamma  SumGamma;
};

std::map<int, TargetParams> nucleonFitParams;

void makeNucleonParams() {
    nucleonFitParams[1000060120] = {  // C12
        // DifMean (Difference Mean)
        {
            {0.0397, 1.5895},  // hA
            {0.0067, 1.4666},  // hN
            {0.1033, 1.3434},  // INCL
            {-0.0921, 1.6805}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.0256, 2.0761},  // hA
            {-0.0241, 1.5651},  // hN
            {-0.0673, 0.9921},  // INCL
            {-0.0403, 1.1312}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.0000, 0.2097, 0.0031, 4.9989},  // hA
            {-1.5993, 1.7594, 0.0000, 0.1032},  // hN
            {-3.8267, 2.6684, 0.4702, 0.7241},  // INCL
            {-4.4605, 0.9880, 0.3776, 1.6196}  // G4
        }
    };
    nucleonFitParams[1000180400] = {  // Ar40
        // DifMean (Difference Mean)
        {
            {-0.0086, 0.4249},  // hA
            {-0.1964, 0.9925},  // hN
            {-0.0992, 0.3686},  // INCL
            {-0.5177, 0.8338}  // G4
        },
        // DifStd (Difference Std)
        {
            {-0.0230, 2.4054},  // hA
            {0.2043, 2.0487},  // hN
            {0.2596, 1.3365},  // INCL
            {0.1272, 1.5463}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.0447, 0.1173, 0.0070, 3.2000},  // hA
            {-4.6356, 0.1916, 0.1374, 2.7061},  // hN
            {-3.9135, 2.1172, 0.2181, 0.7484},  // INCL
            {-3.3876, 1.9572, 0.0907, 0.7590}  // G4
        }
    };
    nucleonFitParams[1000080160] = {  // O16
        // DifMean (Difference Mean)
        {
            {0.0187, 1.5989},  // hA
            {-0.0029, 1.4802},  // hN
            {0.1176, 1.3384},  // INCL
            {-0.1389, 1.8132}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.0011, 2.1667},  // hA
            {-0.0089, 1.6803},  // hN
            {-0.0349, 1.0824},  // INCL
            {-0.0136, 1.1882}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.0000, 0.1696, 0.0115, 0.9472},  // hA
            {-6.7320, 0.1636, 0.2784, 3.4920},  // hN
            {-3.6009, 2.4897, 0.4714, 0.8123},  // INCL
            {-4.2426, 1.0249, 0.3022, 1.6248}  // G4
        }
    };
    nucleonFitParams[1000260560] = {  // Fe56
        // DifMean (Difference Mean)
        {
            {0.0005, 0.7074},  // hA
            {-0.2305, 1.1033},  // hN
            {-0.1406, 1.1194},  // INCL
            {-0.5439, 1.8611}  // G4
        },
        // DifStd (Difference Std)
        {
            {-0.0029, 2.5204},  // hA
            {0.3196, 2.1606},  // hN
            {0.3096, 1.3404},  // INCL
            {0.2359, 1.4816}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.3339, 0.0528, 0.0595, 2.9499},  // hA
            {-4.7378, 0.1965, 0.1073, 2.4599},  // hN
            {-4.0297, 2.8468, 0.1586, 0.5933},  // INCL
            {-3.3271, 3.2049, 0.0000, 0.4905}  // G4
        }
    };
}
void makePionParams() {
    pionFitParams[1000060120] = {  // C12
        // DifMean (Difference Mean)
        {
            {3.2148, 1.4954},  // hA
            {-0.1639, 2.3667},  // hN
            {0.1236, 2.2404},  // INCL
            {-0.2144, 2.1879}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.7429, 2.9893},  // hA
            {0.0169, 1.3000},  // hN
            {-0.0912, 0.9988},  // INCL
            {0.1210, 0.9577}  // G4
        },
        // SumMean (Sum Mean)
        {
            {1.1616, 4.0239},  // hA
            {2.3657, 4.5485},  // hN
            {1.4105, 3.3038},  // INCL
            {0.1020, 2.6253}  // G4
        },
        // SumStd (Sum Std)
        {
            {0.5463, 0.8969},  // hA
            {0.3593, 3.1787},  // hN
            {0.9076, 1.3596},  // INCL
            {1.5146, 2.5651}  // G4
        },
    };
    pionFitParams[1000180400] = {  // Ar40
        // DifMean (Difference Mean)
        {
            {2.1697, 0.9996},  // hA
            {-0.9212, 1.7142},  // hN
            {-0.5002, 1.1292},  // INCL
            {-1.0078, 1.3965}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.5534, 3.0303},  // hA
            {0.4847, 2.1299},  // hN
            {0.7082, 1.2766},  // INCL
            {0.2975, 1.3637}  // G4
        },
        // SumMean (Sum Mean)
        {
            {3.9453, 4.7679},  // hA
            {7.9786, 6.1138},  // hN
            {4.9160, 1.9157},  // INCL
            {8.6519, 1.7847}  // G4
        },
        // SumStd (Sum Std)
        {
            {1.6313, 2.0764},  // hA
            {0.9784, 4.7197},  // hN
            {4.0375, 2.2754},  // INCL
            {1.8597, 2.6864}  // G4
        },
    };
    pionFitParams[1000080160] = {  // O16
        // DifMean (Difference Mean)
        {
            {3.4857, 1.3851},  // hA
            {-0.1127, 2.2500},  // hN
            {0.1788, 2.1465},  // INCL
            {-0.1397, 2.2609}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.8390, 3.1123},  // hA
            {0.1195, 1.4698},  // hN
            {0.0448, 1.0726},  // INCL
            {0.0965, 1.0811}  // G4
        },
        // SumMean (Sum Mean)
        {
            {1.5891, 4.1252},  // hA
            {3.7442, 4.6430},  // hN
            {1.8433, 2.9569},  // INCL
            {2.2825, 3.4394}  // G4
        },
        // SumStd (Sum Std)
        {
            {0.8294, 1.0126},  // hA
            {0.6607, 3.4986},  // hN
            {0.8183, 1.4634},  // INCL
            {1.5140, 1.8598}  // G4
        },
    };
    pionFitParams[1000260560] = {  // Fe56
        // DifMean (Difference Mean)
        {
            {1.7247, 1.0279},  // hA
            {-0.9912, 1.8750},  // hN
            {-0.6421, 1.9355},  // INCL
            {-0.6761, 2.2294}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.6204, 3.1639},  // hA
            {0.8871, 2.1829},  // hN
            {0.7532, 1.2858},  // INCL
            {0.4003, 1.4204}  // G4
        },
        // SumMean (Sum Mean)
        {
            {5.4464, 5.3506},  // hA
            {10.1712, 5.8249},  // hN
            {5.7534, 3.5781},  // INCL
            {9.8912, 1.8829}  // G4
        },
        // SumStd (Sum Std)
        {
            {1.7712, 2.8499},  // hA
            {1.4227, 4.9020},  // hN
            {1.7790, 2.0019},  // INCL
            {2.3325, 2.3918}  // G4
        },
    };
}

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

    if (weight_nom == 0.) return 1.;


    return weight_alt / weight_nom;
}
TH1D* makeHistogram(TTree* t, std::string label, Mode mode, bool useThreshold,int A1, int target, bool compound, bool RW, std::string pstr="pionp",std::string reweightModel="INCL")
{

    int min=0; int max=40;
    if (A1<40){ min=0; max=A1+2;}
    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),
                          Form("hist_%s",label.c_str()),12, -1, 0.2);
    if (!t || t->GetEntries() < 1) return hist;

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

    
    int required_fsi = (mode == kAbs) ? 5 : 6;
    double thresh = useThreshold ? 0.005 : 0.0;

    for (Long64_t j=0; j<t->GetEntries(); j++){
        t->GetEntry(j);
       int mult=0;
       int piMult=0;
       int diff=0;
    double visE=0;

       for(int i=0; i<nh; i++){
            int apdg = abs(pdgh[i]);
            double KE = Eh[i] - mh[i];
            if (apdg==211|| apdg==111){
                piMult++;
            }
            
            if (KE <= thresh) continue;

            if (apdg==2212 || apdg==2112){
                mult++;
                if (apdg==2212) diff++;
                else diff--;
            }
           
            if (apdg > 1000000000){
                int Z = (apdg/10000) % 1000;
                int A = (apdg/10) % 1000;
                int N = A - Z;

                   if (Z > 2) continue;
             if (KE<=thresh*A) continue;

             if (compound==true)    mult += Z + N;
            if (compound==true) diff+=Z-N;
            }
            
        }
        if (required_fsi==6 && (mult<3 || piMult!=0)) continue;
        if (required_fsi==5 && (piMult!=0 || mult<2)) continue;


        for (int i = 0; i < nh; i++) {           // fixed: nf/pdgf/Ef -> nh/pdgh/Eh
            int    apdg = abs(pdgh[i]);
            double KE   = Eh[i] - mh[i];

            if (apdg == 211) {
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

        double Ehad = visE;                       // fixed: removed re-declaration of visE, named clearly
        double Ebias = 0.0;
        if (mode==kKO)
            Ebias = (KEini - Ehad) / KEini;
        else 
            Ebias = (Eini - Ehad) / Eini;
        Ebias=-Ebias;
    double visERW=1;
    if (RW) {
    TFile* f2 = TFile::Open("FSI_KOAbs_Evis_reweight_template.root", "READ");
        TH2D* hist_nom = (TH2D*)f2->Get(Form("hA2018_%s_KEini_vs_Ebias", pstr.c_str()));
    TH2D* hist_alt = (TH2D*)f2->Get(Form("%s_%s_KEini_vs_Ebias",     reweightModel.c_str(), pstr.c_str()));
    double visERW=GetVisEReweight(hist_nom, hist_alt, KEini, -Ebias);

    double rw     = visERW;
    if(rw<0.01) rw=0.01; if (rw>100) rw=100;
            hist->Fill(mult, rw);
        } else {
            hist->Fill(mult, 1);
        }
    }

    return hist;
}
// =======================================================
// Main plotting function
// =======================================================
void plotVisibleEwNuclRW(int probe=211, double ke=0.175, int target=1000180400,
                  bool useThreshold = true, std::string reweightModel="hN2018"){
    makePionParams();
    makeNucleonParams();
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
    std::string particle;
    if (probe == 2212) particle = "protonPlus";
    else if (probe == 2112) particle = "neutron";
    else if (probe == 211) particle = "piPlus";
    else if (probe == 111) particle = "pi0";
    else if (probe == -211) particle = "piMinus";

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
    int compound=1;
    TH1D* h1 = makeHistogram(t1,"20i",mode,useThreshold,A,target,compound,1,particle,reweightModel);
    TH1D* h2 = makeHistogram(t2,"20j",mode,useThreshold,A,target,compound,0);
    TH1D* h3 = makeHistogram(t3,"20k",mode,useThreshold,A,target,compound,0);
    TH1D* h4 = makeHistogram(t4,"20l",mode,useThreshold,A,target,compound,0);

h1->Scale(1.0 / h1->Integral());
h2->Scale(1.0 / h2->Integral());
h3->Scale(1.0 / h3->Integral());
h4->Scale(1.0 / h4->Integral());

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


    h1->GetXaxis()->SetTitle("(E^{vis}-E^{ini})/E^{ini}");
    if (probe>211) h1->GetXaxis()->SetTitle("(E^{vis}-T^{ini})/T^{ini}");

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
            "visESumThreshPionAb" :
            "visESumPionAb";
    if (mode==kKO)
        channelStr = useThreshold ?
            "visESumThreshKO" :
            "visESumKO";
    // Argon Sum parameters




    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%dRW%s.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound,reweightModel.c_str()));



}