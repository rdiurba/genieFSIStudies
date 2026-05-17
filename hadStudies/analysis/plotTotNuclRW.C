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

// ============================================================
// Nucleon: PDG 2212 (p+), -2212 (p-), 2112 (n)
// ============================================================
// Probe PDG: 2212

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
            {0.5746, 1.1866},  // hA
            {0.2722, 1.2417},  // hN
            {0.2678, 1.1802},  // INCL
            {-1.2307, 2.3133}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.3543, 1.8746},  // hA
            {-0.0353, 1.5883},  // hN
            {-0.0449, 1.2859},  // INCL
            {0.2465, 1.3957}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.7849, 0.7495, 0.0000, 1.0000},  // hA
            {-0.7622, 1.0998, 0.0000, 1.0000},  // hN
            {-5.5661, 4.5083, 0.0000, 1.0000},  // INCL
            {-5.8711, 4.7078, 0.0000, 1.0000}  // G4
        }
    };
    nucleonFitParams[1000180400] = {  // Ar40
        // DifMean (Difference Mean)
        {
            {0.0028, 0.4404},  // hA
            {-0.0195, 0.8275},  // hN
            {0.2749, 0.6338},  // INCL
            {-1.6409, 1.9140}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.6351, 1.9862},  // hA
            {0.4676, 1.7612},  // hN
            {0.4292, 1.3969},  // INCL
            {0.4111, 1.5907}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.9028, 0.7306, 0.0000, 1.0000},  // hA
            {-0.6253, 0.7244, 0.0000, 1.0000},  // hN
            {-5.6058, 4.4220, 0.0000, 1.0000},  // INCL
            {-7.7269, 5.6687, 0.0000, 1.0000}  // G4
        }
    };
    nucleonFitParams[1000080160] = {  // O16
        // DifMean (Difference Mean)
        {
            {0.5969, 1.1537},  // hA
            {0.2018, 1.2612},  // hN
            {0.3277, 1.1343},  // INCL
            {-1.3469, 2.3772}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.4093, 1.9114},  // hA
            {0.0217, 1.6573},  // hN
            {0.1314, 1.2863},  // INCL
            {0.1971, 1.4729}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.7918, 0.6834, 0.0000, 1.0000},  // hA
            {-0.8751, 1.0721, 0.0000, 1.0000},  // hN
            {-5.7255, 4.6972, 0.0000, 1.0000},  // INCL
            {-6.4376, 5.0006, 0.0000, 1.0000}  // G4
        }
    };
    nucleonFitParams[1000260560] = {  // Fe56
        // DifMean (Difference Mean)
        {
            {0.0474, 0.6469},  // hA
            {-0.0670, 0.9154},  // hN
            {0.2489, 0.8166},  // INCL
            {-1.7560, 2.1978}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.7882, 1.9893},  // hA
            {0.6754, 1.7866},  // hN
            {0.4306, 1.4387},  // INCL
            {0.6319, 1.5114}  // G4
        },
        // SumGamma (exponential: intercept*exp(KE^exp * slope) + floor)
        {
            {-0.9633, 0.7441, 0.0000, 1.0000},  // hA
            {-0.6125, 0.6577, 0.0000, 1.0000},  // hN
            {-5.0279, 4.0393, 0.0000, 1.0000},  // INCL
            {-6.5398, 4.9023, 0.0000, 1.0000}  // G4
        }
    };
}

// ============================================================
// Pion: PDG 211 (pi+), -211 (pi-), 111 (pi0)
// ============================================================
// Probe PDG: 211

struct PionTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsLinear SumMean;
    ModelParamsLinear SumStd;
};

std::map<int, PionTargetParams> pionFitParams;

void makePionParams() {
    pionFitParams[1000060120] = {  // C12
        // DifMean (Difference Mean)
        {
            {2.9002, 1.7511},  // hA
            {-0.1257, 2.3318},  // hN
            {0.1201, 2.2436},  // INCL
            {-0.2030, 2.1775}  // G4
        },
        // DifStd (Difference Std)
        {
            {-0.1002, 3.7254},  // hA
            {-0.0056, 1.3204},  // hN
            {-0.1042, 1.0109},  // INCL
            {0.1389, 0.9416}  // G4
        },
        // SumMean (Sum Mean)
        {
            {1.0699, 4.1030},  // hA
            {2.0520, 4.8404},  // hN
            {1.2387, 3.4527},  // INCL
            {0.7994, 1.9895}  // G4
        },
        // SumStd (Sum Std)
        {
            {0.6212, 0.8323},  // hA
            {0.3243, 3.2116},  // hN
            {0.9083, 1.3591},  // INCL
            {1.2666, 2.7622}  // G4
        },
    };
    pionFitParams[1000180400] = {  // Ar40
        // DifMean (Difference Mean)
        {
            {2.2903, 0.8905},  // hA
            {-0.6633, 1.4800},  // hN
            {-0.4752, 1.1072},  // INCL
            {-1.0096, 1.3981}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.2987, 3.2604},  // hA
            {0.3182, 2.2811},  // hN
            {0.7081, 1.2767},  // INCL
            {0.2957, 1.3652}  // G4
        },
        // SumMean (Sum Mean)
        {
            {3.7611, 4.9256},  // hA
            {6.6766, 7.1275},  // hN
            {4.9647, 1.8701},  // INCL
            {8.6358, 1.7998}  // G4
        },
        // SumStd (Sum Std)
        {
            {2.1374, 1.6616},  // hA
            {1.2266, 4.5166},  // hN
            {3.7494, 2.4962},  // INCL
            {1.4930, 3.0141}  // G4
        },
    };
    pionFitParams[1000080160] = {  // O16
        // DifMean (Difference Mean)
        {
            {2.8132, 1.9100},  // hA
            {-0.0463, 2.1889},  // hN
            {0.1818, 2.1437},  // INCL
            {-0.1161, 2.2395}  // G4
        },
        // DifStd (Difference Std)
        {
            {-0.1951, 3.9639},  // hA
            {0.0106, 1.5698},  // hN
            {0.0372, 1.0794},  // INCL
            {0.1149, 1.0644}  // G4
        },
        // SumMean (Sum Mean)
        {
            {1.4981, 4.2037},  // hA
            {3.2590, 5.0669},  // hN
            {1.8017, 2.9937},  // INCL
            {2.5676, 3.1655}  // G4
        },
        // SumStd (Sum Std)
        {
            {0.9073, 0.9469},  // hA
            {0.7763, 3.4013},  // hN
            {0.8655, 1.4225},  // INCL
            {1.0681, 2.2491}  // G4
        },
    };
    pionFitParams[1000260560] = {  // Fe56
        // DifMean (Difference Mean)
        {
            {1.8171, 0.9453},  // hA
            {-0.7655, 1.6732},  // hN
            {-0.6367, 1.9307},  // INCL
            {-0.6924, 2.2437}  // G4
        },
        // DifStd (Difference Std)
        {
            {0.3591, 3.3994},  // hA
            {0.6809, 2.3669},  // hN
            {0.7452, 1.2930},  // INCL
            {0.3858, 1.4334}  // G4
        },
        // SumMean (Sum Mean)
        {
            {5.0295, 5.7037},  // hA
            {9.1261, 6.7006},  // hN
            {5.6704, 3.6514},  // INCL
            {10.2044, 1.6007}  // G4
        },
        // SumStd (Sum Std)
        {
            {2.5278, 2.2597},  // hA
            {1.7786, 4.6312},  // hN
            {1.6104, 2.1473},  // INCL
            {2.0127, 2.6716}  // G4
        },
    };
}

// =======================================================
// Unified histogram builder
// =======================================================
TH1D* makeHistogram(TTree* t, std::string label, Mode mode, bool useThreshold,int A1, int target, bool compound, bool RW, std::string reweightModel="INCL++")
{
    int min=0; int max=40;
    if (A1<40){ min=0; max=A1+2;}
    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),
                          Form("hist_%s",label.c_str()),max,min,max);

    if (!t || t->GetEntries() < 1) return hist;

    int pdgf[100], nf, probe_fsi;
    double Ef[100], mh[100];
    int np; int nn; double KEini;
    t->SetBranchAddress("pdgh",&pdgf);
    t->SetBranchAddress("ke",&KEini);
    t->SetBranchAddress("nh",&nf);
    t->SetBranchAddress("Eh",&Ef);
    t->SetBranchAddress("mh",&mh);
    t->SetBranchAddress("np",&np);
    t->SetBranchAddress("nn",&nn);
    t->SetBranchAddress("probe_fsi",&probe_fsi);
    
    int required_fsi = (mode == kAbs) ? 5 : 6;
    double thresh = useThreshold ? 0.005 : 0.0;

    for (Long64_t j=0; j<t->GetEntries(); j++){
        t->GetEntry(j);
        if (probe_fsi!=5 && probe_fsi!=6) continue;
       int mult=0;
       int piMult=0;
       int diff=0;
       for(int i=0; i<nf; i++){
            int apdg = abs(pdgf[i]);
            double KE = Ef[i] - mh[i];
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

        if (useThreshold==0 && compound==false) mult=np+nn;
if (RW) {
    double hAEst=1;
    double AltEst=1;
    if (mode==kAbs){
    const auto& tp = pionFitParams.at(target);

        const LinearParams& altSumMean = (reweightModel=="hN2018") ? tp.SumMean.hN   :
                                         (reweightModel=="Geant4")  ? tp.SumMean.G4   :
                                                                       tp.SumMean.INCL;
        const LinearParams& altSumStd  = (reweightModel=="hN2018") ? tp.SumStd.hN    :
                                         (reweightModel=="Geant4")  ? tp.SumStd.G4    :
                                                                       tp.SumStd.INCL;

        double meanAltSum = altSumMean.slope * KEini + altSumMean.intercept;
        double stdAltSum  = altSumStd.slope  * KEini + altSumStd.intercept;
        double meanhASum  = tp.SumMean.hA.slope * KEini + tp.SumMean.hA.intercept;
        double stdhASum   = tp.SumStd.hA.slope  * KEini + tp.SumStd.hA.intercept;



        const LinearParams& altDifMean = (reweightModel=="hN2018") ? tp.DifMean.hN   :
                                         (reweightModel=="Geant4")  ? tp.DifMean.G4   :
                                                                       tp.DifMean.INCL;
        const LinearParams& altDifStd  = (reweightModel=="hN2018") ? tp.DifStd.hN    :
                                         (reweightModel=="Geant4")  ? tp.DifStd.G4    :
                                                                       tp.DifStd.INCL;

        double meanAltDif = altDifMean.slope * KEini + altDifMean.intercept;
        double stdAltDif  = altDifStd.slope  * KEini + altDifStd.intercept;
        double meanhADif  = tp.DifMean.hA.slope * KEini + tp.DifMean.hA.intercept;
        double stdhADif   = tp.DifStd.hA.slope  * KEini + tp.DifStd.hA.intercept;

        double hAEstSum  = TMath::Gaus(mult, meanhASum,  stdhASum,  1);
        double AltEstSum = TMath::Gaus(mult, meanAltSum, stdAltSum, 1);
        double hAEstDif  = TMath::Gaus(diff, meanhADif,  stdhADif,  1);
        double AltEstDif = TMath::Gaus(diff, meanAltDif, stdAltDif, 1);
        hAEst=hAEstSum*hAEstDif; AltEst=AltEstSum*AltEstDif;
    
    }
    else{


        const auto& tp = nucleonFitParams.at(target);

        const GammaParams& altGamma = (reweightModel=="hN2018") ? tp.SumGamma.hN   :
                                      (reweightModel=="Geant4")  ? tp.SumGamma.G4   :
                                                                   tp.SumGamma.INCL;

        double gammaAlt = altGamma.intercept * TMath::Exp(TMath::Power(KEini, altGamma.exp) * altGamma.slope) + altGamma.floor;
        double gammaHA  = tp.SumGamma.hA.intercept * TMath::Exp(TMath::Power(KEini, tp.SumGamma.hA.exp) * tp.SumGamma.hA.slope) + tp.SumGamma.hA.floor;
        hAEst  = TMath::Exp(-(mult - 3) * gammaHA);
        AltEst = TMath::Exp(-(mult - 3) * gammaAlt);
    }
    double rw     = (hAEst > 0) ? AltEst / hAEst : 1.0;
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
void plotTotNuclRW(int probe=211, double ke=0.175, int target=1000180400,
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
    int compound=0;
    TH1D* h1 = makeHistogram(t1,"20i",mode,useThreshold,A,target,compound,1,reweightModel);
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
    std::string channelStr="nucleonSumThreshPionAbs";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "nucleonSumThreshPionAbs" :
            "nucleonSumPionAbs";
    if (mode==kKO)
        channelStr = useThreshold ?
            "nucleonSumThreshKO" :
            "nucleonSumKO";
    // Argon Sum parameters




    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%dRW%s.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound,reweightModel.c_str()));

    // ---------------------------------------------------
    // Fit
    // ---------------------------------------------------
    double edge1 = -1.0;
    double edge2 = -1.0;
    double edge3 = -1.0;
    double edge4 = -1.0;
    
    TH1* hists[] = {h1, h2, h3, h4};
    double* edges[] = {&edge1, &edge2, &edge3, &edge4};


   
    
    
    TF1* g1 = (mode==kAbs) ? new TF1("g1","gaus",3,max) : new TF1("g1","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g2 = (mode==kAbs) ? new TF1("g2","gaus",3,max) : new TF1("g2","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g3 = (mode==kAbs) ? new TF1("g3","gaus",3,max) : new TF1("g3","[0]*exp(-(x-3)*[1])",3,max);
    TF1* g4 = (mode==kAbs) ? new TF1("Geant4","gaus",3,max) : new TF1("Geant4","[0]*exp(-(x-3)*[1])",3,max);

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
    h1->Fit(g1,"LRQ");
    h2->Fit(g2,"LRQ");
    h3->Fit(g3,"LRQ");
    h4->Fit(g4,"LRQ");

}