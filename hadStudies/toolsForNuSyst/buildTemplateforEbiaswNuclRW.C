#include "TH1D.h"
#include <TCanvas.h>
#include <THStack.h>
#include <TLatex.h>
#include <TMath.h>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"
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

void normalizePerEnergy(TH2D* h)
{
    int nX = h->GetNbinsX();
    int nY = h->GetNbinsY();

    // include X underflow (0) and overflow (nX+1)
    for (int ix = 0; ix <= nX + 1; ++ix) {

        double sum = 0.0;

        // include Y underflow (0) and overflow (nY+1)
        for (int iy = 0; iy <= nY + 1; ++iy)
            sum += h->GetBinContent(ix, iy);

        if (sum == 0.0) continue;

        for (int iy = 0; iy <= nY + 1; ++iy) {
            int bin = h->GetBin(ix, iy);

            double val = h->GetBinContent(bin);
            double err = h->GetBinError(bin);

            h->SetBinContent(bin, val / sum);
            h->SetBinError(bin, err / sum);
        }
    }
}
TH2* getHist(std::string reweightModel, int pdg, std::string particle) { // get histogram from the TTree (e.g. gst/ginuke file)

  const int nummax = 99;


    int eBiasBins=25;
    Float_t vEdges[26];
    vEdges[0]=-1;
      for (int i = 1; i < 26; ++i){
        vEdges[i]=0+0.05*(i-1);
}
        int bins=40;
     Float_t binEdges[41];
      std::vector<double> v;
      for (int i = 0; i < 40; ++i){
        //v.push_back( std::round(0.02 * std::pow(10., i/23.) * 1e4) / 1e4 );
        v.push_back( std::round((0.025 +0.05*i)*1e4) / 1e4 );

}
    for (int i = 1; i < 40; ++i)
    binEdges[i] = 0.5 * (v[i-1] + v[i]);


    binEdges[0]  = v[0]  - (binEdges[1]  - v[0]);
    binEdges[40] = v[39] + (v[39] - binEdges[39]);
    //TH2D* h_KEini_Ebias = new TH2D(Form("%s_%s",simulation.c_str(),particle.c_str()), "Ebias vs KEini; KEini [GeV]; Relative Ebias", bins, binEdges, eBiasBins, vEdges);
    TH2D* h_KEini_Ebias = new TH2D(Form("%s_%s_KEini_vs_Ebias", reweightModel.c_str(), particle.c_str()), Form("%s_%s_KEini_vs_Ebias", reweightModel.c_str(), particle.c_str()), 40, 0, 2, 20, -1, 1);

    TH2D* hist;
    for (int i=1; i<40; i++){
    TFile* file = new TFile(Form("$PSCRATCH/hadronFiles/%s_1000180400_%gGeV_hA2018_100k.ginuke.root",particle.c_str(),v.at(i) ));



  
    TTree* tree = (TTree*)file->Get("ginuke");
    // get tree branches and define histograms
    double Eini;
    tree->SetBranchAddress("e", &Eini);
    double Pini;
    tree->SetBranchAddress("p", &Pini);
    double KEini;
    tree->SetBranchAddress("ke", &KEini);
    int probe_fsi; // 1==nofsi; 2==cex; 3==elas; 4==inel; 5=abs; ...
    tree->SetBranchAddress("probe_fsi", &probe_fsi);
    
    int nh; // number of FS hadrons
    tree->SetBranchAddress("nh", &nh);
    double Eh[nummax];
    tree->SetBranchAddress("Eh", Eh);
    double mh[nummax];
    tree->SetBranchAddress("mh", mh);
    int pdgh[nummax]; // {pi+, pi-, pi0, n, p, gamma, K0, K+, Lambda, ...} == {211, -211, 111, 2112, 2212, 22, 311, 321, 3122, ...}
    tree->SetBranchAddress("pdgh", pdgh);
    
    int np;
    tree->SetBranchAddress("np", &np);
    TH1I* h_np = new TH1I("h_np", "Number of FS protons", 20, 0, 20);
    int nn;
    tree->SetBranchAddress("nn", &nn);
    TH1I* h_nn = new TH1I("h_nn", "Number of FS neutrons", 20, 0, 20);
    int npip;
    tree->SetBranchAddress("npip", &npip);
    TH1I* h_npip = new TH1I("h_npip", "Number of FS #pi^{+}", 5, 0, 5);
    int npim;
    tree->SetBranchAddress("npim", &npim);
    TH1I* h_npim = new TH1I("h_npim", "Number of FS #pi^{-}", 5, 0, 5);
    int npi0;
    tree->SetBranchAddress("npi0", &npi0);
    int probe=pdg;
    // fill histograms
        double thresh=0.005;
        bool compound=0;
    for (int i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
       int    mult   = 0;
        int    piMult = 0;
        int diff =0;

      if (Eh[0]==Eini) continue;
        for (int i = 0; i < nh; i++) {           // fixed: nf/pdgf/Ef -> nh/pdgh/Eh
            int    apdg = abs(pdgh[i]);
            double KE   = Eh[i] - mh[i];

            if (apdg == 211) {

                piMult++;
            }
            if (apdg == 22) {
                //if (KE <= 0.5 * thresh) continue;
                  // fixed: Ef -> Eh

            }
            
            if (apdg == 111) {
                 // fixed: Ef -> Eh
                            piMult++;
            }
            if (KE <= thresh) continue;

            if (apdg == 2212 || apdg == 2112) {
                mult++;
                if (apdg==2212) diff++;
                else diff--;
            }

            if (apdg > 1000000000) {
                int Z = (apdg / 10000) % 1000;
                int A = (apdg / 10)    % 1000;
                int N = A - Z;

                if (Z > 2)          continue;
                if (KE <= thresh*A) continue;


                    mult += Z + N;
                    diff += Z-N;

            }
        }

        int fsi = -1;
        
        if (Eini==Eh[0]) {
            fsi = 1;
        }
        else if ((abs(probe) == 211 || abs(probe) == 111) &&
                 (piMult == 0)) {
            fsi = 5; // Absorption
        }
        else if ((abs(probe) == 2212 || abs(probe) == 2112) &&
                 (piMult == 0) && mult > 2) {  // fixed: missing opening {
            fsi = 6; // Knock-out
        }
        else if ((abs(probe)==211 && piMult >1) || (abs(probe)==2212 && piMult>0)) {
            fsi = 7;
        }

        else if (mult + piMult >= 2) {            
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
            if (((mult== 1 && probe>2100) || ((piMult==1) && probe<212)) && ( fsi != 1))
                fsi = 3;                          // fixed: fsi==3 -> fsi=3
        }


      double maxKE = -999;
      int maxidx = -1;
      double minKE = 999;
      int minidx = -1;
      double Ehad = 0;
      int ngamma = 0;
      int nother = 0;
      for (int idx = 0; idx < nh; ++idx) {
        double tmpKE = Eh[idx] - mh[idx];
        if (tmpKE > maxKE) {
          maxKE = tmpKE;
          maxidx = idx;
        }
        if (tmpKE < minKE) {
          minKE = tmpKE;
          minidx = idx;
        }
        
        if (pdgh[idx]==211||pdgh[idx]==-211||pdgh[idx]==111)
          Ehad += Eh[idx]; // E_pi
          //Ehad += tmpKE; // KE_pi
        else if (pdgh[idx]==2212){
         if (tmpKE>thresh) Ehad += tmpKE;} //KE_p}
        else if (pdgh[idx]==2112)
          ;
        else if (pdgh[idx]==22) {
          ++ngamma;
          Ehad += Eh[idx];
        }
        else if (pdgh[idx] > 1000000000) {
                int Z = (pdgh[idx] / 10000) % 1000;
                int A = (pdgh[idx] / 10)    % 1000;
                int N = A - Z;

                if (Z > 2)          continue;
                if (tmpKE <= thresh*A) continue;


                 if (compound==true) Ehad += tmpKE;
            }
      }
      
      if ( (fsi==5 || fsi==6) && fsi!=7 && KEini>0&&KEini<20) { // selections (e.g. exclusive channel & KE range)

    double hAEst=0;
    double AltEst=0;
    if (pdg<300){
    const auto& tp = pionFitParams.at(1000180400);

        const LinearParams& altSumMean = (reweightModel=="hN2018") ? tp.SumMean.hN   :
                                         (reweightModel=="HG4BertCasc")  ? tp.SumMean.G4   :
                                                                       tp.SumMean.INCL;
        const LinearParams& altSumStd  = (reweightModel=="hN2018") ? tp.SumStd.hN    :
                                         (reweightModel=="HG4BertCasc")  ? tp.SumStd.G4    :
                                                                       tp.SumStd.INCL;

        double meanAltSum = altSumMean.slope * KEini + altSumMean.intercept;
        double stdAltSum  = altSumStd.slope  * KEini + altSumStd.intercept;
        double meanhASum  = tp.SumMean.hA.slope * KEini + tp.SumMean.hA.intercept;
        double stdhASum   = tp.SumStd.hA.slope  * KEini + tp.SumStd.hA.intercept;



        const LinearParams& altDifMean = (reweightModel=="hN2018") ? tp.DifMean.hN   :
                                         (reweightModel=="HG4BertCasc")  ? tp.DifMean.G4   :
                                                                       tp.DifMean.INCL;
        const LinearParams& altDifStd  = (reweightModel=="hN2018") ? tp.DifStd.hN    :
                                         (reweightModel=="HG4BertCasc")  ? tp.DifStd.G4    :
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


        const auto& tp = nucleonFitParams.at(1000180400);

        const GammaParams& altGamma = (reweightModel=="hN2018") ? tp.SumGamma.hN   :
                                      (reweightModel=="HG4BertCasc")  ? tp.SumGamma.G4   :
                                                                   tp.SumGamma.INCL;

            const LinearParams& altDifMean = (reweightModel=="hN2018") ? tp.DifMean.hN   :
                                         (reweightModel=="HG4BertCasc")  ? tp.DifMean.G4   :
                                                                       tp.DifMean.INCL;
        const LinearParams& altDifStd  = (reweightModel=="hN2018") ? tp.DifStd.hN    :
                                         (reweightModel=="HG4BertCasc")  ? tp.DifStd.G4    :
                                                                       tp.DifStd.INCL;

        double gammaAlt = altGamma.intercept * TMath::Exp(TMath::Power(KEini, altGamma.exp) * altGamma.slope) + altGamma.floor;
        double gammaHA  = tp.SumGamma.hA.intercept * TMath::Exp(TMath::Power(KEini, tp.SumGamma.hA.exp) * tp.SumGamma.hA.slope) + tp.SumGamma.hA.floor;
        double hAEstSum  = TMath::Exp(-(mult - 3) * gammaHA);
        double AltEstSum = TMath::Exp(-(mult - 3) * gammaAlt);



        double meanAltDif = altDifMean.slope * KEini + altDifMean.intercept;
        double stdAltDif  = altDifStd.slope  * KEini + altDifStd.intercept;
        double meanhADif  = tp.DifMean.hA.slope * KEini + tp.DifMean.hA.intercept;
        double stdhADif   = tp.DifStd.hA.slope  * KEini + tp.DifStd.hA.intercept;

        double hAEstDif  = TMath::Gaus(diff, meanhADif,  stdhADif,  1);
        double AltEstDif = TMath::Gaus(diff, meanAltDif, stdAltDif, 1);
        hAEst=hAEstSum*hAEstDif; AltEst=AltEstSum*AltEstDif;
    }
    double rw     = (hAEst > 0) ? AltEst / hAEst : 1.0;
    if(rw<0.01) rw=0.01; if (rw>100) rw=100;

        
        double Ebias;
        if (pdg==2212 || pdg==2112)
          Ebias = (KEini - Ehad) / KEini;
        else if (pdg==211 || pdg==111 || pdg==-211)
          Ebias = (Eini - Ehad) / Eini;
        //h_Ebias->Fill(Ebias);
        h_KEini_Ebias->Fill(KEini, Ebias,rw);
      }
    }

    }
    
        normalizePerEnergy(h_KEini_Ebias);
  return h_KEini_Ebias;
}

void buildTemplate(int pdg=211) {

  TH3 *hA2018_hist, *hN2018_hist, *INCL_hist, *G4BC_hist;
  string particle;

    if (pdg == 2212) particle = "protonPlus";
    else if (pdg == 2112) particle = "neutron";
    else if (pdg == 211) particle = "piPlus";
    else if (pdg == 111) particle = "pi0";
    else if (pdg == -211) particle = "piMinus";
    hN2018_hist = (TH3*)getHist(Form("hN2018"), pdg, particle);
    INCL_hist = (TH3*)getHist(Form("HINCL"), pdg, particle);
    G4BC_hist = (TH3*)getHist(Form("HG4BertCasc"), pdg, particle);
  
  
  
  // column normalization
  /*for (int i=1; i<=hA2018_hist->GetNbinsX(); ++i) { // Eini
    double sum1 = hA2018_hist->ProjectionX()->GetBinContent(i);
    double sum2 = hN2018_hist->ProjectionX()->GetBinContent(i);
    double sum3 = INCL_hist->ProjectionX()->GetBinContent(i);
    double sum4 = G4BC_hist->ProjectionX()->GetBinContent(i);
    for (int j=1; j<=hA2018_hist->GetNbinsY(); ++j) { // Ebias
      double tmp = hA2018_hist->GetBinContent(i,j)/sum1;
      hA2018_hist->SetBinContent(i,j, tmp);
      
      tmp = hN2018_hist->GetBinContent(i,j)/sum2;
      hN2018_hist->SetBinContent(i,j, tmp);
      
      tmp = INCL_hist->GetBinContent(i,j)/sum3;
      INCL_hist->SetBinContent(i,j, tmp);
      
      tmp = G4BC_hist->GetBinContent(i,j)/sum4;
      G4BC_hist->SetBinContent(i,j, tmp);
    }
  }
  */
    gStyle->SetOptStat(0);

  TCanvas* cc1 = new TCanvas();


  hN2018_hist->GetZaxis()->SetRangeUser(0,0.5);
  hN2018_hist->Draw("colz");

  cc1->Print(Form("RW2hN2018_%s_Ebias.png",particle.c_str()) );
  INCL_hist->GetZaxis()->SetRangeUser(0,0.5);

  INCL_hist->Draw("colz");
  cc1->Print(Form("RW2INCL_%s_Ebias.png",particle.c_str()) );
  G4BC_hist->GetZaxis()->SetRangeUser(0,0.5);

  G4BC_hist->Draw("colz");
  cc1->Print(Form("RW2G4BC_%s_Ebias.png",particle.c_str()) );
    if (pdg == 2212) particle = "proton";
    else if (pdg == 2112) particle = "neutron";
    else if (pdg == 211) particle = "pionp";
    else if (pdg == 111) particle = "pion0";
    else if (pdg == -211) particle = "pionm";
  TFile *outfile = new TFile("FSI_KOAbs_hA2018rwNucl_templates.root","UPDATE");
  hN2018_hist->Write(Form("rwhN2018_%s_KEini_vs_Ebias", particle.c_str()));
  INCL_hist->Write(Form("rwINCL++_%s_KEini_vs_Ebias", particle.c_str()));
  G4BC_hist->Write(Form("rwGeant4_%s_KEini_vs_Ebias", particle.c_str()));

  outfile->Write();  outfile->Close(); }
void buildTemplateforEbiaswNuclRW(){
    makePionParams();
    makeNucleonParams();
buildTemplate(211);

buildTemplate(2212);



}