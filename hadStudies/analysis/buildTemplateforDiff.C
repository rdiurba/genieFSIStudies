#include "TH1D.h"
#include <TCanvas.h>
#include <THStack.h>
#include <TLatex.h>
#include <TMath.h>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"
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
TH2* getHist(std::string simulation, int pdg, std::string particle) { // get histogram from the TTree (e.g. gst/ginuke file)

  const int nummax = 99;


    int NuclDiffBins=25;
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
    //TH2D* h_KEini_NuclDiff = new TH2D(Form("%s_%s",simulation.c_str(),particle.c_str()), "NuclDiff vs KEini; KEini [GeV]; Relative NuclDiff", bins, binEdges, NuclDiffBins, vEdges);
    TH2D* h_KEini_NuclDiff = new TH2D(Form("%s_%s_NuclDiff", simulation.c_str(), particle.c_str()), Form("%s_%s_NuclDiff", simulation.c_str(), particle.c_str()), 40, 0, 2, 40,-20,20);

    TH2D* hist;
    for (int i=1; i<40; i++){
    TFile* file = new TFile(Form("$PSCRATCH/hadronFiles/%s_1000180400_%gGeV_%s_100k.ginuke.root",particle.c_str(),v.at(i),simulation.c_str() ));



  
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
        bool compound=0;
    for (int i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
       int    mult   = 0; int diff=0;
        int    piMult = 0;
        double thresh=0.005;
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


      
      if ( (fsi==5 || fsi==6) && fsi!=7 && KEini>0&&KEini<20) { // selections (e.g. exclusive channel & KE range)



    
        //h_NuclDiff->Fill(NuclDiff);
        h_KEini_NuclDiff->Fill(KEini, diff);
      }
    }

    }
    
        normalizePerEnergy(h_KEini_NuclDiff);
  return h_KEini_NuclDiff;
}

void buildTemplate(int pdg=211) {

  TH3 *hA2018_hist, *hN2018_hist, *INCL_hist, *G4BC_hist;
  string particle;


    if (pdg == 2212) particle = "protonPlus";
    else if (pdg == 2112) particle = "neutron";
    else if (pdg == 211) particle = "piPlus";
    else if (pdg == 111) particle = "pi0";
    else if (pdg == -211) particle = "piMinus";
    hA2018_hist = (TH3*)getHist(Form("hA2018"), pdg, particle);
    hN2018_hist = (TH3*)getHist(Form("hN2018"), pdg, particle);
    INCL_hist = (TH3*)getHist(Form("HINCL"), pdg, particle);
    G4BC_hist = (TH3*)getHist(Form("HG4BertCasc"), pdg, particle);
  
  
    gStyle->SetOptStat(0);

  TCanvas* cc1 = new TCanvas();
  hA2018_hist->GetZaxis()->SetRangeUser(0,0.5);
  hA2018_hist->Draw("colz");
  cc1->Print(Form("hA2018_%s_NuclDiff.png",particle.c_str()) );

  hN2018_hist->GetZaxis()->SetRangeUser(0,0.5);
  hN2018_hist->Draw("colz");

  cc1->Print(Form("hN2018_%s_NuclDiff.png",particle.c_str()) );
  INCL_hist->GetZaxis()->SetRangeUser(0,0.5);

  INCL_hist->Draw("colz");
  cc1->Print(Form("INCL_%s_NuclDiff.png",particle.c_str()) );
  G4BC_hist->GetZaxis()->SetRangeUser(0,0.5);

  G4BC_hist->Draw("colz");
  cc1->Print(Form("G4BC_%s_NuclDiff.png",particle.c_str()) );

  TFile *outfile = new TFile("FSI_KOAbs_Mult_reweight_template.root","UPDATE");
  hA2018_hist->Write(Form("hA2018_%s_NuclDiff", particle.c_str()));
  hN2018_hist->Write(Form("hN2018_%s_NuclDiff", particle.c_str()));
  INCL_hist->Write(Form("INCL++_%s_NuclDiff", particle.c_str()));
  G4BC_hist->Write(Form("Geant4_%s_NuclDiff", particle.c_str()));

  outfile->Write();  outfile->Close(); }
void buildTemplateforDiff(){
buildTemplate(-211);
buildTemplate(111);
buildTemplate(2212);
buildTemplate(2112);


}