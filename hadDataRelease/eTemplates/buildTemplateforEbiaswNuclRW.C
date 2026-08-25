#include "TH1D.h"
#include <TCanvas.h>
#include <THStack.h>
#include <TLatex.h>
#include <TMath.h>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include "FSIReweight.h"
enum Mode { kKO, kAbs };
#include <map>
#include <variant>
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
TH2* getHist(std::string reweightModel, int pdg, std::string particle, int target) { // get histogram from the TTree (e.g. gst/ginuke file)

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
    TH2D* h_KEini_Ebias = new TH2D(Form("rw%s_%s_%d_KEini_vs_Ebias", reweightModel.c_str(), particle.c_str(),target), Form("rw%s_%s_%d_KEini_vs_Ebias", reweightModel.c_str(), particle.c_str(),target), 40, 0, 2, 20, -1, 1);

    TH2D* hist;
    for (int i=1; i<40; i++){
    TFile* file = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%gGeV_hA2018_100k.ginuke.root",particle.c_str(),target,v.at(i) ));



  
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
        bool compound=1;
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

int max=40;
 double hAEst=1;
    double AltEst=1;

     bool isPion = (probe==211 || probe==-211 || probe==111);

            computeMultDiffWeights(pdg, target, mult, diff, KEini, max,
                                   reweightModel, isPion,
                                   hAEst, AltEst);
 
            double rw = (AltEst / hAEst);
         
        
        double Ebias;
        if (pdg==2212 || pdg==2112)
          Ebias = (KEini - Ehad) / KEini;
        else if (pdg==211 || pdg==111 || pdg==-211)
          Ebias = (Eini - Ehad) / Eini;
        //h_Ebias->Fill(Ebias);
        h_KEini_Ebias->Fill(KEini,  Ebias,rw);
      }
    }
     file->Close();

    }
    
        normalizePerEnergy(h_KEini_Ebias);
  return h_KEini_Ebias;
}

void buildTemplate(int pdg=211, int target=1000180400) {

  TH2 *hA2018_hist, *hN2018_hist, *INCL_hist, *G4BC_hist;
  string particle;

    if (pdg == 2212) particle = "protonPlus";
    else if (pdg == 2112) particle = "neutron";
    else if (pdg == 211) particle = "piPlus";
    else if (pdg == 111) particle = "pi0";
    else if (pdg == -211) particle = "piMinus";
    hN2018_hist = (TH2*)getHist(Form("hN2018"), pdg, particle,target);
    INCL_hist = (TH2*)getHist(Form("HINCL"), pdg, particle,target);
    G4BC_hist = (TH2*)getHist(Form("HG4BertCasc"), pdg, particle,target);
  
  
  
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

  cc1->Print(Form("RW2hN2018_%s_%d_Ebias.png",particle.c_str(),target) );
  INCL_hist->GetZaxis()->SetRangeUser(0,0.5);

  INCL_hist->Draw("colz");
  cc1->Print(Form("RW2INCL_%s_%d_Ebias.png",particle.c_str(),target) );
  G4BC_hist->GetZaxis()->SetRangeUser(0,0.5);

  G4BC_hist->Draw("colz");
  cc1->Print(Form("RW2G4BC_%s_%d_Ebias.png",particle.c_str(),target) );

  hN2018_hist->SetName(Form("rwhN2018_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  INCL_hist->SetName(Form("rwINCL++_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  G4BC_hist->SetName(Form("rwGeant4_%s_%d_KEini_vs_Ebias", particle.c_str(),target));

  hN2018_hist->SetTitle(Form("rwhN2018_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  INCL_hist->SetTitle(Form("rwINCL++_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  G4BC_hist->SetTitle(Form("rwGeant4_%s_%d_KEini_vs_Ebias", particle.c_str(),target));

  TFile *outfile = new TFile("FSI_KOAbs_Evis_reweight_template.root","UPDATE");
  hN2018_hist->Write(Form("rwhN2018_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  INCL_hist->Write(Form("rwINCL++_%s_%d_KEini_vs_Ebias", particle.c_str(),target));
  G4BC_hist->Write(Form("rwGeant4_%s_%d_KEini_vs_Ebias", particle.c_str(),target));

  outfile->Write();  outfile->Close(); }
void buildTemplateforEbiaswNuclRW(){
makeAllParams();
buildTemplate(211,1000060120);
buildTemplate(-211,1000060120);
buildTemplate(111,1000060120);
buildTemplate(2212,1000060120);
buildTemplate(2112,1000060120);
buildTemplate(211,1000180400);
buildTemplate(-211,1000180400);
buildTemplate(111,1000180400);
buildTemplate(2212,1000180400);
buildTemplate(2112,1000180400);

buildTemplate(211,1000260560);
buildTemplate(-211,1000260560);
buildTemplate(111,1000260560);
buildTemplate(2212,1000260560);
buildTemplate(2112,1000260560);


}