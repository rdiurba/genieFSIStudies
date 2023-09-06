#include <iostream>
#include <TH2.h>
#include <TH1.h>
#include <TH1D.h>
#include <chrono>
#include <thread>
#include <vector>
#include <TStyle.h>
#include <TCanvas.h>
#include "TObject.h"
#include "TString.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TMath.h"
#include <string>
#include <vector>
#include <map>
#include "TDecompChol.h"

int
   lowx = 0,
   highx = 15,
   bins = 15;

double threshPi[6] = {0, 0.0067, 0.00895, 0.01335, 0.02490, 0.03380};     // 0 mm, 3 mm, 5 mm, 1 cm, 3 cm, 5 cm
double threshMu[6] = {0, 0.0059, 0.00795, 0.01188, 0.02220, 0.03040};
double threshP[6]  = {0, 0.0150, 0.02060, 0.03065, 0.05655, 0.07575};     // threshholds for pion, muon, proton, K+ meson
double threshK[6]  = {0, 0.0116, 0.01555, 0.02315, 0.04275, 0.05740};     // in GeV

double
   mpi = 0.1396,
   mmu = 0.1057,
   mp  = 0.9383,
   mk  = 0.4937;     // masses of particles in GeV

std::vector<TH1F*> g_mult20i;

TMatrixD covMatrix(std::vector<TH1F*> multihist) {
   
   TMatrixD covmat(bins, bins);
   double covariance, avg1, avg2;
   int histamount = multihist.size();
   
   for (Int_t bin1 = 0; bin1 < bins; bin1++) {
     for (Int_t bin2 = 0; bin2 < bins; bin2++) {
         covariance = 0, avg1 = 0, avg2 = 0;
         for (int i = 0; i < histamount; ++i) {
            avg1 += multihist[i]->GetBinContent(bin1+1);
            avg2 += multihist[i]->GetBinContent(bin2+1);
         }
         avg1 /= histamount; avg2 /= histamount;
         
         for (int i = 0; i < histamount; ++i) {
            covariance += (multihist[i]->GetBinContent(bin1+1) - avg1) * (multihist[i]->GetBinContent(bin2+1) - avg2);
         }
         covariance /= (histamount-1);
         
         covmat(bin1, bin2) = covariance;
      }
   }
   
   return covmat;
};



TH1D* makeHistogram(TTree* t, std::string label,  int numEvents=0, int thresh=0) {                       // define a custom histogram with listed arguments, thresh corresponds to a track length
   
   TH1D *hist_Mult = new TH1D(Form("hist_Mult%s", label.c_str()), Form("hist_Mult%s", label.c_str()), bins, lowx, highx); 
   
   
   int pdgf[100]; int nf; double Ef[100];    // pdgf - particle code, nf - amount of particles in final state, Ef - energy of particles in final state. [100] apparently is expected to be bigger than nf
   int fspl; double El;    // fspl is also a particle code, El is also an energy
   
   if (numEvents==0) numEvents=t->GetEntries();    // t->GetEntries() returns amount of entries, here numEvent ~ 200000
   
   //std::cout << t->GetEntries();
   
   t->SetBranchAddress("pdgf",&pdgf);
   t->SetBranchAddress("nf",&nf);      // there are nfp, nfm, nfpip, nfpim, but not nf. Meaning that "nf" hoards them all?
   t->SetBranchAddress("Ef",&Ef);
   t->SetBranchAddress("fspl",&fspl);
   t->SetBranchAddress("El",&El);
   
   for (Long64_t j=0;j<numEvents;j++) {      // cycle through all events
      
      int mult = 0; t->GetEntry(j);    // GetEntry reads all branches of entry and returns total number of bytes read.  
      
      
      
      for(int i=0; i<nf; i++) {
         if (abs(pdgf[i])==211  && Ef[i]-mpi > threshPi[thresh]) mult++;      // abs() to count antiparticles as well?
         if (abs(pdgf[i])==13   && Ef[i]-mmu > threshMu[thresh]) mult++;
         if (abs(pdgf[i])==2212 && Ef[i]-mp  > threshP[thresh])  mult++;
         if (abs(pdgf[i])==321  && Ef[i]-mk  > threshK[thresh])  mult++;
      }
      
      if (abs(fspl)==13 && El-mmu>threshMu[thresh]) mult++;    // El is just like Ef, but for muons?
      hist_Mult->Fill(mult);
      
   }
   
   return hist_Mult;
   
}
   
TProfile* makeProfile(TTree* t, TTree* t_weights,  std::string label,  int numEvents=0, int thresh=0) {     // create a TProfile fonction
   
   TProfile *profile_Mult = new TProfile(Form("profile_Mult%s", label.c_str()), Form("profile_Mult%s", label.c_str()), bins, lowx, highx, "s");    // arguments are name and title, amount of bins, low x, high x, option "s" stands for "standart deviation"
   
   std::vector<TH1F*> hist_Mult20i_multisimVec;
   std::vector<TH1F*> hist_Nnuc20i_multisimVec;    // these are *vectors* of floating point histograms
   
   for(int i=0; i<100; i++) {
      TH1F* tempPiHist = new TH1F(Form("hist_Npi20i_multisim_%d", i), Form("hist_Npi20i_multisim_%d", i), bins, lowx, highx);
      hist_Mult20i_multisimVec.push_back(tempPiHist);    // 100 histograms are created and put into hist_Mult20i_multisimVec
   }
   
   TH1D *hist_Mult20i = new TH1D(Form("hist_Mult20i%s", label.c_str()), Form("hist_Mult20i%s", label.c_str()), bins, lowx, highx);
   
   int pdgf[100]; int nf; double Ef[100];
   int fspl; double El;
   if (numEvents==0) numEvents=t->GetEntries();
   
   std::vector<int> tmpvector(100,0);
   
   t->SetBranchAddress("pdgf",&pdgf);
   t->SetBranchAddress("nf",&nf);
   t->SetBranchAddress("Ef",&Ef);
   t->SetBranchAddress("fspl",&fspl);
   t->SetBranchAddress("El",&El);
   
   Double_t weights20iP[100];       // weights for proton
   Double_t weights20iPi[100];      // weights for pion
   
   t_weights->SetBranchAddress("tweak_responses_FSI_pi_VariationResponse",&weights20iPi);
   t_weights->SetBranchAddress("tweak_responses_FSI_N_VariationResponse",&weights20iP);
   
   
   
   for (Long64_t j=0;j<numEvents;j++) { 
      
      int mult=0; t->GetEntry(j); t_weights->GetEntry(j);
      
      for (int i=0; i<nf; i++) {
         if (abs(pdgf[i])==211  && Ef[i]-mpi > threshPi[thresh]) mult++;
         if (abs(pdgf[i])==13   && Ef[i]-mmu > threshMu[thresh]) mult++;
         if (abs(pdgf[i])==2212 && Ef[i]-mp  > threshP[thresh])  mult++;
         if (abs(pdgf[i])==321  && Ef[i]-mk  > threshK[thresh])  mult++;
      }
      
      if (abs(fspl)==13 && El-mmu>threshMu[thresh]) mult++;
      hist_Mult20i->Fill(mult);
      
      for (int i=0; i<100; i++){  
         if(std::isnan(weights20iPi[i]) || weights20iPi[i]<0 || std::isnan(weights20iP[i]) || weights20iP[i]>10) continue;    // avoid this entry if it's invalid
         hist_Mult20i_multisimVec.at(i)->Fill(mult, weights20iPi[i] * weights20iP[i]);    // fill histogram with weight
      }
      
   }
   
   for (int i=0; i<100; i++){
    //g_mult20i.push_back(hist_Mult20i_multisimVec.at(i));
  
      hist_Mult20i_multisimVec.at(i)->Scale(1.f/hist_Mult20i->GetEntries());     // normalizing the weight distribution
      
      for (int index=1; index<bins+1; index++){
         
         profile_Mult->Fill(profile_Mult->GetXaxis()->GetBinCenter(index),hist_Mult20i_multisimVec.at(i)->GetBinContent(index));    // for each bin index profile_Mult get filled with bin X, and from 100 histograms
         //std::cout<<profile_Mult->GetBinContent(index)<<","<<hist_Mult20i_multisimVec.at(i)->GetBinContent(index)<<std::endl;
         
      }
   }
   
   for(int i=0; i<100; i++) {
      
      g_mult20i.push_back(hist_Mult20i_multisimVec.at(i));
   }
   
   
   return profile_Mult;
   
}

void plotMultWithCovar(int thresh=0) {
   
   gROOT->LoadMacro("protoDUNEStyle.C");
   gROOT->SetStyle("protoDUNEStyle");
   gROOT->ForceStyle();
   
   gStyle->SetTitleX(0.35);
   gStyle->SetOptFit(111);
   
   
   TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);
   
   std::cout<<"Loading files"<<std::endl;
   //exit(0);
   
   TFile *f_numu          = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/2x2_AR23_20i_00_000_2M.gst.root");
   TFile *f_weights       = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/AR23_20i_2x2_fsiRW_400k.root");
   TFile *f_weights1sigma = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/AR23_20i_2x2_fsiInclRW.root");
   TFile *f_nue           = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/2x2_AR23_20j_00_000_2M.gst.root");
   TFile *f_numu10b       = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/2x2_AR23_20k_00_000_212k.gst.root");
   TFile *f_nue10b        = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/2x2_AR23_20l_00_000_500k.gst.root");
   TFile *f_nu21 = new TFile("/pnfs/dune/persistent/users/rdiurba/2x2_RHC_gst/2x2_G21_11a_00_000_200k.gst.root"); 
   TTree* t_numu    = (TTree*) f_numu->Get("gst");
   TTree* t_weights = (TTree*) f_weights->Get("events");
   TTree* t_1sigma  = (TTree*) f_weights1sigma->Get("events");
   TTree* t_nu20j   = (TTree*) f_nue->Get("gst");
   TTree* t_nu20k   = (TTree*) f_numu10b->Get("gst");
   TTree* t_nu20l   = (TTree*) f_nue10b->Get("gst");
   TTree* t_nu21    = (TTree*) f_nu21->Get("gst");
   std::cout<<"Loaded files"<<std::endl;
   
   TH1D* mult20i = makeHistogram(t_numu,"20i",400000,thresh);
   TH1D* mult20j = makeHistogram(t_nu20j,"20j",2000000,thresh);
   TH1D* mult20k = makeHistogram(t_nu20k,"20k",200000,thresh);
   TH1D* mult20l = makeHistogram(t_nu20l,"20l",500000,thresh);
   TH1D* profile20i = makeProfile(t_numu,t_weights,"20l",400000,thresh);
   TH1D* mult21 = makeHistogram(t_nu21,"G21",200000,thresh);
   //TH1D* mult20i_copy=(TH1D*) mult20i->Clone("mult20i_copy");
   //TH1D* mult20i_allUncert=(TH1D*) mult20i_allUncert->Clone("mult20i_allUncert"); 
   mult20i->Scale(1.f/mult20i->GetEntries());  
  mult20j->Scale(1.f/mult20j->GetEntries());
  mult20l->Scale(1.f/mult20l->GetEntries());
  mult20k->Scale(1.f/mult20k->GetEntries());
    mult21->Scale(1.f/mult21->GetEntries());
    TH1D* mult20i_copy=(TH1D*) mult20i->Clone("mult20i_copy");
 


 
   for (int index=1; index<bins; index++) {
      mult20i->SetBinError(index,profile20i->GetBinError(index));    // get bin errors from profile20i to mult20i
      // std::cout<<mult20i->GetBinContent(index)<<","<<profile20i->GetBinContent(index)<<std::endl;
      profile20i->SetBinError(index,TMath::Sqrt(TMath::Power(profile20i->GetBinError(index),2)+TMath::Power(mult20i_copy->GetBinError(index),2)));


    }
  
   mult20i->SetLineColor(kBlack);
   mult20j->SetLineStyle(2);
   mult20j->SetLineColor(kRed);
   mult20l->SetLineStyle(3);
   mult20l->SetLineColor(kGreen);
   mult20k->SetLineStyle(4);
   mult20k->SetLineColor(kBlue);
   mult21->SetLineStyle(5);
   mult21->SetLineColor(kGreen+1); 
   std::string label="0 mm";
   if (thresh==1) label="3 mm";
   if (thresh==2) label="5 mm";
   if (thresh==3) label="1 cm";
   if (thresh==4) label="3 cm";
   if (thresh==5) label="5 cm";
   
   mult20i->SetTitle(Form("Multiplicity of Particles with Min. Length %s",label.c_str()));
   mult20i->GetXaxis()->SetTitle("Number of Charged Particles");
   mult20i->GetYaxis()->SetTitle("Fraction of Events");
   
   gPad->SetLeftMargin(0.15);
   mult20i->GetYaxis()->SetRangeUser(0,2.0*mult20i->GetMaximum());
   mult20i->GetYaxis()->SetTitleOffset(1.4);
   
   TLegend *lErr = new TLegend(0.4,0.55,0.8,0.85);
   
   lErr->SetTextFont(133);
   lErr->SetTextSize(25);
   lErr->AddEntry(mult20i,"GENIE v3.4 hA","l");
   lErr->AddEntry(mult20j,"GENIE v3.4 hN","l");
   lErr->AddEntry(mult20k,"GENIE v3.4 INCL++","l");
   lErr->AddEntry(mult20l,"GENIE v3.4 Geant4","l");
   
   TLegend *lOnlyExtreme = new TLegend(0.4,0.55,0.8,0.85);
   
   lOnlyExtreme->SetTextFont(133);
   lOnlyExtreme->SetTextSize(25);
   lOnlyExtreme->AddEntry(mult20i,"AR23 hA (stat.+FSI uncert.)","l");
   lOnlyExtreme->AddEntry(mult20j,"AR23 hN","l");
   //lOnlyExtreme->AddEntry(mult20l,"AR23 Geant4","l");
   lOnlyExtreme->AddEntry(mult21,"SUSAv2 hA","l");
   mult20i->GetXaxis()->CenterTitle();
   mult20i->GetYaxis()->CenterTitle();
   mult20i->Draw("HIST");
   mult20j->Draw("HIST SAME");
   mult20k->Draw("HIST SAME");
   mult20l->Draw("HIST SAME");
   
   lErr->Draw("SAME");
   c1->Print(Form("P_multNominal_%d.png",thresh));
   mult20i->Draw("E0 P");
   mult20j->Draw("HIST SAME");
   //mult20l->Draw("HIST SAME");
   mult21->Draw("HIST SAME");
   lOnlyExtreme->Draw("SAME");
   c1->Print(Form("P_multNominalWeights_%d.png",thresh));
   
   
//   cout << g_mult20i[0]->GetBinContent(5) << '\n';
   
   TMatrixD covarianceMatrix = covMatrix(g_mult20i);
   
   
   // Invert the covariance matrix
   TDecompChol invMC(covarianceMatrix);
   bool test = invMC.Decompose();      // check if matrix positive definite
   
   auto inverseCM = invMC.Invert();    // inverse matrix
   auto diagonalMatrix=covarianceMatrix;
   auto covarianceMatrixAll=covarianceMatrix;
   auto diagonalMatrixAll=covarianceMatrix;
   for (int i = 0; i < bins; ++i) {
      for (int j = 0; j < bins; ++j) {
         if (i!=j) diagonalMatrix(i,j)=0;
         else{ diagonalMatrixAll(i,j)=diagonalMatrix(i,j)+TMath::Power(mult20i_copy->GetBinError(i+1),2)+TMath::Power(mult20j->GetBinError(i+1),2);
           covarianceMatrixAll(i,j)=diagonalMatrixAll(i,j);


        }
      }
  }
   auto inverseDiagonal=diagonalMatrix.Invert();
   TDecompChol invCMAll(covarianceMatrixAll);
   auto inverseCMAll=invCMAll.Invert(); 
   
   
   // Calculate the chi2 for any histogram (hN, INCL++, Geant4)
   double delta1, delta2;
   
   double chi2 = 0;
   double chi2Diagonal=0;
   double chi2FromProfile=0;
   double chi2All=0;
   for (int bin1 = 0; bin1 < bins; ++bin1) {
      
      delta1 = profile20i->GetBinContent(bin1+1) - mult20j->GetBinContent(bin1+1);
      
      chi2FromProfile+=(delta1*delta1)/TMath::Power(profile20i->GetBinError(bin1+1),2);
      //std::cout<<"inverse of diag. covar,Inverse, inverse full from diagonal, std from covar, profile covar: "<<inverseDiagonal(bin1,bin1)<<","<<inverseCM(bin1,bin1)<<","<<TMath::Sqrt(1.f/covarianceMatrix(bin1,bin1))<<","<<TMath::Sqrt(covarianceMatrix(bin1,bin1))<<","<<profile20i->GetBinError(bin1+1)<<std::endl;
      for (int bin2 = 0; bin2 < bins; ++bin2) {
         
         delta2 = profile20i->GetBinContent(bin2+1) - mult20j->GetBinContent(bin2+1);
         if (bin1==bin2) chi2Diagonal += delta1 * inverseDiagonal(bin1, bin2) * delta2;
         chi2 += delta1 * inverseCM(bin1, bin2) * delta2;
         chi2All += delta1*inverseCMAll(bin1,bin2)*delta2;     
      }
   }
   
   cout << "chi2 (syst. only, all): = " << chi2 << ","<<chi2All<<'\n';
   cout << "chi2 from diagonal (covar and profile): "<<chi2Diagonal<<","<<chi2FromProfile<<std::endl;
   //double critChi2 = TMath::ChisquareQuantile(0.95, bins);
   //cout << "critChi2 = " << critChi2 << '\n';
   TF1* fitLandau=new TF1("landFit","landau");   
   profile20i->Fit("landFit");
   TH2F *heatmap1 = new TH2F("h1", "Values of the covariance matrix", bins, 0, bins, bins, 0, bins);
   TH2F *heatmap2 = new TH2F("h1", "Values of the inverse covariance matrix", bins, 0, bins, bins, 0, bins);
   
   for (int i = 0; i < bins; ++i) {
       for (int j = 0; j < bins; ++j) {
           heatmap1->Fill(i, j, covarianceMatrix(i,j));
           heatmap2->Fill(i, j, inverseCM(i,j));
       }
   }
   
   // Create a canvas to display the heatmap
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
   gPad->SetRightMargin(0.15);
   heatmap1->Draw("COLZ"); 
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
   gPad->SetRightMargin(0.15);
   heatmap2->Draw("COLZ"); 
   
}
