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
TH1D* makeHistogram(TTree* t, std::string label,  int numEvents=0, int thresh=0){
TH1D *hist_Mult=new TH1D(Form("hist_Mult%s",label.c_str()),Form("hist_Mult%s",label.c_str()),15,0,15);
   std::vector<double> threshPi={0.0,0.0067,0.00895,0.01335,0.0249,0.0338};
   std::vector<double> threshMu={0,0.0059,0.00795,0.01188,0.0222,0.0304};
   std::vector<double> threshP={0,0.015,0.0206,0.03065,0.05655,0.07575};
   std::vector<double> threshK={0,0.0116,0.01555,0.02315,0.04275,0.0574};
   double mk=0.4937; double mp=0.9383; double mpi=0.1396; double mmu=0.1057;
   int pdgf[100]; int nf; double Ef[100];
   int fspl; double El;
   if (numEvents==0) numEvents=t->GetEntries();
   t->SetBranchAddress("pdgf",&pdgf);
   t->SetBranchAddress("nf",&nf);
   t->SetBranchAddress("Ef",&Ef);
   t->SetBranchAddress("fspl",&fspl);
   t->SetBranchAddress("El",&El);
   for (Long64_t j=0;j<numEvents;j++){
   int mult=0; t->GetEntry(j);
   for(int i=0; i<nf; i++){
   if (abs(pdgf[i])==211 && Ef[i]-mpi>threshPi.at(thresh)) mult++;
   if (abs(pdgf[i])==321 && Ef[i]-mk>threshK.at(thresh)) mult++;
   if (abs(pdgf[i])==2212 && Ef[i]-mp>threshP.at(thresh)) mult++;
   if (abs(pdgf[i])==13 && Ef[i]-mmu>threshMu.at(thresh)) mult++;
   }
   if (abs(fspl)==13 && El-mmu>threshMu.at(thresh)) mult++;
   hist_Mult->Fill(mult);
   }
   return hist_Mult;
}
TProfile* makeProfile(TTree* t, TTree* t_weights,  std::string label,  int numEvents=0, int thresh=0){
TProfile *profile_Mult=new TProfile(Form("profile_Mult%s",label.c_str()),Form("profile_Mult%s",label.c_str()),15,0,15,"s");
std::vector<TH1F*> hist_Mult20i_multisimVec;
std::vector<TH1F*> hist_Nnuc20i_multisimVec;
for(int i=0; i<100; i++){
TH1F* tempPiHist=new TH1F(Form("hist_Npi20i_multisim_%d",i),Form("hist_Npi20i_multisim_%d",i),15,0,15);
hist_Mult20i_multisimVec.push_back(tempPiHist);
}
TH1D *hist_Mult20i=new TH1D(Form("hist_Mult20i%s",label.c_str()),Form("hist_Mult20i%s",label.c_str()),15,0,15);
   std::vector<double> threshPi={0.0,0.0067,0.00895,0.01335,0.0249,0.0338};
   std::vector<double> threshMu={0,0.0059,0.00795,0.01188,0.0222,0.0304};
   std::vector<double> threshP={0,0.015,0.0206,0.03065,0.05655,0.07575};
   std::vector<double> threshK={0,0.0116,0.01555,0.02315,0.04275,0.0574};
   double mk=0.4937; double mp=0.9383; double mpi=0.1396; double mmu=0.1057;
   int pdgf[100]; int nf; double Ef[100];
   int fspl; double El;
   if (numEvents==0) numEvents=t->GetEntries();
   t->SetBranchAddress("pdgf",&pdgf);
   t->SetBranchAddress("nf",&nf);
   t->SetBranchAddress("Ef",&Ef);
   t->SetBranchAddress("fspl",&fspl);
   t->SetBranchAddress("El",&El);
   Double_t weights20iP[100];
   Double_t weights20iPi[100];
   t_weights->SetBranchAddress("tweak_responses_FSI_pi_VariationResponse",&weights20iPi);
   t_weights->SetBranchAddress("tweak_responses_FSI_N_VariationResponse",&weights20iP);
   for (Long64_t j=0;j<numEvents;j++){
   int mult=0; t->GetEntry(j); t_weights->GetEntry(j);
   for(int i=0; i<nf; i++){
   if (abs(pdgf[i])==211 && Ef[i]-mpi>threshPi.at(thresh)) mult++;
   if (abs(pdgf[i])==321 && Ef[i]-mk>threshK.at(thresh)) mult++;
   if (abs(pdgf[i])==2212 && Ef[i]-mp>threshP.at(thresh)) mult++;
   if (abs(pdgf[i])==13 && Ef[i]-mmu>threshMu.at(thresh)) mult++;
   }
   if (abs(fspl)==13 && El-mmu>threshMu.at(thresh)) mult++;
    hist_Mult20i->Fill(mult);

    for (int i=0; i<100; i++){  
    if(isnan(weights20iPi[i]) || weights20iPi[i]<0 || isnan(weights20iP[i]) || weights20iP[i]>10) continue;
    //std::cout<<weights20iPi[i]<<std::endl;
    hist_Mult20i_multisimVec.at(i)->Fill(mult,weights20iPi[i]*weights20iP[i]);
	}
   }
     for (int i=0; i<100; i++){
      hist_Mult20i_multisimVec.at(i)->Scale(1.f/hist_Mult20i->GetEntries());
     for (int index=1; index<16; index++){
     profile_Mult->Fill(profile_Mult->GetXaxis()->GetBinCenter(index),hist_Mult20i_multisimVec.at(i)->GetBinContent(index));
      std::cout<<profile_Mult->GetBinContent(index)<<","<<hist_Mult20i_multisimVec.at(i)->GetBinContent(index)<<std::endl;
      }}
   return profile_Mult;
}

void plotMultiplicity(int thresh=0)
{


gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);




  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);


   std::cout<<"Loading files"<<std::endl;
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
 TH1D* mult20i_copy=(TH1D*)mult20i->Copy("mult20i_copy");   
 mult20i->Scale(1.f/mult20i->GetEntries());  
   for (int index=1; index<11; index++){
     mult20i->SetBinError(index,profile20i->GetBinError(index));
   std::cout<<mult20i->GetBinContent(index)<<","<<profile20i->GetBinContent(index)<<std::endl;
      }
   
   mult20j->Scale(1.f/mult20j->GetEntries());
   mult20l->Scale(1.f/mult20l->GetEntries());
   mult20k->Scale(1.f/mult20k->GetEntries());
   mult20i->SetLineColor(kBlack);
   mult20j->SetLineStyle(2);
   mult20j->SetLineColor(kRed);
   mult20l->SetLineStyle(3);
   mult20l->SetLineColor(kGreen);
   mult20k->SetLineStyle(4);
   mult20k->SetLineColor(kBlue);
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
   mult20i->GetYaxis()->SetRangeUser(0,0.4);
   mult20i->GetYaxis()->SetTitleOffset(1.4);
  TLegend *lErr = new TLegend(0.4,0.55,0.8,0.85);
  lErr->SetTextFont(133);
  lErr->SetTextSize(25);
  lErr->AddEntry(mult20i,"GENIE v3.4 AR23 hA","l");
  //lErr->AddEntry(multG21, "GENIE v3.4 G21 hA","l");
  lErr->AddEntry(mult20j,"GENIE v3.4 AR23 hN","l");
  lErr->AddEntry(mult20k,"GENIE v3.4 AR23 INCL++","l");
  lErr->AddEntry(mult20l,"GENIE v3.4 AR23 Geant4","l");
  TLegend *lOnlyExtreme = new TLegend(0.4,0.55,0.8,0.85);
  lOnlyExtreme->SetTextFont(133);
  lOnlyExtreme->SetTextSize(25);
  lOnlyExtreme->AddEntry(mult20i,"hA with FSI Uncert.","l");
  lOnlyExtreme->AddEntry(mult20k,"INCL++","l");
  lOnlyExtreme->AddEntry(mult20l,"Geant4","l");
   mult20i->GetXaxis()->CenterTitle();
   mult20i->GetYaxis()->CenterTitle();
   mult20i->Draw("HIST");
   mult20j->Draw("HIST SAME");
   mult20k->Draw("HIST SAME");
   mult20l->Draw("HIST SAME");
   lErr->Draw("SAME");
   c1->Print(Form("multNominal_%d.png",thresh));
   mult20i->Draw("E0 P");
   mult20k->Draw("HIST SAME");
   mult20l->Draw("HIST SAME");
   lOnlyExtreme->Draw("SAME");
   c1->Print(Form("multNominalWeights_%d.png",thresh));
}
