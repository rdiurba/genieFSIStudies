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
TH1D* makeHistogram(TTree* t, std::string label){
std::cout<<t->GetEntries()<<","<<label<<std::endl;
TH1D *hist_Mult=new TH1D(Form("hist_Mult%s",label.c_str()),Form("hist_Mult%s",label.c_str()),60,0,60);
if (t->GetEntries()<1) return hist_Mult;
 
   double mk=0.4937; double mp=0.9383; double mpi=0.1396; double mmu=0.1057;
   int pdgf[100]; int nf; double Ef[100];
   int fspl; double El; int probe_fsi;
   int numEvents=t->GetEntries();
   t->SetBranchAddress("pdgh",&pdgf);
   t->SetBranchAddress("nh",&nf);
   t->SetBranchAddress("Eh",&Ef);
   t->SetBranchAddress("probe_fsi",&probe_fsi);
   Double_t weights20iP[100];
   Double_t weights20iPi[100];
   //t_weights->SetBranchAddress("tweak_responses_FSI_pi_VariationResponse",&weights20iPi);
   //t_weights->SetBranchAddress("tweak_responses_FSI_N_VariationResponse",&weights20iP);
   for (Long64_t j=0;j<numEvents;j++){
   int mult=0; t->GetEntry(j); // t_weights->GetEntry(j);
   if(probe_fsi!=6) continue;
   for(int i=0; i<nf; i++){

   //if (abs(pdgf[i])==211 && Ef[i]-mpi>threshPi.at(thresh)) mult++;
   //if (abs(pdgf[i])==321 && Ef[i]-mk>threshK.at(thresh)) mult++;
   if (abs(pdgf[i])==2212 && Ef[i]-mp>0.02) mult++;
   if (abs(pdgf[i])==2112 && Ef[i]-mp>0.02) mult++;
   //if (abs(pdgf[i])==13 && Ef[i]-mmu>threshMu.at(thresh)) mult++;
   }
   //if (abs(fspl)==13 && El-mmu>threshMu.at(thresh)) mult++;
    hist_Mult->Fill(mult);
   }



   return hist_Mult;


}



void plotKOSumThresh(int probe=2212, double ke=0.125, int target=1000822080)
{
std::string probeStr="protonPlus";
if (probe==-2212) probeStr="protonMinus";
if (probe==2112) probeStr="neutron";


std::string keStr="0.125";
if (ke==0.5) keStr="0.5";
if (ke==0.25) keStr="0.25";


gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);




  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);


   std::cout<<"Loading files"<<std::endl;
   TFile *f_numu = new TFile(Form("/exp/dune/data/users/rdiurba/genieFSIStudies/%s_%d_%sGeV_hA2018_1M.ginuke.root",probeStr.c_str(),target,keStr.c_str()));
  TFile *f_nue = new TFile(Form("/exp/dune/data/users/rdiurba/genieFSIStudies/%s_%d_%sGeV_hN2018_1M.ginuke.root",probeStr.c_str(),target,keStr.c_str()));
   TTree* t_numu = (TTree*)f_numu->Get("ginuke");
   TTree* t_nu20j = (TTree*)f_nue->Get("ginuke");


   TFile *f_numu10b = new TFile(Form("/exp/dune/data/users/rdiurba/genieFSIStudies/%s_%d_%sGeV_HINCL_1M.ginuke.root",probeStr.c_str(),target,keStr.c_str()));
   TFile *f_nue10b = new TFile(Form("/exp/dune/data/users/rdiurba/genieFSIStudies/%s_%d_%sGeV_HG4BertCasc_1M.ginuke.root",probeStr.c_str(),target,keStr.c_str()));
   TTree* t_nu20k = (TTree*)f_numu10b->Get("ginuke");
   TTree* t_nu20l = (TTree*)f_nue10b->Get("ginuke");
   std::cout<<"Loaded files"<<std::endl;
TH1D* mult20i=makeHistogram(t_numu,"20i");
TH1D* mult20j=makeHistogram(t_nu20j,"20j");
TH1D* mult20k=makeHistogram(t_nu20k,"20k");
TH1D* mult20l=makeHistogram(t_nu20l,"20l");
    mult20i->Scale(1.f/mult20i->GetEntries());  

   
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
   mult20i->SetTitle(Form("Number of Nucleons in Nucleon Knockout"));
   mult20i->GetXaxis()->SetTitle("N_{p}+N_{n} (KE Thres. of 20 MeV)");
   mult20i->GetYaxis()->SetTitle("Fraction of Events");
   gPad->SetLeftMargin(0.15);
   mult20i->GetYaxis()->SetRangeUser(0,0.8);
   mult20i->GetYaxis()->SetTitleOffset(1.4);
  TLegend *lErr = new TLegend(0.3,0.55,0.8,0.85);
  lErr->SetTextFont(133);
  lErr->SetTextSize(25);

  lErr->SetHeader(Form("p^{+} on Carbon with T=%s GeV",keStr.c_str()));
  if (probe==2212){
  if (target==1000822080) lErr->SetHeader(Form("p^{+} on Lead with T=%s GeV",keStr.c_str()));
  if (target==1000260560) lErr->SetHeader(Form("p^{+} on Iron with T=%s GeV",keStr.c_str()));
  }
    if (probe==-2212){
  lErr->SetHeader(Form("p^{-} on Carbon with T=%s GeV",keStr.c_str()));
  if (target==1000822080) lErr->SetHeader(Form("p^{-} on Lead with T=%s GeV",keStr.c_str()));
  if (target==1000260560) lErr->SetHeader(Form("p^{-} on Iron with T=%s GeV",keStr.c_str()));
  }
  if (probe==2112){
  lErr->SetHeader(Form("Neutron on Carbon with T=%s GeV",keStr.c_str()));
  if (target==1000822080) lErr->SetHeader(Form("Neutron on Lead with T=%s GeV",keStr.c_str()));
  if (target==1000260560) lErr->SetHeader(Form("Neutron on Iron with T=%s GeV",keStr.c_str()));
  }
 
  lErr->AddEntry(mult20i,"GENIE v3.4 hA2018","l");
  //lErr->AddEntry(multG21, "GENIE v3.4 SuSAv2 hA","l");
  lErr->AddEntry(mult20j,"GENIE v3.4 hN2018","l");
  lErr->AddEntry(mult20k,"GENIE v3.4 INCL++","l");
  lErr->AddEntry(mult20l,"GENIE v3.4 Geant4","l");
  
   mult20i->GetXaxis()->CenterTitle();
   mult20i->GetYaxis()->CenterTitle();
   mult20i->Draw("HIST");
   mult20j->Draw("HIST SAME");
   mult20k->Draw("HIST SAME");
   mult20l->Draw("HIST SAME");
   lErr->Draw("SAME");
   c1->Print(Form("nucleonThreshSumKO_%s_%s_%d.png",probeStr.c_str(),keStr.c_str(),target));
      c1->Print(Form("nucleonThreshSumKO_%s_%s_%d.pdf",probeStr.c_str(),keStr.c_str(),target));
   TF1* g1=new TF1("g1","expo",0,40);
   TF1* g2=new TF1("g2","expo",0,40);
   TF1* g3=new TF1("g3","expo",0,40);
   TF1* g4=new TF1("g4","expo",0,40);
   mult20i->Fit(g1);
   mult20j->Fit(g2);
   mult20k->Fit(g3);
   mult20l->Fit(g4);

std::ofstream file;
file.open("nucleonKOSumThresh.txt",std::ios_base::app);
file<<probe<<","<<ke<<","<<target<<","<<g1->GetParameter(0)<<","<<g1->GetParError(0)<<","<<g1->GetParameter(1)<<","<<g1->GetParError(1)<<","<<g2->GetParameter(0)<<","<<g2->GetParError(0)<<","<<g2->GetParameter(1)<<","<<g2->GetParError(1)<<","<<g3->GetParameter(0)<<","<<g3->GetParError(0)<<","<<g3->GetParameter(1)<<","<<g3->GetParError(1)<<","<<g4->GetParameter(0)<<","<<g4->GetParError(0)<<","<<g4->GetParameter(1)<<","<<g4->GetParError(1)<<std::endl;
file.clear();



}
