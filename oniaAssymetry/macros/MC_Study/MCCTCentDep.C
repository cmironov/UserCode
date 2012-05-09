//===================== Macro to calculate PrJpsi Eff or Acc x Eff ============================================
//read tree from DiMuonOnia2DPlots.cc
#ifndef __CINT__
#endif
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "FitFunctions.h"
#include "TObjArray.h"
//iSpec = 1: pT, = 2: Rap, = 3: Cent, = 4: Phi, = 5: dPhi, PutWeight = 1 means put weight
// Prompt =1 prompt PrJpsi 
bool IsAccept(Double_t pt, Double_t eta); //you can define acceptance here 
double FindCenWeight(int Bin);//gives you weight according to cent
double FindV2Weight(double v2wg);
void DoEffCor3D(TFile* a, int b, double c, double d, double e, double *f);
void DoEffCor1D(TFile* a, int b, double i, double *c);
void MCCTCentDep()   
{
  for(int iSpec = 0; iSpec < 3; iSpec++){
    int Prompt =1; int PutWeight = 1;
    bool bDefault = true; // true : cowboy only, false : salior
    bool bCowboy = false; // true : cowboy only, false : salior
    bool bSailor = false; // true : cowboy only, false : salior
    double fake_v2 = 0.3;

    // iSpec : choose the condition for default/cowboy/sailor
    if(iSpec == 0) {bDefault = true;bCowboy = false;bSailor = false;}
    if(iSpec == 1) {bDefault = false;bCowboy = true;bSailor = false;}
    if(iSpec == 2) {bDefault = false;bCowboy = false;bSailor = true;}

    char cCond[512];
    int iCond = 0; // the number of the cases, 1 : default, 2 : cowboy, 3 : sailor
    if(bDefault) {sprintf(cCond, "default"); iCond = 0;}
    if(bCowboy)  {sprintf(cCond, "cowboy"); iCond = 1;}
    if(bSailor)  {sprintf(cCond, "sailor"); iCond = 2;}

    // iCat : decide the categories for Centrality (1: 0 - 10, 2: 10 - 20, 3: 20 - 30, 4: 30 - 60 %)
    for(int iCat = 0; iCat < 4; iCat++){
      gROOT->SetStyle("Plain");
      gStyle->SetPalette(1);
      gStyle->SetFrameBorderMode(0);
      gStyle->SetFrameFillColor(0);
      gStyle->SetCanvasColor(0);
      gStyle->SetTitleFillColor(0);
      gStyle->SetStatColor(0);
      gStyle->SetPadBorderSize(0);
      gStyle->SetCanvasBorderSize(0);
      gStyle->SetOptTitle(0); // at least most of the time
      gStyle->SetOptStat(1); // most of the time, sometimes "nemriou" might be useful to display name, 
      gStyle->SetOptFit(1); // set to 1 only if you want to display fit results
      //==================================== Define Histograms====================================================
      char OutTextFile[100]; 
      if(iCat == 0) sprintf(OutTextFile,"MCCT_PrJpsi_Cent_0010_dPhi_%s.tex", cCond);
      if(iCat == 1) sprintf(OutTextFile,"MCCT_PrJpsi_Cent_1020_dPhi_%s.tex", cCond);
      if(iCat == 2) sprintf(OutTextFile,"MCCT_PrJpsi_Cent_2030_dPhi_%s.tex", cCond);
      if(iCat == 3) sprintf(OutTextFile,"MCCT_PrJpsi_Cent_3060_dPhi_%s.tex", cCond);
      ofstream dataFile(Form(OutTextFile));
      TH1D *diMuonsInvMass_Gen = new TH1D("diMuonsInvMass_Gen","diMuonsInvMass_Gen", 100,2.95,3.25);

      TH1D *diMuonsPt_Gen = new TH1D("diMuonsPt_Gen","diMuonsPt_Gen", 100,0,50);
      TH1D *Bin_Gen = new TH1D("Bin_Gen","Bin_Gen", 40,0,40);
      //==============================================Define Acc Eff Stuff here===========================================
      // Pt bin sizes
      // 0-1.5, 1.5-3, 3-4.5, 4.5-6, 6-7.5...

      const int ndPhiBins = 4;
      const int nFiles = 6;
      double dphi_bound[100] = {0};

      dphi_bound[0] = 0.0;
      dphi_bound[1] = TMath::Pi()*2/16;
      dphi_bound[2] = TMath::Pi()*4/16;
      dphi_bound[3] = TMath::Pi()*6/16;
      dphi_bound[4] = TMath::Pi()*8/16;

      //X Axis error on Eff graph 
      double xdphi_bound[ndPhiBins] = {0.0};
      for(int i = 0; i < ndPhiBins; i++){
        xdphi_bound[i] = dphi_bound[i] + (dphi_bound[i+1]-dphi_bound[i])/2;
        cout<<"xdphi_bound["<<i<<"] : "<<xdphi_bound[i]<<endl;
        dataFile<<"xdphi_bound["<<i<<"] : "<<xdphi_bound[i]<<endl;
      }

      double genError, recError;
      double gen_pt[100]={0}, gen_ptError[100]={0}; 
      double rec_pt[100]={0}, rec_ptError[100]={0}; 

      // Histogram 2D Arrays
      TH1D *diMuonsInvMass_GenA[10][1000];
      TH1D *diMuonsInvMass_RecA[10][1000];
      TH1D *diMuonsPt_GenA[10][1000];
      TH1D *diMuonsPt_RecA[10][1000];
      char nameGen[10][500], nameRec[10][500], nameGenPt[10][500], nameRecPt[10][500];
      char namePt_1B[500];//for bkg func
      for (int ifile = 0; ifile <= 6; ifile++) {
        for (Int_t idphi = 0; idphi < ndPhiBins; idphi++) {
          sprintf(nameGen[ifile],"DiMuonMassGen_pt_%d_%d_%d",idphi,ifile,iCat);
          sprintf(nameRec[ifile],"DiMuonMassRec_pt_%d_%d_%d",idphi,ifile,iCat);

          sprintf(nameGenPt[ifile],"DiMuonPtGen_pt_%d_%d_%d",idphi,ifile,iCat);
          sprintf(nameRecPt[ifile],"DiMuonPtRec_pt_%d_%d_%d",idphi,ifile,iCat);

          diMuonsInvMass_GenA[ifile][idphi]= new TH1D(nameGen[ifile],nameGen[ifile],  100,2.95,3.25); //for eff Gen;
          diMuonsInvMass_GenA[ifile][idphi]->Sumw2();
          diMuonsInvMass_GenA[ifile][idphi]->SetMarkerStyle(7);
          diMuonsInvMass_GenA[ifile][idphi]->SetMarkerColor(4);
          diMuonsInvMass_GenA[ifile][idphi]->SetLineColor(4);

          diMuonsInvMass_RecA[ifile][idphi] = new TH1D(nameRec[ifile],nameRec[ifile], 100,2.95,3.25); //for eff Rec;
          diMuonsInvMass_RecA[ifile][idphi]->Sumw2();
          diMuonsInvMass_RecA[ifile][idphi]->SetMarkerStyle(8);
          diMuonsInvMass_RecA[ifile][idphi]->SetMarkerColor(4);
          diMuonsInvMass_RecA[ifile][idphi]->SetLineColor(4);

          diMuonsPt_GenA[ifile][idphi]= new TH1D(nameGenPt[ifile],nameGenPt[ifile],  100,0,40); //for eff Gen;
          diMuonsPt_RecA[ifile][idphi]= new TH1D(nameRecPt[ifile],nameRecPt[ifile],  100,0,40); //for eff Rec;
        }
      }

      //===========================================Input Root File============================================================
      char fileName[10][500];
      //scales for different pT bins
      double scale[10]={0};  

      if(Prompt==1){
        cout<<" prompt weight "<<endl;
        scale[0]=2.35829e-07;
        scale[1]=1.99854e-07;
        scale[2]=4.48263e-08;
        scale[3]=1.01144e-08;
        scale[4]=4.89604e-09;
        scale[5]=2.62102e-09;


      }

      if(PutWeight==0){scale[0]=(1);scale[1]=(1);scale[2]=(1);scale[3]=(1);scale[4]=(1);scale[5]=(1);}

      if(Prompt ==1){
        cout<<"==================Prompt PrJpsi================================================"<<endl;
        sprintf(fileName[0],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0003_total.root");
        sprintf(fileName[1],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0306_total.root");
        sprintf(fileName[2],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0609_total.root");
        sprintf(fileName[3],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0912_total.root");
        sprintf(fileName[4],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1215_total.root");
        sprintf(fileName[5],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1530_total.root");
      }

      TFile *infile;
      TTree *tree;
      TTree *gentree;
      //======================  File loop Starts ============================
      for(int ifile = 0; ifile < nFiles; ifile++){
        infile=new TFile(fileName[ifile],"R");
        tree=(TTree*)infile->Get("SingleMuonTree");
        gentree=(TTree*)infile->Get("SingleGenMuonTree");
        //Event variables
        int eventNb,runNb,lumiBlock, gbin, rbin;
        //Jpsi Variables
        Double_t JpsiMass,JpsiPt,JpsiRap, JpsiCharge;
        Double_t JpsiVprob;
        Double_t JpsiPhi;
        Double_t JpsiEta;
        Double_t JpsiPsi[38];
        Double_t JpsiGenPsi;
        //2.) muon variables RECO                                                                       
        double muPosPx, muPosPy, muPosPz,  muPosEta, muPosPt, muPosP, muPosPhi;
        double muNegPx, muNegPy, muNegPz,  muNegEta, muNegPt, muNegP, muNegPhi;
        //(1).Positive Muon                                     
        double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
        int muPos_found, muPos_pixeLayers, muPos_nValidMuHits,muPos_arbitrated;
        bool muPos_matches,muPos_tracker;
        //(2).Negative Muon                                     
        double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
        int muNeg_found, muNeg_pixeLayers, muNeg_nValidMuHits,muNeg_arbitrated;
        bool muNeg_matches,muNeg_tracker;

        //Gen Level variables
        //Gen PrJpsi Variables
        double GenJpsiMass, GenJpsiPt, GenJpsiRap;
        double GenJpsiPx, GenJpsiPy, GenJpsiPz;
        double GenJpsiPhi;
        double GenJpsiEta;
        double GenJpsiPsi;
        //2.) Gen muon variables 
        double GenmuPosPx, GenmuPosPy, GenmuPosPz,  GenmuPosEta, GenmuPosPt, GenmuPosPhi;
        double GenmuNegPx, GenmuNegPy, GenmuNegPz,  GenmuNegEta, GenmuNegPt, GenmuNegPhi;

        // HLT_HIL1DoubleMu0_NHitQ
        int hbit1;

        //Event variables
        tree->SetBranchAddress("eventNb",&eventNb);
        tree->SetBranchAddress("runNb",&runNb);
        tree->SetBranchAddress("lumiBlock",&lumiBlock);
        tree->SetBranchAddress("hbit1",&hbit1);

        //Jpsi Variables
        tree->SetBranchAddress("JpsiCharge",&JpsiCharge);
        tree->SetBranchAddress("JpsiMass",&JpsiMass);
        tree->SetBranchAddress("JpsiPt",&JpsiPt);
        tree->SetBranchAddress("JpsiPhi",&JpsiPhi);
        tree->SetBranchAddress("JpsiEta",&JpsiEta);
        tree->SetBranchAddress("JpsiPsi",&JpsiPsi);
        tree->SetBranchAddress("JpsiGenPsi",&JpsiGenPsi);
        tree->SetBranchAddress("JpsiRap",&JpsiRap);
        tree->SetBranchAddress("JpsiVprob",&JpsiVprob);
        tree->SetBranchAddress("rbin",&rbin);

        //muon variable
        tree->SetBranchAddress("muPosPx",&muPosPx);
        tree->SetBranchAddress("muPosPy",&muPosPy);
        tree->SetBranchAddress("muPosPz",&muPosPz);
        tree->SetBranchAddress("muPosEta",&muPosEta);
        tree->SetBranchAddress("muPosPhi",&muPosPhi);
        tree->SetBranchAddress("muNegPx", &muNegPx);
        tree->SetBranchAddress("muNegPy", &muNegPy);
        tree->SetBranchAddress("muNegPz", &muNegPz);
        tree->SetBranchAddress("muNegEta", &muNegEta);
        tree->SetBranchAddress("muNegPhi", &muNegPhi);


        //1). Positive Muon
        tree->SetBranchAddress("muPos_nchi2In", &muPos_nchi2In);
        tree->SetBranchAddress("muPos_dxy", &muPos_dxy);
        tree->SetBranchAddress("muPos_dz", &muPos_dz);
        tree->SetBranchAddress("muPos_nchi2Gl", &muPos_nchi2Gl);
        tree->SetBranchAddress("muPos_found", &muPos_found);
        tree->SetBranchAddress("muPos_pixeLayers", &muPos_pixeLayers);
        tree->SetBranchAddress("muPos_nValidMuHits", &muPos_nValidMuHits);
        tree->SetBranchAddress("muPos_matches", &muPos_matches);
        tree->SetBranchAddress("muPos_tracker", &muPos_tracker);
        tree->SetBranchAddress("muPos_arbitrated", &muPos_arbitrated);

        //2). Negative Muon                                                                            
        tree->SetBranchAddress("muNeg_nchi2In", &muNeg_nchi2In);
        tree->SetBranchAddress("muNeg_dxy", &muNeg_dxy);
        tree->SetBranchAddress("muNeg_dz", &muNeg_dz);
        tree->SetBranchAddress("muNeg_nchi2Gl", &muNeg_nchi2Gl);
        tree->SetBranchAddress("muNeg_found", &muNeg_found);
        tree->SetBranchAddress("muNeg_pixeLayers", &muNeg_pixeLayers);
        tree->SetBranchAddress("muNeg_nValidMuHits", &muNeg_nValidMuHits);
        tree->SetBranchAddress("muNeg_matches", &muNeg_matches);
        tree->SetBranchAddress("muNeg_tracker", &muNeg_tracker);
        tree->SetBranchAddress("muNeg_arbitrated", &muNeg_arbitrated);
        //====================================Gen Variables=========================================================
        //Gen Jpsi Variables
        gentree->SetBranchAddress("GenJpsiMass",   &GenJpsiMass);
        gentree->SetBranchAddress("GenJpsiPt",     &GenJpsiPt);
        gentree->SetBranchAddress("GenJpsiPhi",    &GenJpsiPhi);
        gentree->SetBranchAddress("GenJpsiPsi",    &GenJpsiPsi);
        gentree->SetBranchAddress("GenJpsiRap",    &GenJpsiRap);
        gentree->SetBranchAddress("GenJpsiEta",    &GenJpsiEta);
        gentree->SetBranchAddress("GenJpsiPx",     &GenJpsiPx);
        gentree->SetBranchAddress("GenJpsiPy",     &GenJpsiPy);
        gentree->SetBranchAddress("GenJpsiPz",     &GenJpsiPz);
        gentree->SetBranchAddress("gbin",&gbin);
        //muon variable
        gentree->SetBranchAddress("GenmuPosPx",    &GenmuPosPx);
        gentree->SetBranchAddress("GenmuPosPy",    &GenmuPosPy);
        gentree->SetBranchAddress("GenmuPosPz",    &GenmuPosPz);
        gentree->SetBranchAddress("GenmuPosEta",    &GenmuPosEta);
        gentree->SetBranchAddress("GenmuPosPhi",    &GenmuPosPhi);
        gentree->SetBranchAddress("GenmuNegPx",    &GenmuNegPx);
        gentree->SetBranchAddress("GenmuNegPy",    &GenmuNegPy);
        gentree->SetBranchAddress("GenmuNegPz",    &GenmuNegPz);
        gentree->SetBranchAddress("GenmuNegEta",    &GenmuNegEta);
        gentree->SetBranchAddress("GenmuNegPhi",    &GenmuNegPhi);

        //====================================================== Gen tree loop ================================================
        int NAccep=0;
        int nGenEntries=gentree->GetEntries();
        cout<<" Total Entries in GenLevel Tree for pT range: "<<fileName[ifile]<<"  "<<   nGenEntries<< " ==============="<<endl;

        for(int i=0; i< nGenEntries; i++)  {        
          gentree->GetEntry(i);
          //Only printing 
          if(i%100000==0){
            cout<<" processing record "<<i<<"/"<<nGenEntries<<endl;
            //cout<<" Mass "<< GenJpsiMass<< " pT "<< GenJpsiPt << " Y " <<GenJpsiRap<<endl;
          }

          bool GenPosIn=0, GenNegIn=0;
          GenmuPosPt= TMath::Sqrt(GenmuPosPx*GenmuPosPx + GenmuPosPy*GenmuPosPy); 
          GenmuNegPt= TMath::Sqrt(GenmuNegPx*GenmuNegPx + GenmuNegPy*GenmuNegPy); 

          diMuonsInvMass_Gen->Fill(GenJpsiMass);
          diMuonsPt_Gen->Fill(GenJpsiPt);

          if(IsAccept(GenmuPosPt, GenmuPosEta)) {GenPosIn=1;}
          if(IsAccept(GenmuNegPt, GenmuNegEta)) {GenNegIn=1;}

          if((GenPosIn ==1 && GenNegIn==1)) NAccep++;

          Bin_Gen->Fill(gbin);

          int AccJpsi = 0;
          if(GenJpsiPt < 6.5) continue;
          if(GenJpsiPt >= 6.5 && fabs(GenJpsiRap) < 2.4 && (GenPosIn ==1 && GenNegIn==1)) AccJpsi = 1;

          double gdPhi2mu = GenmuPosPhi - GenmuNegPhi;
          while (gdPhi2mu > TMath::Pi()) gdPhi2mu -= 2*TMath::Pi();
          while (gdPhi2mu <= -TMath::Pi()) gdPhi2mu += 2*TMath::Pi();
          double gchkCowboy = 1*gdPhi2mu;

          if(bCowboy && !(gchkCowboy > 0.)) continue;
          if(bSailor && (gchkCowboy > 0.)) continue;

          double GenCenWeight =0, GenWeight =0;
          GenCenWeight=FindCenWeight(gbin);   

          double gJpsidPhi = TMath::Abs(GenJpsiPsi);
          double gmean_dphi = TMath::Pi()/16;
          if (gJpsidPhi>TMath::Pi()/8) gmean_dphi = 3*TMath::Pi()/16;
          if (gJpsidPhi>TMath::Pi()/4) gmean_dphi = 5*TMath::Pi()/16;
          if (gJpsidPhi>3*TMath::Pi()/8) gmean_dphi = 7*TMath::Pi()/16;
          GenWeight=GenCenWeight*scale[ifile]*(2.0/TMath::Pi()*(1+2*fake_v2*cos(2*gJpsidPhi)*TMath::Gaus(GenJpsiPt,15,4,0)));

          if(PutWeight==0) GenWeight=1; 

          for (Int_t idphi = 0; idphi < ndPhiBins; idphi++) {
              if(iCat == 0)  {
                if( AccJpsi == 1 && (gbin >= 0 && gbin < 4) && (GenPosIn==1 && GenNegIn==1) && 
                  (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) && (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 0.0 ) && 
                  (TMath::Abs(GenJpsiPsi)>=dphi_bound[idphi] && TMath::Abs(GenJpsiPsi)<dphi_bound[idphi+1])){diMuonsInvMass_GenA[ifile][idphi]->Fill(GenJpsiMass,GenWeight);}
              }
              if(iCat == 1)  {
                if( AccJpsi == 1 && (gbin >= 4 && gbin < 8) && (GenPosIn==1 && GenNegIn==1) && 
                  (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) && (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 0.0 ) && 
                  (TMath::Abs(GenJpsiPsi)>=dphi_bound[idphi] && TMath::Abs(GenJpsiPsi)<dphi_bound[idphi+1])){diMuonsInvMass_GenA[ifile][idphi]->Fill(GenJpsiMass,GenWeight);}
              }
              if(iCat == 2)  {
                if( AccJpsi == 1 && (gbin >= 8 && gbin < 12) && (GenPosIn==1 && GenNegIn==1) && 
                  (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) && (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 0.0 ) && 
                  (TMath::Abs(GenJpsiPsi)>=dphi_bound[idphi] && TMath::Abs(GenJpsiPsi)<dphi_bound[idphi+1])){diMuonsInvMass_GenA[ifile][idphi]->Fill(GenJpsiMass,GenWeight);}
              }
              if(iCat == 3)  {
                if( AccJpsi == 1 && (gbin >= 12 && gbin < 24) && (GenPosIn==1 && GenNegIn==1) && 
                  (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) && (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 0.0 ) && 
                  (TMath::Abs(GenJpsiPsi)>=dphi_bound[idphi] && TMath::Abs(GenJpsiPsi)<dphi_bound[idphi+1])){diMuonsInvMass_GenA[ifile][idphi]->Fill(GenJpsiMass,GenWeight);}
              }
          }
        }//gen loop end

        cout<<" accepted no "<< NAccep<<endl;
        dataFile<<" accepted no "<< NAccep<<endl;

        //dataFile<<" accepted no "<< NAccep<<endl;
        //   new TCanvas;
        //diMuonsInvMass_Gen->Draw();
        //gPad->Print("plots/diMuonsInvMass_Gen.png");

        new TCanvas;
        diMuonsPt_Gen->Draw();

        //=============== Rec Tree Loop ==============================================================================

        // start to fill up reco corrected by efficiency
        
        char tmp_eff_input[512], tmp_input_histo[512];
        sprintf(tmp_eff_input,"../../EffRoots_New/PrJpsi_HighPt_%s.root", cCond);

        cout<<"tmp_eff_input : "<<tmp_eff_input<<endl;
        TFile *eff_input;
        eff_input=new TFile(tmp_eff_input,"R");
        
        int nRecEntries=tree->GetEntries();
        cout<<" Total Entries in reconstructed Tree for pT range "<<fileName[ifile]<<"  "<<nRecEntries<< "====="<<endl;
        for(int i=0; i<nRecEntries; i++)  {     
          tree->GetEntry(i);
          //Only printing 
          if(i%10000==0){
            cout<<" processing record "<<i<<"/"<<nRecEntries<<endl;
            //cout<<" processing Run  " <<runNb <<" event "<<eventNb<<" lum block "<<lumiBlock<<endl;    
            //cout<<" Mass "<< JpsiMass<< " pT "<< JpsiPt << " Y " <<JpsiRap<<"  "<<JpsiVprob<<" charge "<<JpsiCharge<<" rbin "<<rbin<<endl; 
          }
          bool PosPass=0, NegPass=0, AllCut=0 ,PosIn=0, NegIn=0;
          muPosPt= TMath::Sqrt(muPosPx*muPosPx + muPosPy*muPosPy); 
          muPosP = TMath::Sqrt(muPosPx*muPosPx + muPosPy*muPosPy+ muPosPz*muPosPz); 
          muNegPt= TMath::Sqrt(muNegPx*muNegPx + muNegPy*muNegPy); 
          muNegP = TMath::Sqrt(muNegPx*muNegPx + muNegPy*muNegPy +muNegPz*muNegPz); 

          if(IsAccept(muPosPt, muPosEta)){PosIn=1;}
          if(IsAccept(muNegPt, muNegEta)){NegIn=1;}

          int AccJpsi = 0;
          if(JpsiPt < 6.5) continue;
          if(JpsiPt >= 6.5 && fabs(JpsiRap) < 2.4 && (PosIn == 1 && NegIn == 1)) AccJpsi = 1;

          bool mu_Global = ((muPos_nchi2Gl >=0) && (muNeg_nchi2Gl >=0));
          bool mu_Tracker = ((muPos_tracker==1) && (muNeg_tracker==1));

          if(muPos_found > 10 && muPos_pixeLayers > 0 && muPos_nchi2In < 4.0 && TMath::Abs(muPos_dxy) < 3 && TMath::Abs(muPos_dz) < 15 && muPos_nchi2Gl < 20 && muPos_arbitrated==1 && muPos_tracker==1){PosPass=1;}     

          if(muNeg_found > 10 && muNeg_pixeLayers > 0 && muNeg_nchi2In < 4.0 && TMath::Abs(muNeg_dxy) < 3 && TMath::Abs(muNeg_dz) < 15 && muNeg_nchi2Gl < 20 && muNeg_arbitrated==1 && muNeg_tracker==1){NegPass=1;}


          double dPhi2mu = muPosPhi - muNegPhi;
          while (dPhi2mu > TMath::Pi()) dPhi2mu -= 2*TMath::Pi();
          while (dPhi2mu <= -TMath::Pi()) dPhi2mu += 2*TMath::Pi();
          double chkCowboy = 1*dPhi2mu;

          if(bCowboy && !(chkCowboy > 0.)) continue;
          if(bSailor && (chkCowboy > 0.)) continue;

          //cout<<"Cut checks, muPos_matches : "<<muPos_matches<<", muNeg_matches : "<<muNeg_matches<<", PosIn : "<<PosIn<<", NegIn : "<<NegIn
          //    <<", PosPass : "<<PosPass<<", NegPass : "<<NegPass<<", mu_Global : "<<mu_Global<<", mu_Tracker : "<<mu_Tracker<<endl;
          //if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}
          if(hbit1 == 1 && (muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

          //AllCut = 1;
          //without ID cut
          // if((muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && mu_Global && mu_Tracker){AllCut=1;}

          //without trigger matched
          //if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

          
          double eff_cor[2];
          double RecCenWeight=0,RecWeight=0;
          RecCenWeight=FindCenWeight(rbin);   
          
          double JpsidPhi = 0.0;
          JpsidPhi = TMath::Abs(JpsiGenPsi);

          double rmean_dphi = TMath::Pi()/16;
          if (JpsidPhi>TMath::Pi()/8) rmean_dphi = 3*TMath::Pi()/16;
          if (JpsidPhi>TMath::Pi()/4) rmean_dphi = 5*TMath::Pi()/16;
          if (JpsidPhi>3*TMath::Pi()/8) rmean_dphi = 7*TMath::Pi()/16;
          RecWeight=RecCenWeight*scale[ifile]*(2.0/TMath::Pi()*(1+2*fake_v2*cos(2*JpsidPhi)*TMath::Gaus(JpsiPt,15,4,0)));

          if(PutWeight==0)RecWeight=1;

          //Eff loop for reco
          for (Int_t idphi = 0; idphi < ndPhiBins; idphi++) {
            if((JpsiCharge == 0) && (JpsiVprob > 0.01)) {     
              if(iCat == 0) {
                if( AccJpsi == 1 && (AllCut==1) && (rbin >= 0 && rbin < 4) &&  
                    (JpsiPt>=6.5 && JpsiPt<40.0) && (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 0.0) && 
                    (TMath::Abs(JpsiGenPsi) > dphi_bound[idphi] && TMath::Abs(JpsiGenPsi) <=dphi_bound[idphi+1])){
                  DoEffCor3D(eff_input, iCond, rbin, JpsiPt, fabs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][idphi]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));
                }
              }
              if(iCat == 1) {
                if( AccJpsi == 1 && (AllCut==1) && (rbin >= 4 && rbin < 8) && 
                    (JpsiPt>=6.5 && JpsiPt<40.0) && (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 0.0) && 
                    (TMath::Abs(JpsiGenPsi) > dphi_bound[idphi] && TMath::Abs(JpsiGenPsi) <=dphi_bound[idphi+1])){
                  DoEffCor3D(eff_input, iCond, rbin, JpsiPt, fabs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][idphi]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
              }
              if(iCat == 2) {
                if(AccJpsi == 1 && (AllCut==1) && (rbin >= 8 && rbin < 12) && 
                  (JpsiPt>=6.5 && JpsiPt<40.0) && (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 0.0) && 
                  (TMath::Abs(JpsiGenPsi) > dphi_bound[idphi] && TMath::Abs(JpsiGenPsi) <=dphi_bound[idphi+1])){
                  DoEffCor3D(eff_input, iCond, rbin, JpsiPt, fabs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][idphi]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
              }
              if(iCat == 3) {
                if(AccJpsi == 1 && (AllCut==1) && (rbin >= 12 && rbin < 24) &&
                  (JpsiPt>=6.5 && JpsiPt<40.0) && (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 0.0) && 
                  (TMath::Abs(JpsiGenPsi) > dphi_bound[idphi] && TMath::Abs(JpsiGenPsi) <=dphi_bound[idphi+1])){
                  DoEffCor3D(eff_input, iCond, rbin, JpsiPt, fabs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][idphi]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
              }
            }
          }
        }//rec tree loop ends
      }  // file loop ends

      ///////////////////////////////////////////////////////////////////
      cout<< " adding "<<endl;
      TH1D *diMuonsInvMass_RecA1[100];
      TH1D *diMuonsInvMass_GenA1[100];
      TH1D *diMuonsPt_GenA1[100];
      TH1D *diMuonsPt_RecA1[100];
      TF1 *backfun_1;

      for(Int_t idphi = 0; idphi < ndPhiBins; idphi++){
        diMuonsInvMass_RecA1[idphi] = diMuonsInvMass_RecA[0][idphi];
        diMuonsInvMass_GenA1[idphi] = diMuonsInvMass_GenA[0][idphi];
        diMuonsPt_GenA1[idphi] = diMuonsPt_GenA[0][idphi];
        diMuonsPt_RecA1[idphi] = diMuonsPt_RecA[0][idphi];

        for (int ifile = 1; ifile < nFiles; ifile++) {
          diMuonsInvMass_RecA1[idphi]->Add(diMuonsInvMass_RecA[ifile][idphi]);
          diMuonsInvMass_GenA1[idphi]->Add(diMuonsInvMass_GenA[ifile][idphi]);     
          diMuonsPt_GenA1[idphi]->Add(diMuonsPt_GenA[ifile][idphi]); 
          diMuonsPt_RecA1[idphi]->Add(diMuonsPt_RecA[ifile][idphi]); 
        }
      }
      //===========================Fitting=================================================================================//
      // Fit ranges
      double mass_low, mass_high;
      double MassJpsi, WeidthJpsi;
      // Fit Function crystall ball
      // Jpsi Settings
      MassJpsi = 3.096; 
      WeidthJpsi = 0.028;
      mass_low = 2.945; 
      mass_high = 3.24;  // Fit ranges

      TF1 *GAUSPOL = new TF1("GAUSPOL",CrystalBall,2.4,3.8,6);//2.4,3.8,6);
      GAUSPOL->SetParNames("Yield (J/#psi)","BinWidth","Mean","Sigma","#alpha","n");
      GAUSPOL->SetParameter(2, MassJpsi);
      GAUSPOL->SetParameter(3, WeidthJpsi);
      //GAUSPOL->SetParLimits(3, 0.1*WeidthJpsi,2.0*WeidthJpsi);
      GAUSPOL->SetParameter(4, 1.2);
      GAUSPOL->SetParameter(5, 20.0);
      GAUSPOL->SetLineWidth(2.0);
      GAUSPOL->SetLineColor(2);


      //=====================Loop for eff========================================================================================//
      //define stuff here for error on weighted samples
      double GenNo[100]={0};
      double Eff[100]={0};
      double GenError[100]={0};
      double RecError[100]={0};

      for (Int_t idphi = 0; idphi < ndPhiBins; idphi++) {

        gen_pt[idphi] = diMuonsInvMass_GenA1[idphi]->IntegralAndError(1, 100, genError);
        gen_ptError[idphi]= genError;
        cout<<" gen_pt["<<idphi<<"] "<< gen_pt[idphi] <<" error   "<<  gen_ptError[idphi]<<endl;

        // cout<<" *********************** "<<diMuonsInvMass_RecA1[idphi]->GetMaximum()<<endl;
        //giving inetial value for crystall ball fourth parameter 
        GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[idphi]->Integral(0,50));
        GAUSPOL->FixParameter(1, diMuonsInvMass_RecA1[idphi]->GetBinWidth(1));
        //GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[idphi]->GetMaximum());
        //new TCanvas;
        //diMuonsInvMass_RecA1[idphi]->Draw();

        new TCanvas;
        diMuonsInvMass_RecA1[idphi]->Fit("GAUSPOL","EMRQ", "", mass_low, mass_high); // Jpsi
        //diMuonsInvMass_RecA1[idphi]->Fit("GAUSPOL","LLMERQ", "", 8.5,10.5); // Jpsi
        //diMuonsInvMass_RecA1[idphi]->Fit("GAUSPOL","LLMER", "", mass_low, mass_high); // Jpsi
        diMuonsInvMass_RecA1[idphi]->DrawCopy("EPLsame");

        double JpsiMass = GAUSPOL->GetParameter(2);
        double JpsiWidth = GAUSPOL->GetParameter(3);

        double JpsiYield = GAUSPOL->GetParameter(0); 
        Double_t JpsiYieldError = GAUSPOL->GetParError(0); 

        //cout<<JpsiYieldError<<"*****************"<<endl;

        //if(TMath::IsNan(JpsiYieldError)=1) {JpsiYieldError=TMath::Sqrt(JpsiYield);}

        double par[20];
        GAUSPOL->GetParameters(par);
        sprintf(namePt_1B,"pt_1B_%d",idphi);

        backfun_1 = new TF1(namePt_1B, Pol2, mass_low, mass_high, 3);
        backfun_1->SetParameters(&par[3]);

        double MassLow=(JpsiMass-3*JpsiWidth);
        double MassHigh=(JpsiMass+3*JpsiWidth);


        int binlow =diMuonsInvMass_RecA1[idphi]->GetXaxis()->FindBin(MassLow);
        int binhi =diMuonsInvMass_RecA1[idphi]->GetXaxis()->FindBin(MassHigh);

        //double binwidth=diMuonsInvMass_RecA1[idphi]->GetBinWidth(1);
        //yield by function 
        //rec_pt[idphi] = JpsiYield;
        //rec_ptError[idphi]= JpsiYieldError;

        //yield by histogram integral
        binlow = 1;
        binhi = 100;
        /*
        for(int i = 0; i < diMuonsInvMass_RecA1[idphi]->GetEntries(); i++){
          cout<<" contents : "<<diMuonsInvMass_RecA1[idphi]->GetBinContent(i)<<endl;
        }
        cout<<" 1 : "<<diMuonsInvMass_RecA1[idphi]->GetEntries()<<endl;
        */
        rec_pt[idphi] = diMuonsInvMass_RecA1[idphi]->IntegralAndError(1, 100,recError);
        rec_ptError[idphi]= recError;
        cout<<" rec_pt["<<idphi<<"] "<< rec_pt[idphi] <<" error   "<<  rec_ptError[idphi]<<endl;

      }

      TFile *outfile;
      char tmp_output[512];
      if(iCat == 0) sprintf(tmp_output,"MCCT_PrJpsi_Cent_0010_dPhi_%s.root",cCond);
      if(iCat == 1) sprintf(tmp_output,"MCCT_PrJpsi_Cent_1020_dPhi_%s.root",cCond);
      if(iCat == 2) sprintf(tmp_output,"MCCT_PrJpsi_Cent_2030_dPhi_%s.root",cCond);
      if(iCat == 3) sprintf(tmp_output,"MCCT_PrJpsi_Cent_3060_dPhi_%s.root",cCond);
      outfile =new TFile(tmp_output,"Recreate");

      double gsum = 0.0;
      double gen_pt_Norm[ndPhiBins];
      double gen_ptError_Norm[ndPhiBins];
      for(int i = 0; i < ndPhiBins; i++){
        gsum += gen_pt[i];
      }
      gsum = gsum*(dphi_bound[1]-dphi_bound[0]);

      for(int i = 0; i < ndPhiBins; i++){
        gen_pt_Norm[i] = gen_pt[i]/gsum;
        gen_ptError_Norm[i] = gen_ptError[i]/gsum;
      }


      double sum = 0.0;
      double rec_pt_Norm[ndPhiBins];
      double rec_ptError_Norm[ndPhiBins];
      for(int i = 0; i < ndPhiBins; i++){
        sum += rec_pt[i];
      }
      sum = sum*(dphi_bound[1]-dphi_bound[0]);

      for(int i = 0; i < ndPhiBins; i++){
        rec_pt_Norm[i] = rec_pt[i]/sum;
        rec_ptError_Norm[i] = rec_ptError[i]/sum;
      }


      double mom_err[ndPhiBins] = {0.0};
      TGraphErrors *Jpsi_Reco_Norm = new TGraphErrors(ndPhiBins, xdphi_bound, rec_pt_Norm, mom_err, rec_ptError_Norm);
      TGraphErrors *Jpsi_Gen_Norm = new TGraphErrors(ndPhiBins, xdphi_bound, gen_pt_Norm, mom_err, gen_ptError_Norm);
      TH1F *hJpsi_Gen = new TH1F("hJpsi_Gen","hJpsi_Gen",ndPhiBins,0,TMath::Pi()/2);
      TH1F *hJpsi_Reco = new TH1F("hJpsi_Reco","hJpsi_Reco",ndPhiBins,0,TMath::Pi()/2);
      hJpsi_Gen->Sumw2();
      hJpsi_Reco->Sumw2();
      for(int idphi = 0; idphi < ndPhiBins; idphi++){
        hJpsi_Gen->SetBinContent(idphi+1,gen_pt[idphi]);
        hJpsi_Gen->SetBinError(idphi+1,gen_ptError[idphi]);
        hJpsi_Reco->SetBinContent(idphi+1,rec_pt[idphi]);
        hJpsi_Reco->SetBinError(idphi+1,rec_ptError[idphi]);
      }

      if(iCat == 0){
        hJpsi_Gen->SetName("hGen_Jpsi_Cent_0010");
        hJpsi_Reco->SetName("hReco_Jpsi_Cent_0010");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Cent_0010");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Cent_0010");
      }
      if(iCat == 1){
        hJpsi_Gen->SetName("hGen_Jpsi_Cent_1020");
        hJpsi_Reco->SetName("hReco_Jpsi_Cent_1020");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Cent_1020");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Cent_1020");
      }
      if(iCat == 2){
        hJpsi_Gen->SetName("hGen_Jpsi_Cent_2030");
        hJpsi_Reco->SetName("hReco_Jpsi_Cent_2030");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Cent_2030");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Cent_2030");
      }
      if(iCat == 3){
        hJpsi_Gen->SetName("hGen_Jpsi_Cent_3060");
        hJpsi_Reco->SetName("hReco_Jpsi_Cent_3060");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Cent_3060");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Cent_3060");
      }
      hJpsi_Gen->Write();
      hJpsi_Reco->Write();
      Jpsi_Reco_Norm->Write();
      Jpsi_Gen_Norm->Write();

      outfile->Write();
      outfile->Close();
    }
  }
}


bool IsAccept(Double_t pt, Double_t eta)
{

  //return(1);
  //return (fabs(eta) < 2.4); 

  return (fabs(eta) < 2.4 &&
      (    ( fabs(eta) < 1.0 && pt >= 3.4 ) ||
           (  1.0 <= fabs(eta) && fabs(eta) < 1.5 && pt >= 5.8-2.4*fabs(eta) ) ||
           (  1.5 <= fabs(eta) && pt >= 3.3667 - 7.0/9.0*fabs(eta)) ));
}


double FindCenWeight(int Bin)
{
  double NCollArray[40]={
    1747.8600, 1567.5300, 1388.3900, 1231.7700, 1098.2000, 980.4390, 861.6090, 766.0420, 676.5150, 593.4730,
    521.9120, 456.5420, 398.5460, 346.6470, 299.3050, 258.3440, 221.2160, 188.6770, 158.9860, 134.7000,
    112.5470, 93.4537, 77.9314, 63.5031, 52.0469, 42.3542, 33.9204, 27.3163, 21.8028, 17.2037,
    13.5881, 10.6538, 8.3555, 6.4089, 5.1334, 3.7322, 3.0663, 2.4193, 2.1190, 1.7695
  };
  return(NCollArray[Bin]);
}

double FindV2Weight(double v2wg)
{
  int v2wgt = 0;
  if(v2wg > 0 && v2wg < TMath::Pi()/8) v2wgt = 0;
  if(v2wg > TMath::Pi()/8 && v2wg < 2*TMath::Pi()/8) v2wgt = 1;
  if(v2wg > 2*TMath::Pi()/8 && v2wg < 3*TMath::Pi()/8) v2wgt = 2;
  if(v2wg > 3*TMath::Pi()/8 && v2wg < 4*TMath::Pi()/8) v2wgt = 3;
  double v2Array[4] = {
    1.3631, 1.1486, 0.8429, 0.6274
  };
  return(v2Array[v2wgt]);
}

void DoEffCor1D(TFile *a, int b, double i, double *c){

  /*
  const char *cSp[3] = {"Cent","Pt","Rap"};
  char tmp[512];
  if(a == 0) sprintf(tmp,"../EffRoots/PrJpsi_HighPt_%s_default_Bit1.root",cSp[b]);
  if(a == 1) sprintf(tmp,"../EffRoots/PrJpsi_HighPt_%s_cowboy_Bit1.root",cSp[b]);
  if(a == 2) sprintf(tmp,"../EffRoots/PrJpsi_HighPt_%s_sailor_Bit1.root",cSp[b]);
  //cout<<"file name : "<<tmp<<endl;
  //TFile *f = new TFile(tmp,"READ");
  */
  TH1F *heff;
  if(b==0) heff = (TH1F*)a->Get("hEff_Cents");
  if(b==1) heff = (TH1F*)a->Get("hEff_Pts");
  if(b==2) heff = (TH1F*)a->Get("hEff_Raps");
  if(!heff) cout<<"this is empty !!!"<<endl;

  if(i > 30) i = 29.9;
  int sbin = heff->FindBin(i);
  //cout<<" x value : "<<i<<", sbin : "<<sbin<<endl;
  c[0] = heff->GetBinContent(sbin);
  c[1] = heff->GetBinError(sbin);
  //cout<<" eff : "<<c[0]<<endl;
  if(c[0] == 0) {
    cout<<" eff : 0 "<<endl;
    cout<<" pt : "<<i<<endl;
    c[0] = 1.0;
  }
  //c[0] = 0.3;

}

void DoEffCor3D(TFile* a, int b, double c, double d, double e, double *f){
  char tmp[512];
  if(b == 0) sprintf(tmp,"eff_default");
  if(b == 1) sprintf(tmp,"eff_cowboy");
  if(b == 2) sprintf(tmp,"eff_sailor");
  TH3D *hEff = (TH3D*)a->Get(tmp);
  int sbin = 0;
  sbin = hEff->FindBin(c, d, e);
  f[0] = hEff->GetBinContent(sbin);
  f[1] = hEff->GetBinError(sbin);
  if(f[0] == 0) {
    cout<<"b : "<<b<<", c : "<<c<<", d : "<<d<<", e : "<<e<<endl;
    cout<<"sbin : "<<sbin<<", eff : "<<f[0]<<", err : "<<f[1]<<endl;
  }
  //return f;
}

