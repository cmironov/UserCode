//===================== Macro to calculate PrJpsi Eff or Acc x Eff ============================================
//read tree from DiMuonOnia2DPlots.cc
#ifndef __CINT__
#endif
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
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
// Prompt =1 prompt PrJpsi 
bool IsAccept(Double_t pt, Double_t eta); //you can define acceptance here 
double FindCenWeight(int Bin);//gives you weight according to cent
void Eff3DJpsiLowPt()   
{
  int Prompt =1; int PutWeight = 1;

  for(int iCat = 0; iCat < 3; iCat++){
  // How to use 
  // Decide default/cowboy/sailor
  // Example: if you want to run the case of cowboy
  // bool bDefault = false;
  // bool bCowboy = true;
  // bool bSailor = flase;

  bool bDefault = true; // true : default, false : sailor or cowboy
  bool bCowboy  = false; // true : cowboy only, false : salior
  bool bSailor  = false; // true : cowboy only, false : salior

  if(iCat == 0) {bDefault = true;bCowboy = false;bSailor = false;} // default
  if(iCat == 1) {bDefault = false;bCowboy = true;bSailor = false;} // cowboy
  if(iCat == 2) {bDefault = false;bCowboy = false;bSailor = true;} // sailor

  char cCd[512];
  if(bDefault) sprintf(cCd, "default");
  if(bCowboy)  sprintf(cCd, "cowboy");
  if(bSailor)  sprintf(cCd, "sailor");

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
  gStyle->SetOptStat(0); // most of the time, sometimes "nemriou" might be useful to display name, 
  //number of entries, mean, rms, integral, overflow and underflow
  gStyle->SetOptFit(0); // set to 1 only if you want to display fit results
  //==================================== Define Histograms====================================================

  //==============================================Define Acc Eff Stuff here===========================================
  // Pt bin sizes
  // 0-1.5, 1.5-3, 3-4.5, 4.5-6, 6-7.5...

  const int nCentBins = 5;
  const int nPtBins = 3;
  const int nRapBins = 2;
  const int ndPhiBins = 4;
  const int nFiles = 6;
  double ct_bound[nCentBins+1] = {0, 4, 8, 12, 24, 40};
  double xct_bound[nCentBins] = {0.0};
  double pt_bound[nPtBins+1] = {0.0, 3.0, 4.5, 6.5};
  double xpt_bound[nPtBins] = {0.0};
  double rap_bound[nRapBins+1] = {0.0, 1.6, 2.4};
  double xrap_bound[nRapBins] = {0.0};
  double dphi_bound[ndPhiBins+1] = {0.0, TMath::Pi()/8, TMath::Pi()/4, 3*TMath::Pi()/8, TMath::Pi()/2};
  double xdphi_bound[ndPhiBins] = {0.0};
  char OutTextFile[100];
  sprintf(OutTextFile,"eff_%s.tex", cCd);
  ofstream dataFile(Form(OutTextFile));

  char tmp_start[512];
  sprintf(tmp_start,"%%%% Getting Efficiency starts, Category : %s !!!!! %%%%%", cCd);  
  cout<< tmp_start << endl;
  dataFile<< tmp_start << endl;

  // x, y, z - axis 
  dataFile<<""<<endl;
  dataFile<<"xaxis of Cent"<<endl;
  for(int i = 0; i < nCentBins; i++){
    xct_bound[i] = ct_bound[i] + (ct_bound[i+1]-ct_bound[i])/2;
    cout<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
    dataFile<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
  }
  dataFile<<""<<endl;
  dataFile<<"xaxis of pT"<<endl;
  for(int i = 0; i < nPtBins; i++){
    xpt_bound[i] = pt_bound[i] + (pt_bound[i+1]-pt_bound[i])/2;
    cout<<"xpt_bound["<<i<<"] : "<<xpt_bound[i]<<endl;
    dataFile<<"xpt_bound["<<i<<"] : "<<xpt_bound[i]<<endl;
  }
  dataFile<<""<<endl;
  dataFile<<"xaxis of rap"<<endl;
  for(int i = 0; i < nRapBins; i++){
    xrap_bound[i] = rap_bound[i] + (rap_bound[i+1]-rap_bound[i])/2;
    cout<<"xrap_bound["<<i<<"] : "<<xrap_bound[i]<<endl;
    dataFile<<"xrap_bound["<<i<<"] : "<<xrap_bound[i]<<endl;
  }
  dataFile<<""<<endl;
  dataFile<<"xaxis of dphi"<<endl;
  for(int i = 0; i < ndPhiBins; i++){
    xdphi_bound[i] = dphi_bound[i] + (dphi_bound[i+1]-dphi_bound[i])/2;
    cout<<"xdphi_bound["<<i<<"] : "<<xdphi_bound[i]<<endl;
    dataFile<<"xdphi_bound["<<i<<"] : "<<xdphi_bound[i]<<endl;
  }


  TH1F *hTempMass = new TH1F("hTempMass","",100, 2.95, 3.25);
  //TH1F *hTempMass = new TH1F("hTempMass","",100, 2.98, 3.16);
  //TH1F *hTempMass = new TH1F("hTempMass","",100, 2.6, 3.5);
  TH1F *hGenDiMuonf[nFiles][nCentBins][nPtBins][nRapBins];
  TH1F *hRecoDiMuonf[nFiles][nCentBins][nPtBins][nRapBins];
  TH1F *hGenDiMuon[nCentBins][nPtBins][nRapBins];
  TH1F *hRecoDiMuon[nCentBins][nPtBins][nRapBins];
  double genNo[nCentBins][nPtBins][nRapBins];
  double genErr[nCentBins][nPtBins][nRapBins];
  double recoNo[nCentBins][nPtBins][nRapBins];
  double recoErr[nCentBins][nPtBins][nRapBins];
  double eff[nCentBins][nPtBins][nRapBins];
  double effErr[nCentBins][nPtBins][nRapBins];
  for(int fl = 0; fl < nFiles; fl++){
    for(int i = 0; i < nCentBins; i++){
      for(int j = 0; j < nPtBins; j++){
        for(int k = 0; k < nRapBins; k++){
          hGenDiMuonf[fl][i][j][k] = (TH1F*)hTempMass->Clone();
          hRecoDiMuonf[fl][i][j][k] = (TH1F*)hTempMass->Clone();
          hGenDiMuon[i][j][k] = (TH1F*)hTempMass->Clone();
          hRecoDiMuon[i][j][k] = (TH1F*)hTempMass->Clone();
          hGenDiMuonf[fl][i][j][k]->Sumw2();
          hRecoDiMuonf[fl][i][j][k]->Sumw2();
          hGenDiMuon[i][j][k]->Sumw2();
          hRecoDiMuon[i][j][k]->Sumw2();
        }
      }
    }
  }

  char fileName[10][512];
  double scale[6];
  scale[0]=2.35829e-07;
  scale[1]=1.99854e-07;
  scale[2]=4.48263e-08;
  scale[3]=1.01144e-08;
  scale[4]=4.89604e-09;
  scale[5]=2.62102e-09;

  if(PutWeight==0){scale[0]=(1);scale[1]=(1);scale[2]=(1);scale[3]=(1);scale[4]=(1);scale[5]=(1);}

  // loop for pT
  cout<<"==================Prompt PrJpsi================================================"<<endl;
  sprintf(fileName[0],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0003_total.root");
  sprintf(fileName[1],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0306_total.root");
  sprintf(fileName[2],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0609_total.root");
  sprintf(fileName[3],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0912_total.root");
  sprintf(fileName[4],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1215_total.root");
  sprintf(fileName[5],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1530_total.root");  

  TFile *infile;
  TTree *tree;
  TTree *gentree;

  for(int ifile = 0; ifile < nFiles; ifile++){
    infile=new TFile(fileName[ifile],"R");
    tree=(TTree*)infile->Get("SingleMuonTree");
    gentree=(TTree*)infile->Get("SingleGenMuonTree");
    //Event variables
    int eventNb,runNb,lumiBlock, gbin, rbin;
    int hbit1;
    //Jpsi Variables
    Double_t JpsiMass,JpsiPt,JpsiRap, JpsiCharge;
    Double_t JpsiVprob;
    Double_t JpsiPhi;
    Double_t JpsiEta;
    Double_t JpsiPsi[38];
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
    gentree->SetBranchAddress("gbin",          &gbin);
    //muon variable
    gentree->SetBranchAddress("GenmuPosPx",    &GenmuPosPx);
    gentree->SetBranchAddress("GenmuPosPy",    &GenmuPosPy);
    gentree->SetBranchAddress("GenmuPosPz",    &GenmuPosPz);
    gentree->SetBranchAddress("GenmuPosEta",   &GenmuPosEta);
    gentree->SetBranchAddress("GenmuPosPhi",   &GenmuPosPhi);
    gentree->SetBranchAddress("GenmuNegPx",    &GenmuNegPx);
    gentree->SetBranchAddress("GenmuNegPy",    &GenmuNegPy);
    gentree->SetBranchAddress("GenmuNegPz",    &GenmuNegPz);
    gentree->SetBranchAddress("GenmuNegEta",   &GenmuNegEta);
    gentree->SetBranchAddress("GenmuNegPhi",   &GenmuNegPhi);

    //====================================================== Gen tree loop ================================================
    int NAccep=0;
    int nGenEntries=gentree->GetEntries();
    cout<<" Total Entries in GenLevel Tree for pT range: "<<fileName[ifile]<<"  "<<   nGenEntries<< " ==============="<<endl;

    for(int i=0; i< nGenEntries; i++)  {        
      //cout<<"i : "<<i<<endl;
      if(!(gentree->GetEntry(i))) continue;
      //cout<<" gentree ("<<i<<")"<<endl;
      //Only printing 
      if(i%10000==0){
        cout<<" processing record "<<i<<"/"<<nGenEntries<<endl;
        //cout<<" Mass "<< GenJpsiMass<< " pT "<< GenJpsiPt << " Y " <<GenJpsiRap<<endl;
      }

      bool GenPosIn=0, GenNegIn=0;
      GenmuPosPt= TMath::Sqrt(GenmuPosPx*GenmuPosPx + GenmuPosPy*GenmuPosPy); 
      GenmuNegPt= TMath::Sqrt(GenmuNegPx*GenmuNegPx + GenmuNegPy*GenmuNegPy); 

      if(IsAccept(GenmuPosPt, GenmuPosEta)) {GenPosIn=1;}
      if(IsAccept(GenmuNegPt, GenmuNegEta)) {GenNegIn=1;}

      int AccJpsi = 0;
      if(GenJpsiPt >= 6.5) continue;
      if((GenJpsiPt >= 3.0 && GenJpsiPt < 6.5 && fabs(GenJpsiRap) > 1.6 && fabs(GenJpsiRap) <= 2.4 && GenPosIn == 1 && GenNegIn == 1)) {AccJpsi = 1;}

      double gdPhi2mu = GenmuPosPhi - GenmuNegPhi;
      while (gdPhi2mu > TMath::Pi()) gdPhi2mu -= 2*TMath::Pi();
      while (gdPhi2mu <= -TMath::Pi()) gdPhi2mu += 2*TMath::Pi();
      double gchkCowboy = 1*gdPhi2mu;

      if(bCowboy) if(!(gchkCowboy > 0.)) continue;
      if(bSailor) if((gchkCowboy > 0.)) continue;

      if((GenPosIn ==1 && GenNegIn==1)) NAccep++;

      double GenCenWeight =0, GenWeight =0;
      GenCenWeight=FindCenWeight(gbin);   
      GenWeight=GenCenWeight*scale[ifile];

      if(PutWeight==0) GenWeight=1; 

      //adding pT of all pt bins to see diss is cont

      //cout<<"1. GenJpsiPt : "<<GenJpsiPt<<", GenJpsiEta : "<<GenJpsiEta<<", GenJpsiRap : "<<GenJpsiRap<<", |GenJpsiPsi| : "<<TMath::Abs(GenJpsiPsi)<<endl;
      //cout<<"1. GenPosIn : "<<GenPosIn<<", GenNegIn : "<<GenNegIn<<endl;
      for(int ict = 0; ict < nCentBins; ict++){
        for(int ipt = 0; ipt < nPtBins; ipt++){
          for(int irap = 0; irap < nRapBins; irap++){
            if( (AccJpsi==1) && 
                (fabs(GenJpsiRap) >= rap_bound[irap]  && fabs(GenJpsiRap) < rap_bound[irap+1] ) && 
                (GenJpsiPt >= pt_bound[ipt] && GenJpsiPt < pt_bound[ipt+1]) &&
                (gbin >= ct_bound[ict] && gbin < ct_bound[ict+1])) {
              hGenDiMuonf[ifile][ict][ipt][irap]->Fill(GenJpsiMass,GenWeight);
            }
          }
        }
      }
    }//gen loop end

    //cout<<" accepted no "<< NAccep<<endl;

    //=============== Rec Tree Loop ==============================================================================

    int nRecEntries=tree->GetEntries();
    cout<<"Total Entries in reconstructed Tree for pT range "<<fileName[ifile]<<"  "<<nRecEntries<< "====="<<endl;
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

      if(JpsiPt >= 6.5) continue;

      if(IsAccept(muPosPt, muPosEta)){PosIn=1;}
      if(IsAccept(muNegPt, muNegEta)){NegIn=1;}

      double dPhi2mu = muPosPhi - muNegPhi;
      while (dPhi2mu > TMath::Pi()) dPhi2mu -= 2*TMath::Pi();
      while (dPhi2mu <= -TMath::Pi()) dPhi2mu += 2*TMath::Pi();

      double chkCowboy = 1*dPhi2mu;

      if(bCowboy) if(!(chkCowboy > 0.)) continue;
      if(bSailor) if((chkCowboy > 0.)) continue;


      bool mu_Global = ((muPos_nchi2Gl >=0) && (muNeg_nchi2Gl >=0));
      bool mu_Tracker = ((muPos_tracker==1) && (muNeg_tracker==1));

      if(muPos_found > 10 && muPos_pixeLayers > 0 && muPos_nchi2In < 4.0 && TMath::Abs(muPos_dxy) < 3 && TMath::Abs(muPos_dz) < 15 && muPos_nchi2Gl < 20 && muPos_arbitrated==1 && muPos_tracker==1){PosPass=1;}     

      if(muNeg_found > 10 && muNeg_pixeLayers > 0 && muNeg_nchi2In < 4.0 && TMath::Abs(muNeg_dxy) < 3 && TMath::Abs(muNeg_dz) < 15 && muNeg_nchi2Gl < 20 && muNeg_arbitrated==1 && muNeg_tracker==1){NegPass=1;}

      if(hbit1 == 1 && (muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

      // AllCut = 1;
      // without ID cut
      // if((muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && mu_Global && mu_Tracker){AllCut=1;}

      // without trigger matched
      // if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

      int AccJpsi = 0;
      if(JpsiPt >= 3.0 && JpsiPt < 6.5 && fabs(JpsiRap) > 1.6 && fabs(JpsiRap) <= 2.4 && (PosIn==1 && NegIn==1)) AccJpsi = 1;

      double RecCenWeight=0,RecWeight=0;
      RecCenWeight=FindCenWeight(rbin);   
      RecWeight=RecCenWeight*scale[ifile];

      if(PutWeight==0)RecWeight=1;

      //Eff loop for reco
      if((JpsiCharge == 0) && (JpsiVprob > 0.01)) {     
        for(int ict = 0; ict < nCentBins; ict++){
          for(int ipt = 0; ipt < nPtBins; ipt++){
            for(int irap = 0; irap < nRapBins; irap++){
              if(AccJpsi == 1 && (AllCut==1) && 
                  (JpsiPt >= pt_bound[ipt] && JpsiPt < pt_bound[ipt+1]) && 
                  (fabs(JpsiRap) >= rap_bound[irap] && fabs(JpsiRap) < rap_bound[irap+1]) &&  
                  (rbin >= ct_bound[ict] && rbin < ct_bound[ict+1])){
                hRecoDiMuonf[ifile][ict][ipt][irap]->Fill(JpsiMass,RecWeight);
              }
            }
          }
        }
      }
    } //rec tree loop ends
  }  // file loop ends

  //======================  File loop Starts ============================

  ///////////////////////////////////////////////////////////////////
  cout<< " adding "<<endl;
  TFile *outfile;
  char tmp_output[512];
  sprintf(tmp_output,"PrJpsi_LowPt_%s.root", cCd);
  outfile =new TFile(tmp_output, "Recreate");
  TCanvas *c1 = new TCanvas();
  char gtmp[512], gtmp1[512];
  char rtmp[512], rtmp1[512];

  for(int ict = 0; ict < nCentBins; ict++){
    for(int ipt = 0; ipt < nPtBins; ipt++){
      for(int irap = 0; irap < nRapBins; irap++){
        for(int ifile = 0; ifile < nFiles; ifile++){
          hGenDiMuon[ict][ipt][irap]->Add(hGenDiMuonf[ifile][ict][ipt][irap],1);
          hRecoDiMuon[ict][ipt][irap]->Add(hRecoDiMuonf[ifile][ict][ipt][irap],1);
        }
        sprintf(gtmp,"hGenDiMuon_Cent_%1.f_%1.f_pt_%0.1f_%0.1f_rap_%0.1f_%0.1f",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],rap_bound[irap],rap_bound[irap+1]);
        sprintf(gtmp1,"plots/Gen/hGenDiMuon_Cent_%1.f_%1.f_pt_%0.1f_%0.1f_rap_%0.1f_%0.1f_%s.png",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],rap_bound[irap],rap_bound[irap+1],cCd);
        sprintf(rtmp,"hRecDiMuon_Cent_%1.f_%1.f_pt_%0.1f_%0.1f_rap_%0.1f_%0.1f",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],rap_bound[irap],rap_bound[irap+1]);
        sprintf(rtmp1,"plots/Rec/hRecDiMuon_Cent_%1.f_%1.f_pt_%0.1f_%0.1f_rap_%0.1f_%0.1f_%s.png",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],rap_bound[irap],rap_bound[irap+1],cCd);
        hGenDiMuon[ict][ipt][irap]->SetName(gtmp);
        hGenDiMuon[ict][ipt][irap]->Write();
        //hGenDiMuon[ipt][irap][idphi]->Draw();
        //c1->SaveAs(gtmp1);
        hRecoDiMuon[ict][ipt][irap]->SetName(rtmp);
        hRecoDiMuon[ict][ipt][irap]->Write();
        //hRecoDiMuon[ipt][irap][idphi]->Draw();
        //c1->SaveAs(rtmp1);
      }
    }
  }

  cout<<"Starts to calculate efficiency"<<endl;
  dataFile<<""<<endl;
  //=====================Loop for eff========================================================================================//
  //define stuff here for error on weighted samples
  TH3F *eff3D = new TH3F("eff3D","Effiicency;Centrality;p_{T} [GeV/c];rapidity", nCentBins, ct_bound, nPtBins, pt_bound, nRapBins, rap_bound);
  TH1F *eff1DCentPt[nCentBins][nPtBins];
  TH1F *eff1D_tmp;
  eff1D_tmp = new TH1F("eff1D_tmp","Effiicency;rapidity", nRapBins, rap_bound);
  for(int i = 0; i < nCentBins; i++){
    for(int j = 0; j < nPtBins; j++){
      eff1DCentPt[i][j] = (TH1F*)eff1D_tmp->Clone();
    }
  } 


  for(int ict = 0; ict < nCentBins; ict++){
    for(int ipt = 0; ipt < nPtBins; ipt++){
      for(int irap = 0; irap < nRapBins; irap++){
        genNo[ict][ipt][irap] = hGenDiMuon[ict][ipt][irap]->IntegralAndError(1, 100, genErr[ict][ipt][irap]); 
        recoNo[ict][ipt][irap] = hRecoDiMuon[ict][ipt][irap]->IntegralAndError(1, 100, recoErr[ict][ipt][irap]);

        //calculate Eff         
        if(genNo[ict][ipt][irap] == 0 || recoNo[ict][ipt][irap] == 0) {
          eff[ict][ipt][irap] = 0;
          effErr[ict][ipt][irap] = 0;
        }else{
          eff[ict][ipt][irap] = recoNo[ict][ipt][irap]/genNo[ict][ipt][irap]; 

          double tmpGenNo = genNo[ict][ipt][irap];
          double tmpGenErr = genErr[ict][ipt][irap];
          double tmpRecNo = recoNo[ict][ipt][irap];
          double tmpRecErr = recoErr[ict][ipt][irap];
          double tmpEff = eff[ict][ipt][irap];
          double tmpEffErr = 0.0;

          //error    
          double tmp_err_s1_1 = (tmpEff * tmpEff)/(tmpGenNo * tmpGenNo);
          double tmp_err_s1_2 = (tmpRecErr * tmpRecErr);
          double tmp_err_cat_s1 = tmp_err_s1_1 * tmp_err_s1_2;


          double tmp_err_s2_1 = ( (1 - tmpEff)*(1 - tmpEff) ) / (tmpGenNo * tmpGenNo);
          double tmp_err_s2_2 = TMath::Abs(( tmpGenErr*tmpGenErr ) - ( tmpRecErr * tmpRecErr));
          double tmp_err_cat_s2 = tmp_err_s2_1 * tmp_err_s2_2;
          tmpEffErr = sqrt( tmp_err_cat_s1 + tmp_err_cat_s2 );

          effErr[ict][ipt][irap] = tmpEffErr;

          //error without weight
          dataFile<<" Bin ["<<ict<<"]["<<ipt<<"]["<<irap<<"] - "<< " Reco Jpsi : "<< tmpRecNo  <<", Gen Jpsi : "<< tmpGenNo <<endl;
          dataFile<<" Eff ["<<ict<<"]["<<ipt<<"]["<<irap<<"] - "<< tmpEff <<" Error "<< tmpEffErr <<endl;
          cout<<" Bin ["<<ict<<"]["<<ipt<<"]["<<irap<<"] - "<< " Reco Jpsi : "<< tmpRecNo  <<", Gen Jpsi : "<< tmpGenNo <<endl;
          cout<<" Eff ["<<ict<<"]["<<ipt<<"]["<<irap<<"] - "<< tmpEff <<" Error "<< tmpEffErr <<endl;
        }

        int sbin = eff3D->FindBin(xct_bound[ict],xpt_bound[ipt],xrap_bound[irap]);
        eff3D->SetBinContent(sbin,eff[ict][ipt][irap]);
        eff3D->SetBinError(sbin,effErr[ict][ipt][irap]);

        int sbin3 = eff1DCentPt[ict][ipt]->FindBin(xrap_bound[irap]);
        //cout<<"sbin3 : "<<sbin3<<" eff : "<<eff[ict][ipt][irap]<<endl;
        eff1DCentPt[ict][ipt]->SetBinContent(sbin3,eff[ict][ipt][irap]);
        eff1DCentPt[ict][ipt]->SetBinError(sbin3,effErr[ict][ipt][irap]);
        //cout<<"Test Eff : "<<eff1DCentPt[ict][ipt]->GetBinContent(sbin3)<<endl;
      }
    }
  }


  char tmp_eff3D[512], tmp_eff3D_png[512];
  sprintf(tmp_eff3D,"eff_%s",cCd); sprintf(tmp_eff3D_png,"eff_%s.png",cCd);
  eff3D->SetName(tmp_eff3D);
  eff3D->Draw();
  eff3D->Write();
  //c1->SaveAs(tmp_eff3D_png);

  for(int ict = 0; ict < nCentBins; ict++){
    for(int ipt = 0; ipt < nPtBins; ipt++){
      char tmp_eff1D[512], tmp_eff1D_png[512];
      sprintf(tmp_eff1D,"eff1D_Cent_%1.f_%1.f_Pt_%0.1f_%0.1f_%s",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],cCd); sprintf(tmp_eff1D_png,"eff_1D_Cent_%1.f_%1.f_Pt_%0.1f_%0.1f_%s.png",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],cCd);
      eff1DCentPt[ict][ipt]->SetName(tmp_eff1D);
      eff1DCentPt[ict][ipt]->Draw("P");
      eff1DCentPt[ict][ipt]->Write();
      //c1->SaveAs(tmp_eff1D_png);
    }
  }


  //c1->SaveAs("eff_3D_default_etHFm.png");
  outfile->Write();
  //outfile->Close();
  dataFile<<""<<endl;
  dataFile.close();
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

