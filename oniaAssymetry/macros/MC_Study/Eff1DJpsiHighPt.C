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
void Eff1DJpsiHighPt()   
{
  int Prompt =1; int PutWeight = 1;

  for(int iCat = 0; iCat < 3; iCat++){
    for(int iTrg = 0; iTrg < 3; iTrg++){
    // How to use 
    // Decide default/cowboy/sailor
    // Example: if you want to run the case of cowboy
    // bool bDefault = false;
    // bool bCowboy = true;
    // bool bSailor = flase;

    bool bDefault = true; // true : default, false : sailor or cowboy
    bool bCowboy  = false; // true : cowboy only, false : salior
    bool bSailor  = false; // true : cowboy only, false : salior

    // Trigger Selection; You should take only one case among bAllTrg/bNominal/bHBit1
    bool bAllTrg = false;
    bool bNominal = true;
    bool bHBit1 = false;

    if(iTrg == 0) {bAllTrg = true; bNominal = false; bHBit1 = false;}
    if(iTrg == 1) {bAllTrg = false; bNominal = true; bHBit1 = false;}
    if(iTrg == 2) {bAllTrg = false; bNominal = false; bHBit1 = true;}

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

    const char *cTrg[3] = {"AllTrg", "Nominal", "Bit1"};
    TFile *outfile;
    char tmp_output[512];
    if(bAllTrg) sprintf(tmp_output,"PrJpsi_HighPt_%s_%s.root", cCd, cTrg[0]);
    if(bNominal) sprintf(tmp_output,"PrJpsi_HighPt_%s_%s.root", cCd, cTrg[1]);
    if(bHBit1) sprintf(tmp_output,"PrJpsi_HighPt_%s_%s.root", cCd, cTrg[2]);
    outfile =new TFile(tmp_output, "Recreate");

    for(int iSpec = 0; iSpec < 3; iSpec++){

      const int nCentBins = 5;
      const int nPtBins = 7;
      const int nRapBins = 7;
      const int ndPhiBins = 4;
      const int nFiles = 6;
      double ct_bound[nCentBins+1] = {0, 4, 8, 12, 24, 40};
      double xct_bound[nCentBins] = {0.0};
      double pt_bound[nPtBins+1] = {0.0, 4.0, 7.0, 11.0, 14.0, 18.0, 25.0, 30.0};
      //double pt_bound[nPtBins+1] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 30.0};
      //double pt_bound[nPtBins+1] = {3.0, 4.5, 6.5, 8.0, 10.0, 40.0};
      double xpt_bound[nPtBins] = {0.0};
      double rap_bound[nRapBins+1] = {0.0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.4};
      //double rap_bound[nRapBins+1] = {-2.4, -1.6, -1.2, -0.8, 0.8, 1.2, 1.6, 2.4}; 
      double xrap_bound[nRapBins] = {0.0};
      double dphi_bound[ndPhiBins+1] = {0.0, TMath::Pi()/8, TMath::Pi()/4, 3*TMath::Pi()/8, TMath::Pi()/2};
      double xdphi_bound[ndPhiBins] = {0.0};

      const char *cSp[3] = {"Cents","Pts","Raps"};
      char OutTextFile[100];
      sprintf(OutTextFile,"eff_%s_%s.tex", cSp[iSpec], cCd);
      ofstream dataFile(Form(OutTextFile));

      char tmp_start[512];
      sprintf(tmp_start,"%%%% Getting Efficiency starts, Category : %s !!!!! %%%%%", cCd);  
      cout<< tmp_start << endl;
      dataFile<< tmp_start << endl;

      // x, y, z - axis 
      dataFile<<""<<endl;
      dataFile<<"xaxis of Cent"<<endl;
      for(int i = 0; i < nCentBins; i++){
        if(i == (nCentBins-1)){
          xct_bound[i] = ct_bound[i-4] + (ct_bound[i-1]-ct_bound[i-4])/2;
          cout<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
          dataFile<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
        }else{
          xct_bound[i] = ct_bound[i] + (ct_bound[i+1]-ct_bound[i])/2;
          cout<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
          dataFile<<"xct_bound["<<i<<"] : "<<xct_bound[i]<<endl;
        }
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

      int nBins_tmp;
      if(iSpec == 0) { nBins_tmp = nCentBins; }
      if(iSpec == 1) { nBins_tmp = nPtBins; }
      if(iSpec == 2) { nBins_tmp = nRapBins; }
      const int nBins = nBins_tmp;

      TH1F *hTempMass = new TH1F("hTempMass","",100, 2.6, 3.5);
      TH1F *hGenDiMuonf[nFiles][nBins];
      TH1F *hRecoDiMuonf[nFiles][nBins];
      TH1F *hGenDiMuon[nBins];
      TH1F *hRecoDiMuon[nBins];
      double genNo[nBins];
      double genErr[nBins];
      double recoNo[nBins];
      double recoErr[nBins];
      double eff[nBins];
      double effErr[nBins];
      for(int fl = 0; fl < nFiles; fl++){
        for(int i = 0; i < nBins; i++){
          hGenDiMuonf[fl][i] = (TH1F*)hTempMass->Clone();
          hRecoDiMuonf[fl][i] = (TH1F*)hTempMass->Clone();
          hGenDiMuon[i] = (TH1F*)hTempMass->Clone();
          hRecoDiMuon[i] = (TH1F*)hTempMass->Clone();
          hGenDiMuonf[fl][i]->Sumw2();
          hRecoDiMuonf[fl][i]->Sumw2();
          hGenDiMuon[i]->Sumw2();
          hRecoDiMuon[i]->Sumw2();
        }
      }

      char fileName[10][512];
      //scales for different pT bins
      //double scale[5]={1.99236e-07, 4.56079e-08, 9.70856e-09, 4.75729e-09, 2.64798e-09};
      //2.16547e-07 1.99236e-07 4.56079e-08 9.70856e-09 4.75729e-09 2.64798e-09
      double scale[6]={2.16547e-07, 1.99236e-07, 4.56079e-08, 9.70856e-09, 4.75729e-09, 2.64798e-09};
      if(PutWeight==0){scale[0]=(1);scale[1]=(1);scale[2]=(1);scale[3]=(1);scale[4]=(1);scale[5]=(1);}

      // loop for pT
      cout<<"==================Prompt PrJpsi================================================"<<endl;
      /*
         sprintf(fileName[0],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0306_total.root");
         sprintf(fileName[1],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0609_total.root");
         sprintf(fileName[2],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0912_total.root");
         sprintf(fileName[3],"../RootFiles/DiMuonTTree_PromptJpsi_Pt1215_total.root");
         sprintf(fileName[4],"../RootFiles/DiMuonTTree_PromptJpsi_Pt1530_total.root");
         */

      sprintf(fileName[0],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0003_total.root");
      sprintf(fileName[1],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0306_total.root");
      sprintf(fileName[2],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0609_total.root");
      sprintf(fileName[3],"../RootFiles/DiMuonTTree_PromptJpsi_Pt0912_total.root");
      sprintf(fileName[4],"../RootFiles/DiMuonTTree_PromptJpsi_Pt1215_total.root");
      sprintf(fileName[5],"../RootFiles/DiMuonTTree_PromptJpsi_Pt1530_total.root");


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
          if(GenJpsiPt < 6.5) continue;
          if((GenJpsiPt >= 6.5 && GenJpsiPt <= 40.0 && fabs(GenJpsiRap) <= 2.4 && GenPosIn == 1 && GenNegIn == 1)) {AccJpsi = 1;}

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
          double vars = 0.0, bin1 = 0.0, bin2 = 0.0;
          if(iSpec == 0) vars = gbin;
          if(iSpec == 1) vars = GenJpsiPt;
          if(iSpec == 2) vars = GenJpsiRap;
          for(int i = 0; i < nBins; i++){
            if(iSpec == 0){
              bin1 = ct_bound[i]; bin2 = ct_bound[i+1];
              if( (AccJpsi==1) && (vars >= bin1 && vars < bin2) && GenJpsiPt >= 20.0 && GenJpsiPt <= 30.0 && fabs(GenJpsiRap) <= 1.2) {
                hGenDiMuonf[ifile][i]->Fill(GenJpsiMass,GenWeight);
              }
            }
            if(iSpec == 1){
              bin1 = pt_bound[i]; bin2 = pt_bound[i+1];
              if( (AccJpsi==1) && (vars >= bin1 && vars < bin2) && gbin >= 0 && gbin <= 40 && fabs(GenJpsiRap) <= 1.2) {
                hGenDiMuonf[ifile][i]->Fill(GenJpsiMass,GenWeight);
              }
            }
            if(iSpec == 2){
              bin1 = rap_bound[i]; bin2 = rap_bound[i+1];
              if( (AccJpsi==1) && (vars >= bin1 && vars < bin2) && gbin >= 0 && gbin <= 40 && GenJpsiPt >= 20.0 && GenJpsiPt <= 30.0) {
                hGenDiMuonf[ifile][i]->Fill(GenJpsiMass,GenWeight);
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


          if(IsAccept(muPosPt, muPosEta)){PosIn=1;}
          if(IsAccept(muNegPt, muNegEta)){NegIn=1;}

          double dPhi2mu = muPosPhi - muNegPhi;
          while (dPhi2mu > TMath::Pi()) dPhi2mu -= 2*TMath::Pi();
          while (dPhi2mu <= -TMath::Pi()) dPhi2mu += 2*TMath::Pi();

          double chkCowboy = 1*dPhi2mu;

          if(bCowboy) if(!(chkCowboy > 0.)) continue;
          if(bSailor) if((chkCowboy > 0.)) continue;

          if(JpsiPt < 6.5) continue;

          bool mu_Global = ((muPos_nchi2Gl >=0) && (muNeg_nchi2Gl >=0));
          bool mu_Tracker = ((muPos_tracker==1) && (muNeg_tracker==1));

          if(muPos_found > 10 && muPos_pixeLayers > 0 && muPos_nchi2In < 4.0 && TMath::Abs(muPos_dxy) < 3 && TMath::Abs(muPos_dz) < 15 && muPos_nchi2Gl < 20 && muPos_arbitrated==1 && muPos_tracker==1){PosPass=1;}     

          if(muNeg_found > 10 && muNeg_pixeLayers > 0 && muNeg_nchi2In < 4.0 && TMath::Abs(muNeg_dxy) < 3 && TMath::Abs(muNeg_dz) < 15 && muNeg_nchi2Gl < 20 && muNeg_arbitrated==1 && muNeg_tracker==1){NegPass=1;}

          if(bAllTrg){
            if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}
          }
          if(bNominal){
            if((muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}
          }
          if(bHBit1){
            if(hbit1 == 1 && (muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}
          }


          // AllCut = 1;
          // without ID cut
          // if((muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && mu_Global && mu_Tracker){AllCut=1;}

          // without trigger matched
          // if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

          double RecCenWeight=0,RecWeight=0;
          RecCenWeight=FindCenWeight(rbin);   
          RecWeight=RecCenWeight*scale[ifile];

          if(PutWeight==0)RecWeight=1;

          //Eff loop for reco
          double vars = 0.0, bin1 = 0.0, bin2 = 0.0;
          if(iSpec == 0) vars = rbin;
          if(iSpec == 1) vars = JpsiPt;
          if(iSpec == 2) vars = JpsiRap;
          if((JpsiCharge == 0) && (JpsiVprob > 0.01)) {     
            for(int i = 0; i < nBins; i++){
              if(iSpec == 0){
                bin1 = ct_bound[i]; bin2 = ct_bound[i+1];
                if( (AllCut == 1) && (vars >= bin1 && vars < bin2) && JpsiPt >= 20.0 && JpsiPt <= 30.0 && fabs(JpsiRap) <= 1.2) {
                  hRecoDiMuonf[ifile][i]->Fill(JpsiMass,RecWeight);
                }
              }
              if(iSpec == 1){
                bin1 = pt_bound[i]; bin2 = pt_bound[i+1];
                if( (AllCut == 1) && (vars >= bin1 && vars < bin2) && rbin >= 0 && rbin <= 40 && fabs(JpsiRap) <= 1.2) {
                  hRecoDiMuonf[ifile][i]->Fill(JpsiMass,RecWeight);
                }
              }
              if(iSpec == 2){ 
                bin1 = rap_bound[i]; bin2 = rap_bound[i+1];
                if( (AllCut == 1) && (vars >= bin1 && vars < bin2) && rbin >= 0 && rbin <= 40 && JpsiPt >= 20.0 && JpsiPt <= 30.0) {
                  hRecoDiMuonf[ifile][i]->Fill(JpsiMass,RecWeight);
                }
              }
            }
          }
        } //rec tree loop ends
      }  // file loop ends

      //======================  File loop Starts ============================

      ///////////////////////////////////////////////////////////////////
      cout<< " adding "<<endl;
      TCanvas *c1 = new TCanvas();
      char gtmp[512], gtmp1[512];
      char rtmp[512], rtmp1[512];

      for(int i = 0; i < nBins; i++){
        for(int ifile = 0; ifile < nFiles; ifile++){
          hGenDiMuon[i]->Add(hGenDiMuonf[ifile][i],1);
          hRecoDiMuon[i]->Add(hRecoDiMuonf[ifile][i],1);
        }
        if(iSpec == 0) sprintf(gtmp,"hGenDiMuon_Cent_%1.f_%1.f",ct_bound[i],ct_bound[i+1]);
        if(iSpec == 1) sprintf(gtmp,"hGenDiMuon_Pt_%0.1f_%0.1f",pt_bound[i],pt_bound[i+1]);
        if(iSpec == 2) sprintf(gtmp,"hGenDiMuon_Rap_%0.1f_%0.1f",rap_bound[i],rap_bound[i+1]);
        if(iSpec == 0) sprintf(gtmp1,"plots/Gen/hGenDiMuon_Cent_%1.f_%1.f_%s.png",ct_bound[i],ct_bound[i+1],cCd);
        if(iSpec == 1) sprintf(gtmp1,"plots/Gen/hGenDiMuon_Pt_%0.1f_%0.1f_%s.png",pt_bound[i],pt_bound[i+1],cCd);
        if(iSpec == 2) sprintf(gtmp1,"plots/Gen/hGenDiMuon_Rap_%0.1f_%0.1f_%s.png",rap_bound[i],rap_bound[i+1],cCd);
        if(iSpec == 0) sprintf(rtmp,"hRecDiMuon_Cent_%1.f_%1.f",ct_bound[i],ct_bound[i+1]);
        if(iSpec == 1) sprintf(rtmp,"hRecDiMuon_Pt_%0.1f_%0.1f",pt_bound[i],pt_bound[i+1]);
        if(iSpec == 2) sprintf(rtmp,"hRecDiMuon_Rap_%0.1f_%0.1f",rap_bound[i],rap_bound[i+1]);
        if(iSpec == 0) sprintf(rtmp1,"plots/Rec/hRecDiMuon_Cent_%1.f_%1.f_%s.png",ct_bound[i],ct_bound[i+1],cCd);
        if(iSpec == 1) sprintf(rtmp1,"plots/Rec/hRecDiMuon_Pt_%0.1f_%0.1f_%s.png",pt_bound[i],pt_bound[i+1],cCd);
        if(iSpec == 2) sprintf(rtmp1,"plots/Rec/hRecDiMuon_Rap_%0.1f_%0.1f_%s.png",rap_bound[i],rap_bound[i+1],cCd);
        hGenDiMuon[i]->SetName(gtmp);
        //hGenDiMuon[i]->Write();
        //hGenDiMuon[ipt][irap][idphi]->Draw();
        //c1->SaveAs(gtmp1);
        hRecoDiMuon[i]->SetName(rtmp);
        //hRecoDiMuon[i]->Write();
        //hRecoDiMuon[ipt][irap][idphi]->Draw();
        //c1->SaveAs(rtmp1);
      }

      cout<<"Starts to calculate efficiency"<<endl;
      dataFile<<""<<endl;
      //=====================Loop for eff========================================================================================//
      //define stuff here for error on weighted samples
      for(int i = 0; i < nBins; i++){
        genNo[i] = hGenDiMuon[i]->IntegralAndError(1, 100, genErr[i]); 
        recoNo[i] = hRecoDiMuon[i]->IntegralAndError(1, 100, recoErr[i]);

        //calculate Eff         
        if(genNo[i] == 0 || recoNo[i] == 0) {
          eff[i] = 0;
          effErr[i] = 0;
        }else{
          eff[i] = recoNo[i]/genNo[i]; 

          double tmpGenNo = genNo[i];
          double tmpGenErr = genErr[i];
          double tmpRecNo = recoNo[i];
          double tmpRecErr = recoErr[i];
          double tmpEff = eff[i];
          double tmpEffErr = 0.0;

          //error    
          double tmp_err_s1_1 = (tmpEff * tmpEff)/(tmpGenNo * tmpGenNo);
          double tmp_err_s1_2 = (tmpRecErr * tmpRecErr);
          double tmp_err_cat_s1 = tmp_err_s1_1 * tmp_err_s1_2;


          double tmp_err_s2_1 = ( (1 - tmpEff)*(1 - tmpEff) ) / (tmpGenNo * tmpGenNo);
          double tmp_err_s2_2 = TMath::Abs(( tmpGenErr*tmpGenErr ) - ( tmpRecErr * tmpRecErr));
          double tmp_err_cat_s2 = tmp_err_s2_1 * tmp_err_s2_2;
          tmpEffErr = sqrt( tmp_err_cat_s1 + tmp_err_cat_s2 );

          effErr[i] = tmpEffErr;

          //error without weight
          dataFile<<" Bin ["<<i<<"] - "<< " Reco Jpsi : "<< tmpRecNo  <<", Gen Jpsi : "<< tmpGenNo <<endl;
          dataFile<<" Eff ["<<i<<"] - "<< tmpEff <<" Error "<< tmpEffErr <<endl;
          cout<<" Bin ["<<i<<"] - "<< " Reco Jpsi : "<< tmpRecNo  <<", Gen Jpsi : "<< tmpGenNo <<endl;
          cout<<" Eff ["<<i<<"] - "<< tmpEff <<" Error "<< tmpEffErr <<endl;
        }
      }
      TH1F *hEff = new TH1F();
      int nEffBins;
      if(iSpec == 0) {hEff = new TH1F("hEff","hEff;Centrality (%);Efficiency",nCentBins,ct_bound); nEffBins = nCentBins;}
      if(iSpec == 1) {hEff = new TH1F("hEff","hEff;p_{T} GeV/c;Efficiency",nPtBins,pt_bound); nEffBins = nPtBins;}
      if(iSpec == 2) {hEff = new TH1F("hEff","hEff;y;Efficiency",nRapBins,rap_bound); nEffBins = nRapBins;}
      
      for(int i = 0; i < nEffBins; i++){
        hEff->SetBinContent(i+1,eff[i]);
        hEff->SetBinError(i+1,effErr[i]);
        cout<<"Trying to measure eff vs "<<cSp[iSpec]<<endl;
        cout<<"Eff : "<<eff[i]<<", err : "<<effErr[i]<<endl;
        dataFile<<"Eff : "<<eff[i]<<", err : "<<effErr[i]<<endl;
      }

      outfile->cd();
      hEff->Draw();
      char tmp_histo[512];
      sprintf(tmp_histo,"hEff_%s",cSp[iSpec]);
      hEff->SetName(tmp_histo);
      hEff->Write();
      //c1->SaveAs("eff_3D_default_etHFm.png");
      //outfile->Close();
      dataFile<<""<<endl;
      dataFile.close();
    }
    outfile->Write();
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

