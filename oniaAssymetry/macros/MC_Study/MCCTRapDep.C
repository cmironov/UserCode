//===================== Macro to calculate PrJpsi Eff or Acc x Eff ============================================
//read tree from DiMuonOnia2DPlots.cc
#ifndef __CINT__
#endif
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
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
//iSpec = 1: pT, = 2: Rap, = 3: Cent, = 4: Phi, = 5: dPhi, PutWeight = 1 means put weight
// Prompt =1 prompt PrJpsi 
bool IsAccept(Double_t pt, Double_t eta); //you can define acceptance here 
double FindCenWeight(int Bin);//gives you weight according to cent
double FindV2Weight(double v2wg);
void DoEffCor(TFile* a, int b, double c, double d, double e, double *f);
//void DoEffCor(TFile* a, char *b, double c, double d, double e, double *f);
void MCCTRapDep()   
{
  for(int iSpec = 0; iSpec < 3; iSpec++){
    int Prompt =1; int PutWeight = 1;
    bool bDefault = true; // true : cowboy only, false : salior
    bool bCowboy = false; // true : cowboy only, false : salior
    bool bSailor = false; // true : cowboy only, false : salior
    double fake_v2 = 0.3;

    if(iSpec == 0) {bDefault = true;bCowboy = false;bSailor = false;}
    if(iSpec == 1) {bDefault = false;bCowboy = true;bSailor = false;}
    if(iSpec == 2) {bDefault = false;bCowboy = false;bSailor = true;}

    char cCd[512];
    int iCd = 0;
    if(bDefault) {sprintf(cCd, "default"); iCd = 0;}
    if(bCowboy)  {sprintf(cCd, "cowboy"); iCd = 1;}
    if(bSailor)  {sprintf(cCd, "sailor"); iCd = 2;}

    double eff[10] = {0.0};
    double effErr[10] = {0.0};

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
      //gStyle->SetOptStat(1); // most of the time, sometimes "nemriou" might be useful to display name, 
      gStyle->SetOptStat("emri");
      //number of entries, mean, rms, integral, overflow and underflow
      gStyle->SetOptFit(1); // set to 1 only if you want to display fit results
      //==================================== Define Histograms====================================================
      char OutTextFile[100]; 
      if(iCat == 0) sprintf(OutTextFile,"MCCT_PrJpsi_Raps_0012_dPhi_%s.tex", cCd);
      if(iCat == 1) sprintf(OutTextFile,"MCCT_PrJpsi_Raps_1216_dPhi_%s.tex", cCd);
      if(iCat == 2) sprintf(OutTextFile,"MCCT_PrJpsi_Raps_1624H_dPhi_%s.tex", cCd);
      if(iCat == 3) sprintf(OutTextFile,"MCCT_PrJpsi_Raps_1624L_dPhi_%s.tex", cCd);
      ofstream dataFile(Form(OutTextFile));
      TH1D *diMuonsInvMass_Gen = new TH1D("diMuonsInvMass_Gen","diMuonsInvMass_Gen", 100,2.6,3.5);

      TH1D *diMuonsPt_Gen = new TH1D("diMuonsPt_Gen","diMuonsPt_Gen", 100,0,50);
      TH1D *Bin_Gen = new TH1D("Bin_Gen","Bin_Gen", 40,0,40);
      //==============================================Define Acc Eff Stuff here===========================================
      // Pt bin sizes
      // 0-1.5, 1.5-3, 3-4.5, 4.5-6, 6-7.5...

      int Nptbin = 4;
      double pt_bound[100] = {0};

      pt_bound[0] = 0.0;
      pt_bound[1] = TMath::Pi()*2/16;
      pt_bound[2] = TMath::Pi()*4/16;
      pt_bound[3] = TMath::Pi()*6/16;
      pt_bound[4] = TMath::Pi()*8/16;

      TH1D *EffJpsiGen = new TH1D("EffJpsiGen","",Nptbin,pt_bound[0],pt_bound[Nptbin]);
      TH1D *EffJpsiReco = new TH1D("EffJpsiReco","",Nptbin,pt_bound[0],pt_bound[Nptbin]);
      EffJpsiGen->Sumw2();
      EffJpsiReco->Sumw2();


      //X Axis error on Eff graph 
      double PT[100], DelPT[100], mom_err[100];
      for (Int_t ih = 0; ih < Nptbin; ih++) {
        PT[ih] = (pt_bound[ih] + pt_bound[ih+1])/2.0;
        DelPT[ih] = pt_bound[ih+1] - pt_bound[ih];
        mom_err[ih] = DelPT[ih]/2.0;
        cout<<PT[ih]<<" X axis graph "<<DelPT[ih]<<endl;
      }

      double genError, recError;
      double gen_pt[100]={0}, gen_ptError[100]={0}; 
      double rec_pt[100]={0}, rec_ptError[100]={0}; 
      double Eff_cat_1[100]={0},Err_Eff_cat_1[100]={0};  

      // Histogram 2D Arrays
      TH1D *diMuonsInvMass_GenA[10][1000];
      TH1D *diMuonsInvMass_RecA[10][1000];
      TH1D *diMuonsPt_GenA[10][1000];
      TH1D *diMuonsPt_RecA[10][1000];
      char nameGen[10][500], nameRec[10][500], nameGenPt[10][500], nameRecPt[10][500];
      char namePt_1B[500];//for bkg func
      for (int ifile = 0; ifile <= 6; ifile++) {
        for (Int_t ih = 0; ih < Nptbin; ih++) {
          sprintf(nameGen[ifile],"DiMuonMassGen_pt_%d_%d_%d",ifile,iCat,ih);
          sprintf(nameRec[ifile],"DiMuonMassRec_pt_%d_%d_%d",ifile,iCat,ih);

          sprintf(nameGenPt[ifile],"DiMuonPtGen_pt_%d_%d_%d",ifile,iCat,ih);
          sprintf(nameRecPt[ifile],"DiMuonPtRec_pt_%d_%d_%d",ifile,iCat,ih);

          diMuonsInvMass_GenA[ifile][ih]= new TH1D(nameGen[ifile],nameGen[ifile],  100,2.6,3.5); //for eff Gen;
          diMuonsInvMass_GenA[ifile][ih]->Sumw2();
          diMuonsInvMass_GenA[ifile][ih]->SetMarkerStyle(7);
          diMuonsInvMass_GenA[ifile][ih]->SetMarkerColor(4);
          diMuonsInvMass_GenA[ifile][ih]->SetLineColor(4);

          diMuonsInvMass_RecA[ifile][ih] = new TH1D(nameRec[ifile],nameRec[ifile], 100,2.6,3.5); //for eff Rec;
          diMuonsInvMass_RecA[ifile][ih]->Sumw2();
          diMuonsInvMass_RecA[ifile][ih]->SetMarkerStyle(8);
          diMuonsInvMass_RecA[ifile][ih]->SetMarkerColor(4);
          diMuonsInvMass_RecA[ifile][ih]->SetLineColor(4);

          diMuonsPt_GenA[ifile][ih]= new TH1D(nameGenPt[ifile],nameGenPt[ifile],  100,0,40); //for eff Gen;
          diMuonsPt_RecA[ifile][ih]= new TH1D(nameRecPt[ifile],nameRecPt[ifile],  100,0,40); //for eff Rec;
        }
      }

      //===========================================Input Root File============================================================
      char fileName[10][500];
      //scales for different pT bins
      double scale[10]={0};  


      //2.06586e-07 1.99236e-07 4.44677e-08 9.49795e-09 4.7002e-09 2.57384e-09
      //2.06586e-07 1.99236e-07 4.44677e-08 9.06712e-09 4.7002e-09 3.2173e-08
      //2.21184e-07 2.14494e-07 4.71056e-08 1.06043e-08 4.75729e-09 2.66996e-09
      //2.16547e-07 1.99236e-07 4.56079e-08 9.70856e-09 4.75729e-09 2.64798e-09
      scale[0]=2.35829e-07;
      scale[1]=1.99854e-07;
      scale[2]=4.48263e-08;
      scale[3]=1.01144e-08;
      scale[4]=4.89604e-09;
      scale[5]=2.62102e-09;

      if(PutWeight==0){scale[0]=(1);scale[1]=(1);scale[2]=(1);scale[3]=(1);scale[4]=(1);scale[5]=(1);}

      if(Prompt ==1){
        cout<<"==================Prompt PrJpsi================================================"<<endl;
        sprintf(fileName[0],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0003_total.root");
        sprintf(fileName[1],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0306_total.root");
        sprintf(fileName[2],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0609_total.root");
        sprintf(fileName[3],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt0912_total.root");
        sprintf(fileName[4],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1215_total.root");
        sprintf(fileName[5],"/Users/donghomoon/Analysis/2011HIRData/2011JpsiV2/20111126_Jpsi_Psi/EffCorrection/JpsiEff_New_0329_NewBins_3DEff/0430_4Deff_lowPt_With_dPhi_JpsiGenPsi_HighPt/gRpRootFiles/DiMuonTTree_PromptJpsi_Pt1530_total.root");

        //sprintf(fileName[6],"../RootFiles/DiMuonTTree_PromptJpsi_Pt30XX_total.root");
      }

      TFile *infile;
      TTree *tree;
      TTree *gentree;

      //======================  File loop Starts ============================
      for(int ifile =0; ifile<=5; ifile++){
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
        double muPosPx, muPosPy, muPosPz,  muPosEta, muPosPt,muPosP, muPosPhi;
        double muNegPx, muNegPy, muNegPz,  muNegEta, muNegPt,muNegP, muNegPhi;
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

        // HLTrigger
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
        //tree->SetBranchAddress("JpsiPsi",&JpsiPsi[38]);
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
        gentree->SetBranchAddress("GenJpsiEta",    &GenJpsiEta);
        gentree->SetBranchAddress("GenJpsiPsi",    &GenJpsiPsi);
        gentree->SetBranchAddress("GenJpsiRap",    &GenJpsiRap);
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

          double gdPhi2mu = GenmuPosPhi - GenmuNegPhi;
          while (gdPhi2mu > TMath::Pi()) gdPhi2mu -= 2*TMath::Pi();
          while (gdPhi2mu <= -TMath::Pi()) gdPhi2mu += 2*TMath::Pi();
          double gchkCowboy = 1*gdPhi2mu;

          if(bCowboy && !(gchkCowboy > 0.)) continue;
          if(bSailor && (gchkCowboy > 0.)) continue;

          if((GenPosIn ==1 && GenNegIn==1)) NAccep++;

          Bin_Gen->Fill(gbin);

          double GenCenWeight =0, GenWeight =0;
          GenCenWeight=FindCenWeight(gbin);
          
          double gJpsidPhi = TMath::Abs(GenJpsiPsi);
          double gmean_dphi = TMath::Pi()/16;
          if (gJpsidPhi>TMath::Pi()/8) gmean_dphi = 3*TMath::Pi()/16;
          if (gJpsidPhi>TMath::Pi()/4) gmean_dphi = 5*TMath::Pi()/16;
          if (gJpsidPhi>3*TMath::Pi()/8) gmean_dphi = 7*TMath::Pi()/16;
          GenWeight=GenCenWeight*scale[ifile]*(2.0/TMath::Pi()*(1+2*fake_v2*cos(2*gmean_dphi)*TMath::Gaus(GenJpsiPt,15,4,0)));

          if(PutWeight==0) GenWeight=1; 


          for (Int_t ih = 0; ih < Nptbin; ih++) {

            //adding pT of all pt bins to see diss is cont

            if(iCat == 0) {
              if(  (gbin >= 4 && gbin < 24) && (GenPosIn==1 && GenNegIn==1) && (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) &&  (TMath::Abs(GenJpsiRap)<1.2 && TMath::Abs(GenJpsiRap) >= 0.0 ) && (TMath::Abs(GenJpsiPsi)>=pt_bound[ih] && TMath::Abs(GenJpsiPsi)<pt_bound[ih+1])){diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass,GenWeight);}
            }
            if(iCat == 1) {
              if(  (gbin >= 4 && gbin < 24) && (GenPosIn==1 && GenNegIn==1) && (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) &&  (TMath::Abs(GenJpsiRap)<1.6 && TMath::Abs(GenJpsiRap) >= 1.2 ) && (TMath::Abs(GenJpsiPsi)>=pt_bound[ih] && TMath::Abs(GenJpsiPsi)<pt_bound[ih+1])){diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass,GenWeight);}
            }
            if(iCat == 2) {
              if(  (gbin >= 4 && gbin < 24) && (GenPosIn==1 && GenNegIn==1) && (GenJpsiPt >= 6.5 && GenJpsiPt < 40.0) &&  (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 1.6 ) && (TMath::Abs(GenJpsiPsi)>=pt_bound[ih] && TMath::Abs(GenJpsiPsi)<pt_bound[ih+1])){diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass,GenWeight);}
            }
            if(iCat == 3) {
              if(  (gbin >= 4 && gbin < 24) && (GenPosIn==1 && GenNegIn==1) && (GenJpsiPt >= 3.0 && GenJpsiPt < 6.5) &&  (TMath::Abs(GenJpsiRap)<2.4 && TMath::Abs(GenJpsiRap) >= 1.6 ) && (TMath::Abs(GenJpsiPsi)>=pt_bound[ih] && TMath::Abs(GenJpsiPsi)<pt_bound[ih+1])){diMuonsInvMass_GenA[ifile][ih]->Fill(GenJpsiMass,GenWeight);}
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
        sprintf(tmp_eff_input,"../../EffRoots/PrJpsi_HighPt_%s.root", cCd);
        sprintf(tmp_input_histo,"eff_%s", cCd);

        TFile *eff_input;
        eff_input=new TFile(tmp_eff_input,"R");

        int nRecEntries=tree->GetEntries();
        cout<<"Total Entries in reconstructed Tree for pT range "<<fileName[ifile]<<"  "<<nRecEntries<< "====="<<endl;
        for(int i=0; i<nRecEntries; i++)  {     
          tree->GetEntry(i);
          //if(bCowboy && !(gchkCowboy > 0.)) {cout<<"This is not Cowboy from Gen"<<endl; continue;}
          //if(bSailor && (gchkCowboy > 0.)) {cout<<"This is not Sailor from Gen"<<endl; continue;}
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

          if(hbit1 == 0) continue;

          if(IsAccept(muPosPt, muPosEta)){PosIn=1;}
          if(IsAccept(muNegPt, muNegEta)){NegIn=1;}

          bool mu_Global = ((muPos_nchi2Gl >=0) && (muNeg_nchi2Gl >=0));
          bool mu_Tracker = ((muPos_tracker==1) && (muNeg_tracker==1));

          double dPhi2mu = muPosPhi - muNegPhi;
          while (dPhi2mu > TMath::Pi()) dPhi2mu -= 2*TMath::Pi();
          while (dPhi2mu <= -TMath::Pi()) dPhi2mu += 2*TMath::Pi();
          double chkCowboy = 1*dPhi2mu;

          if(bCowboy && !(chkCowboy > 0.)) continue;
          if(bSailor && (chkCowboy > 0.)) continue;

          if(muPos_found > 10 && muPos_pixeLayers > 0 && muPos_nchi2In < 4.0 && TMath::Abs(muPos_dxy) < 3 && TMath::Abs(muPos_dz) < 15 && muPos_nchi2Gl < 20  &&
              muPos_arbitrated==1 && muPos_tracker==1){PosPass=1;}
          if(muNeg_found >10 && muNeg_pixeLayers >0 && muNeg_nchi2In <4.0 && TMath::Abs(muNeg_dxy) < 3 && TMath::Abs(muNeg_dz) < 15 && muNeg_nchi2Gl < 20 &&
              muNeg_arbitrated==1 && muNeg_tracker==1){NegPass=1;}

          //cout<<"Cut checks, muPos_matches : "<<muPos_matches<<", muNeg_matches : "<<muNeg_matches<<", PosIn : "<<PosIn<<", NegIn : "<<NegIn
          //    <<", PosPass : "<<PosPass<<", NegPass : "<<NegPass<<", mu_Global : "<<mu_Global<<", mu_Tracker : "<<mu_Tracker<<endl;
          //if((PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

          if((muPos_matches==1 && muNeg_matches==1) && (PosIn==1 && NegIn==1) && (PosPass==1 && NegPass==1)&& mu_Global && mu_Tracker){AllCut=1;}

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
          RecWeight=RecCenWeight*scale[ifile]*(2.0/TMath::Pi()*(1+2*fake_v2*cos(2*rmean_dphi)*TMath::Gaus(JpsiPt,15,4,0)));


          if(PutWeight==0)RecWeight=1;


          //Eff loop for reco
          for (Int_t ih = 0; ih < Nptbin; ih++) {
            if((JpsiCharge == 0) && (JpsiVprob > 0.01)) {     
                if(iCat == 0) {
                  if((AllCut==1) && (rbin >= 4 && rbin < 24) && (JpsiPt>=6.5 && JpsiPt<40.0) && 
                      (TMath::Abs(JpsiRap) < 1.2 && TMath::Abs(JpsiRap) >= 0.0) && 
                      (TMath::Abs(JpsiGenPsi) > pt_bound[ih] && TMath::Abs(JpsiGenPsi) <=pt_bound[ih+1])){
                  DoEffCor(eff_input, iCd, rbin, JpsiPt, JpsiRap, eff_cor);
                  //DoEffCor(eff_input, tmp_input_histo, TMath::Abs(JpsiGenPsi), JpsiPt, TMath::Abs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
                }
                if(iCat == 1) {
                  if((AllCut==1) && (rbin >= 4 && rbin < 24) && (JpsiPt>=6.5 && JpsiPt<40.0) && 
                      (TMath::Abs(JpsiRap) < 1.6 && TMath::Abs(JpsiRap) >= 1.2) && 
                      (TMath::Abs(JpsiGenPsi) > pt_bound[ih] && TMath::Abs(JpsiGenPsi) <=pt_bound[ih+1])){
                  DoEffCor(eff_input, iCd, rbin, JpsiPt, JpsiRap, eff_cor);
                  //DoEffCor(eff_input, tmp_input_histo, TMath::Abs(JpsiGenPsi), JpsiPt, TMath::Abs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
                }
                if(iCat == 2) {
                  if((AllCut==1) && (rbin >= 4 && rbin < 24) && (JpsiPt>=6.5 && JpsiPt<40.0) && 
                      (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 1.6) && 
                      (TMath::Abs(JpsiGenPsi) > pt_bound[ih] && TMath::Abs(JpsiGenPsi) <=pt_bound[ih+1])){
                  DoEffCor(eff_input, iCd, rbin, JpsiPt, JpsiRap, eff_cor);
                  //DoEffCor(eff_input, tmp_input_histo, TMath::Abs(JpsiGenPsi), JpsiPt, TMath::Abs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
                }
                if(iCat == 3) {
                  if((AllCut==1) && (rbin >= 4 && rbin < 24) && (JpsiPt>=3.0 && JpsiPt<6.5) && 
                      (TMath::Abs(JpsiRap) < 2.4 && TMath::Abs(JpsiRap) >= 1.6) && 
                      (TMath::Abs(JpsiGenPsi) > pt_bound[ih] && TMath::Abs(JpsiGenPsi) <=pt_bound[ih+1])){
                  DoEffCor(eff_input, iCd, rbin, JpsiPt, JpsiRap, eff_cor);
                  //DoEffCor(eff_input, tmp_input_histo, TMath::Abs(JpsiGenPsi), JpsiPt, TMath::Abs(JpsiRap), eff_cor);
                  diMuonsInvMass_RecA[ifile][ih]->Fill(JpsiMass,RecWeight*((double)1.0/eff_cor[0]));}
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

      for(Int_t ih = 0; ih < Nptbin; ih++){
        diMuonsInvMass_RecA1[ih] = diMuonsInvMass_RecA[0][ih];
        diMuonsInvMass_GenA1[ih] = diMuonsInvMass_GenA[0][ih];
        diMuonsPt_GenA1[ih] = diMuonsPt_GenA[0][ih];
        diMuonsPt_RecA1[ih] = diMuonsPt_RecA[0][ih];

        for (int ifile = 1; ifile <= 6; ifile++) {
          diMuonsInvMass_RecA1[ih]->Add(diMuonsInvMass_RecA[ifile][ih]);
          diMuonsInvMass_GenA1[ih]->Add(diMuonsInvMass_GenA[ifile][ih]);     
          diMuonsPt_GenA1[ih]->Add(diMuonsPt_GenA[ifile][ih]); 
          diMuonsPt_RecA1[ih]->Add(diMuonsPt_RecA[ifile][ih]); 
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
      double errEff_cat_S1[100]={0};
      double errEff_cat_S2[100]={0};
      double errEff_cat_S1_1[100]={0},errEff_cat_S1_2[100]={0};
      double errEff_cat_S2_1[100]={0},errEff_cat_S2_2[100]={0};
      char PlotName[500],PlotName1[500];
      char GPlotName[500],GPlotName1[500];

      for (Int_t ih = 0; ih < Nptbin; ih++) {

        gen_pt[ih] = diMuonsInvMass_GenA1[ih]->IntegralAndError(1, 100, genError);
        gen_ptError[ih]= genError;
        cout<<" gen_pt[ih] "<< gen_pt[ih] <<" error   "<<  gen_ptError[ih]<<endl;

        sprintf(PlotName,"plots/DiMuonMass_Raps_%d_PsiBin_%d_%s.png",iCat, ih, cCd);
        sprintf(PlotName1,"plots/DiMuonMass_Raps_%d_PsiBin_%d_%s.pdf",iCat, ih, cCd);


        // cout<<" *********************** "<<diMuonsInvMass_RecA1[ih]->GetMaximum()<<endl;
        //giving inetial value for crystall ball fourth parameter 
        diMuonsInvMass_RecA1[ih]->Rebin(2);
        GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[ih]->Integral(0,50));
        GAUSPOL->FixParameter(1, diMuonsInvMass_RecA1[ih]->GetBinWidth(1));
        //GAUSPOL->SetParameter(0, diMuonsInvMass_RecA1[ih]->GetMaximum());
        //new TCanvas;
        //diMuonsInvMass_RecA1[ih]->Draw();

        new TCanvas;
        diMuonsInvMass_RecA1[ih]->Fit("GAUSPOL","EMRQ", "", mass_low, mass_high); // Jpsi
        //diMuonsInvMass_RecA1[ih]->Fit("GAUSPOL","LLMERQ", "", 8.5,10.5); // Jpsi
        //diMuonsInvMass_RecA1[ih]->Fit("GAUSPOL","LLMER", "", mass_low, mass_high); // Jpsi
        diMuonsInvMass_RecA1[ih]->DrawCopy("EPLsame");

        // new TCanvas;
        //diMuonsInvMass_RecA1[ih]->Fit("GAUSPOL_1","LLMER", "", mass_low, mass_high);
        //diMuonsInvMass_RecA1[ih]->DrawCopy("EPLsame");
        //gPad->Print(PlotName);
        //gPad->Print(PlotName1);

        //  cout << GAUSPOL_1->GetChisquare()<<endl;
        //for(int i=0;i<=100;i++) {cout<<i<<"  "<<diMuonsInvMass_RecA1[ih]->GetBinContent(i)<<endl;}
        //return;
        //double JpsiMass = GAUSPOL_1->GetParameter(2);
        //double JpsiWidth = GAUSPOL_1->GetParameter(3);
        //double JpsiYield = GAUSPOL_1->GetParameter(4); 


        double JpsiMass = GAUSPOL->GetParameter(2);
        double JpsiWidth = GAUSPOL->GetParameter(3);

        double JpsiYield = GAUSPOL->GetParameter(0); 
        Double_t JpsiYieldError = GAUSPOL->GetParError(0); 

        //cout<<JpsiYieldError<<"*****************"<<endl;

        //if(TMath::IsNan(JpsiYieldError)=1) {JpsiYieldError=TMath::Sqrt(JpsiYield);}

        double par[20];
        GAUSPOL->GetParameters(par);
        sprintf(namePt_1B,"pt_1B_%d",ih);

        backfun_1 = new TF1(namePt_1B, Pol2, mass_low, mass_high, 3);
        backfun_1->SetParameters(&par[3]);

        double MassLow=(JpsiMass-3*JpsiWidth);
        double MassHigh=(JpsiMass+3*JpsiWidth);


        int binlow =diMuonsInvMass_RecA1[ih]->GetXaxis()->FindBin(MassLow);
        int binhi =diMuonsInvMass_RecA1[ih]->GetXaxis()->FindBin(MassHigh);

        //double binwidth=diMuonsInvMass_RecA1[ih]->GetBinWidth(1);
        //yield by function 
        //rec_pt[ih] = JpsiYield;
        //rec_ptError[ih]= JpsiYieldError;

        //yield by histogram integral
        rec_pt[ih] = diMuonsInvMass_RecA1[ih]->IntegralAndError(binlow, binhi,recError);
        rec_ptError[ih]= recError;

        //calculate Eff         
        Eff_cat_1[ih] = rec_pt[ih]/gen_pt[ih]; 

        //calculate error on eff
        GenNo[ih]=gen_pt[ih];
        Eff[ih]= Eff_cat_1[ih];
        GenError[ih]=gen_ptError[ih];
        RecError[ih]=rec_ptError[ih];

        //error    
        errEff_cat_S1_1[ih]= ( (Eff[ih] * Eff[ih]) /(GenNo[ih] * GenNo[ih]) );
        errEff_cat_S1_2[ih]= (RecError[ih] * RecError[ih]);  
        errEff_cat_S1[ih]= (errEff_cat_S1_1[ih] * errEff_cat_S1_2[ih]);


        errEff_cat_S2_1[ih]= ( (1 - Eff[ih])* (1 - Eff[ih]) ) / ( GenNo[ih] * GenNo[ih]);
        errEff_cat_S2_2[ih]= (GenError[ih] * GenError[ih] ) - ( RecError[ih] * RecError[ih] );  
        errEff_cat_S2[ih]=errEff_cat_S2_1[ih]*errEff_cat_S2_2[ih];
        Err_Eff_cat_1[ih]=sqrt(errEff_cat_S1[ih] + errEff_cat_S2[ih]);

        //error without weight
        if(PutWeight==0){Err_Eff_cat_1[ih]= Eff_cat_1[ih]*TMath::Sqrt(gen_ptError[ih]*gen_ptError[ih]/(gen_pt[ih]*gen_pt[ih]) + rec_ptError[ih]*rec_ptError[ih]/(rec_pt[ih]* rec_pt[ih]));}
        cout<<"=================================== This is bin "<<ih<<"================================================="<<endl;
        cout<<"Jpsi Yield by integral of histo:  "<< diMuonsInvMass_RecA1[ih]->IntegralAndError(binlow, binhi,recError) <<"  error "<< rec_ptError[ih]<<endl; 
        cout<<"JpsiYield by CB yield determ:     "<< JpsiYield << " JpsiWidth "<< JpsiWidth<<" JpsiMass "<<JpsiMass <<" error "<< TMath::Sqrt((JpsiYield))<<endl;
        cout<<"Jpsi Yield by Function integral:  "<< GAUSPOL->Integral(MassLow,MassHigh)<<endl;
        cout<<" Bin "<<ih<< " Reco Jpsi  "<<  rec_pt[ih] <<" Gen Jpsi "<<gen_pt[ih]<<endl;
        cout<<" Eff "<< Eff_cat_1[ih]<<" Error "<<Err_Eff_cat_1[ih]<<endl;

        char Spectra[100];
        sprintf(Spectra,"Psi");

        dataFile<<"--------------------------------------------------------------------------------------------------------------------"<<endl;
        dataFile<<"Gen Jpsi "<< gen_pt[ih] <<" gen error "<<gen_ptError[ih]<<" rec Jpsi "<<rec_pt[ih] <<" rec error "<<rec_ptError[ih]<<"  "<<Spectra<<"["<<pt_bound[ih]<<" - "<<pt_bound[ih+1]<<"]"<<endl;
        dataFile<<"Jpsi eff "<< Eff_cat_1[ih]<<" error "<<Err_Eff_cat_1[ih]<<"  "<<Spectra<<"["<<pt_bound[ih]<<" - "<<pt_bound[ih+1]<<"]"<<endl;


        EffJpsiGen->SetBinContent(ih+1,gen_pt[ih]);
        EffJpsiGen->SetBinError(ih+1,gen_ptError[ih]);
        EffJpsiReco->SetBinContent(ih+1,rec_pt[ih]);
        EffJpsiReco->SetBinError(ih+1,rec_ptError[ih]);


        backfun_1->SetLineColor(4);
        backfun_1->SetLineWidth(1);
        //backfun_1->Draw("same");

        sprintf(GPlotName,"plots/DiMuonGenMass_Raps_%d_PsiBin_%d_%s.png",iCat,ih,cCd);
        sprintf(GPlotName1,"plots/DiMuonGenMass_Raps_%d_PsiBin_%d_%s.pdf",iCat,ih,cCd);

        //Drawing histo
        new TCanvas;
        diMuonsInvMass_GenA1[ih]->Draw("EPL");
        //gPad->Print(GPlotName);
        //gPad->Print(GPlotName1);
        //  if (iSpec==1){ new TCanvas; diMuonsPt_GenA1[ih]->Draw(); new TCanvas; diMuonsPt_RecA1[ih]->Draw();}
      }

      dataFile<<"Efficiency"<<"\t"<<"Error"<<endl;
      dataFile<<eff[0]<<"\t"<<effErr[0]<<"\t"<<eff[1]<<"\t"<<effErr[1]<<"\t"<<eff[2]<<"\t"<<effErr[2]<<"\t"<<eff[3]<<"\t"<<effErr[3]<<"\t"<<endl;
      dataFile.close();

      TFile *outfile;
      char tmp_output[512];
      if(iCat == 0) sprintf(tmp_output,"MCCT_PrJpsi_Raps_0012_dPhi_%s.root",cCd);
      if(iCat == 1) sprintf(tmp_output,"MCCT_PrJpsi_Raps_1216_dPhi_%s.root",cCd);
      if(iCat == 2) sprintf(tmp_output,"MCCT_PrJpsi_Raps_1624H_dPhi_%s.root",cCd);
      if(iCat == 3) sprintf(tmp_output,"MCCT_PrJpsi_Raps_1624L_dPhi_%s.root",cCd);
      outfile =new TFile(tmp_output,"Recreate");

      double gsum = 0.0;
      double gen_pt_Norm[4];
      double gen_ptError_Norm[4];
      for(int i = 0; i < Nptbin; i++){
        gsum += gen_pt[i];
      }
      gsum = gsum*(pt_bound[1]-pt_bound[0]);

      for(int i = 0; i < Nptbin; i++){
        gen_pt_Norm[i] = gen_pt[i]/gsum;
        gen_ptError_Norm[i] = gen_ptError[i]/gsum;
      }



      double sum = 0.0;
      double rec_pt_Norm[4];
      double rec_ptError_Norm[4];
      for(int i = 0; i < Nptbin; i++){
        sum += rec_pt[i];
      }
      sum = sum*(pt_bound[1]-pt_bound[0]);

      for(int i = 0; i < Nptbin; i++){
        rec_pt_Norm[i] = rec_pt[i]/sum;
        rec_ptError_Norm[i] = rec_ptError[i]/sum;
      }

      TH1F *hJpsi_Gen = new TH1F("hJpsi_Gen","hJpsi_Gen",4,0,TMath::Pi()/2);
      TH1F *hJpsi_Reco = new TH1F("hJpsi_Reco","hJpsi_Reco",4,0,TMath::Pi()/2);
      hJpsi_Gen->Sumw2();
      hJpsi_Reco->Sumw2();
      for(int i = 0; i < 4; i++){
        hJpsi_Gen->SetBinContent(i+1,gen_pt[i]);
        hJpsi_Gen->SetBinError(i+1,gen_ptError[i]);
        hJpsi_Reco->SetBinContent(i+1,rec_pt[i]);
        hJpsi_Reco->SetBinError(i+1,rec_ptError[i]);
      }

      TGraphErrors *Eff_Jpsi = new TGraphErrors(Nptbin, PT, Eff_cat_1, mom_err,Err_Eff_cat_1);
      TGraphErrors *Jpsi_Gen = new TGraphErrors(Nptbin, PT, gen_pt, mom_err,gen_ptError);
      TGraphErrors *Jpsi_Reco = new TGraphErrors(Nptbin, PT, rec_pt, mom_err,rec_ptError);
      TGraphErrors *Jpsi_Reco_Norm = new TGraphErrors(Nptbin, PT, rec_pt_Norm, mom_err, rec_ptError_Norm);
      TGraphErrors *Jpsi_Gen_Norm = new TGraphErrors(Nptbin, PT, gen_pt_Norm, mom_err, gen_ptError_Norm);
      Eff_Jpsi->SetMarkerStyle(20);
      Eff_Jpsi->SetMarkerColor(kBlack);
      Eff_Jpsi->SetMarkerSize(1.2);
      Eff_Jpsi->GetYaxis()->SetTitle("Reconstruction Efficiency");
      Eff_Jpsi->GetXaxis()->CenterTitle(true);;
      Eff_Jpsi->GetXaxis()->SetTitle("#phi^{J/#psi} - #Psi_{RP}");
      //if(iSpec==5) Eff_Jpsi->GetXaxis()->SetTitle("#phi^{J/#psi} - #Psi_{RP}"); // J/psi
      Eff_Jpsi->GetYaxis()->SetRangeUser(0,1.0);
      if(iCat == 0){
        Eff_Jpsi->SetName("Eff_Jpsi_Raps_0012");
        Jpsi_Gen->SetName("Gen_Jpsi_Raps_0012");
        Jpsi_Reco->SetName("Reco_Jpsi_Raps_0012");
        hJpsi_Gen->SetName("hGen_Jpsi_Raps_0012");
        hJpsi_Reco->SetName("hReco_Jpsi_Raps_0012");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Raps_0012");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Raps_0012");
      }
      if(iCat == 1){
        Eff_Jpsi->SetName("Eff_Jpsi_Raps_1216");
        Jpsi_Gen->SetName("Gen_Jpsi_Raps_1216");
        Jpsi_Reco->SetName("Reco_Jpsi_Raps_1216");
        hJpsi_Gen->SetName("hGen_Jpsi_Raps_1216");
        hJpsi_Reco->SetName("hReco_Jpsi_Raps_1216");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Raps_1216");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Raps_1216");
      }
      if(iCat == 2){
        Eff_Jpsi->SetName("Eff_Jpsi_Raps_1624H");
        Jpsi_Gen->SetName("Gen_Jpsi_Raps_1624H");
        Jpsi_Reco->SetName("Reco_Jpsi_Raps_1624H");
        hJpsi_Gen->SetName("hGen_Jpsi_Raps_1624H");
        hJpsi_Reco->SetName("hReco_Jpsi_Raps_1624H");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Raps_1624H");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Raps_1624H");
      }
      if(iCat == 3){
        Eff_Jpsi->SetName("Eff_Jpsi_Raps_1624L");
        Jpsi_Gen->SetName("Gen_Jpsi_Raps_1624L");
        Jpsi_Reco->SetName("Reco_Jpsi_Raps_1624L");
        hJpsi_Gen->SetName("hGen_Jpsi_Raps_1624L");
        hJpsi_Reco->SetName("hReco_Jpsi_Raps_1624L");
        Jpsi_Reco_Norm->SetName("nReco_Jpsi_Raps_1624L");
        Jpsi_Gen_Norm->SetName("nGen_Jpsi_Raps_1624L");
      }
      Eff_Jpsi->Write();
      Jpsi_Gen->Write();
      Jpsi_Reco->Write();
      hJpsi_Gen->Write();
      hJpsi_Reco->Write();
      Jpsi_Reco_Norm->Write();
      Jpsi_Gen_Norm->Write();


      double v1, v2;
      v1 = 0.0; v2 = TMath::Pi()/2;
      TF1 *stline = new TF1("stline","[0]",v1,v2);
      stline->SetLineColor(9);
      stline->SetLineWidth(2);

      //TLegend *legend_GP = new TLegend( 0.45,0.85,0.80,0.95);
      TLegend *legend_GP = new TLegend(0.20,0.80,0.54,0.90);
      legend_GP->SetBorderSize(0);
      legend_GP->SetFillStyle(0);
      legend_GP->SetFillColor(0);
      legend_GP->SetTextSize(0.032);
      //legend_GP->AddEntry(Eff_Jpsi,"PythiaEvtGen + HydjetBass", "PL");
      TLatex *lt1 = new TLatex();
      lt1->SetNDC();
      lt1->SetTextAlign(12);
      lt1->SetTextFont(42);
      lt1->SetTextSize(0.04);

      gStyle->SetOptFit(0);
      new TCanvas;
      Eff_Jpsi->Draw("APZ");
      Eff_Jpsi->Fit("stline","ER","",v1,v2);
      //Eff_Jpsi->Fit("stline","ER","",0,TMath::Pi());
      //Eff_Jpsi->Fit("stline","ER","",-TMath::Pi(),TMath::Pi());
      //Eff_Jpsi->Fit("stline","ER","",-TMath::Pi()/2,TMath::Pi()/2);
      double effNo = stline->GetParameter(0);
      double effErrNo = stline->GetParError(0);
      cout<<" efficiency : "<<effNo<<endl;
      char tmp[512], tmp1[512];
      sprintf(tmp,"PythiaEvtGen + HydjetBass : %.2f #pm %.2f %%", effNo*100, effErrNo*100);
      legend_GP->AddEntry(Eff_Jpsi,tmp, "PL");
      legend_GP->Draw("Same");

      if(iCat == 0) sprintf(tmp1,"0.0 < |y| < 1.2");
      if(iCat == 1) sprintf(tmp1,"1.2 < |y| < 1.6");
      if(iCat == 2) sprintf(tmp1,"1.6 < |y| < 2.4");
      if(iCat == 3) sprintf(tmp1,"1.6 < |y| < 2.4");
      lt1->DrawLatex(0.21,0.91,tmp1);

      //Eff_Jpsi->Write();

      char tmp_pdf[512], tmp_png[512];
      if(iCat == 0) {sprintf(tmp_pdf,"plots/Eff_Jpsi_Raps_0012_dPhi_%s.pdf",cCd);sprintf(tmp_png,"plots/Eff_Jpsi_Raps_0012_dPhi_%s.png",cCd);}
      if(iCat == 1) {sprintf(tmp_pdf,"plots/Eff_Jpsi_Raps_1216_dPhi_%s.pdf",cCd);sprintf(tmp_png,"plots/Eff_Jpsi_Raps_1216_dPhi_%s.png",cCd);}
      if(iCat == 2) {sprintf(tmp_pdf,"plots/Eff_Jpsi_Raps_1624H_dPhi_%s.pdf",cCd);sprintf(tmp_png,"plots/Eff_Jpsi_Raps_1624H_dPhi_%s.png",cCd);}
      if(iCat == 3) {sprintf(tmp_pdf,"plots/Eff_Jpsi_Raps_1624L_dPhi_%s.pdf",cCd);sprintf(tmp_png,"plots/Eff_Jpsi_Raps_1624L_dPhi_%s.png",cCd);}
      gPad->Print(tmp_pdf);
      gPad->Print(tmp_png);

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

// a : eff file, c : centrality, d : pT, e : rapidity, f[] : array of eff & err
void DoEffCor(TFile* a, int b, double c, double d, double e, double *f){
  char tmp[512];
  double bct[6] = {0, 4, 8, 12, 24, 40};
  double bpt[6] = {3.0, 4.5, 6.5, 8.0, 10.0, 40.0};
  int ict = 0;
  int ipt = 0;
  if(c >= 4 && c < 8) ict = 0; 
  if(c >= 8 && c < 12) ict = 1; 
  if(c >= 12 && c < 24) ict = 2; 

  if(d >= 3.0 && d < 4.5) ipt = 0;
  if(d >= 4.5 && d < 6.5) ipt = 1;
  if(d >= 6.5 && d < 8.0) ipt = 2;
  if(d >= 8.0 && d < 10.0) ipt = 3;
  if(d >= 10.0 && d < 40.0) ipt = 4;

  char icCd[512];
  if(b == 0) sprintf(icCd,"default");
  if(b == 1) sprintf(icCd,"cowboy");
  if(b == 2) sprintf(icCd,"sailor");

  sprintf(tmp,"eff1D_Cent_%1.f_%1.f_Pt_%0.1f_%0.1f_%s",bct[1],bct[4],bpt[ipt],bpt[ipt+1], icCd);
  TH1D *hEff = (TH1D*)a->Get(tmp);
  int sbin = 0;
  sbin = hEff->FindBin(e);
  f[0] = hEff->GetBinContent(sbin);
  f[1] = hEff->GetBinError(sbin);
  //cout<<"sbin : "<<sbin<<", eff : "<<f[0]<<", err : "<<f[1]<<endl;
  //return f;
}

/*
void DoEffCor(TFile* a, char *b, double c, double d, double e, double *f){
  TH3D *hEff = (TH3D*)a->Get(b);
  int sbin = 0;
  sbin = hEff->FindBin(c, d, e);
  f[0] = hEff->GetBinContent(sbin);
  f[1] = hEff->GetBinError(sbin);
  //cout<<"sbin : "<<sbin<<", eff : "<<f[0]<<", err : "<<f[1]<<endl;
  //return f;
}
*/

