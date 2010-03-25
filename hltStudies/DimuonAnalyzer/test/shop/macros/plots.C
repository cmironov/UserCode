#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>

// miscellaneous  
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation/test/headers/utilities.h"
//#include "/Users/eusmartass/Software/DijetCorrelation/test/headers/utilities.h"
#endif


//________________________________________________________________
void plots(const char* pDimu="jpsi",
	   const char* pData="mix",
	   const char* pMuon="GLB",
	   const char* pHisto="Efficiency",
	   const char* pB="central",
	   const char* pEta="barrel",
	   const char* pInFilePath="/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation/test/lxjobs/shop",
	   const char* pOutFilePath="figs",
	   Int_t bodytype=2,
	   Bool_t bSavePlots=false
	   )
{
  gROOT->Macro("/afs/cern.ch/user/m/mironov/utilities/setStyle.C+");
  // gStyle->SetPalette(1);
  bool dofit =true;
  bool printIntegratedFit = true;

  char szBuf[256];
 
  sprintf(szBuf,"%s",pOutFilePath);
  TString outFileLocation(szBuf);
  cout<<"outfileLocation "<<outFileLocation<<endl;

  // in file name
  sprintf(szBuf,"%s/%s_%s_offline_%s_%s_%s_hists.root",pInFilePath,pDimu,pData,pB,pMuon,pEta);
  TString fileName(szBuf);
  cout << "In file ..." << fileName <<endl;
  TFile inFile(fileName);
  
 // set the drawing
  TH1 *phPtAxis = new TH1F("phPtAxis",";p_{T};Efficiency",1,0.,20.);
  phPtAxis->SetDirectory(0);
  phPtAxis->SetMinimum(0.);
  phPtAxis->SetMaximum(1.2);  
   
  TH1 *phYAxis = new TH1F("phYAxis",";y;Efficiency",1,-3.5,3.5);
  phYAxis->SetDirectory(0);
  phYAxis->SetMinimum(0.);
  phYAxis->SetMaximum(1.2);  
       
  TH1 *phEtaAxis = new TH1F("phEtaAxis",";#eta;Efficiency",1,-3.5,3.5);
  phEtaAxis->SetDirectory(0);
  phEtaAxis->SetMaximum(1.2);
  phEtaAxis->SetMinimum(0.);
 
  TH1 *phPtPairAxis  = new TH1F("phPtPairAxis",";p_{T} [GeV/c]; Efficiency",1,0.,20.);
  phPtPairAxis->SetDirectory(0);
  phPtPairAxis->SetMaximum(1.2);
  phPtPairAxis->SetMinimum(0.);

 
  //----- histogram names
  // directory of the histogram 
  sprintf(szBuf,"plot/phPtEta_trk");
  if(bodytype==2) sprintf(szBuf,"plot/phPtY_pair");
  TString baseHistName(szBuf);
  TString histUpName; 
  TString histDownName; 

  sprintf(szBuf,"%s",pHisto);
  TString with(szBuf);
  TString hist;
  if ( with.CompareTo("Efficiency") == 0) 
    {
      histUpName = baseHistName+"SimReco";
      histDownName = baseHistName+"Sim";
    }
  if ( with.CompareTo("Fake")       == 0) 
    {
      histUpName = baseHistName+"RecoFake";
      histDownName = baseHistName+"Reco";
      
      phEtaAxis->GetYaxis()->SetTitle("Fake rate");
      phPtAxis->GetYaxis()->SetTitle("Fake rate");
      phEtaAxis->SetMaximum(1.);
      phPtAxis->SetMaximum(1.);
      phYAxis->GetYaxis()->SetTitle("Fake rate");
      phPtPairAxis->GetYaxis()->SetTitle("Fake rate");
      phPtPairAxis->SetMaximum(1.);
      phYAxis->SetMaximum(1.);
      
    }
  cout<<"Histograms name are : "<< histUpName <<"\t " << histDownName<< endl;


  //---------------------------------------------------------------------------------
  TH2 *phUp = (TH2 *)(inFile.Get(histUpName)->Clone("phUp"));
  phUp->SetDirectory(0);
	
  TH2 *phDown = (TH2 *)inFile.Get(histDownName)->Clone("phDown");
  phDown->SetDirectory(0);
   

  TGraphAsymmErrors* phTheta = new TGraphAsymmErrors(); 
  TGraphAsymmErrors* phPt    = new TGraphAsymmErrors(); 

  phPt    = (TGraphAsymmErrors*)divideBayesian(phUp,phDown,1);
  phTheta = (TGraphAsymmErrors*)divideBayesian(phUp,phDown,2);

  //-----------

  Double_t ndown = phDown->GetEntries();
  Double_t nup   = phUp->GetEntries();
  Double_t ndown_err = sqrt(ndown);
  Double_t nup_err   = sqrt(nup);
  cout << "=====" << ndown << " +- err " << ndown_err << endl;
  cout << "=====" << nup << " +- err " << nup_err << endl;
     
  double ratio     = nup/ndown;
  double ratio_err = ratio * sqrt( pow(nup_err/nup,2) + pow(ndown_err/ndown,2) );
  //-----------
  TLatex lx;
  lx.SetTextColor(2);
  lx.SetTextSize(0.03);

  sprintf(szBuf,"pc_pt%s_%s_%s_offline_%s_%s_%s_%d",pHisto,pDimu,pData,pB,pMuon,pEta,bodytype);
  TCanvas *pc_pt = new TCanvas(szBuf,szBuf);
  pc_pt->cd();

  Float_t fitmax = 20;
  if(bodytype==1) {phPtAxis->Draw(); fitmax=60;}
  else phPtPairAxis->Draw();
 
  if(dofit && bodytype==2) phPt->Fit("pol0","RQ","",0.,fitmax);
  if(dofit && bodytype==1) phPt->Fit("pol0","RQ","",5.,fitmax);

  phPt->Draw("Psame");
  if(printIntegratedFit)
   {
     cout << "=====ratio: " <<ratio << " #pm " << ratio_err << endl;
     sprintf(szBuf, "Integrated eff: %.0f%%#pm%.0f",ratio*100,ratio_err*100);
     lx.SetTextColor(9);
     lx.DrawLatex(0.,1.,szBuf);
   }

  if(bodytype==1) sprintf(szBuf,"muon: %s_%s_offline_%s_%s_%s",pDimu,pData,pB,pMuon,pEta);
  else sprintf(szBuf,"dimuon: %s_%s_offline_%s_%s_%s",pDimu,pData,pB,pMuon,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(0.,0.1,szBuf);
  pc_pt->Update();

  //___________________________________________________________________
  //eta plots
  sprintf(szBuf,"pc_eta%s_%s_%s_offline_%s_%s_%s_%d",pHisto,pDimu,pData,pB,pMuon,pEta,bodytype);
  TCanvas *pc_eta = new TCanvas(szBuf,szBuf);
  pc_eta->cd();

  if(bodytype==1) {phEtaAxis->Draw();}
  else phYAxis->Draw();

  
  sprintf(szBuf,"%s",pEta);
  TString where(szBuf);
  if(dofit){
    if ( where.CompareTo("barrel") == 0) {phTheta->Fit("pol0","RQ","",-0.8,0.8);}
    else {phTheta->Fit("pol0","RQ","",-2.4,2.4);}
  }
  phTheta->Draw("Psame");
    
  if(printIntegratedFit)
    {
     cout << "=====ratio: " <<ratio << " #pm " << ratio_err << endl;
     sprintf(szBuf, "Integrated eff: %.0f%%#pm%.0f",ratio*100,ratio_err*100);
     lx.SetTextColor(9);
     lx.DrawLatex(0.,1.,szBuf);
   }

  if(bodytype==1) sprintf(szBuf,"muon: %s_%s_offline_%s_%s_%s",pDimu,pData,pB,pMuon,pEta);
  else sprintf(szBuf,"dimuon: %s_%s_offline_%s_%s_%s",pDimu,pData,pB,pMuon,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(-2.,0.1,szBuf);
  pc_pt->Update();

  //________________________________________________________________
  if(bSavePlots)
    {
      sprintf(szBuf,"%s/%s",pOutFilePath,pc_pt->GetTitle());
      TString outFileBase(szBuf);
      TString outFileGif = outFileBase+".gif";
      pc_pt->Print(outFileGif.Data(),"gifLandscape");
     
   
      sprintf(szBuf,"%s/%s",pOutFilePath,pc_eta->GetTitle());
      TString outFileBase1(szBuf);
      outFileGif =  outFileBase1+".gif";
      pc_eta->Print(outFileGif.Data(),"gifLandscape");
    }
}
 
