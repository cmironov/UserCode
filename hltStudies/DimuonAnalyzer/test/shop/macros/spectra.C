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
void spectra(const char* pDimu="z0",
	   const char* pData="mix",
	   const char* pB="minbias",
	   const char* pEta="full",
	   const char* pInFilePath="/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/Correlation/DijetCorrelation/test/lxjobs/shop",
	   const char* pOutFilePath="figs",
	   Int_t bodytype=2,
	   Bool_t bSavePlots=false
	   )
{
  gROOT->Macro("/afs/cern.ch/user/m/mironov/utilities/setStyle.C+");
  // gStyle->SetPalette(1);
  char szBuf[256];


  sprintf(szBuf,"%s",pOutFilePath);
  TString outFileLocation(szBuf);
  cout<<"outfileLocation "<<outFileLocation<<endl;

  // in file name
  sprintf(szBuf,"%s/%s_%s_offline_%s_STA_%s_hists.root",pInFilePath,pDimu,pData,pB,pEta);
  TString fileName1(szBuf);
  cout << "STA file ..." << fileName1 <<endl;
  TFile inFileSta(fileName1);

  sprintf(szBuf,"%s/%s_%s_offline_%s_GLB_%s_hists.root",pInFilePath,pDimu,pData,pB,pEta);
  TString fileName2(szBuf);
  cout << "GLB file ..." << fileName2 <<endl;
  TFile inFileGlb(fileName2);
 
  //----- histograms
  // directory of the histogram 
  sprintf(szBuf,"plot/phPtEta_trk");
  if(bodytype==2) sprintf(szBuf,"plot/phPtY_pair");
  TString baseHistName(szBuf);
  TString histSimName  = baseHistName+"Sim"; 
  TString histRecoName = baseHistName+"Reco"; 
  TString histFakeName = baseHistName+"RecoFake"; 

  cout<<"Histograms name are : "<< histSimName <<"\t"<< histRecoName<<"\t"<<histFakeName<<endl;

  // sta hist
  TH2 *phSimSta  = (TH2 *)(inFileSta.Get(histSimName)->Clone("phSimSta"));
  phSimSta->SetDirectory(0);
  TH1 *phPtSimSta = (TH1D*)phSimSta->ProjectionX();
  phPtSimSta->SetDirectory(0);
  TH1 *phThSimSta = (TH1D*)phSimSta->ProjectionY();
  phThSimSta->SetDirectory(0);

  TH2 *phRecoSta = (TH2 *)inFileSta.Get(histRecoName)->Clone("phRecoSta");
  phRecoSta->SetDirectory(0);
  TH1 *phPtRecoSta = (TH1D*)phRecoSta->ProjectionX();
  phPtRecoSta->SetDirectory(0);
  TH1 *phThRecoSta = (TH1D*)phRecoSta->ProjectionY();
  phThRecoSta->SetDirectory(0);

   
  TH2 *phFakeSta = (TH2 *)inFileSta.Get(histFakeName)->Clone("phFakeSta");
  phFakeSta->SetDirectory(0);
  TH1 *phPtFakeSta = (TH1D*)phFakeSta->ProjectionX();
  phPtFakeSta->SetDirectory(0);
  TH1 *phThFakeSta = (TH1D*)phFakeSta->ProjectionY();
  phThFakeSta->SetDirectory(0);

  // glb hists
  TH2 *phSimGlb = (TH2 *)(inFileGlb.Get(histSimName)->Clone("phSimGlb"));
  phSimGlb->SetDirectory(0);
  TH1 *phPtSimGlb = (TH1D*)phSimGlb->ProjectionX();
  phPtSimGlb->SetDirectory(0);
  TH1 *phThSimGlb = (TH1D*)phSimGlb->ProjectionY();
  phThSimGlb->SetDirectory(0);

  TH2 *phRecoGlb = (TH2 *)inFileGlb.Get(histRecoName)->Clone("phRecoGlb");
  phRecoGlb->SetDirectory(0);
  TH1 *phPtRecoGlb = (TH1D*)phRecoGlb->ProjectionX();
  phPtRecoGlb->SetDirectory(0);
  TH1 *phThRecoGlb = (TH1D*)phRecoGlb->ProjectionY();
  phThRecoGlb->SetDirectory(0);
   
  TH2 *phFakeGlb = (TH2 *)inFileGlb.Get(histFakeName)->Clone("phFakeGlb");
  phFakeGlb->SetDirectory(0);
  TH1 *phPtFakeGlb = (TH1D*)phFakeGlb->ProjectionX();
  phPtFakeGlb->SetDirectory(0);
  TH1 *phThFakeGlb = (TH1D*)phFakeGlb->ProjectionY();
  phThFakeGlb->SetDirectory(0);


 // ------------ axis
  TH1 *phPtAxis = new TH1F("phPtAxis",";p_{T};entries",1,0.,30.);
  phPtAxis->SetDirectory(0);
  phPtAxis->SetMinimum(1.);
  phPtAxis->SetMaximum(600); 

  TH1 *phPtAxisZ = new TH1F("phPtAxisZ",";p_{T};entries",1,0.,100.);
  phPtAxisZ->SetDirectory(0);
  phPtAxisZ->SetMinimum(1.);
  phPtAxisZ->SetMaximum(600); 
   
  TH1 *phYAxis = new TH1F("phYAxis",";y;entries",1,-3.5,3.5);
  phYAxis->SetDirectory(0);
  phYAxis->SetMinimum(1.);
  phYAxis->SetMaximum(400);  
       
  TH1 *phEtaAxis = new TH1F("phEtaAxis",";#eta;entries",1,-3.5,3.5);
  phEtaAxis->SetDirectory(0);
  phEtaAxis->SetMinimum(1.);
  phEtaAxis->SetMaximum(400.);

  TH1 *phPtPairAxis  = new TH1F("phPtPairAxis",";p_{T} [GeV/c]; entries",1,0.,20.);
  phPtPairAxis->SetDirectory(0);
  phPtPairAxis->SetMaximum(400.);
  phPtPairAxis->SetMinimum(1);

 sprintf(szBuf,"%s",pData);
 TString what(szBuf);
 
 if( what.CompareTo("sgn")==0 )
   { 
     phPtAxis->SetMaximum(1000); 
     phPtAxisZ->SetMaximum(1000);
     phYAxis->SetMaximum(700); 
     phEtaAxis->SetMaximum(1000.);
     phPtPairAxis->SetMaximum(1000);
   }


  //-------------  set marker type
  //pt
  phPtSimSta->SetMarkerStyle(20);
  phPtRecoSta->SetMarkerStyle(21);
  phPtFakeSta->SetMarkerStyle(22);

  phPtSimGlb->SetMarkerStyle(20);
  phPtRecoGlb->SetMarkerStyle(21); 
  phPtFakeGlb->SetMarkerStyle(22); 

  //eta/y
  phThSimSta->SetMarkerStyle(20);
  phThRecoSta->SetMarkerStyle(21);
  phThFakeSta->SetMarkerStyle(22);

  phThSimGlb->SetMarkerStyle(20);
  phThRecoGlb->SetMarkerStyle(21); 
  phThFakeGlb->SetMarkerStyle(22);

  // color
  // phPtSimSta->SetMarkerColor(20);
  phPtRecoSta->SetMarkerColor(9);
  phPtFakeSta->SetMarkerColor(2);

  // phPtSimGlb->SetMarkerColor(20);
  phPtRecoGlb->SetMarkerColor(9); 
  phPtFakeGlb->SetMarkerColor(2); 

  //eta/y
  // phThSimSta->SetMarkerColor(20);
  phThRecoSta->SetMarkerColor(9);
  phThFakeSta->SetMarkerColor(2);

  //  phThSimGlb->SetMarkerColor(20);
  phThRecoGlb->SetMarkerColor(9); 
  phThFakeGlb->SetMarkerColor(2);

  TLatex lx;
  lx.SetTextColor(2);
  lx.SetTextSize(0.03);
  //---------
  TLegend *pl=new TLegend(0.8,0.8,0.95,0.95);
  pl->SetFillColor(10);
  pl->SetBorderSize(0);
  //--------
  sprintf(szBuf,"pc_ptSpectra_%s_%s_offline_%s_%s_%d",pDimu,pData,pB,pEta,bodytype);
  TCanvas *pc_pt = new TCanvas(szBuf,szBuf,900,500);
  pc_pt->Divide(2,1);
  pc_pt->cd(1);
  //sta
  if(bodytype==1) 
    {
      sprintf(szBuf,"%s",pDimu);
      TString with(szBuf);
      cout<<"###### "<< with<<endl;
      if( with.CompareTo("z0")==0 ){ phPtAxisZ->Draw();}
      else {phPtAxis->Draw();}
    }
  else phPtPairAxis->Draw();
  // gPad->SetLogy();

  phPtSimSta->Rebin();
  phPtRecoSta->Rebin();
  phPtFakeSta->Rebin();

  phPtSimSta->Draw("Psame");
  phPtRecoSta->Draw("Psame");
  phPtFakeSta->Draw("Psame");

  if(bodytype==1) sprintf(szBuf,"STA muon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  else sprintf(szBuf,"STA dimuon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(1.,560.,szBuf);
 
  pl->AddEntry(phPtSimSta,"sim","P");
  pl->AddEntry(phPtRecoSta,"reco","P");
  pl->AddEntry(phPtFakeSta,"fake","P");
  pl->Draw("same");
  //glb
  pc_pt->cd(2);
 
  if(bodytype==1) 
    {
      sprintf(szBuf,"%s",pDimu);
      TString with(szBuf);
      if(with.CompareTo("z0")==0) {phPtAxisZ->Draw();}
      else {phPtAxis->Draw();}
    }
  else phPtPairAxis->Draw();

  phPtSimGlb->Rebin();
  phPtRecoGlb->Rebin();
  phPtFakeGlb->Rebin();

  phPtSimGlb->Draw("Psame");
  phPtRecoGlb->Draw("Psame");
  phPtFakeGlb->Draw("Psame");

  pl->Draw("same");

  if(bodytype==1) sprintf(szBuf,"GLB muon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  else sprintf(szBuf,"GLB dimuon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(1.,560.,szBuf);
  pc_pt->Update();

  //___________________________________________________________________
  TLegend *pl2=new TLegend(0.8,0.8,0.95,0.95);
  pl2->SetFillColor(10);
  pl2->SetBorderSize(0);
  //eta plots
  sprintf(szBuf,"pc_etaSpectra_%s_%s_offline_%s_%s_%d",pDimu,pData,pB,pEta,bodytype);
  TCanvas *pc_eta = new TCanvas(szBuf,szBuf,900,500);
  pc_eta->Divide(2,1);
  pc_eta->cd(1);

  //sta
  if(bodytype==1) {phEtaAxis->Draw();}
  else phYAxis->Draw();

  phThSimSta->Draw("Psame");
  phThRecoSta->Draw("Psame");
  phThFakeSta->Draw("Psame");

  pl2->AddEntry(phThSimSta,"sim","P");
  pl2->AddEntry(phThRecoSta,"reco","P");
  pl2->AddEntry(phThFakeSta,"fake","P");
  pl2->Draw("same");

  
  if(bodytype==1) sprintf(szBuf,"STA muon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  else sprintf(szBuf,"STA dimuon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(-3.,370.,szBuf);
  pc_eta->Update();
  //glb
  pc_eta->cd(2);
 
  if(bodytype==1) {phEtaAxis->Draw();}
  else phYAxis->Draw();

  phThSimGlb->Draw("Psame");
  phThRecoGlb->Draw("Psame");
  phThFakeGlb->Draw("Psame");

  pl2->Draw("same");

  if(bodytype==1) sprintf(szBuf,"GLB muon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  else sprintf(szBuf,"GLB dimuon: %s_%s_offline_%s_%s",pDimu,pData,pB,pEta);
  lx.SetTextColor(2);
  lx.DrawLatex(-3.,370,szBuf);
  pc_eta->Update();


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
 
