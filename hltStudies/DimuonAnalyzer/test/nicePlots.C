#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>
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
#include <TGraphAsymmErrors.h>

// miscellaneous  
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>

#endif

//________________________________________________________________
void nicePlots()
{
  gROOT->Macro("setStyle.C+");

  // get the files
  //reco
  char szBuf[256];
  bool bSavePlots = false;
  
  const char* inputFileName     = "dimuonGenEff.root";
  const char* inputFileLocation = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/UserCode/DimuonAnalyzer/test";
  const char* outputFileLocation = "/afs/cern.ch/user/m/mironov/scratch0/CMSSW_3_4_0/src/UserCode/DimuonAnalyzer/test/figs";
 
  sprintf(szBuf,"%s/%s",inputFileLocation,inputFileName);
  TFile inputFileRoot(szBuf);

  // get the histogram 

  //pt
  TH1 *phPtGen = (TH1 *)inputFileRoot.Get("demo/phPt_genDimuon")->Clone("phPtGen");
  phPtGen->SetDirectory(0);
  TH1 *phPtReco = (TH1 *)inputFileRoot.Get("demo/phPt_recoDimuon")->Clone("phPtReco");
  phPtReco->SetDirectory(0);

  //y
  TH1 *phYGen = (TH1 *)inputFileRoot.Get("demo/phY_genDimuon")->Clone("phYGen");
  phYGen->SetDirectory(0);
  TH1 *phYReco = (TH1 *)inputFileRoot.Get("demo/phY_recoDimuon")->Clone("phYReco");
  phYReco->SetDirectory(0);

  // calculate efficiency 
  TGraphAsymmErrors* phPt_dimuEff = new TGraphAsymmErrors();
  phPt_dimuEff->BayesDivide(phPtReco,phPtGen,"");
  TGraphAsymmErrors* phY_dimuEff = new TGraphAsymmErrors(); 
  phY_dimuEff->BayesDivide(phYReco,phYGen,"");


 // set drawing axis
  TH1 *phYAxis = new TH1F("phYAxis",";y;reco/gen",50,-2.5,2.5);
  phYAxis->SetDirectory(0);
  phYAxis->SetMaximum(2.);
 
  TH1 *phPtAxis  = new TH1F("phPtAxis",";p_{T} [GeV/c]; reco/gen",200,0.,100.);
  phPtAxis->SetDirectory(0);
  phPtAxis->SetMaximum(2.);

 
  // colors:
  phPt_dimuEff->SetMarkerStyle(22);
  phY_dimuEff->SetMarkerStyle(22);
  phPt_dimuEff->SetMarkerColor(6);
  phY_dimuEff->SetMarkerColor(6);

  sprintf(szBuf,"pc_pt_%s",inputFileName);
  TCanvas *pc_pt = new TCanvas(szBuf,szBuf);
  pc_pt->cd();
  phPt_dimuEff->Fit("pol0","+","",20.,100.);
  phPtAxis->Draw();
  phPt_dimuEff->Draw("Psame");
  pc_pt->Update();
    
  sprintf(szBuf,"pc_y_%s",inputFileName);
  TCanvas *pc_y = new TCanvas(szBuf,szBuf);
  pc_y->cd();
  phY_dimuEff->Fit("pol0","+","",-2.5,2.5);
  phYAxis->Draw();
  phY_dimuEff->Draw("Psame");
  pc_y->Update();
  
 
  if(bSavePlots)
    {
      sprintf(szBuf,"%s/%s",outputFileLocation,pc_pt->GetTitle());
      TString outFileBase(szBuf);
      TString outFileGif = outFileBase+".gif";
      pc_pt->Print(outFileGif.Data(),"gifLandscape");
     
      sprintf(szBuf,"%s/%s",outputFileLocation,pc_y->GetTitle());
      TString outFileBase1(szBuf);
      outFileGif =  outFileBase1+".gif";
      pc_y->Print(outFileGif.Data(),"gifLandscape");
  
    }
 
}
  
