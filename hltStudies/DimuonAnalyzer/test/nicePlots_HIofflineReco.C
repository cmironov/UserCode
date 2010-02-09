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
void nicePlots_HIofflineReco()
{
  //  gROOT->Macro("setStyle.C+");
  gROOT->ProcessLine(".x setStyle.C");
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.07,"");
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleW(1);
  

  char szBuf[256];
  // booleans
  bool fit_eff =false;

  bool bSavePlots = true;
  bool is_allPlots =true;
  bool is_PtYMass = true;
  bool is_event_info = true;

  bool is_upsilon = true;
  bool is_signal_only = false;
  const char* inputFileName     = "DimuonAnalyzer_Upsilon_embedding_500events";//dimuonGenEff_Z_500evts";
  const char* inputFileLocation = "root";//"/afs/cern.ch/user/s/silvest//public/html/reco_muons/reconstruction_workshop/";
  const char* outputFileLocation = "figs";

  int marker_color = 9;
  int marker_size = 1;
  int marker_style = 20;
  if(is_signal_only&&is_upsilon) 
    {
      int maxXMass2 = 20;
      int maxXMass1 = 10;
      int minXMass1 = 9;
      int minXMass2 = 0;
      int maxMass1 = 2010;
      int maxMass2 = 500;
      int nbinMass1 =50;
      int nbinMass2 =50;
      int maxPt = 300;
      int maxY =300;
      int maxNtrk =4;
    }
      else   if(is_signal_only&&!is_upsilon) 
	{
	  int maxXMass = 140;
	  int minXMass = 0;
	  int maxMass = 140;
	  int nbinMass1 =100;
	  int nbinMass2 =100;
	  int maxPt = 150;
	  int maxY =120;
	  int maxNtrk =4;
	}
    
      else   if(!is_signal_only&&is_upsilon) 
    	{
	  int minXMass1 = 9;
	  double maxXMass1 = 9.5;
	  int nbinMass1 =3;
	  int maxMass1 = 600;

	  int minXMass2 = 0;
	  int maxXMass2 = 20;
	  int nbinMass2 =20;
	  int maxMass2 = 200;
	  int maxPt = 70;
	  int maxY =70;
	  int maxNtrk =2200;
	}

      else
	{
	  int maxXMass2 = 140;
	  int maxXMass1 = 140;
	  int minXMass1 = 0;
	  int minXMass2 = 0;
	  int maxMass1 = 100;
	  int maxMass2 = 100;
	  int nbinMass1 =100;
	  int nbinMass2 =100;

	  int maxPt = 70;
	  int maxY =70;
	  int maxNtrk =2200;
	}	  
    
      

  // get the files
  //reco
 
  sprintf(szBuf,"%s/%s.root",inputFileLocation,inputFileName);
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

  // mass
  TH1 *phMGen = (TH1 *)inputFileRoot.Get("demo/phM_genDimuon")->Clone("phMGen");
  phMGen->SetDirectory(0);
  TH1 *phMReco = (TH1 *)inputFileRoot.Get("demo/phM_recoDimuon")->Clone("phMReco");
  phMReco->SetDirectory(0);



  phPtReco->Rebin(4);
  phPtGen->Rebin(4);
  phYGen->Rebin(5);
  phYReco->Rebin(5);
  

  // calculate efficiency 
  TGraphAsymmErrors* phPt_dimuEff = new TGraphAsymmErrors();
  phPt_dimuEff->BayesDivide(phPtReco,phPtGen,"");
  cout << "phPt_dimuEff done " <<endl;
  TGraphAsymmErrors* phY_dimuEff = new TGraphAsymmErrors(); 
  phY_dimuEff->BayesDivide(phYReco,phYGen,"");


 // set drawing axis
  TH1 *phYAxis = new TH1F("phYAxis",";y;efficiency(reco/gen)",25,-2.5,2.5);
  phYAxis->SetDirectory(0);
  phYAxis->SetMaximum(1);
 
  TH1 *phPtAxis  = new TH1F("phPtAxis",";p_{T} [GeV/c];efficiency(reco/gen)",20,0.,20.);
  phPtAxis->SetDirectory(0);
  phPtAxis->SetMaximum(1);

 
  // colors:
  phPt_dimuEff->SetMarkerStyle(marker_style);
  phPt_dimuEff->SetLineColor(marker_color); 
  phPt_dimuEff->SetLineWidth(marker_size);

  phY_dimuEff->SetMarkerStyle(marker_style);
  phY_dimuEff->SetLineColor(marker_color); 
  phY_dimuEff->SetLineWidth(marker_size);

   phPt_dimuEff->SetMarkerColor(marker_color);
  phY_dimuEff->SetMarkerColor(marker_color);

  sprintf(szBuf,"%s_eff_pt",inputFileName);
  TCanvas *pc_pt = new TCanvas(szBuf,szBuf);
  pc_pt->cd();
  if(fit_eff)phPt_dimuEff->Fit("pol0","+","",20.,100.);
  phPtAxis->Draw();
  phPt_dimuEff->Draw("Psame");
  pc_pt->Update();
    
  sprintf(szBuf,"%s_eff_y",inputFileName);
  TCanvas *pc_y = new TCanvas(szBuf,szBuf);
  pc_y->cd();
  if(fit_eff)phY_dimuEff->Fit("pol0","+","",-2.5,2.5);
  phYAxis->Draw();
  phY_dimuEff->Draw("Psame");
  pc_y->Update();
  
  if(is_PtYMass){
  // pt y and mass plots
  sprintf(szBuf,"%s_pt_y_mass",inputFileName);
  TCanvas *c1 = new TCanvas(szBuf,szBuf);
  c1->Divide(2,3);
 
  c1->cd(1);
  TH1 *phPtAxis2  = new TH1F("phPtAxis2",";p_{T} [GeV/c]",40,0.,20.);
  phPtAxis2->SetDirectory(0);
  phPtAxis2->SetMinimum(0);
  phPtAxis2->SetMaximum(maxPt);
  phPtAxis2->SetTitle("GEN");
  phPtAxis2->Draw();
  phPtGen->SetMarkerStyle(marker_style);
  phPtGen->SetMarkerColor(marker_color);
  phPtGen->SetMarkerSize(marker_size);
  phPtGen->SetLineColor(marker_color); 
  phPtGen->SetLineWidth(marker_size);
  phPtGen->Draw("same");
  c1->Update();
  c1->cd(2);
  TH1 *phPtAxis2  = new TH1F("phPtAxis2",";p_{T} [GeV/c]",40,0.,20.);
  // phPtAxis2->SetDirectory(0);
  phPtAxis2->SetMinimum(0);
  phPtAxis2->SetMaximum(maxPt);

  phPtAxis2->SetTitle("RECO");
  phPtAxis2->Draw();
  phPtReco->SetMarkerStyle(marker_style);
  phPtReco->SetMarkerColor(marker_color);
  phPtReco->SetMarkerSize(marker_size);
  phPtReco->SetLineColor(marker_color); 
  phPtReco->SetLineWidth(marker_size);
  phPtReco->Draw("same");
  c1->Update();

  c1->cd(3);
  TH1 *phYAxis2 = new TH1F("phYAxis2",";y",30,-2.5,2.5);
  phYAxis2->SetDirectory(0);
  phYAxis2->SetMinimum(0);
  phYAxis2->SetMaximum(maxY);
  phYAxis2->SetTitle("GEN");
  phYAxis2->Draw();
  phYGen->SetTitle("GEN");
  phYGen->SetMarkerStyle(marker_style);
  phYGen->SetMarkerColor(marker_color);
  phYGen->SetMarkerSize(marker_size);
  phYGen->SetLineColor(marker_color); 
  phYGen->SetLineWidth(marker_size);
  phYGen->Draw("same");
  c1->Update();
  c1->cd(4);
  TH1 *phYAxis2 = new TH1F("phYAxis2",";y",30,-2.5,2.5);
  phYAxis2->SetDirectory(0);
  phYAxis2->SetMinimum(0);
  phYAxis2->SetMaximum(maxY);
  phYAxis2->SetTitle("RECO");
  phYAxis2->Draw();
  phYReco->SetMarkerStyle(marker_style);
  phYReco->SetMarkerColor(marker_color);
  phYReco->SetMarkerSize(marker_size);
  phYReco->SetLineColor(marker_color); 
  phYReco->SetLineWidth(marker_size);
  phYReco->Draw("sames");
  c1->Update();

  c1->cd(5);
  TH1 *phMAxis2  = new TH1F("phMAxis2",";mass [GeV/c^{2}]",nbinMass1,minXMass1,maxXMass1);
  phMAxis2->SetDirectory(0);
  phMAxis2->SetMaximum(maxMass1);
  phMAxis2->SetTitle("GEN");
  phMAxis2->Draw();
  phMGen->SetMarkerStyle(marker_style);
  phMGen->SetMarkerColor(marker_color);
  phMGen->SetMarkerSize(marker_size);
  phMGen->SetLineColor(marker_color); 
  phMGen->SetLineWidth(marker_size);
  phMGen->Draw("Psame");
  c1->Update();
  c1->cd(6);
  TH1 *phMAxis2  = new TH1F("phMAxis2",";mass [GeV/c^{2}]",nbinMass2,minXMass2,maxXMass2);
  phMAxis2->SetDirectory(0);
  phMAxis2->SetTitle("RECO");
  phMAxis2->SetMaximum(maxMass2);
  phMAxis2->Draw();
  phMReco->SetMarkerStyle(marker_style);
  phMReco->SetMarkerColor(marker_color);
  phMReco->SetMarkerSize(marker_size);
  phMReco->SetLineColor(marker_color); 
  phMReco->SetLineWidth(marker_size);
  phMReco->Draw("Psame");
  c1->Update();
  }
  if(is_allPlots){
  sprintf(szBuf,"%s_allplots",inputFileName);
  TCanvas *c2 = new TCanvas(szBuf,szBuf, 900, 300);
  c2->Divide(3,1);
  // pt y gen Dimuon All
  TH2 *phPtYGenAll = (TH2 *)inputFileRoot.Get("demo/phPtY_genDimuonAll")->Clone("phPtYGenAll");
  phPtYGenAll->SetDirectory(0);
  TH2 *phPtYAxis  = new TH2F("phPtYAxis",";p_{T} [GeV/c];y",40,0.,20.,50,-2.5,2.5);
  phPtYAxis->SetDirectory(0);
  c2->cd(1);
  phPtYAxis->SetTitle("GEN dimuon All");
  phPtYAxis->Draw();
  phPtYGenAll->Draw("Psame");
  c2->Update();

  // pt y roecp Dimuon All
  TH2 *phPtYRecoAll = (TH2 *)inputFileRoot.Get("demo/phPtY_recoDimuonAll")->Clone("phPtYRecoAll");
  phPtYRecoAll->SetDirectory(0);
  c2->cd(2);
  phPtYAxis->SetTitle("RECO dimuon All");
  phPtYAxis->Draw();
  phPtYRecoAll->Draw("Psame");
  c2->Update();
  // pt y roecp Dimuon All
  TH2 *phPtEtaRecoTrack = (TH2 *)inputFileRoot.Get("demo/phPtEta_recoTrack")->Clone("phPtEtaRecoTrack");
  phPtEtaRecoTrack->SetDirectory(0);
  TH2 *phPtEtaAxis  = new TH2F("phPtYAxis",";p_{T} [GeV/c];#eta",40,0.,20.,50,-2.5,2.5);
  phPtEtaAxis->SetDirectory(0);
  c2->cd(3);
  phPtEtaAxis->SetTitle("RECO tracks All");
  phPtEtaAxis->Draw();
  phPtEtaRecoTrack->Draw("Psame");
  c2->Update();
  }
  if(is_event_info)
  {
    // ntracks, nmuon, b ntrk:nmu:b

    // retreive ntuple
    TNtuple *nEventInfo = (TNtuple *)inputFileRoot.Get("demo/pnEventInfo")->Clone("nEventInfo");
    
    // canvas
    sprintf(szBuf,"%s_EventInfo",inputFileName);
    TCanvas *c3 = new TCanvas(szBuf,szBuf, 900, 300);
    c3->Divide(3,1);

    c3->cd(1);
    TH1F *h_ntrk= new TH1F("h_ntrk",";number of tracks per event",50,0,maxNtrk);
    nEventInfo->Project("h_ntrk","ntrk");
    h_ntrk->SetMarkerStyle(marker_style);
    h_ntrk->SetMarkerColor(marker_color);
    h_ntrk->SetMarkerSize(marker_size);
    h_ntrk->SetLineColor(marker_color); 
    h_ntrk->SetLineWidth(marker_size);
    h_ntrk->SetNdivisions( 5, "X" );
    h_ntrk->Draw();
    c3->Update();
    
    c3->cd(2);
    TH1 *h_nmuon  = new TH1F("h_nmuon",";number of muon per event",10,0.,5.);
    nEventInfo->Project("h_nmuon","nmu");
    h_nmuon->SetMarkerStyle(marker_style);
    h_nmuon->SetMarkerColor(marker_color);
    h_nmuon->SetMarkerSize(marker_size);
    h_nmuon->SetLineColor(marker_color); 
    h_nmuon->SetLineWidth(marker_size);
    h_nmuon->SetNdivisions( 5, "X" );
    h_nmuon->Draw();
    c3->Update();
 
    c3->cd(3);
    TH1 *h_b  = new TH1F("h_b",";impact parameter (b)",40,0.,18.);
    nEventInfo->Project("h_b","b");
    h_b->SetMarkerStyle(marker_style);
    h_b->SetMarkerColor(marker_color);
    h_b->SetMarkerSize(marker_size);
    h_b->SetLineColor(marker_color); 
    h_b->SetLineWidth(marker_size);
    h_b->SetNdivisions( 5, "X" );
    h_b->Draw();
    c3->Update();
  }

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
  
      sprintf(szBuf,"%s/%s",outputFileLocation,c1->GetTitle());
      TString outFileBase1(szBuf);
      outFileGif =  outFileBase1+".gif";
      if(is_PtYMass)c1->Print(outFileGif.Data(),"gifLandscape");

      sprintf(szBuf,"%s/%s",outputFileLocation,c2->GetTitle());
      TString outFileBase1(szBuf);
      outFileGif =  outFileBase1+".gif";
      if(is_allPlots)     c2->Print(outFileGif.Data(),"gifLandscape");

      sprintf(szBuf,"%s/%s",outputFileLocation,c3->GetTitle());
      TString outFileBase1(szBuf);
      outFileGif =  outFileBase1+".gif";
      if(is_event_info)     c3->Print(outFileGif.Data(),"gifLandscape");
    }
 
}
  
