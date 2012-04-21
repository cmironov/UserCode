#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TInterpreter.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>

#endif

void v2SummaryPlots_cent_EPCorr()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  //gROOT->Macro("./rootlogon.C");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","Pr","NBkg"};
  const char* legend[5]   = {"","Inclusive J/#psi","Non-prompt J/#psi", "Prompt J/#psi","Background"};
  int choseSignal         = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominal_bit1_cent_bkg"};


  // options
  bool bSavePlots      = true;
 
  bool bAddBkg         = true;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bAddTrigger     = false; 
  bool bAddAutoCor     = false;
  bool bDoCowboySailor = false;

  bool bDoEffCorrection= false;

  bool bDoCS               = false;
  bool bDoEffCorrection_cs = false;
  bool bDoMuHits_cs        = false;
  bool bDoBarrel_cs        = false;
  bool bDoVtx_cs           = false;
  bool doBkg_cs            = false;

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  double ptIntegrated[2]    = {6.5, 40};  

  const int ncentbins1                    = 4;
  double ncoll1[ncentbins1]                = {381.3,   329.4,  224.3,   89.9};
  // centrlaity bins: 0-10, 10-20, 20-30, 30-60
  const int ncentbins                    = 4;
  double ncoll[ncentbins]                = {355.4, 261.4178, 187.1470, 89.9};
  double ncoll_err[ncentbins]            = { 0.,     0.,       0.,      0.};

 
  TFile *f0 = new TFile(Form("./a1_corrV2.root"));
  if(!f0->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}

  TFile *f1 = new TFile(Form("./extrastud/a1_corrV2.root"));
  if(!f1->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
  
  char histname[200];
  //inclusive
  sprintf(histname,"final_centDependence_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_centbins = (TGraphErrors*)f0->Get(histname);
  if (!pg_jpsi_centbins) { cout << "cannot load nominal_NSig case." << endl; return;}
  // prompt
  sprintf(histname,"final_centDependence_NPr_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_centbins = (TGraphErrors*)f0->Get(histname);
  // non-prompt
  sprintf(histname,"final_centDependence_NNp_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_centbins = (TGraphErrors*)f0->Get(histname);

  // systematic studies
  // bkg
  sprintf(histname,"default_bit1_NBkg_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkg_centbins = (TGraphErrors*)f1->Get(histname);

  sprintf(histname,"zVtxLT10_bit1_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiVtx10_centbins = (TGraphErrors*)f1->Get(histname);
  
  // flat no flat
  sprintf(histname,"noFlat_bit1_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiNoFlat_centbins = (TGraphErrors*)f0->Get(histname);

  // enhanced autocorr
  sprintf(histname,"autoCorr_bit1_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiAutoCorr_centbins = (TGraphErrors*)f0->Get(histname);
 
  // with efficiency correction
  sprintf(histname,"default_bit1_weight_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiEffCorr_centbins = (TGraphErrors*)f0->Get(histname);

  // all triggers
  sprintf(histname,"default_bit1_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_centbins = (TGraphErrors*)f1->Get(histname);

  // ######################## sailors and cowboys
  sprintf(histname,"default_cowboy_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_centbins = (TGraphErrors*)f1->Get(histname);

  // with corrections
  sprintf(histname,"default_cowboy_weight_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_weight_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_weight_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_weight_centbins = (TGraphErrors*)f1->Get(histname);


  // mu valid hits
  sprintf(histname,"nMuValHits12_cowboy_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"nMuValHits12_sailor_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins = (TGraphErrors*)f1->Get(histname);


  // zvtx
  sprintf(histname,"zVtxLT10_cowboy_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"zVtxLT10_sailor_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins = (TGraphErrors*)f1->Get(histname);

  // midrapidity
 sprintf(histname,"singleMuLTeta1.2_cowboy_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"singleMuLTeta1.2_sailor_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins = (TGraphErrors*)f1->Get(histname);

  // cowboy bkg
    sprintf(histname,"default_cowboy_NBkg_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_onlyCow_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_NBkg_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_noCow_centbins = (TGraphErrors*)f1->Get(histname);

  //_____________________________________________________________________________
  // centrality dependence
  // N_{part} axis
  TCanvas *pcCentrality = new TCanvas("pcCentrality","pcCentrality");
  pcCentrality->cd();
  TH1F *phPadCentrality = new TH1F("phPadCentrality",";N_{part};v_{2};",400,0,400);
  phPadCentrality->GetXaxis()->SetLabelSize(20);
  phPadCentrality->GetXaxis()->SetLabelFont(43);
  phPadCentrality->GetXaxis()->SetTitleSize(30);
  phPadCentrality->GetXaxis()->SetTitleFont(43);
  phPadCentrality->GetXaxis()->SetTitleOffset(1.1);
  phPadCentrality->GetXaxis()->CenterTitle();

  phPadCentrality->GetYaxis()->SetLabelSize(20);
  phPadCentrality->GetYaxis()->SetLabelFont(43);
  phPadCentrality->GetYaxis()->SetTitleSize(32);
  phPadCentrality->GetYaxis()->SetTitleFont(43);
  phPadCentrality->GetYaxis()->SetTitleOffset(1.1);
  phPadCentrality->GetYaxis()->CenterTitle();

  phPadCentrality->SetMaximum(0.25);
  phPadCentrality->SetMinimum(-0.1);
  // phPadCentrality->SetMaximum(1);
  //phPadCentrality->SetMinimum(-1);
  if(bAddBkg)
    {
      phPadCentrality->SetMaximum(0.4);
    }
  phPadCentrality->Draw();
  
  switch(choseSignal){
  case 1:
    pg_jpsi_centbins->SetMarkerStyle(20);
    pg_jpsi_centbins->SetMarkerSize(1.8);
    pg_jpsi_centbins->SetMarkerColor(602);
    pg_jpsi_centbins->SetLineColor(602);
    if(!bDoCS)pg_jpsi_centbins->Draw("[P]");
    break;
  case 2:
    pg_pr_jpsi_centbins->SetMarkerStyle(20);
    pg_pr_jpsi_centbins->SetMarkerSize(1.8);
    pg_pr_jpsi_centbins->SetMarkerColor(602);
    pg_pr_jpsi_centbins->SetLineColor(602);
    pg_pr_jpsi_centbins->Draw("[P]");
    break;
  case 3:
    pg_npr_jpsi_centbins->SetMarkerStyle(20);
    pg_npr_jpsi_centbins->SetMarkerSize(1.8);
    pg_npr_jpsi_centbins->SetMarkerColor(602);
    pg_npr_jpsi_centbins->SetLineColor(602);
    pg_npr_jpsi_centbins->Draw("[P]");
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }

  if(bAddBkg)
    {
      pg_bkg_centbins->SetMarkerStyle(24);
      pg_bkg_centbins->SetMarkerSize(1.8);
      pg_bkg_centbins->SetMarkerColor(kBlue);
      pg_bkg_centbins->SetLineColor(kBlue);
      pg_bkg_centbins->Draw("[P]");
      TLegend *legCent = new TLegend(0.7,0.6,0.9,0.7);
      legCent->SetFillColor(0);
      legCent->SetBorderSize(0);
      legCent->SetTextSize(0.03);
      legCent->AddEntry(pg_jpsi_centbins,"Signal","P");
      legCent->AddEntry(pg_bkg_centbins,"Bkg","P");
      legCent->Draw("same");
    }

 if(bAddVtx)
    {
      pg_jpsiVtx10_centbins->SetMarkerStyle(24);
      pg_jpsiVtx10_centbins->SetMarkerSize(1.8);
      pg_jpsiVtx10_centbins->SetMarkerColor(kBlue);
      pg_jpsiVtx10_centbins->SetLineColor(kBlue);
      pg_jpsiVtx10_centbins->Draw("[P]");
      
      TLegend *legVtx = new TLegend(0.6,0.55,0.9,0.7);
      legVtx->SetFillColor(0);
      legVtx->SetBorderSize(0);
      legVtx->SetTextSize(0.03);
      legVtx->AddEntry(pg_jpsi_centbins,"Default: |z_{vtx}|<25 cm","P");
      legVtx->AddEntry(pg_jpsiVtx10_centbins,"|z_{vtx}|<10 cm","P");
      legVtx->Draw("same");
    }

  
  if(bAddNoFlat)
    {
      pg_jpsiNoFlat_centbins->SetMarkerStyle(24);
      pg_jpsiNoFlat_centbins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_centbins->SetMarkerColor(kBlue);
      pg_jpsiNoFlat_centbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.6,0.55,0.9,0.7);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      legEp->AddEntry(pg_jpsi_centbins,"Default: w/ flattening","P");
      legEp->AddEntry(pg_jpsiNoFlat_centbins,"No flattening","P");
      
      legEp->Draw("same");
    }
  if(bAddAutoCor)
    {
      pg_jpsiAutoCorr_centbins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_centbins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_centbins->SetMarkerColor(kBlue);
      pg_jpsiAutoCorr_centbins->Draw("p");
      
      TLegend *legAuto = new TLegend(0.7,0.55,0.9,0.7);
      legAuto->SetFillColor(0);
      legAuto->SetBorderSize(0);
      legAuto->SetTextSize(0.03);
           
      legAuto->AddEntry(pg_jpsi_centbins,"+- && -+","P");
      legAuto->AddEntry(pg_jpsiAutoCorr_centbins,"++ && --","P");
      
      legAuto->Draw("same");
    }

  if(bDoCowboySailor)
    {
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_centbins->Draw("[P]");
      pg_jpsiL1NHitTrig_noCow_centbins->Draw("[P]");
      
      TLegend *legTrig = new TLegend(0.2,0.55,0.57,0.75);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      legTrig->AddEntry(pg_jpsi_centbins,"Default: sum","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
      legTrig->Draw("same");
    }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(bDoCS) 
    {
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerStyle(24);

      pg_jpsiL1NHitTrig_onlyCow_weight_centbins->SetMarkerStyle(21);
      pg_jpsiL1NHitTrig_noCow_weight_centbins->SetMarkerStyle(21);
	
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_centbins->Draw("[P]");
      pg_jpsiL1NHitTrig_noCow_centbins->Draw("[P]");

      if(bDoEffCorrection_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_weight_centbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_centbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_weight_centbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_centbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_weight_centbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_weight_centbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_weight_centbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_weight_centbins->Draw("[P]");
	  
	  TLegend *legEffCS = new TLegend(0.2,0.55,0.4,0.75);
	  legEffCS->SetFillColor(0);
	  legEffCS->SetBorderSize(0);
	  legEffCS->SetTextSize(0.03);
	  
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_weight_centbins,"Cowboys+Eff","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_weight_centbins,"Sailors+Eff","P");
	  
	  legEffCS->Draw("same");

	}

      if(bDoMuHits_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins->Draw("[P]");
	  
	  TLegend *legHitsCS = new TLegend(0.2,0.55,0.4,0.75);
	  legHitsCS->SetFillColor(0);
	  legHitsCS->SetBorderSize(0);
	  legHitsCS->SetTextSize(0.03);
	  
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins,"Cowboys+12 muValidHits","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins,"Sailors+12 muValidHits","P");
	  
	  legHitsCS->Draw("same");

	}

      if(bDoBarrel_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins->Draw("[P]");
	  
	  TLegend *legBarrelCS = new TLegend(0.2,0.55,0.4,0.75);
	  legBarrelCS->SetFillColor(0);
	  legBarrelCS->SetBorderSize(0);
	  legBarrelCS->SetTextSize(0.03);
	  
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins,"Cowboys |#eta^{#mu}|<1.2","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins,"Sailorsin |#eta^{#mu}|<1.2","P");
	  
	  legBarrelCS->Draw("same");
	}

      if(bDoVtx_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins->Draw("[P]");
	  
	  TLegend *legBarrelCS = new TLegend(0.2,0.55,0.4,0.75);
	  legBarrelCS->SetFillColor(0);
	  legBarrelCS->SetBorderSize(0);
	  legBarrelCS->SetTextSize(0.03);
	  
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys |z_{vtx}|<25 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins,"Cowboys |z_{vtx}|<10 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors |z_{vtx}|<25 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins,"Sailors |z_{vtx}|<10 cm","P");
	  
	  legBarrelCS->Draw("same");
	}
      if(doBkg_cs)
	{
	  pg_bkgL1NHitTrig_onlyCow_centbins->SetMarkerColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_centbins->SetMarkerColor(kMagenta);
	  pg_bkgL1NHitTrig_onlyCow_centbins->SetLineColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_centbins->SetLineColor(kMagenta);
	      
	  pg_bkgL1NHitTrig_onlyCow_centbins->SetMarkerSize(1.8);
	  pg_bkgL1NHitTrig_noCow_centbins->SetMarkerSize(1.8);

	  pg_bkgL1NHitTrig_onlyCow_centbins->Draw("[P]");
	  pg_bkgL1NHitTrig_noCow_centbins->Draw("[P]");
	  
	  TLegend *legBkgCS = new TLegend(0.2,0.55,0.4,0.75);
	  legBkgCS->SetFillColor(0);
	  legBkgCS->SetBorderSize(0);
	  legBkgCS->SetTextSize(0.03);
	  
	  legBkgCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
	  legBkgCS->AddEntry(pg_bkgL1NHitTrig_onlyCow_centbins,"Bkg Cowboys","P");
	  legBkgCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
	  legBkgCS->AddEntry(pg_bkgL1NHitTrig_noCow_centbins,"Bkg Sailors","P");
	  
	  legBkgCS->Draw("same");
	}
    }// do cowboy-sailor systematics

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(bDoEffCorrection)
    {
      pg_jpsiEffCorr_centbins->SetMarkerStyle(24);
      pg_jpsiEffCorr_centbins->SetMarkerSize(1.8);
      pg_jpsiEffCorr_centbins->SetMarkerColor(kBlue);
      pg_jpsiEffCorr_centbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.6,0.55,0.9,0.7);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      legEp->AddEntry(pg_jpsi_centbins,"Default: w/o eff correction","P");
      legEp->AddEntry(pg_jpsiEffCorr_centbins,"W/ eff correction","P");
      
      legEp->Draw("same");

    }
  //_______ stuff to write
  TLatex *tex_c1 = new TLatex(0.53,0.92,"CMS Preliminary");
  tex_c1->SetNDC();
  tex_c1->SetTextAlign(13);
  tex_c1->SetTextFont(43);
  tex_c1->SetTextSize(25);
  tex_c1->SetLineWidth(1);
  tex_c1->Draw();

  TLatex *tex_c2 = new TLatex(0.53,0.86,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  tex_c2->SetNDC();
  tex_c2->SetTextAlign(13);
  tex_c2->SetTextFont(43);
  tex_c2->SetTextSize(25);
  tex_c2->SetLineWidth(2);
  tex_c2->Draw();

  TLatex *tex_c3 = new TLatex(0.53,0.80,"L_{int} = 150 #mub^{-1}");
  tex_c3->SetNDC();
  tex_c3->SetTextAlign(13);
  tex_c3->SetTextFont(43);
  tex_c3->SetTextSize(25);
  tex_c3->SetLineWidth(2);
  tex_c3->Draw();

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",rapIntegrated[1]));       // rapidity
  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[0], ptIntegrated[1])); 
  //  lt1->DrawLatex(0.18,0.71,Form("%s",eventPlane[1]));

  if(!bAddTrigger)
    {
      lt1->DrawLatex(0.25, 0.2, Form("30-60%%"));
      lt1->DrawLatex(0.45, 0.2, Form("20-30%%"));
      lt1->DrawLatex(0.6, 0.2, Form("10-20%%"));
      lt1->DrawLatex(0.8, 0.2, Form("0-10%%"));
    }
  else
    {
      lt1->DrawLatex(0.25, 0.34, Form("30-60%%"));
      lt1->DrawLatex(0.45, 0.34, Form("20-30%%"));
      lt1->DrawLatex(0.6, 0.34, Form("10-20%%"));
      lt1->DrawLatex(0.8, 0.34, Form("0-10%%"));
    }

  if(bSavePlots)
    {
      gSystem->mkdir(Form("./figs/png/"),kTRUE);
      gSystem->mkdir(Form("./figs/pdf/"),kTRUE);
      pcCentrality->SaveAs(Form("figs/png/%s_%s_v2_centBins.png",chosenSignal,outputName[0]));
      pcCentrality->SaveAs(Form("figs/pdf/%s_%s_v2_centBins.pdf",chosenSignal,outputName[0]));
    }
 
}
