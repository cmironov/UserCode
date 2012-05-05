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

  const char* signal[5]   = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]   = {"","Inclusive J/#psi", "Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal         = 3; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominal_bit1_cent_"};


  // options
  bool bSavePlots      = true;
  bool bDoCowboySailor = false;
  bool bDoFitting      = false;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bDoEffCorrection= false;
  bool bDo6Bins        = false;
  
  bool bAddBkg         = false;
  bool bAddAutoCor     = false;
  bool bAddTrigger     = false;

  bool bDoSummary      = false;

  bool bDoCS               = false; // this always on for the bellow options
  bool bDoSegMatchGT1      = true;
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

 
  TFile *f0;
  TFile *f1;
  if(choseSignal==1)
    {
      f0 = new TFile(Form("./a1_corrV2.root"));
      if(!f0->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
      f1 = new TFile(Form("./extrastud/a1_corrV2.root"));
      if(!f1->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
    }
  else
    {
      f0 = new TFile(Form("./a1_corrV2_pnp.root"));
      if(!f0->IsOpen()) { cout << "cannot open a1_corrV2_pnp.root" << endl; return;}
      // f1 = new TFile(Form("./extrastud/a1_corrV2_pnp.root"));
      // if(!f1->IsOpen()) { cout << "cannot open a1_corrV2_pnp.root" << endl; return;}
    }
  f1 = new TFile(Form("./extrastud/a1_corrV2.root"));

  char histname[200];
  //inclusive
  sprintf(histname,"final_centDependence_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_centbins = (TGraphErrors*)f0->Get(histname);

 sprintf(histname,"final_centDependence_statErr_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
 TGraphErrors *pg_jpsi_centbins_statErr = (TGraphErrors*)f0->Get(histname);

  // if (!pg_jpsi_centbins) { cout << "cannot load nominal_NSig case." << endl; return;}
  // prompt
  sprintf(histname,"final_centDependence_NPr_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_centbins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"final_centDependence_statErr_NPr_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_centbins_statErr = (TGraphErrors*)f0->Get(histname);

  // non-prompt
  sprintf(histname,"final_centDependence_NNp_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_centbins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_centDependence_statErr_NNp_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_centbins_statErr = (TGraphErrors*)f0->Get(histname);


  // systematic studies
  // fitting
  sprintf(histname,"default_constrained_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_constrained_centbins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_polFunct_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_polFunct_centbins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"default_signalCB3WN_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_signalCB3WN_centbins = (TGraphErrors*)f0->Get(histname);

  //------ special extra for p-np, 
  sprintf(histname,"default_bit1_1GaussResol_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_1GaussResol_centbins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_bit1_ResolFixToPRMC_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_ResolFixToPRMC_centbins = (TGraphErrors*)f0->Get(histname);


   //zvtx
  sprintf(histname,"zVtxLT10_bit1_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiVtx10_centbins = (TGraphErrors*)f0->Get(histname);
  
  // flat no flat
  sprintf(histname,"noFlat_bit1_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiNoFlat_centbins = (TGraphErrors*)f0->Get(histname);
 
  // with efficiency correction
  sprintf(histname,"default_bit1_weight_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiEffCorr_centbins = (TGraphErrors*)f0->Get(histname);

  // all triggers
  sprintf(histname,"default_bit1_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_centbins = (TGraphErrors*)f1->Get(histname);

  // //// extra studies
  // 6bins : compare among themselves the v2 obtianed with stat errors only
   TFile *f2 = new TFile(Form("./a1_corrV2_6bin.root"));
   sprintf(histname,"default_bit1_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
   TGraphErrors *pg_jpsi6bin_centbins = (TGraphErrors*)f2->Get(histname);
  TGraphErrors *pg_jpsi4bin_centbins = (TGraphErrors*)f0->Get(histname);
 //  if (!pg_jpsi4bin_centbins) { cout << "cannot load 4bincase." << endl; return;}
//   if (!pg_jpsi6bin_centbins) { cout << "cannot load 6bincase." << endl; return;}

 // enhanced autocorr
  sprintf(histname,"autoCorr_bit1_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiAutoCorr_centbins = (TGraphErrors*)f1->Get(histname);

  // bkg
  sprintf(histname,"default_bit1_NBkg_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkg_centbins = (TGraphErrors*)f1->Get(histname);
  //  if (!pg_bkg_centbins) { cout << "cannot load nominal_Bkg case." << endl; return;}
 


 // ######################## sailors and cowboys
  sprintf(histname,"default_cowboy_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_centbins = (TGraphErrors*)f1->Get(histname);
  // if (!pg_jpsiL1NHitTrig_onlyCow_centbins) { cout << "cannot load baseline cowboy." << endl; return;}
  sprintf(histname,"default_sailor_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_centbins = (TGraphErrors*)f1->Get(histname);
  // if (!pg_jpsiL1NHitTrig_noCow_centbins) { cout << "cannot load baseline sailor and signal case." << endl; return;}
  
  // with corrections
  sprintf(histname,"default_cowboy_weight_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_weight_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_weight_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_weight_centbins = (TGraphErrors*)f1->Get(histname);


  // mu valid hits
  sprintf(histname,"nMuValHits12_cowboy_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"nMuValHits12_sailor_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_nMuValHits12_centbins = (TGraphErrors*)f1->Get(histname);


  // zvtx
  sprintf(histname,"zVtxLT10_cowboy_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"zVtxLT10_sailor_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_zVtxLT10_centbins = (TGraphErrors*)f1->Get(histname);

  // midrapidity
 sprintf(histname,"singleMuLTeta1.2_cowboy_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"singleMuLTeta1.2_sailor_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_centbins = (TGraphErrors*)f1->Get(histname);

  // nSegments match
  sprintf(histname,"segMatchGT1_cowboy_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"segMatchGT1_sailor_%s_rap%.1f-%.1f_pT%.1f-%.1f",chosenSignal,rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins = (TGraphErrors*)f1->Get(histname);
  // if (!pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins) { cout << "cannot load cowboy." << endl; return;}
  // if (!pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins) { cout << "cannot load sailor." << endl; return;}

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

  if(bDoSummary)
    {
      phPadCentrality->SetMaximum(0.3);
      phPadCentrality->SetMinimum(-0.2);
    }
  phPadCentrality->Draw();
  
  switch(choseSignal){
  case 1:
    //pg_jpsi_centbins->SetFillStyle(0);
    pg_jpsi_centbins->SetFillColor(kBlue-10);
  
    pg_jpsi_centbins_statErr->SetMarkerStyle(20);
    pg_jpsi_centbins_statErr->SetMarkerSize(1.8);
    pg_jpsi_centbins_statErr->SetMarkerColor(kBlue+2);
    pg_jpsi_centbins_statErr->SetLineColor(kBlue+2);
    if(!bDoCS && !bDo6Bins)
      {
	pg_jpsi_centbins->Draw("2");
	pg_jpsi_centbins_statErr->Draw("[P]");
      }
    break;
  case 2:
    pg_pr_jpsi_centbins->SetFillColor(kRed-10);

    pg_pr_jpsi_centbins_statErr->SetMarkerStyle(21);
    pg_pr_jpsi_centbins_statErr->SetMarkerSize(1.8);
    pg_pr_jpsi_centbins_statErr->SetMarkerColor(kRed+2);
    pg_pr_jpsi_centbins_statErr->SetLineColor(kRed+2);
    if(!bDoCS && !bDo6Bins)
      {
	pg_pr_jpsi_centbins->Draw("2");
	pg_pr_jpsi_centbins_statErr->Draw("[P]");
      }
    break;
  case 3:
    pg_npr_jpsi_centbins->SetFillColor(kOrange-9);

    pg_npr_jpsi_centbins_statErr->SetMarkerStyle(33);
    pg_npr_jpsi_centbins_statErr->SetMarkerSize(1.8);
    pg_npr_jpsi_centbins_statErr->SetMarkerColor(kOrange+2);
    pg_npr_jpsi_centbins_statErr->SetLineColor(kOrange+2);
    if(!bDoCS && !bDo6Bins)
      {
	pg_npr_jpsi_centbins->Draw("2");
	pg_npr_jpsi_centbins_statErr->Draw("[P]");
      }
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }
  cout<<"####### GATE !"<<endl;
  if(bDo6Bins)
    {
      pg_jpsi6bin_centbins->SetMarkerStyle(24);
      pg_jpsi6bin_centbins->SetMarkerSize(1.8);
      pg_jpsi6bin_centbins->SetMarkerColor(kBlue+2);
      pg_jpsi6bin_centbins->SetLineColor(kBlue+2);
      
      pg_jpsi4bin_centbins->SetMarkerStyle(20);
      pg_jpsi4bin_centbins->SetMarkerSize(1.8);
      pg_jpsi4bin_centbins->SetMarkerColor(kBlue+2);
      pg_jpsi4bin_centbins->SetLineColor(kBlue+2);

      pg_jpsi4bin_centbins->Draw("[P]");
      pg_jpsi6bin_centbins->Draw("[P]");

      TLegend *leg6Bin = new TLegend(0.2,0.6,0.5,0.7);
      leg6Bin->SetFillColor(0);
      leg6Bin->SetBorderSize(0);
      leg6Bin->SetTextSize(0.03);
      leg6Bin->AddEntry(pg_jpsi4bin_centbins,"Default:4 d#phi bins","P");
      leg6Bin->AddEntry(pg_jpsi6bin_centbins,"6 d#phi bins","P");
      leg6Bin->Draw("same");
    }

  if(bDoFitting)
    {
      pg_constrained_centbins->SetMarkerStyle(24);
      pg_polFunct_centbins->SetMarkerStyle(25);
      pg_signalCB3WN_centbins->SetMarkerStyle(27);

      pg_constrained_centbins->SetMarkerColor(kBlue);
      pg_polFunct_centbins->SetMarkerColor(kBlue);
      pg_signalCB3WN_centbins->SetMarkerColor(kBlue);

      pg_constrained_centbins->SetLineColor(kBlue);
      pg_polFunct_centbins->SetLineColor(kBlue);
      pg_signalCB3WN_centbins->SetLineColor(kBlue);

      pg_constrained_centbins->SetMarkerSize(1.8);
      pg_polFunct_centbins->SetMarkerSize(1.8);
      pg_signalCB3WN_centbins->SetMarkerSize(1.8);

      pg_constrained_centbins->Draw("[P]");
      pg_polFunct_centbins->Draw("[P]");
      pg_signalCB3WN_centbins->Draw("[P]");

      TLegend *legFit = new TLegend(0.15,0.6,0.5,0.75);
      legFit->SetFillColor(0);
      legFit->SetBorderSize(0);
      legFit->SetTextSize(0.03);
      legFit->AddEntry(pg_jpsi_centbins,"Nominal","P");
      legFit->AddEntry(pg_constrained_centbins,"Constrained","P");
      legFit->AddEntry(pg_signalCB3WN_centbins,"Sgn: C-B only","P");
      legFit->AddEntry(pg_polFunct_centbins,"Bkg:straight line ","P");


      if(choseSignal!=1)
	{
	   // for p=np only
	  pg_1GaussResol_centbins->SetMarkerStyle(28);
	  pg_ResolFixToPRMC_centbins->SetMarkerStyle(30);
	  
	  pg_1GaussResol_centbins->SetMarkerColor(kBlue);
	  pg_ResolFixToPRMC_centbins->SetMarkerColor(kBlue);
	  pg_1GaussResol_centbins->SetLineColor(kBlue);
	  pg_ResolFixToPRMC_centbins->SetLineColor(kBlue);
	  
	  pg_1GaussResol_centbins->SetMarkerSize(1.8);
	  pg_ResolFixToPRMC_centbins->SetMarkerSize(1.8);
	  
	  pg_1GaussResol_centbins->Draw("[P]");
	  pg_ResolFixToPRMC_centbins->Draw("[P]");
	  legFit->AddEntry(pg_1GaussResol_centbins,"Sgn Res: 1 Gauss","P");
	  legFit->AddEntry(pg_ResolFixToPRMC_centbins,"Sgn Res: fixed to MC ","P");
	}
      legFit->Draw("same");
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
      
      TLegend *legVtx = new TLegend(0.2,0.6,0.5,0.75);
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
      
      TLegend *legEp = new TLegend(0.2,0.6,0.5,0.75);
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
      
      TLegend *legTrig = new TLegend(0.2,0.62,0.5,0.76);
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

      if(bDoSegMatchGT1)
	{
	  pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins->Draw("[P]");
	  
	  TLegend *legEffCS = new TLegend(0.2,0.55,0.4,0.75);
	  legEffCS->SetFillColor(0);
	  legEffCS->SetBorderSize(0);
	  legEffCS->SetTextSize(0.03);
	  
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"Cowboys","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_segMatchGT1_centbins,"Cowboys: nSegMatched > 1 && 1dilepton/event","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"Sailors","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_segMatchGT1_centbins,"Sailors: nSegMatched > 1 && 1dilepton/event","P");
	  
	  legEffCS->Draw("same");
	}

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
      
      TLegend *legEp = new TLegend(0.2,0.6,0.5,0.7);
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
  cout<<"####### GATE 222!"<<endl;

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
      if(bDoCS)
	{
	  pcCentrality->SaveAs(Form("figs/png/%s_%s_v2cs_centBins.png",chosenSignal,outputName[0]));
	  pcCentrality->SaveAs(Form("figs/pdf/%s_%s_v2cs_centBins.pdf",chosenSignal,outputName[0]));
	}
      else
	{
	  pcCentrality->SaveAs(Form("figs/png/%s_%s_v2_centBins.png",chosenSignal,outputName[0]));
	  pcCentrality->SaveAs(Form("figs/pdf/%s_%s_v2_centBins.pdf",chosenSignal,outputName[0]));
	}

    }
 

  //****************************** RESULTS PLOTS!!!!
//   if(bDoSummary)
//     {
//       TCanvas *pcSummary = new TCanvas("pcSummary","pcSummary",1200,400);
//       pcSummary->Divide(3,1);
      
//       TH1F *phPadSummary = new TH1F("phPadSummary",";N_{part};v_{2};",400,0,400);

//       phPadSummary->GetXaxis()->SetLabelSize(20);
//       phPadSummary->GetXaxis()->SetLabelFont(43);
//       phPadSummary->GetXaxis()->SetTitleSize(30);
//       phPadSummary->GetXaxis()->SetTitleFont(43);
//       phPadSummary->GetXaxis()->SetTitleOffset(1.1);
//       phPadSummary->GetXaxis()->CenterTitle();
      
//       phPadSummary->GetYaxis()->SetLabelSize(20);
//       phPadSummary->GetYaxis()->SetLabelFont(43);
//       phPadSummary->GetYaxis()->SetTitleSize(32);
//       phPadSummary->GetYaxis()->SetTitleFont(43);
//       phPadSummary->GetYaxis()->SetTitleOffset(1.1);
//       phPadSummary->GetYaxis()->CenterTitle();
      
//       phPadSummary->SetMaximum(0.5);
//       phPadSummary->SetMinimum(-0.3);
   

//       pcSummary->cd(1);
//       phPadSummary->Draw();

//       pcSummary->cd(2);
//       phPadSummary->Draw();

//       pcSummary->cd(3);
//       phPadSummary->Draw();

//     }

}
