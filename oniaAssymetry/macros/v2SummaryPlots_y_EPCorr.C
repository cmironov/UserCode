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


void v2SummaryPlots_y_EPCorr()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  // gROOT->Macro("./rootlogon.C");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]   = {"","Inclusive J/#psi", "Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal         = 3; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominal_bit1_y_6bin"};

  // options
  bool bSavePlots      = true;
  bool bDoLowpt        = false;

  bool bDoSummary      = true;
 
  bool bDoCowboySailor = false;
  bool bDoFitting      = false;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bDoEffCorrection= false;
  bool bDo6Bins        = false;


  bool bAddBkg         = false;
  bool bAddTrigger     = false; 
  bool bAddAutoCor     = false;

  bool bDoCS               = false;
  bool bDoEffCorrection_cs = false;
  bool bDoMuHits_cs        = false;
  bool bDoBarrel_cs        = false;
  bool bDoVtx_cs           = false;
  bool doBkg_cs            = false;

  int centIntegrated[2]    = {10, 60}; 
  double ptIntegrated[3]   = {3.0, 6.5, 40.0}; 
 
   TFile *f0;
  TFile *f1;
  if(!bDoSummary)
    {
      if(choseSignal==1)
	{
	  f0 = new TFile(Form("./a3_corrV2.root"));
	  if(!f0->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
	  f1 = new TFile(Form("./extrastud/a3_corrV2.root"));
	  if(!f1->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
	}
      else
	{
	  f0 = new TFile(Form("./a3_corrV2_pnp.root"));
	  if(!f0->IsOpen()) { cout << "cannot open a3_corrV2_pnp.root" << endl; return;}
	  // f1 = new TFile(Form("./extrastud/a1_corrV2_pnp.root"));
	  // if(!f1->IsOpen()) { cout << "cannot open a1_corrV2_pnp.root" << endl; return;}
	}
      f1 = new TFile(Form("./extrastud/a3_corrV2.root"));
    }
  else
    {
      f0 = new TFile(Form("./a3_corrV2.root"));
      if(!f0->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
      f1 = new TFile(Form("./a3_corrV2_pnp.root"));
      if(!f1->IsOpen()) { cout << "cannot open a3_corrV2_pnp.root" << endl; return;}
    }

  // pt_rap bins, high_pt
  const int nybins_highpt = 3;
  double ybins_highpt_center[nybins_highpt]          = {0.6, 1.4, 2.0};
  double ybins_highpt_center_err[nybins_highpt]      = {0.6, 0.2, 0.4};
  double ybins_highpt_center_errGhost[nybins_highpt] = {0.,0.,0.};

  // sgn
  double jpsi_highpt_ptybins[nybins_highpt]          = {0.05, 0.05, 0.05};
  double jpsi_highpt_ptybins_err[nybins_highpt]      = {0.017, 0.026, 0.024};
  double jpsi_highpt_ptybins_errGhost[nybins_highpt] = {0.,0.,0 };

  // finer bins 
  const int nyfinebins_highpt = 5;
  double yfinebins_highpt_center[nyfinebins_highpt]          = {0.3, 0.9, 1.4, 1.8, 2.2};
  double yfinebins_highpt_center_err[nyfinebins_highpt]      = {0.3, 0.3, 0.2, 0.2, 0.2};
  double yfinebins_highpt_center_errGhost[nyfinebins_highpt] = {0.,  0.,  0.,  0.,  0.};

  double jpsiCorr_highpt_ptyfinebins[nyfinebins_highpt]      = {0.0703,
								0.0738,
								0.0800,
								0.0685,
								0.0590};
  
  double jpsiCorr_highpt_ptyfinebins_err[nyfinebins_highpt]  = {0.0255,
								0.0231,
								0.0256,
								0.0280,
								0.0486};

  char histname[200];
  // correlations in y:
  TGraphErrors *pg_jpsiCorr_highpt_ptyfinebins    = new TGraphErrors(nyfinebins_highpt,yfinebins_highpt_center,jpsiCorr_highpt_ptyfinebins,yfinebins_highpt_center_err,jpsiCorr_highpt_ptyfinebins_err);


  //// High pT
  //inclusive
  sprintf(histname,"final_yDependence_highpt_NSig_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
 
  sprintf(histname,"final_yDependence_highpt_statErr_NSig_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi_highpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);

 
 // prompt
  sprintf(histname,"final_yDependence_highpt_NPr_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_pr_jpsi_highpt_ptybins;
  if(bDoSummary) pg_pr_jpsi_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  else pg_pr_jpsi_highpt_ptybins = (TGraphErrors*)f0->Get(histname);


   sprintf(histname,"final_yDependence_highpt_statErr_NPr_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_pr_jpsi_highpt_ptybins_statErr;
  if(bDoSummary) pg_pr_jpsi_highpt_ptybins_statErr = (TGraphErrors*)f1->Get(histname);
  else pg_pr_jpsi_highpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);
 
  // non-prompt
  sprintf(histname,"final_yDependence_highpt_NNp_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_npr_jpsi_highpt_ptybins;
  if(bDoSummary)pg_npr_jpsi_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  else pg_npr_jpsi_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"final_yDependence_highpt_statErr_NNp_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_npr_jpsi_highpt_ptybins_statErr;
  if(bDoSummary)pg_npr_jpsi_highpt_ptybins_statErr = (TGraphErrors*)f1->Get(histname);
  else pg_npr_jpsi_highpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);

  // systematic studies
  // 6bins : compare among themselves the v2 obtianed with stat errors only
  TFile *f2 = new TFile(Form("./a3_corrV2_6bin.root"));
  sprintf(histname,"default_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi6bin_ptybins = (TGraphErrors*)f2->Get(histname);
  TGraphErrors *pg_jpsi4bin_ptybins = (TGraphErrors*)f0->Get(histname);
  // if (!pg_jpsi4bin_ptybins) { cout << "cannot load 4bincase." << endl; return;}
  // if (!pg_jpsi6bin_ptybins) { cout << "cannot load 6bincase." << endl; return;}

  // fits
  sprintf(histname,"default_constrained_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi_constrained_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_polFunct_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi_polFunct_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_signalCB3WN_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsi_signalCB3WN_highpt_ptybins = (TGraphErrors*)f0->Get(histname);

  //------ special extra for p-np, 
  sprintf(histname,"default_bit1_1GaussResol_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_1GaussResol_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_bit1_ResolFixToPRMC_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_ResolFixToPRMC_highpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // zvtx  
  sprintf(histname,"zVtxLT10_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiVtx10_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  
  // flat no flat
  sprintf(histname,"noFlat_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiNoFlat_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  
 
  // with efficiency correction
  sprintf(histname,"default_bit1_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiEffCorr_highpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // all triggers
  sprintf(histname,"default_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_highpt_ptybins = (TGraphErrors*)f1->Get(histname);

  /// extra studies
  // bkg
  sprintf(histname,"default_bit1_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_bkg_highpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // enhanced autocorr
  sprintf(histname,"autoCorr_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiAutoCorr_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  // ######################## sailors and cowboys
  sprintf(histname,"default_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_highpt_ptybins = (TGraphErrors*)f0->Get(histname);




  // with corrections
  sprintf(histname,"default_cowboy_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins = (TGraphErrors*)f1->Get(histname);


  // mu valid hits
  sprintf(histname,"nMuValHits12_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"nMuValHits12_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // zvtx
  sprintf(histname,"zVtxLT10_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"zVtxLT10_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // midrapidity
  sprintf(histname,"singleMuLTeta1.2_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"singleMuLTeta1.2_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // cowboy bkg
  sprintf(histname,"default_cowboy_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_bkgL1NHitTrig_onlyCow_highpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[1],ptIntegrated[2]);
  TGraphErrors *pg_bkgL1NHitTrig_noCow_highpt_ptybins   = (TGraphErrors*)f1->Get(histname);

  /////////////////////////////////// Low pt
  sprintf(histname,"final_yDependence_lowpt_NSig_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_yDependence_lowpt_statErr_NSig_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_lowpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);
  
  // prompt
  sprintf(histname,"final_yDependence_lowpt_NPr_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_lowpt_ptybins;
  if(bDoSummary) pg_pr_jpsi_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  else pg_pr_jpsi_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_yDependence_lowpt_statErr_NPr_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_lowpt_ptybins_statErr;
  if(bDoSummary) pg_pr_jpsi_lowpt_ptybins_statErr = (TGraphErrors*)f1->Get(histname);
  else pg_pr_jpsi_lowpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);
  
  // non-prompt
  sprintf(histname,"final_yDependence_lowpt_NNp_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_lowpt_ptybins;
  if(bDoSummary)pg_npr_jpsi_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  else pg_npr_jpsi_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_yDependence_lowpt_statErr_NNp_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_lowpt_ptybins_statErr;
  if(bDoSummary)pg_npr_jpsi_lowpt_ptybins_statErr = (TGraphErrors*)f1->Get(histname);
  else pg_npr_jpsi_lowpt_ptybins_statErr = (TGraphErrors*)f0->Get(histname);

  // systematic studies
  //fits
    sprintf(histname,"default_constrained_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_constrained_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_polFunct_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_polFunct_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_signalCB3WN_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_signalCB3WN_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  //------ special extra for p-np, 
  sprintf(histname,"default_bit1_1GaussResol_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_1GaussResol_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_bit1_ResolFixToPRMC_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_ResolFixToPRMC_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);



  // bkg
  sprintf(histname,"default_bit1_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkg_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);
  // if (!pg_bkg_lowpt_ptybins) { cout << "cannot load highpt hist" << endl; return;}


  //vtx
  sprintf(histname,"zVtxLT10_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiVtx10_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  
  // flat no flat
  sprintf(histname,"noFlat_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiNoFlat_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // enhanced autocorr
  sprintf(histname,"autoCorr_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiAutoCorr_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
 
  // with efficiency correction
  sprintf(histname,"default_bit1_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiEffCorr_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // all triggers
  sprintf(histname,"default_bit1_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // ######################## sailors and cowboys
  sprintf(histname,"default_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_lowpt_ptybins = (TGraphErrors*)f0->Get(histname);

  // with corrections
  sprintf(histname,"default_cowboy_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_weight_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // mu valid hits
  sprintf(histname,"nMuValHits12_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"nMuValHits12_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // zvtx
  sprintf(histname,"zVtxLT10_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"zVtxLT10_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // midrapidity
  sprintf(histname,"singleMuLTeta1.2_cowboy_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"singleMuLTeta1.2_sailor_%s_cent%d-%d_pt%.1f-%.1f",chosenSignal,centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);

  // cowboy bkg
    sprintf(histname,"default_cowboy_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_NBkg_cent%d-%d_pt%.1f-%.1f",centIntegrated[0],centIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_noCow_lowpt_ptybins   = (TGraphErrors*)f1->Get(histname);
 
  //*************************************************************************************
  // rapidity dependence
  TCanvas *pcY = new TCanvas("pcY","pcY");
  pcY->cd();
  TH1F *pcPadY = new TH1F("pcPadY",";|y|;v_{2};",100,0,2.4);
  pcPadY->GetXaxis()->SetLabelSize(20);
  pcPadY->GetXaxis()->SetLabelFont(43);
  pcPadY->GetXaxis()->SetTitleSize(27);
  pcPadY->GetXaxis()->SetTitleFont(43);
  pcPadY->GetXaxis()->SetTitleOffset(1.2);
  pcPadY->GetXaxis()->CenterTitle();

  pcPadY->GetYaxis()->SetLabelSize(20);
  pcPadY->GetYaxis()->SetLabelFont(43);
  pcPadY->GetYaxis()->SetTitleSize(32);
  pcPadY->GetYaxis()->SetTitleFont(43);
  pcPadY->GetYaxis()->SetTitleOffset(1.1);
  pcPadY->GetYaxis()->CenterTitle();

  pcPadY->SetMaximum(0.25);
  pcPadY->SetMinimum(-0.1);
  //pcPadY->SetMaximum(0.18);
  //pcPadY->SetMinimum(-0.05);

  if(bAddBkg)
  {
    pcPadY->SetMaximum(0.4);
  }

  if(bDoSummary)
    {
      pcPadY->SetMaximum(0.3);
      pcPadY->SetMinimum(-0.2);
    }
  pcPadY->Draw();
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
  //-----------------------


//   pg_jpsi_highpt_ptybinsGhost->SetMarkerStyle(20);
//   pg_jpsi_highpt_ptybinsGhost->SetMarkerSize(2);
//   pg_jpsi_highpt_ptybinsGhost->SetMarkerColor(kRed+2);
//   pg_jpsi_highpt_ptybinsGhost->SetLineColor(kRed+2);
  //  pg_jpsi_highpt_ptybinsGhost->Draw("P");

 switch(choseSignal){
  case 1:
    pg_jpsi_highpt_ptybins->SetFillColor(kBlue-10);

    pg_jpsi_highpt_ptybins_statErr->SetMarkerStyle(20);
    pg_jpsi_highpt_ptybins_statErr->SetMarkerSize(2.0);
    pg_jpsi_highpt_ptybins_statErr->SetMarkerColor(kBlue+2);
    pg_jpsi_highpt_ptybins_statErr->SetLineColor(kBlue+2);
    if(!bDoCS && !bDoLowpt && !bDo6Bins)
      {
	pg_jpsi_highpt_ptybins->Draw("2");
	pg_jpsi_highpt_ptybins_statErr->Draw("[P]");
      }

    pg_jpsi_lowpt_ptybins->SetFillColor(kBlue-8);

    pg_jpsi_lowpt_ptybins_statErr->SetMarkerStyle(24);
    pg_jpsi_lowpt_ptybins_statErr->SetMarkerSize(2.0);
    pg_jpsi_lowpt_ptybins_statErr->SetMarkerColor(kBlue+2);
    if(!bDoCS && bDoLowpt)
      {
	pg_jpsi_lowpt_ptybins->Draw("2");
	pg_jpsi_lowpt_ptybins_statErr->Draw("[P]");
      }
    if(bDoSummary && !bDoCS && !bDo6Bins) 
      {
	pg_jpsi_lowpt_ptybins->Draw("2");
	pg_jpsi_lowpt_ptybins_statErr->Draw("[P]");
      }
    break;
  case 2:
    pg_pr_jpsi_highpt_ptybins->SetFillColor(kRed-10);
    pg_pr_jpsi_lowpt_ptybins->SetFillColor(kRed-8);
    pg_pr_jpsi_highpt_ptybins_statErr->SetMarkerStyle(21);
    pg_pr_jpsi_highpt_ptybins_statErr->SetMarkerSize(1.8);
    pg_pr_jpsi_highpt_ptybins_statErr->SetMarkerColor(kRed+2);
    pg_pr_jpsi_highpt_ptybins_statErr->SetLineColor(kRed+2);
    if(!bDoCS && !bDoLowpt && !bDo6Bins)
      {
	pg_pr_jpsi_highpt_ptybins->Draw("2");
	pg_pr_jpsi_highpt_ptybins_statErr->Draw("[P]");

      }
    pg_pr_jpsi_lowpt_ptybins_statErr->SetMarkerStyle(25);
    pg_pr_jpsi_lowpt_ptybins_statErr->SetMarkerSize(1.8);
    pg_pr_jpsi_lowpt_ptybins_statErr->SetMarkerColor(kRed+2);
    if(!bDoCS && bDoLowpt)
      {
	pg_pr_jpsi_lowpt_ptybins->Draw("2");
	pg_pr_jpsi_lowpt_ptybins_statErr->Draw("[P]");

      }
    if(bDoSummary && !bDoCS ) 
      {
	pg_pr_jpsi_lowpt_ptybins->Draw("2");
	pg_pr_jpsi_lowpt_ptybins_statErr->Draw("[P]");
      }

    break;
  case 3:
    pg_npr_jpsi_highpt_ptybins->SetFillColor(kOrange-9);
    pg_npr_jpsi_lowpt_ptybins->SetFillColor(kOrange-8);
    pg_npr_jpsi_highpt_ptybins_statErr->SetMarkerStyle(33);
    pg_npr_jpsi_highpt_ptybins_statErr->SetMarkerSize(1.8);
    pg_npr_jpsi_highpt_ptybins_statErr->SetMarkerColor(kOrange+2);
    pg_npr_jpsi_highpt_ptybins_statErr->SetLineColor(kOrange+2);
    if(!bDoCS && !bDoLowpt && !bDo6Bins)
      {
	pg_npr_jpsi_highpt_ptybins->Draw("2");
	pg_npr_jpsi_highpt_ptybins_statErr->Draw("[P]");
      }
    pg_npr_jpsi_lowpt_ptybins_statErr->SetMarkerStyle(27);
    pg_npr_jpsi_lowpt_ptybins_statErr->SetMarkerSize(1.8);
    pg_npr_jpsi_lowpt_ptybins_statErr->SetMarkerColor(kOrange+2);
    if(!bDoCS &&bDoLowpt) 
      {
	pg_npr_jpsi_lowpt_ptybins->Draw("2");
	pg_npr_jpsi_lowpt_ptybins_statErr->Draw("[P]");

      }
    if(bDoSummary && !bDoCS) 
      {
	pg_npr_jpsi_lowpt_ptybins->Draw("2");
	pg_npr_jpsi_lowpt_ptybins_statErr->Draw("[P]");
      }
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }

 if(bDo6Bins)
    {
      pg_jpsi6bin_ptybins->SetMarkerStyle(24);
      pg_jpsi6bin_ptybins->SetMarkerSize(1.8);
      pg_jpsi6bin_ptybins->SetMarkerColor(kBlue+2);
      pg_jpsi6bin_ptybins->SetLineColor(kBlue+2);
      
      pg_jpsi4bin_ptybins->SetMarkerStyle(20);
      pg_jpsi4bin_ptybins->SetMarkerSize(1.8);
      pg_jpsi4bin_ptybins->SetMarkerColor(kBlue+2);
      pg_jpsi4bin_ptybins->SetLineColor(kBlue+2);
      
      pg_jpsi4bin_ptybins->Draw("[P]");
      pg_jpsi6bin_ptybins->Draw("[P]");
      
      TLegend *leg6Bin = new TLegend(0.2,0.6,0.5,0.7);
      leg6Bin->SetFillColor(0);
      leg6Bin->SetBorderSize(0);
      leg6Bin->SetTextSize(0.03);
      leg6Bin->AddEntry(pg_jpsi4bin_ptybins,"Default:4 d#phi bins","P");
      leg6Bin->AddEntry(pg_jpsi6bin_ptybins,"6 d#phi bins","P");
      leg6Bin->Draw("same");
    }


 if(bDoFitting)
   {
     if(!bDoLowpt)
       {
	 pg_jpsi_constrained_highpt_ptybins->SetMarkerStyle(24);
	 pg_jpsi_polFunct_highpt_ptybins->SetMarkerStyle(25);
	 pg_jpsi_signalCB3WN_highpt_ptybins->SetMarkerStyle(27);
	 
	 pg_jpsi_constrained_highpt_ptybins->SetMarkerColor(kBlue);
	 pg_jpsi_polFunct_highpt_ptybins->SetMarkerColor(kBlue);
	 pg_jpsi_signalCB3WN_highpt_ptybins->SetMarkerColor(kBlue);
	 
	 pg_jpsi_constrained_highpt_ptybins->SetLineColor(kBlue);
	 pg_jpsi_polFunct_highpt_ptybins->SetLineColor(kBlue);
	 pg_jpsi_signalCB3WN_highpt_ptybins->SetLineColor(kBlue);
	 
	 pg_jpsi_constrained_highpt_ptybins->SetMarkerSize(1.8);
	 pg_jpsi_polFunct_highpt_ptybins->SetMarkerSize(1.8);
	 pg_jpsi_signalCB3WN_highpt_ptybins->SetMarkerSize(1.8);
	 
	 pg_jpsi_constrained_highpt_ptybins->Draw("[P]");
	 pg_jpsi_polFunct_highpt_ptybins->Draw("[P]");
	 pg_jpsi_signalCB3WN_highpt_ptybins->Draw("[P]");
	 
	 TLegend *legFit = new TLegend(0.15,0.6,0.5,0.75);
	 legFit->SetFillColor(0);
	 legFit->SetBorderSize(0);
	 legFit->SetTextSize(0.03);
	 legFit->AddEntry(pg_jpsi_highpt_ptybins,"Nominal","P");
	 legFit->AddEntry(pg_jpsi_constrained_highpt_ptybins,"Constrained","P");
	 legFit->AddEntry(pg_jpsi_signalCB3WN_highpt_ptybins,"Signal: C-B only","P");
	 legFit->AddEntry(pg_jpsi_polFunct_highpt_ptybins,"Bkg:straight line ","P");


	 if(choseSignal!=1)
	{
	  // for p=np only
	  pg_1GaussResol_highpt_ptybins->SetMarkerStyle(28);
	  pg_ResolFixToPRMC_highpt_ptybins->SetMarkerStyle(30);
	  
	  pg_1GaussResol_highpt_ptybins->SetMarkerColor(kBlue);
	  pg_ResolFixToPRMC_highpt_ptybins->SetMarkerColor(kBlue);
	  pg_1GaussResol_highpt_ptybins->SetLineColor(kBlue);
	  pg_ResolFixToPRMC_highpt_ptybins->SetLineColor(kBlue);
	  
	  pg_1GaussResol_highpt_ptybins->SetMarkerSize(1.8);
	  pg_ResolFixToPRMC_highpt_ptybins->SetMarkerSize(1.8);
	  
	  pg_1GaussResol_highpt_ptybins->Draw("[P]");
	  pg_ResolFixToPRMC_highpt_ptybins->Draw("[P]");
	  legFit->AddEntry(pg_1GaussResol_highpt_ptybins,"Sgn Res: 1 Gauss","P");
	  legFit->AddEntry(pg_ResolFixToPRMC_highpt_ptybins,"Sgn Res: fixed to MC ","P");
	}
	 
	 legFit->Draw("same");
       }// not lopt
     else
       {
	 pg_jpsi_constrained_lowpt_ptybins->SetMarkerStyle(24);
	 pg_jpsi_polFunct_lowpt_ptybins->SetMarkerStyle(25);
	 pg_jpsi_signalCB3WN_lowpt_ptybins->SetMarkerStyle(27);
	 
	 pg_jpsi_constrained_lowpt_ptybins->SetMarkerColor(kBlue);
	 pg_jpsi_polFunct_lowpt_ptybins->SetMarkerColor(kBlue);
	 pg_jpsi_signalCB3WN_lowpt_ptybins->SetMarkerColor(kBlue);
	 
	 pg_jpsi_constrained_lowpt_ptybins->SetLineColor(kBlue);
	 pg_jpsi_polFunct_lowpt_ptybins->SetLineColor(kBlue);
	 pg_jpsi_signalCB3WN_lowpt_ptybins->SetLineColor(kBlue);
	 
	 pg_jpsi_constrained_lowpt_ptybins->SetMarkerSize(1.8);
	 pg_jpsi_polFunct_lowpt_ptybins->SetMarkerSize(1.8);
	 pg_jpsi_signalCB3WN_lowpt_ptybins->SetMarkerSize(1.8);
	 
	 pg_jpsi_constrained_lowpt_ptybins->Draw("[P]");
	 pg_jpsi_polFunct_lowpt_ptybins->Draw("[P]");
	 pg_jpsi_signalCB3WN_lowpt_ptybins->Draw("[P]");
	 
	 TLegend *legFit = new TLegend(0.15,0.6,0.5,0.75);
	 legFit->SetFillColor(0);
	 legFit->SetBorderSize(0);
	 legFit->SetTextSize(0.03);
	 legFit->AddEntry(pg_jpsi_lowpt_ptybins,"Nominal","P");
	 legFit->AddEntry(pg_jpsi_constrained_lowpt_ptybins,"Constrained","P");
	 legFit->AddEntry(pg_jpsi_signalCB3WN_lowpt_ptybins,"Signal: C-B only","P");
	 legFit->AddEntry(pg_jpsi_polFunct_lowpt_ptybins,"Bkg:straight line ","P");


	 if(choseSignal!=1)
	{
	  // for p=np only
	  pg_1GaussResol_lowpt_ptybins->SetMarkerStyle(28);
	  pg_ResolFixToPRMC_lowpt_ptybins->SetMarkerStyle(30);
	  
	  pg_1GaussResol_lowpt_ptybins->SetMarkerColor(kBlue);
	  pg_ResolFixToPRMC_lowpt_ptybins->SetMarkerColor(kBlue);
	  pg_1GaussResol_lowpt_ptybins->SetLineColor(kBlue);
	  pg_ResolFixToPRMC_lowpt_ptybins->SetLineColor(kBlue);
	  
	  pg_1GaussResol_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_ResolFixToPRMC_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  pg_1GaussResol_lowpt_ptybins->Draw("[P]");
	  pg_ResolFixToPRMC_lowpt_ptybins->Draw("[P]");
	  legFit->AddEntry(pg_1GaussResol_lowpt_ptybins,"Sgn Res: 1 Gauss","P");
	  legFit->AddEntry(pg_ResolFixToPRMC_lowpt_ptybins,"Sgn Res: fixed to MC ","P");
	}
	 
	 legFit->Draw("same");
       }
   
   }

  if(bAddBkg)
  {
    pg_bkg_highpt_ptybins->SetMarkerStyle(24);
    pg_bkg_highpt_ptybins->SetMarkerSize(1.8);
    pg_bkg_highpt_ptybins->SetMarkerColor(kRed+2);
    if(!bDoLowpt)pg_bkg_highpt_ptybins->Draw("[P]");

    pg_bkg_lowpt_ptybins->SetMarkerStyle(25);
    pg_bkg_lowpt_ptybins->SetMarkerSize(1.8);
    pg_bkg_lowpt_ptybins->SetMarkerColor(kBlue+2);
    pg_bkg_lowpt_ptybins->SetLineColor(kBlue+2);
    if(bDoLowpt)pg_bkg_lowpt_ptybins->Draw("[P]");

    TLegend *legY = new TLegend(0.6,0.17,0.9,0.3);
    legY->SetFillColor(0);
    legY->SetBorderSize(0);
    legY->SetTextSize(0.03);
    if(bDoLowpt)
      {
	legY->AddEntry(pg_jpsi_lowpt_ptybins,"Signal","P");
	legY->AddEntry(pg_bkg_lowpt_ptybins,"Bkg","P");
      }
    else 
      {
	legY->AddEntry(pg_jpsi_highpt_ptybins,"Signal","P");
	legY->AddEntry(pg_bkg_highpt_ptybins,"Bkg","P");
      }
    legY->Draw("same");
  }

 if(bAddVtx)
    {
      pg_jpsiVtx10_highpt_ptybins->SetMarkerStyle(24);
      pg_jpsiVtx10_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiVtx10_highpt_ptybins->SetMarkerColor(kRed+2);
      if(!bDoLowpt)pg_jpsiVtx10_highpt_ptybins->Draw("[P]");
      

      pg_jpsiVtx10_lowpt_ptybins->SetMarkerStyle(25);
      pg_jpsiVtx10_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiVtx10_lowpt_ptybins->SetMarkerColor(kBlue+2);
      pg_jpsiVtx10_lowpt_ptybins->SetLineColor(kBlue+2);
      cout<<"crushed here 222"<<endl;
      if(bDoLowpt)pg_jpsiVtx10_lowpt_ptybins->Draw("[P]");

      
      cout<<"crushed here"<<endl;
      TLegend *legVtx = new TLegend(0.2,0.25,0.5,0.35);
      legVtx->SetFillColor(0);
      legVtx->SetBorderSize(0);
      legVtx->SetTextSize(0.03);
      if(bDoLowpt)
	{
	  legVtx->AddEntry(pg_jpsi_lowpt_ptybins,"Default: |z_{vtx}|<25 cm","P");
	  legVtx->AddEntry(pg_jpsiVtx10_lowpt_ptybins,"|z_{vtx}|<10 cm","P");
	}
      else
	{
	  legVtx->AddEntry(pg_jpsi_highpt_ptybins,"Default: |z_{vtx}|<25 cm","P");
	  legVtx->AddEntry(pg_jpsiVtx10_highpt_ptybins,"|z_{vtx}|<10 cm","P");
	}

      legVtx->Draw("same");
    }

  
  if(bAddNoFlat)
    {
      pg_jpsiNoFlat_highpt_ptybins->SetMarkerStyle(24);
      pg_jpsiNoFlat_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_highpt_ptybins->SetMarkerColor(kRed+2);
      if(!bDoLowpt)pg_jpsiNoFlat_highpt_ptybins->Draw("p");
   
      pg_jpsiNoFlat_lowpt_ptybins->SetMarkerStyle(24);
      pg_jpsiNoFlat_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_lowpt_ptybins->SetMarkerColor(kBlue+2);
      if(bDoLowpt)pg_jpsiNoFlat_lowpt_ptybins->Draw("p");
      
      
      TLegend *legEp = new TLegend(0.6,0.17,0.9,0.3);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      if(bDoLowpt)
	{
	  legEp->AddEntry(pg_jpsi_lowpt_ptybins,"Default: w/ flattening","P");
	  legEp->AddEntry(pg_jpsiNoFlat_lowpt_ptybins,"No flattening","P");
	}
      else
	{
	  legEp->AddEntry(pg_jpsi_highpt_ptybins,"Default: w/ flattening","P");
	  legEp->AddEntry(pg_jpsiNoFlat_highpt_ptybins,"No flattening","P");
	}

      legEp->Draw("same");
    }


  if(bAddAutoCor)
    {
      pg_jpsiAutoCorr_highpt_ptybins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_highpt_ptybins->SetMarkerColor(kRed+2);
      if(!bDoLowpt)pg_jpsiAutoCorr_highpt_ptybins->Draw("p");

      pg_jpsiAutoCorr_lowpt_ptybins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_lowpt_ptybins->SetMarkerColor(kBlue+2);
      if(bDoLowpt)pg_jpsiAutoCorr_lowpt_ptybins->Draw("p");
      
      TLegend *legAuto = new TLegend(0.6,0.17,0.9,0.3);
      legAuto->SetFillColor(0);
      legAuto->SetBorderSize(0);
      legAuto->SetTextSize(0.03);
           
      if(bDoLowpt)
	{
	  legAuto->AddEntry(pg_jpsi_lowpt_ptybins,"+- && -+","P");
	  legAuto->AddEntry(pg_jpsiAutoCorr_lowpt_ptybins,"++ && --","P");
	}
      else
	{
	  legAuto->AddEntry(pg_jpsi_highpt_ptybins,"+- && -+","P");
	  legAuto->AddEntry(pg_jpsiAutoCorr_highpt_ptybins,"++ && --","P");
	}
      legAuto->Draw("same");
    }

  if(bDoCowboySailor)
    {
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerSize(1.8);
      //
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerSize(1.8);

      TLegend *legTrig = new TLegend(0.7,0.14,0.9,0.26);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      

      if(bDoLowpt)
	{
	  legTrig->AddEntry(pg_jpsi_lowpt_ptybins,"Default: sum","P");
	  pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->Draw("[P]");
	  legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors","P");
	  legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys","P");
	}
      else
	{
	  legTrig->AddEntry(pg_jpsi_highpt_ptybins,"Default: sum","P");
	  pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_highpt_ptybins->Draw("[P]");
	  legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors","P");
	  legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys","P");
	}
	  legTrig->Draw("same");
    }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(bDoCS) 
    {
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerStyle(24);

      pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins->SetMarkerStyle(21);
      pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins->SetMarkerStyle(21);
	
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_highpt_ptybins->SetMarkerSize(1.8);

      //
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerStyle(24);

      pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins->SetMarkerStyle(21);
      pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins->SetMarkerStyle(21);
	
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->SetMarkerSize(1.8);

      if(bDoLowpt)
	{
	  pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_lowpt_ptybins->Draw("[P]");
	}
      else
	{
	  pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_highpt_ptybins->Draw("[P]");
	}
     
 if(bDoEffCorrection_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  
	  TLegend *legEffCS = new TLegend(0.6,0.17,0.9,0.3);
	  legEffCS->SetFillColor(0);
	  legEffCS->SetBorderSize(0);
	  legEffCS->SetTextSize(0.03);
	  
	  if(bDoLowpt)
	    {
	      pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins->Draw("[P]");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_weight_lowpt_ptybins,"Cowboys+Eff","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_weight_lowpt_ptybins,"Sailors+Eff","P");
	    }
	  else
	    {
	      pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins->Draw("[P]");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_weight_highpt_ptybins,"Cowboys+Eff","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors","P");
	      legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_weight_highpt_ptybins,"Sailors+Eff","P");
	    }
	  legEffCS->Draw("same");

	}

      if(bDoMuHits_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  
	  TLegend *legHitsCS = new TLegend(0.6,0.17,0.9,0.3);
	  legHitsCS->SetFillColor(0);
	  legHitsCS->SetBorderSize(0);
	  legHitsCS->SetTextSize(0.03);
	  
	  if(bDoLowpt)
	    {
	      pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins->Draw("[P]");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_lowpt_ptybins,"Cowboys+12 muValidHits","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_nMuValHits12_lowpt_ptybins,"Sailors+12 muValidHits","P");
	    }
	  else
	    {
	      pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins->Draw("[P]");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_highpt_ptybins,"Cowboys+12 muValidHits","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors","P");
	      legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_nMuValHits12_highpt_ptybins,"Sailors+12 muValidHits","P");
	    }
	  legHitsCS->Draw("same");
	}

      if(bDoBarrel_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  
	  TLegend *legBarrelCS = new TLegend(0.6,0.17,0.9,0.3);
	  legBarrelCS->SetFillColor(0);
	  legBarrelCS->SetBorderSize(0);
	  legBarrelCS->SetTextSize(0.03);
	  
	  if(bDoLowpt)
	    {
	      pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins->Draw("[P]");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_lowpt_ptybins,"Cowboys |#eta^{#mu}|<1.2","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_lowpt_ptybins,"Sailors |#eta^{#mu}|<1.2","P");
	    }
	  else
	    {
	      pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins->Draw("[P]");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_highpt_ptybins,"Cowboys |#eta^{#mu}|<1.2","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors","P");
	      legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_highpt_ptybins,"Sailors |#eta^{#mu}|<1.2","P");
	    }

	  legBarrelCS->Draw("same");
	}

      if(bDoVtx_cs)
	{
	  TLegend *legVtxCS = new TLegend(0.6,0.17,0.9,0.3);
	  legVtxCS->SetFillColor(0);
	  legVtxCS->SetBorderSize(0);
	  legVtxCS->SetTextSize(0.03);
	 	  
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  if(bDoLowpt)
	    {
	      pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins->Draw("[P]");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys |z_{vtx}|<25 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_lowpt_ptybins,"Cowboys |z_{vtx}|<10 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors |z_{vtx}|<25 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_noCow_zVtxLT10_lowpt_ptybins,"Sailors |z_{vtx}|<10 cm","P");
	    }
	  else
	    {
	      pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins->Draw("[P]");
	      pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins->Draw("[P]");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys |z_{vtx}|<25 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_highpt_ptybins,"Cowboys |z_{vtx}|<10 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors |z_{vtx}|<25 cm","P");
	      legVtxCS->AddEntry(pg_jpsiL1NHitTrig_noCow_zVtxLT10_highpt_ptybins,"Sailors |z_{vtx}|<10 cm","P");
	    }
	  legVtxCS->Draw("same");

	}


      if(doBkg_cs)
	{
	  pcPadY->SetMaximum(0.5);

	  TLegend *legBkgCS = new TLegend(0.2,0.57,0.4,0.75);
	  legBkgCS->SetFillColor(0);
	  legBkgCS->SetBorderSize(0);
	  legBkgCS->SetTextSize(0.03);
		
	  pg_bkgL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_highpt_ptybins->SetMarkerColor(kMagenta);
	  pg_bkgL1NHitTrig_onlyCow_highpt_ptybins->SetLineColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_highpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_bkgL1NHitTrig_onlyCow_highpt_ptybins->SetMarkerSize(1.8);
	  pg_bkgL1NHitTrig_noCow_highpt_ptybins->SetMarkerSize(1.8);

	  pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_lowpt_ptybins->SetMarkerColor(kMagenta);
	  pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins->SetLineColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_lowpt_ptybins->SetLineColor(kMagenta);
	      
	  pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins->SetMarkerSize(1.8);
	  pg_bkgL1NHitTrig_noCow_lowpt_ptybins->SetMarkerSize(1.8);
	  
	  if(bDoLowpt)
	    {
	      pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins->Draw("[P]");
	      pg_bkgL1NHitTrig_noCow_lowpt_ptybins->Draw("[P]");
	      legBkgCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_lowpt_ptybins,"Cowboys","P");
	      legBkgCS->AddEntry(pg_bkgL1NHitTrig_onlyCow_lowpt_ptybins,"Bkg Cowboys","P");
	      legBkgCS->AddEntry(pg_jpsiL1NHitTrig_noCow_lowpt_ptybins,"Sailors","P");
	      legBkgCS->AddEntry(pg_bkgL1NHitTrig_noCow_lowpt_ptybins,"Bkg Sailors","P");
	    }
	  else
	    {
	      pg_bkgL1NHitTrig_onlyCow_highpt_ptybins->Draw("[P]");
	      pg_bkgL1NHitTrig_noCow_highpt_ptybins->Draw("[P]");
	      legBkgCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_highpt_ptybins,"Cowboys","P");
	      legBkgCS->AddEntry(pg_bkgL1NHitTrig_onlyCow_highpt_ptybins,"Bkg Cowboys","P");
	      legBkgCS->AddEntry(pg_jpsiL1NHitTrig_noCow_highpt_ptybins,"Sailors","P");
	      legBkgCS->AddEntry(pg_bkgL1NHitTrig_noCow_highpt_ptybins,"Bkg Sailors","P");
	    }
	  
	  legBkgCS->Draw("same");
	}
    }// do cowboy-sailor systematics

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(bDoEffCorrection)
    {

      pg_jpsiEffCorr_highpt_ptybins->SetMarkerStyle(27);
      pg_jpsiEffCorr_highpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiEffCorr_highpt_ptybins->SetMarkerColor(kRed+2);
      if(!bDoLowpt)pg_jpsiEffCorr_highpt_ptybins->Draw("p");

      pg_jpsiEffCorr_lowpt_ptybins->SetMarkerStyle(27);
      pg_jpsiEffCorr_lowpt_ptybins->SetMarkerSize(1.8);
      pg_jpsiEffCorr_lowpt_ptybins->SetMarkerColor(kBlue+2);
      if(bDoLowpt)pg_jpsiEffCorr_lowpt_ptybins->Draw("p");
      
      TLegend *legAuto = new TLegend(0.6,0.17,0.9,0.3);
      legAuto->SetFillColor(0);
      legAuto->SetBorderSize(0);
      legAuto->SetTextSize(0.03);
           
      if(bDoLowpt)
	{
	  legAuto->AddEntry(pg_jpsi_lowpt_ptybins,"Default: w/o eff correction","P");
	  legAuto->AddEntry(pg_jpsiEffCorr_lowpt_ptybins,"W/ eff correction","P");
	}
      else
	{
	  legAuto->AddEntry(pg_jpsi_highpt_ptybins,"Default: w/o eff correction","P");
	  legAuto->AddEntry(pg_jpsiEffCorr_highpt_ptybins,"W/ eff correction","P");
	}
      legAuto->Draw("same");

      // some sanity check on what's in these plost
       double y,x2;
       int bin1 = pg_jpsiEffCorr_highpt_ptybins->GetPoint(0,x2,y);
       cout<<"bin1 "<< x2<<"\t"<<y<<endl;

       int bin2 = pg_jpsiEffCorr_highpt_ptybins->GetPoint(1,x2,y);
       cout<<"bin2 "<< x2<<"\t"<<y<<endl;


         int bin3 = pg_jpsiEffCorr_highpt_ptybins->GetPoint(2,x2,y);
       cout<<"bin3 "<< x2<<"\t"<<y<<endl;

    //    int bin1 = pg_jpsi_highpt_ptybins->GetPoint(0,x2,y);
//        cout<<"bin1 "<< x2<<"\t"<<y<<endl;

//        int bin2 = pg_jpsi_highpt_ptybins->GetPoint(1,x2,y);
//        cout<<"bin2 "<< x2<<"\t"<<y<<endl;


//        int bin3 = pg_jpsi_highpt_ptybins->GetPoint(2,x2,y);
//        cout<<"bin3 "<< x2<<"\t"<<y<<endl;
    }

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("Cent. 10 - 60 %%")); 

  TLegend *leg1 = new TLegend(0.2,0.15,0.5,0.25);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);

  
  switch(choseSignal){
  case 1:
    if(!bDoLowpt) leg1->AddEntry(pg_jpsi_highpt_ptybins_statErr,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2]),"P");
    else leg1->AddEntry(pg_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    
    if(bDoSummary)leg1->AddEntry(pg_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    leg1->Draw("same");
    break;
  case 2:
    if(!bDoLowpt) leg1->AddEntry(pg_pr_jpsi_highpt_ptybins_statErr,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2]),"P");
    else leg1->AddEntry(pg_pr_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    
    if(bDoSummary)leg1->AddEntry(pg_pr_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    leg1->Draw("same");
    break;
  case 3:
    if(!bDoLowpt) leg1->AddEntry(pg_npr_jpsi_highpt_ptybins_statErr,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2]),"P");
    else leg1->AddEntry(pg_npr_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    
    if(bDoSummary)leg1->AddEntry(pg_npr_jpsi_lowpt_ptybins_statErr,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    leg1->Draw("same");
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }


 


  if(bSavePlots)
  {
    if(bDoCS)
	{
	  pcY->SaveAs(Form("figs/png/%s_%s_v2cs_ptyBins.png",chosenSignal,outputName[0]));
	  pcY->SaveAs(Form("figs/pdf/%s_%s_v2cs_ptyBins.pdf",chosenSignal,outputName[0]));
	}
      else
	{
	  pcY->SaveAs(Form("figs/png/%s_%s_v2_ptyBins.png",chosenSignal,outputName[0]));
	  pcY->SaveAs(Form("figs/pdf/%s_%s_v2_ptyBins.pdf",chosenSignal,outputName[0]));
	}
  }
  //________________________________________________________________________________________

}
