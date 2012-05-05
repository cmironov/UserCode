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


void v2SummaryPlots_pt_EPCorr()
{
  //  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gROOT->Macro("./rootlogon.C");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]   = {"","Inclusive J/#psi", "Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal         = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominal_bit1_pt"};
   
  // options
  bool bSavePlots      = true;
  bool bDoSummary = false;
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

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  int centIntegrated[2]   = {10, 60}; 
  
  TFile *f0;
  TFile *f1;
  if(choseSignal==1)
    {
      f0 = new TFile(Form("./a2_corrV2.root"));
      if(!f0->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
      f1 = new TFile(Form("./extrastud/a2_corrV2.root"));
      if(!f1->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
    }
  else
    {
      f0 = new TFile(Form("./a2_corrV2_pnp.root"));
      if(!f0->IsOpen()) { cout << "cannot open a1_corrV2_pnp.root" << endl; return;}
      // f1 = new TFile(Form("./extrastud/a1_corrV2_pnp.root"));
      // if(!f1->IsOpen()) { cout << "cannot open a1_corrV2_pnp.root" << endl; return;}
    }
  f1 = new TFile(Form("./extrastud/a2_corrV2.root"));

  char histname[200];
  //inclusive
  sprintf(histname,"final_ptDependence_NSig_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsi_ptbins = (TGraphErrors*)f0->Get(histname);

 sprintf(histname,"final_ptDependence_statErr_NSig_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsi_ptbins_statErr = (TGraphErrors*)f0->Get(histname);
  
  // prompt
  sprintf(histname,"final_ptDependence_NPr_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_ptbins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_ptDependence_statErr_NPr_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_ptbins_statErr = (TGraphErrors*)f0->Get(histname);

  // non-prompt
  sprintf(histname,"final_ptDependence_NNp_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_ptbins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"final_ptDependence_statErr_NNp_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_ptbins_statErr = (TGraphErrors*)f0->Get(histname);

  // systematic studies
  const int nptbins                      = 3;
  double ptbins_center[nptbins]          = {7.3, 9.0, 13.4}; // <pt> values in the AN
  double ptbins_center_err[nptbins]      = {0., 0., 0.};
  double ptbins_center_errGhost[nptbins] = {0., 0., 0.};
  double jpsi_ptbins[nptbins]            = {0.075, 0.1, 0.043};
  double jpsi_ptbins_err[nptbins]        = {0.022, 0.022, 0.020};
  double jpsi_ptbins_errGhost[nptbins]   = {0.,     0.,     0., };

  TGraphErrors *pg_jpsi_ptbinsGhost = new TGraphErrors(nptbins,ptbins_center,jpsi_ptbins,ptbins_center_err,jpsi_ptbins_errGhost);

  // ### fitting
  sprintf(histname,"default_constrained_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_constrained_ptbins = (TGraphErrors*)f0->Get(histname);

  sprintf(histname,"default_polFunct_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_polFunct_ptbins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_signalCB3WN_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_signalCB3WN_ptbins = (TGraphErrors*)f0->Get(histname);

    //------ special extra for p-np, 
  sprintf(histname,"default_bit1_1GaussResol_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_1GaussResol_ptbins = (TGraphErrors*)f0->Get(histname);
  
  sprintf(histname,"default_bit1_ResolFixToPRMC_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_ResolFixToPRMC_ptbins = (TGraphErrors*)f0->Get(histname);


  // bkf
  sprintf(histname,"default_bit1_NBkg_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_bkg_ptbins = (TGraphErrors*)f0->Get(histname);

  //vtx
  sprintf(histname,"zVtxLT10_bit1_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiVtx10_ptbins = (TGraphErrors*)f0->Get(histname);
  
  // flat no flat
  sprintf(histname,"noFlat_bit1_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiNoFlat_ptbins = (TGraphErrors*)f0->Get(histname);

 // with efficiency correction
  sprintf(histname,"default_bit1_weight_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiEffCorr_ptbins = (TGraphErrors*)f0->Get(histname);

  // extra studies
  // 6bins : compare among themselves the v2 obtianed with stat errors only
  TFile *f2 = new TFile(Form("./a2_corrV2_6bin.root"));
   sprintf(histname,"default_bit1_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
   TGraphErrors *pg_jpsi6bin_ptbins = (TGraphErrors*)f2->Get(histname);
  TGraphErrors *pg_jpsi4bin_ptbins = (TGraphErrors*)f0->Get(histname);
  // if (!pg_jpsi4bin_ptbins) { cout << "cannot load 4bincase." << endl; return;}
  // if (!pg_jpsi6bin_ptbins) { cout << "cannot load 6bincase." << endl; return;}


  // enhanced autocorr
  sprintf(histname,"autoCorr_bit1_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiAutoCorr_ptbins = (TGraphErrors*)f1->Get(histname);
 
 

  // all triggers
  sprintf(histname,"default_bit1_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_ptbins = (TGraphErrors*)f1->Get(histname);

  // ######################## sailors and cowboys
  sprintf(histname,"default_cowboy_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_ptbins = (TGraphErrors*)f0->Get(histname);
  sprintf(histname,"default_sailor_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_ptbins = (TGraphErrors*)f0->Get(histname);

  // with corrections
  sprintf(histname,"default_cowboy_weight_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_weight_ptbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_weight_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_weight_ptbins = (TGraphErrors*)f1->Get(histname);


  // mu valid hits
  sprintf(histname,"nMuValHits12_cowboy_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"nMuValHits12_sailor_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins = (TGraphErrors*)f1->Get(histname);


  // zvtx
  sprintf(histname,"zVtxLT10_cowboy_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"zVtxLT10_sailor_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins = (TGraphErrors*)f1->Get(histname);

  // midrapidity
 sprintf(histname,"singleMuLTeta1.2_cowboy_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"singleMuLTeta1.2_sailor_%s_rap%.1f-%.1f_cent%d-%d",chosenSignal,rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins = (TGraphErrors*)f1->Get(histname);

  // cowboy bkg
    sprintf(histname,"default_cowboy_NBkg_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_onlyCow_ptbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"default_sailor_NBkg_rap%.1f-%.1f_cent%d-%d",rapIntegrated[0],rapIntegrated[1],centIntegrated[0],centIntegrated[1]);
  TGraphErrors *pg_bkgL1NHitTrig_noCow_ptbins = (TGraphErrors*)f1->Get(histname);


  //_____________________________________________________________________________
  // pt dependence
  TCanvas *pcPt = new TCanvas("pcPt","pcPt");
  pcPt->cd();
  TH1F *pcPadPT = new TH1F("pcPadPT",";p_{T} (GeV/c);v_{2};",100,0,40);
  pcPadPT->GetXaxis()->SetLabelSize(20);
  pcPadPT->GetXaxis()->SetLabelFont(43);
  pcPadPT->GetXaxis()->SetTitleSize(27);
  pcPadPT->GetXaxis()->SetTitleFont(43);
  pcPadPT->GetXaxis()->SetTitleOffset(1.2);
  pcPadPT->GetXaxis()->CenterTitle();

  pcPadPT->GetYaxis()->SetLabelSize(20);
  pcPadPT->GetYaxis()->SetLabelFont(43);
  pcPadPT->GetYaxis()->SetTitleSize(32);
  pcPadPT->GetYaxis()->SetTitleFont(43);
  pcPadPT->GetYaxis()->SetTitleOffset(1.1);
  pcPadPT->GetYaxis()->CenterTitle();

  pcPadPT->SetMaximum(0.25);
  pcPadPT->SetMinimum(-0.1);
  if(bAddBkg)
    {
      pcPadPT->SetMaximum(0.4);
    }

  if(bDoSummary)
    {
      pcPadPT->SetMaximum(0.3);
      pcPadPT->SetMinimum(-0.2);
    }
  pcPadPT->Draw();

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
  //________________________________________ 

  pg_jpsi_ptbinsGhost->SetMarkerStyle(20);
  pg_jpsi_ptbinsGhost->SetMarkerSize(1.8);
  pg_jpsi_ptbinsGhost->SetMarkerColor(602);
  // pg_jpsi_ptbinsGhost->Draw("P");

  switch(choseSignal){
  case 1:
    //pg_jpsi_ptbins->SetFillStyle(0);
    pg_jpsi_ptbins->SetFillColor(kBlue-10);
  
    pg_jpsi_ptbins_statErr->SetMarkerStyle(20);
    pg_jpsi_ptbins_statErr->SetMarkerSize(1.8);
    pg_jpsi_ptbins_statErr->SetMarkerColor(kBlue+2);
    pg_jpsi_ptbins_statErr->SetLineColor(kBlue+2);
    if(!bDoCS&& !bDo6Bins)
      {
	pg_jpsi_ptbins->Draw("2");
	pg_jpsi_ptbins_statErr->Draw("[P]");
      }
    break;
  case 2:
    pg_pr_jpsi_ptbins->SetFillColor(kRed-10);

    pg_pr_jpsi_ptbins_statErr->SetMarkerStyle(21);
    pg_pr_jpsi_ptbins_statErr->SetMarkerSize(1.8);
    pg_pr_jpsi_ptbins_statErr->SetMarkerColor(kRed+2);
    pg_pr_jpsi_ptbins_statErr->SetLineColor(kRed+2);
    if(!bDoCS && !bDo6Bins)
      {
	pg_pr_jpsi_ptbins->Draw("2");
	pg_pr_jpsi_ptbins_statErr->Draw("[P]");
      }
    break;
  case 3:
    pg_npr_jpsi_ptbins->SetFillColor(kOrange-10);

    pg_npr_jpsi_ptbins_statErr->SetMarkerStyle(33);
    pg_npr_jpsi_ptbins_statErr->SetMarkerSize(1.8);
    pg_npr_jpsi_ptbins_statErr->SetMarkerColor(kOrange+2);
    pg_npr_jpsi_ptbins_statErr->SetLineColor(kOrange+2);
    if(!bDoCS && !bDo6Bins)
      {
	pg_npr_jpsi_ptbins->Draw("2");
	pg_npr_jpsi_ptbins_statErr->Draw("[P]");
      }
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }
  
  if(bDo6Bins)
    {
      pg_jpsi6bin_ptbins->SetMarkerStyle(24);
      pg_jpsi6bin_ptbins->SetMarkerSize(1.8);
      pg_jpsi6bin_ptbins->SetMarkerColor(kBlue+2);
      pg_jpsi6bin_ptbins->SetLineColor(kBlue+2);
      
      pg_jpsi4bin_ptbins->SetMarkerStyle(20);
      pg_jpsi4bin_ptbins->SetMarkerSize(1.8);
      pg_jpsi4bin_ptbins->SetMarkerColor(kBlue+2);
      pg_jpsi4bin_ptbins->SetLineColor(kBlue+2);
      
      pg_jpsi4bin_ptbins->Draw("[P]");
      pg_jpsi6bin_ptbins->Draw("[P]");
      
      TLegend *leg6Bin = new TLegend(0.5,0.6,0.9,0.7);
      leg6Bin->SetFillColor(0);
      leg6Bin->SetBorderSize(0);
      leg6Bin->SetTextSize(0.03);
      leg6Bin->AddEntry(pg_jpsi4bin_ptbins,"Default:4 d#phi bins","P");
      leg6Bin->AddEntry(pg_jpsi6bin_ptbins,"6 d#phi bins","P");
      leg6Bin->Draw("same");
    }


  if(bDoFitting)
    {
      pg_constrained_ptbins->SetMarkerStyle(24);
      pg_polFunct_ptbins->SetMarkerStyle(25);
      pg_signalCB3WN_ptbins->SetMarkerStyle(27);

      pg_constrained_ptbins->SetMarkerColor(kBlue);
      pg_polFunct_ptbins->SetMarkerColor(kBlue);
      pg_signalCB3WN_ptbins->SetMarkerColor(kBlue);

      pg_constrained_ptbins->SetLineColor(kBlue);
      pg_polFunct_ptbins->SetLineColor(kBlue);
      pg_signalCB3WN_ptbins->SetLineColor(kBlue);

      pg_constrained_ptbins->SetMarkerSize(1.8);
      pg_polFunct_ptbins->SetMarkerSize(1.8);
      pg_signalCB3WN_ptbins->SetMarkerSize(1.9);

      pg_constrained_ptbins->Draw("[P]");
      pg_polFunct_ptbins->Draw("[P]");
      pg_signalCB3WN_ptbins->Draw("[P]");
     
      TLegend *legFit = new TLegend(0.6,0.55,0.9,0.7);
      legFit->SetFillColor(0);
      legFit->SetBorderSize(0);
      legFit->SetTextSize(0.03);
      legFit->AddEntry(pg_jpsi_ptbins,"Nominal","P");
      legFit->AddEntry(pg_constrained_ptbins,"Constrained","P");
      legFit->AddEntry(pg_signalCB3WN_ptbins,"Sgn: C-B only","P");
      legFit->AddEntry(pg_polFunct_ptbins,"Bkg:straight line ","P");

      if(choseSignal!=1)
	{
	  // for p=np only
	  pg_1GaussResol_ptbins->SetMarkerStyle(28);
	  pg_ResolFixToPRMC_ptbins->SetMarkerStyle(30);
	  
	  pg_1GaussResol_ptbins->SetMarkerColor(kBlue);
	  pg_ResolFixToPRMC_ptbins->SetMarkerColor(kBlue);
	  pg_1GaussResol_ptbins->SetLineColor(kBlue);
	  pg_ResolFixToPRMC_ptbins->SetLineColor(kBlue);
	  
	  pg_1GaussResol_ptbins->SetMarkerSize(1.8);
	  pg_ResolFixToPRMC_ptbins->SetMarkerSize(1.8);
	  
	  pg_1GaussResol_ptbins->Draw("[P]");
	  pg_ResolFixToPRMC_ptbins->Draw("[P]");
	  legFit->AddEntry(pg_1GaussResol_ptbins,"Sgn Res: 1 Gauss","P");
	  legFit->AddEntry(pg_ResolFixToPRMC_ptbins,"Sgn Res: fixed to MC ","P");
	}

      legFit->Draw("same");
    }

  if(bAddBkg)
    {
      pg_bkg_ptbins->SetMarkerStyle(24);
      pg_bkg_ptbins->SetMarkerSize(1.8);
      pg_bkg_ptbins->SetMarkerColor(kBlue);
      pg_bkg_ptbins->SetLineColor(kBlue);
      pg_bkg_ptbins->Draw("[P]");
      TLegend *legCent = new TLegend(0.6,0.4,0.9,0.6);
      legCent->SetFillColor(0);
      legCent->SetBorderSize(0);
      legCent->SetTextSize(0.03);
      legCent->AddEntry(pg_jpsi_ptbins,"Signal","P");
      legCent->AddEntry(pg_bkg_ptbins,"Bkg","P");
      legCent->Draw("same");
    }

 if(bAddVtx)
    {
      pg_jpsiVtx10_ptbins->SetMarkerStyle(24);
      pg_jpsiVtx10_ptbins->SetMarkerSize(1.8);
      pg_jpsiVtx10_ptbins->SetMarkerColor(kBlue);
      pg_jpsiVtx10_ptbins->SetLineColor(kBlue);
      pg_jpsiVtx10_ptbins->Draw("[P]");
      
      TLegend *legVtx = new TLegend(0.6,0.4,0.9,0.6);
      legVtx->SetFillColor(0);
      legVtx->SetBorderSize(0);
      legVtx->SetTextSize(0.03);
      legVtx->AddEntry(pg_jpsi_ptbins,"Default: |z_{vtx}|<25 cm","P");
      legVtx->AddEntry(pg_jpsiVtx10_ptbins,"|z_{vtx}|<10 cm","P");
      legVtx->Draw("same");
    }

  
  if(bAddNoFlat)
    {
      pg_jpsiNoFlat_ptbins->SetMarkerStyle(24);
      pg_jpsiNoFlat_ptbins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_ptbins->SetMarkerColor(kBlue);
      pg_jpsiNoFlat_ptbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.6,0.4,0.9,0.6);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      legEp->AddEntry(pg_jpsi_ptbins,"Default: w/ flattening","P");
      legEp->AddEntry(pg_jpsiNoFlat_ptbins,"No flattening","P");
      
      legEp->Draw("same");
    }
  if(bAddAutoCor)
    {
      pg_jpsiAutoCorr_ptbins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_ptbins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_ptbins->SetMarkerColor(kBlue);
      pg_jpsiAutoCorr_ptbins->Draw("p");
      
      TLegend *legAuto = new TLegend(0.6,0.4,0.9,0.6);
      legAuto->SetFillColor(0);
      legAuto->SetBorderSize(0);
      legAuto->SetTextSize(0.03);
           
      legAuto->AddEntry(pg_jpsi_ptbins,"+- && -+","P");
      legAuto->AddEntry(pg_jpsiAutoCorr_ptbins,"++ && --","P");
      
      legAuto->Draw("same");
    }

  if(bDoCowboySailor)
    {
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_ptbins->Draw("[P]");
      pg_jpsiL1NHitTrig_noCow_ptbins->Draw("[P]");
      
      TLegend *legTrig = new TLegend(0.6,0.4,0.9,0.6);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      legTrig->AddEntry(pg_jpsi_ptbins,"Default: sum","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys","P");
      legTrig->Draw("same");
    }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(bDoCS) 
    {
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerStyle(24);

      pg_jpsiL1NHitTrig_onlyCow_weight_ptbins->SetMarkerStyle(21);
      pg_jpsiL1NHitTrig_noCow_weight_ptbins->SetMarkerStyle(21);
	
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_ptbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_ptbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_ptbins->Draw("[P]");
      pg_jpsiL1NHitTrig_noCow_ptbins->Draw("[P]");

      if(bDoEffCorrection_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_weight_ptbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_ptbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_weight_ptbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_weight_ptbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_weight_ptbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_weight_ptbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_weight_ptbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_weight_ptbins->Draw("[P]");
	  
	  TLegend *legEffCS = new TLegend(0.6,0.4,0.9,0.6);
	  legEffCS->SetFillColor(0);
	  legEffCS->SetBorderSize(0);
	  legEffCS->SetTextSize(0.03);
	  
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_weight_ptbins,"Cowboys+Eff","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors","P");
	  legEffCS->AddEntry(pg_jpsiL1NHitTrig_noCow_weight_ptbins,"Sailors+Eff","P");
	  
	  legEffCS->Draw("same");

	}

      if(bDoMuHits_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins->Draw("[P]");
	  
	  TLegend *legHitsCS = new TLegend(0.6,0.4,0.9,0.6);
	  legHitsCS->SetFillColor(0);
	  legHitsCS->SetBorderSize(0);
	  legHitsCS->SetTextSize(0.03);
	  
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_nMuValHits12_ptbins,"Cowboys+12 muValidHits","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors","P");
	  legHitsCS->AddEntry(pg_jpsiL1NHitTrig_noCow_nMuValHits12_ptbins,"Sailors+12 muValidHits","P");
	  
	  legHitsCS->Draw("same");

	}

      if(bDoBarrel_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins->Draw("[P]");
	  
	  TLegend *legBarrelCS = new TLegend(0.6,0.4,0.9,0.6);
	  legBarrelCS->SetFillColor(0);
	  legBarrelCS->SetBorderSize(0);
	  legBarrelCS->SetTextSize(0.03);
	  
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_singleMuLTeta12_ptbins,"Cowboys |#eta^{#mu}|<1.2","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_singleMuLTeta12_ptbins,"Sailorsin |#eta^{#mu}|<1.2","P");
	  
	  legBarrelCS->Draw("same");
	}

      if(bDoVtx_cs)
	{
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins->SetMarkerColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins->SetMarkerColor(kMagenta);
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins->SetLineColor(kGreen+2);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins->SetLineColor(kMagenta);
	      
	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins->SetMarkerSize(1.8);
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins->SetMarkerSize(1.8);

	  pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins->Draw("[P]");
	  pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins->Draw("[P]");
	  
	  TLegend *legBarrelCS = new TLegend(0.6,0.4,0.9,0.6);
	  legBarrelCS->SetFillColor(0);
	  legBarrelCS->SetBorderSize(0);
	  legBarrelCS->SetTextSize(0.03);
	  
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys |z_{vtx}|<25 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_zVtxLT10_ptbins,"Cowboys |z_{vtx}|<10 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors |z_{vtx}|<25 cm","P");
	  legBarrelCS->AddEntry(pg_jpsiL1NHitTrig_noCow_zVtxLT10_ptbins,"Sailors |z_{vtx}|<10 cm","P");
	  
	  legBarrelCS->Draw("same");
	}
      if(doBkg_cs)
	{
	  pg_bkgL1NHitTrig_onlyCow_ptbins->SetMarkerColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_ptbins->SetMarkerColor(kMagenta);
	  pg_bkgL1NHitTrig_onlyCow_ptbins->SetLineColor(kGreen+2);
	  pg_bkgL1NHitTrig_noCow_ptbins->SetLineColor(kMagenta);
	      
	  pg_bkgL1NHitTrig_onlyCow_ptbins->SetMarkerSize(1.8);
	  pg_bkgL1NHitTrig_noCow_ptbins->SetMarkerSize(1.8);

	  pg_bkgL1NHitTrig_onlyCow_ptbins->Draw("[P]");
	  pg_bkgL1NHitTrig_noCow_ptbins->Draw("[P]");
	  
	  TLegend *legBkgCS = new TLegend(0.6,0.4,0.9,0.6);
	  legBkgCS->SetFillColor(0);
	  legBkgCS->SetBorderSize(0);
	  legBkgCS->SetTextSize(0.03);
	  
	  legBkgCS->AddEntry(pg_jpsiL1NHitTrig_onlyCow_ptbins,"Cowboys","P");
	  legBkgCS->AddEntry(pg_bkgL1NHitTrig_onlyCow_ptbins,"Bkg Cowboys","P");
	  legBkgCS->AddEntry(pg_jpsiL1NHitTrig_noCow_ptbins,"Sailors","P");
	  legBkgCS->AddEntry(pg_bkgL1NHitTrig_noCow_ptbins,"Bkg Sailors","P");
	  
	  legBkgCS->Draw("same");
	}
    }// do cowboy-sailor systematics

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(bDoEffCorrection)
    {
      pg_jpsiEffCorr_ptbins->SetMarkerStyle(24);
      pg_jpsiEffCorr_ptbins->SetMarkerSize(1.8);
      pg_jpsiEffCorr_ptbins->SetMarkerColor(kBlue);
      pg_jpsiEffCorr_ptbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.6,0.55,0.9,0.7);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      legEp->AddEntry(pg_jpsi_ptbins,"Default: w/o eff correction","P");
      legEp->AddEntry(pg_jpsiEffCorr_ptbins,"W/ eff correction","P");
      
      legEp->Draw("same");

    }

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",rapIntegrated[1]));       // rapidity
  lt1->DrawLatex(0.18,0.77,Form("Cent. 10 - 60 %%")); 
  //  lt1->DrawLatex(0.18,0.71,Form("%s",eventPlane[1]));
  
  if(bSavePlots)
    {
      gSystem->mkdir(Form("./figs/png/"),kTRUE);
      gSystem->mkdir(Form("./figs/pdf/"),kTRUE);
      

       if(bDoCS)
	{
	  pcPt->SaveAs(Form("figs/png/%s_%s_v2cs_ptBins.png",chosenSignal,outputName[0]));
	  pcPt->SaveAs(Form("figs/pdf/%s_%s_v2cs_ptBins.pdf",chosenSignal,outputName[0]));
	}
      else
	{
	  pcPt->SaveAs(Form("figs/png/%s_%s_v2_ptBins.png",chosenSignal,outputName[0]));
	  pcPt->SaveAs(Form("figs/pdf/%s_%s_v2_ptBins.pdf",chosenSignal,outputName[0]));
	}

    }

}

