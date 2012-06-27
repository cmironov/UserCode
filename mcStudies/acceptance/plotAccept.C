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
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
//#include <THIDebug.h>
#include <TStyle.h>
#include <TInterpreter.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>
#endif

void plotAccept()
{
  gROOT->Macro("/afs/cern.ch/user/m/mironov/utilities/setStyle.C+"); 

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.08);

  int bin2 = 25;
  int rebin = 8;
   
  TFile *f1 = new TFile("acc_z0.root");
  TFile *f2 = new TFile("acc_z0photon.root");
  TFile *f3 = new TFile("acc_z0jet.root");
  TFile *f4 = new TFile("acc_z0photonjet.root");
  TFile *f5 = new TFile("acc_powheg.root");

  TH1D *h1n = (TH1D*)f1->Get("NeuDimuPt");
  TH1D *h2n = (TH1D*)f2->Get("NeuDimuPt");
  TH1D *h3n = (TH1D*)f3->Get("NeuDimuPt");
  TH1D *h4n = (TH1D*)f4->Get("NeuDimuPt");
  TH1D *h5n = (TH1D*)f5->Get("NeuDimuPt");

  TH1D *h1d = (TH1D*)f1->Get("DenDimuPt");
  TH1D *h2d = (TH1D*)f2->Get("DenDimuPt");
  TH1D *h3d = (TH1D*)f3->Get("DenDimuPt");
  TH1D *h4d = (TH1D*)f4->Get("DenDimuPt");
  TH1D *h5d = (TH1D*)f5->Get("DenDimuPt");

  TH1D *h1d2 = (TH1D*)f1->Get("Den2DimuPt");
  TH1D *h2d2 = (TH1D*)f2->Get("Den2DimuPt");
  TH1D *h3d2 = (TH1D*)f3->Get("Den2DimuPt");
  TH1D *h4d2 = (TH1D*)f4->Get("Den2DimuPt");
  TH1D *h5d2 = (TH1D*)f5->Get("Den2DimuPt");

    //inv mass plots
  TH1D *h1M = (TH1D*)f1->Get("invmass_histo");
  TH1D *h2M = (TH1D*)f2->Get("invmass_histo");
  TH1D *h3M = (TH1D*)f3->Get("invmass_histo");
  TH1D *h4M = (TH1D*)f4->Get("invmass_histo");
  TH1D *h5M = (TH1D*)f5->Get("invmass_histo");
  
  TH1D *h1nY = (TH1D*)f1->Get("NeuDimuY");
  TH1D *h2nY = (TH1D*)f2->Get("NeuDimuY");
  TH1D *h3nY = (TH1D*)f3->Get("NeuDimuY");
  TH1D *h4nY = (TH1D*)f4->Get("NeuDimuY");
  TH1D *h5nY = (TH1D*)f5->Get("NeuDimuY");
   
  TH1D *h1dY = (TH1D*)f1->Get("DenDimuY");
  TH1D *h2dY = (TH1D*)f2->Get("DenDimuY");
  TH1D *h3dY = (TH1D*)f3->Get("DenDimuY");
  TH1D *h4dY = (TH1D*)f4->Get("DenDimuY");
  TH1D *h5dY = (TH1D*)f5->Get("DenDimuY");  

  TH1D *h1d2Y = (TH1D*)f1->Get("Den2DimuY");
  TH1D *h2d2Y = (TH1D*)f2->Get("Den2DimuY");
  TH1D *h3d2Y = (TH1D*)f3->Get("Den2DimuY");
  TH1D *h4d2Y = (TH1D*)f4->Get("Den2DimuY");
  TH1D *h5d2Y = (TH1D*)f5->Get("Den2DimuY");
  
  TH1D *h1aPt = (TH1D*)f1->Get("AllDimuPt");
  TH1D *h2aPt = (TH1D*)f2->Get("AllDimuPt");
  TH1D *h3aPt = (TH1D*)f3->Get("AllDimuPt");
  TH1D *h4aPt = (TH1D*)f4->Get("AllDimuPt");
  TH1D *h5aPt = (TH1D*)f5->Get("AllDimuPt");
  
  
  TH1D *h1aY = (TH1D*)f1->Get("AllDimuY");
  TH1D *h2aY = (TH1D*)f2->Get("AllDimuY");
  TH1D *h3aY = (TH1D*)f3->Get("AllDimuY");
  TH1D *h4aY = (TH1D*)f4->Get("AllDimuY");
  TH1D *h5aY = (TH1D*)f5->Get("AllDimuY");

    ///// Invariant mass histo
  h1M->SetLineColor(kRed);
  h2M->SetLineColor(kBlue);
  h3M->SetLineColor(kGreen+2);
  h4M->SetLineColor(kOrange+1);

  h1n->SetLineColor(kRed);
  h2n->SetLineColor(kBlue);
  h3n->SetLineColor(kGreen+2);
  h4n->SetLineColor(kOrange+1);

  h1aPt->SetLineColor(kRed);
  h2aPt->SetLineColor(kBlue);
  h3aPt->SetLineColor(kGreen+2);
  h4aPt->SetLineColor(kOrange+1);

  h1aY->SetLineColor(kRed);
  h2aY->SetLineColor(kBlue);
  h3aY->SetLineColor(kGreen+2);
  h4aY->SetLineColor(kOrange+1);

  // nominator
 //  h1n->SetLineColor(kViolet+6);
//   h2n->SetLineColor(kViolet+6);
//   h3n->SetLineColor(kViolet+6);
//   h4n->SetLineColor(kViolet+6);
//   h5n->SetLineColor(kViolet+6);

  h1nY->SetLineColor(kViolet+6);
  h2nY->SetLineColor(kViolet+6);
  h3nY->SetLineColor(kViolet+6);
  h4nY->SetLineColor(kViolet+6);
  h5nY->SetLineColor(kViolet+6);

  // denominator
  h1d->SetLineColor(kAzure+6);
  h2d->SetLineColor(kAzure+6);
  h3d->SetLineColor(kAzure+6);
  h4d->SetLineColor(kAzure+6);
  h5d->SetLineColor(kAzure+6);

  h1dY->SetLineColor(kAzure+6);
  h2dY->SetLineColor(kAzure+6);
  h3dY->SetLineColor(kAzure+6);
  h4dY->SetLineColor(kAzure+6);
  h5dY->SetLineColor(kAzure+6);


  TH1 *phMassAx = new TH1D("phMassAx",";M [GeV/c^{2}];",1,0,150);
  phMassAx->SetMaximum(1e+4);
  phMassAx->SetMinimum(1);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLogy();
  phMassAx->Draw();
  h1M->Draw("same");
  h2M->Draw("same");
  h3M->Draw("same");
  h4M->Draw("same");
  h5M->Draw("same"); 
  TLegend * leg1 = new TLegend(0.7,0.70,0.9,0.95);
  leg1->SetFillColor(0);
 
  leg1->AddEntry(h1M,Form("Z"),"pl");
  leg1->AddEntry(h2M,Form("Z/#gamma^{*}"),"pl");
  leg1->AddEntry(h3M,Form("Z+jet"),"pl");
  leg1->AddEntry(h4M,Form("Z/#gamma^{*}+jet"),"pl");
  leg1->AddEntry(h5M,Form("POWHEG"),"pl");
  leg1->Draw("same");
  lt1->DrawLatex(0.5,0.2,"All generated");
  c1->SaveAs("invMass.gif");

  /*
    ////// Y
  TH1 *phYAx = new TH1D("phYAx",";y;",1,-3.0,3.0);
  phYAx->SetMaximum(250);
  phYAx->SetMinimum(1);
  TCanvas *c1y = new TCanvas("c1y","c1y",1200,350);
  c1y->Divide(5,1);
  c1y->cd(1);
  phYAx->Draw();
  gPad->SetLogy();

  h1aY->Draw("same");
  h1dY->Draw("same");
  h1nY->Draw("same");

  lt1->DrawLatex(0.5,0.2,"Z0");

  c1y->cd(2);
  phYAx->Draw();
  gPad->SetLogy();
  h2aY->Draw("same");
  h2dY->Draw("same");
  h2nY->Draw("same");
  
  lt1->DrawLatex(0.5,0.2,"Z0/#gamma^{*}");
  c1y->cd(3);
  phYAx->Draw();
  gPad->SetLogy();
  h3aY->Draw("same");
  h3dY->Draw("same");
  h3nY->Draw("same");
  lt1->DrawLatex(0.5,0.2,"Z0+jet");

  c1y->cd(4);
  phYAx->Draw();
  gPad->SetLogy();
  h4aY->Draw("same");
  h4dY->Draw("same");
  h4nY->Draw("same");

  lt1->DrawLatex(0.5,0.2,"Z0/#gamma^{*}+jet");

  c1y->cd(5);
  phYAx->Draw();
  gPad->SetLogy();
  h5aY->Draw("same");
  h5dY->Draw("same");
  h5nY->Draw("same");

  lt1->DrawLatex(0.5,0.2,"POWHEG");

  c1y->SaveAs("y.gif");

  */
//


  ////// Pt
  TH1 *phPtAx = new TH1D("phPtAx",";p_{T}[GeV/c^{2}];",1,0,60);
  phPtAx->SetMaximum(1e+3);
  phPtAx->SetMinimum(1);


  TCanvas *c11 = new TCanvas("c11","c11");
  phPtAx->Draw();
  gPad->SetLogy();
  h1n->Draw("same");
  h2n->Draw("same");
  h3n->Draw("same");
  h4n->Draw("same");
  h5n->Draw("same"); 
  TLegend * leg11 = new TLegend(0.7,0.70,0.9,0.95);
  leg11->SetFillColor(0);
 
  leg11->AddEntry(h1n,Form("Z"),"pl");
  leg11->AddEntry(h2n,Form("Z/#gamma^{*}"),"pl");
  leg11->AddEntry(h3n,Form("Z+jet"),"pl");
  leg11->AddEntry(h4n,Form("Z/#gamma^{*}+jet"),"pl");
  leg11->AddEntry(h5n,Form("POWHEG"),"pl");
  leg11->Draw("same");

  lt1->DrawLatex(0.5,0.2,"In acceptance");

  c11->SaveAs("acc_ptAll.gif");

  /*

  TCanvas *c1pt = new TCanvas("c1pt","c1pt",1200,350);
  c1pt->Divide(5,1);
  c1pt->cd(1);
  phPtAx->Draw();
  gPad->SetLogy();
  h1aPt->Draw("same");
  h1d->Draw("same");
  h1n->Draw("same");

  lt1->DrawLatex(0.5,0.2,"Z0");

  c1pt->cd(2);
  phPtAx->Draw();
  gPad->SetLogy();
  h2aPt->Draw("same");
  h2d->Draw("same");
  h2n->Draw("same");
  
  lt1->DrawLatex(0.5,0.2,"Z0/#gamma^{*}");
  c1pt->cd(3);
  phPtAx->Draw();
  gPad->SetLogy();
  h3aPt->Draw("same");
  h3d->Draw("same");
  h3n->Draw("same");
  lt1->DrawLatex(0.5,0.2,"Z0+jet");

  c1pt->cd(4);
  phPtAx->Draw();
  gPad->SetLogy();
  h4aPt->Draw("same");
  h4d->Draw("same");
  h4n->Draw("same");

  lt1->DrawLatex(0.5,0.2,"Z0/#gamma^{*}+jet");


  c1pt->cd(5);
  gPad->SetLogy(0);
  phPtAx->Draw();
  gPad->SetLogy();
  h5aPt->Draw("same");
  h5d->Draw("same");
  h5n->Draw("same");

  lt1->DrawLatex(0.5,0.2,"POWHEG");

  c1pt->SaveAs("pt.gif");

 


  */

//     ///// Comparision accpetance
/*
  h1n->Rebin(rebin);
    h1d->Rebin(rebin);
    h2n->Rebin(rebin);
    h2d->Rebin(rebin);
    h3n->Rebin(rebin);
    h3d->Rebin(rebin);
    h4n->Rebin(rebin);
    h4d->Rebin(rebin);
    h5n->Rebin(rebin);
    h5d->Rebin(rebin);

    h1n->Sumw2();
    h1d->Sumw2();
    h2n->Sumw2();
    h2d->Sumw2();
    h3n->Sumw2();
    h3d->Sumw2();
    h4n->Sumw2();
    h4d->Sumw2();
    h5n->Sumw2();
    h5d->Sumw2();


    TH1D *h1 = new TH1D("h1","h1",bin2,0,200);
    h1->Divide(h1n,h1d,1,1,"B");
    TH1D *h2 = new TH1D("h2","h2",bin2,0,200);
    h2->Divide(h2n,h2d,1,1,"B");
    TH1D *h3 = new TH1D("h3","h3",bin2,0,200);
    h3->Divide(h3n,h3d,1,1,"B");
    TH1D *h4 = new TH1D("h4","h4",bin2,0,200);
    h4->Divide(h4n,h4d,1,1,"B");
    TH1D *h5 = new TH1D("h5","h5",bin2,0,200);
    h5->Divide(h5n,h5d,1,1,"B");

    // calculate the integrated value:
    double err_1n,err_1d;
    double yield_1n = h1n->IntegralAndError(0,60,err_1n);
    double yield_1d = h1d->IntegralAndError(0,60,err_1d);
    double nh1 =  yield_1n/ yield_1d*100;
   

    double err_2n,err_2d;
    double yield_2n = h2n->IntegralAndError(0,60,err_2n);
    double yield_2d = h2d->IntegralAndError(0,60,err_2d);
    double nh2 =  yield_2n/ yield_2d*100;


    double err_3n,err_3d;
    double yield_3n = h3n->IntegralAndError(0,60,err_3n);
    double yield_3d = h3d->IntegralAndError(0,60,err_3d);
    double nh3 =  yield_3n/ yield_3d*100;


    double err_4n,err_4d;
    double yield_4n = h4n->IntegralAndError(0,60,err_4n);
    double yield_4d = h4d->IntegralAndError(0,60,err_4d);
    double nh4 =  yield_4n/ yield_4d*100;


    double err_5n,err_5d;
    double yield_5n = h5n->IntegralAndError(0,60,err_5n);
    double yield_5d = h5d->IntegralAndError(0,60,err_5d);
    double nh5 =  yield_5n/ yield_5d*100;

    // err_acceptance_pt[2][bin] = acceptance_pt[2][bin]* sqrt( pow(err_rec/rec,2)+ pow(err_gen/gen,2));      
    
  //  double nh1 = (h1->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh2 = (h2->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh3 = (h3->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh4 = (h4->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh5 = (h5->GetFunction("pol0")->GetParameter(0)) * 100;


    ////// 
    TH1 *phPtRatioAx = new TH1D("phPtRatioAx",";p_{T}[GeV/c^{2}];",1,0,60);
    phPtRatioAx->SetMaximum(1.2);
    phPtRatioAx->SetMinimum(0.5);

    TCanvas *c2 = new TCanvas("c2","c2");
    c2->cd();
    gStyle->SetOptFit(0);
    phPtRatioAx->Draw();

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kGreen+2);
    h4->SetLineColor(kOrange+2);

    h1->SetMarkerColor(kRed);
    h2->SetMarkerColor(kBlue);
    h3->SetMarkerColor(kGreen+2);
    h4->SetMarkerColor(kOrange+2);

    h1->SetMarkerStyle(20);
    h2->SetMarkerStyle(25);
    h3->SetMarkerStyle(20);
    h4->SetMarkerStyle(22);
    h5->SetMarkerStyle(30);

    h1->SetMarkerSize(2);
    h2->SetMarkerSize(2);
    h3->SetMarkerSize(2);
    h4->SetMarkerSize(2);
    h5->SetMarkerSize(2);

    h1->Draw("same");
    //  h2->Draw("same");
    h3->Draw("same");
    h4->Draw("same");
    h5->Draw("same");

   //  double nh1 = (h1->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh2 = (h2->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh3 = (h3->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh4 = (h4->GetFunction("pol0")->GetParameter(0)) * 100;
//     double nh5 = (h5->GetFunction("pol0")->GetParameter(0)) * 100;



    TLegend * leg2 = new TLegend(0.2,0.70,0.5,0.95);
    leg2->SetFillColor(0);
   
    leg2->AddEntry(h1,Form("Z: %2.1f %%",nh1),"pl");
    // leg2->AddEntry(h2,Form("Z/#gamma^{*}:%2.1f %%",nh2),"pl");
    leg2->AddEntry(h3,Form("Z+jet:%2.1f %%",nh3),"pl");
    leg2->AddEntry(h4,Form("Z/#gamma^{*}+jet:%2.1f %%",nh4),"pl");

    leg2->AddEntry(h5,Form("POWHEG: %2.1f %%",nh5),"pl");
    leg2->Draw("same");


    c2->SaveAs("acceptance_pythia.gif");
    
*/


    //ratio
  //   TCanvas *ca1 = new TCanvas("ca1","ca1");
//     ca1->cd();
//     TH1D *h12 = new TH1D("h12","h12",bin2,0,200);
//     TH1D *h13 = new TH1D("h13","h13",bin2,0,200);

//     h12->Divide(h1,h2,1,1,"b");
//     h13->Divide(h1,h3,1,1,"b");

//     h12->GetXaxis()->SetTitle("p_{T} [GeV/c]");
//     h12->GetYaxis()->SetTitle("#alpha ratio");
//     h12->GetXaxis()->SetRangeUser(0,59);
//     h12->GetYaxis()->SetRangeUser(0.6,1.2);
//     h13->GetXaxis()->SetTitle("p_{T} [GeV/c]");
//     h13->GetYaxis()->SetTitle("#alpha ratio");
//     h13->GetXaxis()->SetRangeUser(0,59);
//     h13->GetYaxis()->SetRangeUser(0.6,1.2);
//     h12->SetLineColor(kBlue);
//     h13->SetLineColor(kGreen);

//     h12->SetMarkerStyle(25);
//     h13->SetMarkerStyle(20);
//     h12->SetMarkerSize(2);
//     h13->SetMarkerSize(2);
//     h12->SetMarkerColor(kBlue);
//     h13->SetMarkerColor(kGreen);


//     h12->Fit("pol0");
//     h13->Fit("pol0");

//     h12->GetFunction("pol0")->SetLineColor(kBlue+1);
//     h13->GetFunction("pol0")->SetLineColor(kGreen+1);

//     h12->Draw("");
//     h13->Draw("same");

//     double nh12 = (h12->GetFunction("pol0")->GetParameter(0));
//     double nh13 = (h13->GetFunction("pol0")->GetParameter(0));

//     char legr[60];
//     TLegend * leg2 = new TLegend(0.52,0.20,0.92,0.33);
//     leg2->SetTextSize(0.02832861);
//     leg2->SetFillColor(0);
//     sprintf(legr,"PYTHIA / MC@NLO           : %1.2f",nh12);
//     leg2->AddEntry(h12,legr,"pl");
//     sprintf(legr,"CTEQ6L1 / MRST2004LO : %1.2f",nh13);
//     leg2->AddEntry(h13,legr,"pl");
//     leg2->Draw("same");

//     ca1->SaveAs("accept_ratio.gif");

}

