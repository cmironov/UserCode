#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"

#include "TFile.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"

#endif
void cmsRAA()
{
  
  bool bDo2PadZoo = false;
  
  //----------------------------------------- charged hadrons
  double ptBins_h[27]={0};
  double ptError_h[27]={0};
  
  double raa_h[27];
  double raaStat_h[27];
  double raaSyst_h[27];
 
  double ptSystXlow_h[27];
  double ptSystXhigh_h[27];

  // reading numbers
  ifstream inData_h; 
  TString inFile("raa_h05.dat");
  inData_h.open(inFile);
  if(inData_h.fail()) {
    cerr << "unable to open file raa_h05.dat for reading" << endl;
    exit(1);
  }
  Int_t j=0;
  Double_t xx_low, xx_high,raa,raa_syst, raa_stat;
  while(!inData_h.eof())
    { 
      inData_h >> xx_low >> xx_high >> raa >> raa_syst >> raa_stat;
      ptBins_h[j]   = xx_low+ (xx_high-xx_low)/2;
      raa_h[j]      = raa;
      raaStat_h[j]  = raa_syst;
      raaSyst_h[j]  = raa_stat;
      
      ptSystXlow_h[j] = (xx_high-xx_low)/2;
      ptSystXhigh_h[j]= (xx_high-xx_low)/2;
      // cout<<"pT"<<ptBins_h[j]<<"\t raa= "<<raa_h[j]<<"\t syst= "<<raaStat_h[j]<<"\t stat ="<<raaSyst_h[j]<<"\t x_low= "<<ptSystXlow_h[j]<<"\t high= "<<ptSystXhigh_h[j]<<endl;
      j++;
    }     
  inData_h.close();
  // done reading numebrs
  
  TGraphErrors *pgRaa_h          = new TGraphErrors(27, ptBins_h, raa_h, ptError_h, raaStat_h);
  TGraphAsymmErrors *pgRaaSyst_h = new TGraphAsymmErrors(27, ptBins_h, raa_h, ptSystXlow_h,ptSystXhigh_h,raaSyst_h,raaSyst_h);
  pgRaa_h->SetName("pgRaa_h");
  pgRaa_h->SetMarkerStyle(24);
  pgRaa_h->SetMarkerSize(1.);
  
  //systm error
  pgRaaSyst_h->SetName("pgRaaSyst_h");
  pgRaaSyst_h->SetFillColor(TColor::GetColor("#33ccff"));

  //----------------------------------------- jet RAA
  double ptBins_jet[]   = {105,115,125,135,145,155,165,175,190,220,270}; 
  double ptError_jet[]  = {0.0, 0.0, 0.0, 0.0, 0.0,0,0,0,0,0,0};
  
  double raa_jet[]      = {0.47,   0.47,   0.46,  0.44,   0.43,  0.47,  0.52,  0.51,  0.44,  0.42,  0.58}; 
  double raaStat_jet[]  = {0.0076, 0.0098, 0.013, 0.016,  0.021, 0.028, 0.038, 0.045, 0.035, 0.039, 0.098};
  double raaSyst_jet[]  = {0.07, 0.07, 0.08, 0.09, 0.08, 0.08, 0.08, 0.08, 0.07, 0.07, 0.10 };
 
  double ptSystXlow_jet[]      = {5, 5, 5, 5, 5,5,5,5,10,20,30};
  double ptSystXhigh_jet[]     = {5, 5, 5, 5, 5,5,5,5,10,20,30};
  
  TGraphErrors *pgRaa_jet          = new TGraphErrors(11, ptBins_jet, raa_jet, ptError_jet, raaStat_jet);
  TGraphAsymmErrors *pgRaaSyst_jet = new TGraphAsymmErrors(11, ptBins_jet, raa_jet, ptSystXlow_jet,ptSystXhigh_jet,raaSyst_jet,raaSyst_jet);
  pgRaa_jet->SetName("pgRaa_jet");
  pgRaa_jet->SetMarkerStyle(20);
  pgRaa_jet->SetMarkerSize(1.);
  
  //systm error
  pgRaaSyst_jet->SetName("pgRaaSyst_jet");
  pgRaaSyst_jet->SetFillColor(TColor::GetColor("#00FFFF"));

 //----------------------------------------- b-jet RAA
  double ptBins_bjet[]   = {108.3}; 
  double ptError_bjet[]  = {0.0};
  
  double raa_bjet[]      = {0.478}; 
  double raaStat_bjet[]  = {0.188};
  double raaSyst_bjet[]  = {0.413};
 
  double ptSystXlow_bjet[]      = {8.3};
  double ptSystXhigh_bjet[]     = {11.7};
  
  TGraphErrors *pgRaa_bjet          = new TGraphErrors(1, ptBins_bjet, raa_bjet, ptError_bjet, raaStat_bjet);
  TGraphAsymmErrors *pgRaaSyst_bjet = new TGraphAsymmErrors(1, ptBins_bjet, raa_bjet, ptSystXlow_bjet,ptSystXhigh_bjet,raaSyst_bjet,raaSyst_bjet);
  pgRaa_bjet->SetName("pgRaa_bjet");
  pgRaa_bjet->SetMarkerStyle(21);
  pgRaa_bjet->SetMarkerSize(1.);
  
  //systm error
  pgRaaSyst_bjet->SetName("pgRaaSyst_bjet");
  pgRaaSyst_bjet->SetFillColor(TColor::GetColor("#FFBF00"));


  //----------------------------------------- photon
  double ptBins_photon[]   = {22.26, 27.30, 34.35, 44.45, 61.72}; 
  double ptError_photon[]  = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  double raa_photon[]      = {1.03, 0.84, 1.37, 0.98, 0.99}; 
  double raaStat_photon[]  = {0.12 ,0.14, 0.22, 0.24, 0.31 };
  double raaSyst_photon[]  = {0.29, 0.27,  0.25, 0.22, 0.21};
 
  double ptSystXlow_photon[]      = {2.5, 2.5, 5, 5, 15};
  double ptSystXhigh_photon[]     = {2.5, 2.5, 5, 5, 15};
  
  TGraphErrors *pgRaa_photon          = new TGraphErrors(5, ptBins_photon, raa_photon, ptError_photon, raaStat_photon);
  TGraphAsymmErrors *pgRaaSyst_photon = new TGraphAsymmErrors(5, ptBins_photon, raa_photon, ptSystXlow_photon,ptSystXhigh_photon,raaSyst_photon,raaSyst_photon);
  pgRaa_photon->SetName("pgRaa_photon");
  pgRaa_photon->SetMarkerStyle(20);
  pgRaa_photon->SetMarkerSize(1.);
  //systm error
  pgRaaSyst_photon->SetName("pgRaaSyst_photon");
  pgRaaSyst_photon->SetFillColor(TColor::GetColor("#ffff00"));
  
  //----------------------------------------- Z 
   double ptBins_z[]   = {91.19}; 
   double ptError_z[]  = {0.0};
   
   double raa_z[]      = {0.95}; 
   double raaStat_z[]  = {0.034};
   double raaSyst_z[]  = {0.13};
   double ptSystXlow_z[]      = {2.5};
   double ptSystXhigh_z[]     = {2.5};

   TGraphErrors *pgRaa_z          = new TGraphErrors(1, ptBins_z, raa_z, ptError_z, raaStat_z);
   TGraphAsymmErrors *pgRaaSyst_z = new TGraphAsymmErrors(1, ptBins_z, raa_z, ptSystXlow_z,ptSystXhigh_z,raaSyst_z,raaSyst_z);
   pgRaa_z->SetName("pgRaa_z");   pgRaa_z->SetMarkerStyle(21);
   pgRaa_z->SetMarkerSize(1.);
  
   //systm error
   pgRaaSyst_z->SetName("pgRaaSyst_z");
   pgRaaSyst_z->SetFillColor(TColor::GetColor("#ff8888"));

 
   // ----------------------------------------- W raa
   double ptBins_w[]   = {80.38}; 
   double ptError_w[]  = {0.0};
   
   double raa_w[]      = {1.04}; 
   double raaStat_w[]  = {0.07};
   double raaSyst_w[]  = {0.12};
   double ptSystXlow_w[]      = {2.5};
   double ptSystXhigh_w[]     = {2.5};
 
   TGraphErrors *pgRaa_w          = new TGraphErrors(1, ptBins_w, raa_w, ptError_w, raaStat_w);
   TGraphAsymmErrors *pgRaaSyst_w = new TGraphAsymmErrors(1, ptBins_w, raa_w, ptSystXlow_w,ptSystXhigh_w,raaSyst_w,raaSyst_w);
   pgRaa_w->SetName("pgRaa_w");
   pgRaa_w->SetMarkerStyle(29);
   pgRaa_w->SetMarkerSize(2.);
  
   //systm error
   pgRaaSyst_w->SetName("pgRaaSyst_w");
   pgRaaSyst_w->SetFillColor(TColor::GetColor("#ff88ff"));


    // ----------------------------------------- 2012 non-prompt Jpsi
   double ptBins_npjpsi[]   = {7.31,8.97,11.32,16.52}; 
   double ptError_npjpsi[]  = {0.0,0.0,0.0,0.0};
   
   double raa_npjpsi[]      = {0.52,0.43,0.43,0.34};  
   double raaStat_npjpsi[]  = {0.12,0.08,0.09,0.07};
   double raaSyst_npjpsi[]  = {0.06,0.05,0.05,0.04};

   double ptSystXlow_npjpsi[]      = {0.81, 0.97, 1.32, 3.52};
   double ptSystXhigh_npjpsi[]     = {0.69, 1.03, 1.68, 13.48};
 
   TGraphErrors *pgRaa_npjpsi          = new TGraphErrors(4, ptBins_npjpsi, raa_npjpsi, ptError_npjpsi, raaStat_npjpsi);
   TGraphAsymmErrors *pgRaaSyst_npjpsi = new TGraphAsymmErrors(4, ptBins_npjpsi, raa_npjpsi, ptSystXlow_npjpsi,ptSystXhigh_npjpsi,raaSyst_npjpsi,raaSyst_npjpsi);
   pgRaa_npjpsi->SetName("pgRaa_npjpsi");
   pgRaa_npjpsi->SetMarkerStyle(25);
   pgRaa_npjpsi->SetMarkerSize(1.);
   pgRaa_npjpsi->SetMarkerColor(1);
  
   //systm error
   pgRaaSyst_npjpsi->SetName("pgRaaSyst_npjpsi");
   pgRaaSyst_npjpsi->SetFillColor(TColor::GetColor("#ee7711"));
  
  
   //---------------------------------------------------------
   // lumi uncert.
   TBox *box = new TBox(0.9,0.9568966,2,1.043103);
   box->SetFillColor(kGray+1);
   box->SetFillStyle(1001);
   

   TBox *box_lin = new TBox(0.9,0.9568966,10,1.043103);
   box_lin->SetFillColor(kGray+1);
   box_lin->SetFillStyle(1001);
  
   //------------------------ single pannel:
  
  TCanvas *pc = new TCanvas("pc","pc");
  TH1 *phAxis = new TH1D("phAxis",";p_{T} [GeV];R_{AA}",1,0.9,300);
  phAxis->GetYaxis()->CenterTitle(true);
  phAxis->SetMinimum(0.);
  phAxis->SetMaximum(2.5);
  gPad->SetLogx();
  phAxis->Draw("");
  
// drawing order:
// A)  stat boxes:
   pgRaaSyst_photon->Draw("2");
   pgRaaSyst_npjpsi->Draw("2");
   pgRaaSyst_h->Draw("2");
   pgRaaSyst_z->Draw("2");
   pgRaaSyst_w->Draw("2");
   pgRaaSyst_jet->Draw("2");
   pgRaaSyst_bjet->Draw("2");
   //nominal+stat
   // photon
   pgRaa_photon->Draw("P z");
   pgRaa_z->Draw("P z");
   pgRaa_w->Draw("P z");
   pgRaa_h->Draw("P z");
   pgRaa_npjpsi->Draw("P z");
   pgRaa_jet->Draw("P z");
   pgRaa_bjet->Draw("P z");

   TLine *line = new TLine(0.9,1,300,1);
   line->SetLineStyle(7);
   line->Draw();
   box->Draw();
// // legends and writings:

   TLegend *leg = new TLegend(0.07,0.87,0.46,0.95,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.028);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
  //  TLegendEntry *entry=leg->AddEntry("hEtSIEIECorrected","CMS PRELIMINARY","");
//    entry->SetLineWidth(1);
//    entry=leg->AddEntry("hEtSIEIECorrected","PbPb #sqrt{s_{NN}} = 2.76 TeV","");
//    entry=leg->AddEntry("hEtSIEIECorrected","","");
//    entry=leg->AddEntry("hEtSIEIECorrected","#int L dt = 7-150 #mub^{-1}","");

   TLegendEntry *entry=leg->AddEntry("hEtSIEIECorrected","CMS PRELIMINARY  PbPb #sqrt{s_{NN}} = 2.76 TeV  #int L dt = 7-150 #mub^{-1}","");
   entry->SetLineWidth(1);
  
  

   // TAA
  // TAA only
   TLegend *leg_taa = new TLegend(0.15,0.5,0.5,0.55,NULL,"brNDC");
   leg_taa->SetBorderSize(0);
   leg_taa->SetTextFont(62);
   leg_taa->SetTextSize(0.025);
   leg_taa->SetLineColor(1);
   leg_taa->SetLineStyle(1);
   leg_taa->SetLineWidth(1);
   leg_taa->SetFillColor(19);
   leg_taa->SetFillStyle(0);
   TLegendEntry *entry_taa= leg_taa->AddEntry("general","T_{AA} uncertainty","f");
   entry_taa->SetFillColor(kGray+1);
   entry_taa->SetFillStyle(1001);
   entry_taa->SetLineColor(0);
   entry_taa->SetLineWidth(1);
   entry_taa->SetMarkerStyle(21);
   entry_taa->SetMarkerSize(1);

   // EWQ legend
   TLegend *leg_ewq = new TLegend(0.35,0.75,0.8,0.85,NULL,"brNDC");
   leg_ewq->SetBorderSize(0);
   leg_ewq->SetTextFont(62);
   leg_ewq->SetTextSize(0.025);
   leg_ewq->SetLineColor(1);
   leg_ewq->SetLineStyle(1);
   leg_ewq->SetLineWidth(1);
   leg_ewq->SetFillColor(19);
   leg_ewq->SetFillStyle(0);
   TLegendEntry *entry_ewq=leg_ewq->AddEntry("raaz","Z  (0-100%) p_{T}^{#mu} > 20 GeV/c^{^{ }}    |y|<2","lpf");
   entry_ewq->SetFillColor(TColor::GetColor("#ff8888"));
   entry_ewq->SetFillStyle(1001);
   entry_ewq->SetLineColor(1);
   entry_ewq->SetLineWidth(1);
   entry_ewq->SetMarkerStyle(21);
   entry_ewq->SetMarkerSize(1.);
  
   entry_ewq=leg_ewq->AddEntry("raaw","W  (0-100%)  p_{T}^{#mu} > 25 GeV/c^{^{ }} |#eta^{#mu}|<2.1","lpf");
   entry_ewq->SetFillColor(TColor::GetColor("#ff88ff"));
   entry_ewq->SetFillStyle(1001);
   entry_ewq->SetLineColor(1);
   entry_ewq->SetLineWidth(1);
   entry_ewq->SetMarkerStyle(29);
   entry_ewq->SetMarkerSize(1.5);
  
   entry_ewq=leg_ewq->AddEntry("raaphoton","Isolated photon  (0-10%)      |#eta|<1.44","lpf");
   entry_ewq->SetFillColor(TColor::GetColor("#ffff00"));
   entry_ewq->SetFillStyle(1001);
   entry_ewq->SetLineColor(1);
   entry_ewq->SetLineWidth(1);
   entry_ewq->SetMarkerStyle(20);
   entry_ewq->SetMarkerSize(1);
   

   // b and h
   TLegend *leg_q = new TLegend(0.15,0.65,0.47,0.74,NULL,"brNDC");
   leg_q->SetBorderSize(0);
   leg_q->SetTextFont(62);
   leg_q->SetTextSize(0.025);
   leg_q->SetLineColor(1);
   leg_q->SetLineStyle(1);
   leg_q->SetLineWidth(1);
   leg_q->SetFillColor(19);
   leg_q->SetFillStyle(0);


   TLegendEntry *entry_q=leg_q->AddEntry("raah","Charged particles^{^{ }} (0-5%)   |#eta|<1","lpf");
   entry_q->SetFillColor(TColor::GetColor("#33ccff"));
   entry_q->SetFillStyle(1001);
   entry_q->SetLineColor(1);
   entry_q->SetLineStyle(1);
   entry_q->SetLineWidth(1);
   entry_q->SetMarkerStyle(24);
   entry_q->SetMarkerSize(1);

   entry_q=leg_q->AddEntry("raab","b-quarks (0-100%)^{ }             |#eta|<2.4","lpf");
   entry_q->SetFillColor(TColor::GetColor("#ee7711"));
   entry_q->SetFillStyle(1001);
   entry_q->SetLineColor(1);
   entry_q->SetLineStyle(1);
   entry_q->SetLineWidth(1);
   entry_q->SetMarkerStyle(25);
   entry_q->SetMarkerSize(1);
   entry_q=leg_q->AddEntry("raabjpsi", "(via secondary J/#psi) ","");
   
     // jet legend
   TLegend *leg_jet = new TLegend(0.57,0.69,0.95,0.74,NULL,"brNDC");
   leg_jet->SetBorderSize(0);
   leg_jet->SetTextFont(62);
   leg_jet->SetTextSize(0.025);
   leg_jet->SetLineColor(1);
   leg_jet->SetLineStyle(1);
   leg_jet->SetLineWidth(1);
   leg_jet->SetFillColor(19);
   leg_jet->SetFillStyle(0);

   TLegendEntry *entry_jet=leg_jet->AddEntry("raaq","q/g-jet^{^{ }} (0-5%)     |#eta|<2","lpf");
   entry_jet->SetFillColor(TColor::GetColor("#00FFFF"));
   entry_jet->SetFillStyle(1001);
   entry_jet->SetLineColor(1);
   entry_jet->SetLineStyle(1);
   entry_jet->SetLineWidth(1);
   entry_jet->SetMarkerStyle(20);
   entry_jet->SetMarkerSize(1);

   entry_jet=leg_jet->AddEntry("raaheavy","b-jet (0-100%)^{ }     |#eta|<2","lpf");
   entry_jet->SetFillColor(TColor::GetColor("#FFBF00"));
   entry_jet->SetFillStyle(1001);
   entry_jet->SetLineColor(1);
   entry_jet->SetLineStyle(1);
   entry_jet->SetLineWidth(1);
   entry_jet->SetMarkerStyle(21);
   entry_jet->SetMarkerSize(1);
  
   leg->Draw();
   leg_taa->Draw();
   leg_ewq->Draw();
   leg_q->Draw();
   leg_jet->Draw();
   gPad->RedrawAxis();
   
   // save the work
   pc->SaveAs("raaZoo_cms_log.pdf");
   pc->SaveAs("raaZoo_cms_log.png");




   //--------------------------------------------------------------------------------
   //-------------------------- 2TPad
   TH1 *phAxis_single = new TH1D("phAxis_single",";p_{T} [GeV/c];R_{AA}",1,0.9,100);
   phAxis_single->GetYaxis()->CenterTitle(true);
   phAxis_single->SetMinimum(0.);
   phAxis_single->SetMaximum(2.);
   TH1 *phAxis_jet = new TH1D("phAxis_jet",";p_{T} [GeV/c];",1,101,300);
   phAxis_jet->GetYaxis()->CenterTitle(true);
   phAxis_jet->SetMinimum(0.);
   phAxis_jet->SetMaximum(2.);

   // redefine legends:
   // TAA only
   TLegend *leg2_taa = new TLegend(0.15,0.85,0.5,0.9,NULL,"brNDC");
   leg2_taa->SetBorderSize(0);
   leg2_taa->SetTextFont(62);
   leg2_taa->SetTextSize(0.025);
   leg2_taa->SetLineColor(1);
   leg2_taa->SetLineStyle(1);
   leg2_taa->SetLineWidth(1);
   leg2_taa->SetFillColor(19);
   leg2_taa->SetFillStyle(0);
   TLegendEntry *entry2_taa= leg2_taa->AddEntry("general","T_{AA} uncertainty","f");
   entry2_taa->SetFillColor(kGray+1);
   entry2_taa->SetFillStyle(1001);
   entry2_taa->SetLineColor(0);
   entry2_taa->SetLineWidth(1);
   entry2_taa->SetMarkerStyle(21);
   entry2_taa->SetMarkerSize(1);

   TLegend *leg2 = new TLegend(0.03,0.7,0.46,0.93,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextFont(62);
   leg2->SetTextSize(0.04);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(19);
   leg2->SetFillStyle(0);
   TLegendEntry *entry2=leg2->AddEntry("general","CMS PRELIMINARY","");
   entry2->SetLineWidth(1);
   entry2=leg2->AddEntry("general","PbPb #sqrt{s_{NN}} = 2.76 TeV","");
   entry2=leg2->AddEntry("general","","");
   entry2=leg2->AddEntry("general","#int L dt = 7-150 #mub^{-1}","");
  
      // jet leg2end
   TLegend *leg2_jet = new TLegend(0.45,0.55,0.9,0.65,NULL,"brNDC");
   leg2_jet->SetBorderSize(0);
   leg2_jet->SetTextFont(62);
   leg2_jet->SetTextSize(0.03);
   leg2_jet->SetLineColor(1);
   leg2_jet->SetLineStyle(1);
   leg2_jet->SetLineWidth(1);
   leg2_jet->SetFillColor(19);
   leg2_jet->SetFillStyle(0);

   TLegendEntry *entry2_jet=leg2_jet->AddEntry("raaq","q/g-jet^{^{ }} (0-5%)     |#eta|<2","lpf");
   entry2_jet->SetFillColor(TColor::GetColor("#00FFFF"));
   entry2_jet->SetFillStyle(1001);
   entry2_jet->SetLineColor(1);
   entry2_jet->SetLineStyle(1);
   entry2_jet->SetLineWidth(1);
   entry2_jet->SetMarkerStyle(20);
   entry2_jet->SetMarkerSize(1);

   entry2_jet=leg2_jet->AddEntry("raaheavy","b-jet (0-100%)^{ }     |#eta|<2","lpf");
   entry2_jet->SetFillColor(TColor::GetColor("#FFBF00"));
   entry2_jet->SetFillStyle(1001);
   entry2_jet->SetLineColor(1);
   entry2_jet->SetLineStyle(1);
   entry2_jet->SetLineWidth(1);
   entry2_jet->SetMarkerStyle(21);
   entry2_jet->SetMarkerSize(1);


   // EWQ legend
   TLegend *leg2_ewq = new TLegend(0.47,0.8,0.9,0.95,NULL,"brNDC");
   leg2_ewq->SetBorderSize(0);
   leg2_ewq->SetTextFont(62);
   leg2_ewq->SetTextSize(0.027);
   leg2_ewq->SetLineColor(1);
   leg2_ewq->SetLineStyle(1);
   leg2_ewq->SetLineWidth(1);
   leg2_ewq->SetFillColor(19);
   leg2_ewq->SetFillStyle(0);
   TLegendEntry *entry2_ewq=leg2_ewq->AddEntry("raaz","Z  (0-100%) p_{T}^{#mu} > 20 GeV/c^{^{ }}    |y|<2","lpf");
   entry2_ewq->SetFillColor(TColor::GetColor("#ff8888"));
   entry2_ewq->SetFillStyle(1001);
   entry2_ewq->SetLineColor(1);
   entry2_ewq->SetLineWidth(1);
   entry2_ewq->SetMarkerStyle(21);
   entry2_ewq->SetMarkerSize(1.);
  
   entry2_ewq=leg2_ewq->AddEntry("raaw","W  (0-100%)  p_{T}^{#mu} > 25 GeV/c^{^{ }}  |#eta^{#mu}|<2.1","lpf");
   entry2_ewq->SetFillColor(TColor::GetColor("#ff88ff"));
   entry2_ewq->SetFillStyle(1001);
   entry2_ewq->SetLineColor(1);
   entry2_ewq->SetLineWidth(1);
   entry2_ewq->SetMarkerStyle(29);
   entry2_ewq->SetMarkerSize(1.5);
  
   entry2_ewq=leg2_ewq->AddEntry("raaphoton","Isolated photon  (0-10%)      |#eta|<1.44","lpf");
   entry2_ewq->SetFillColor(TColor::GetColor("#ffff00"));
   entry2_ewq->SetFillStyle(1001);
   entry2_ewq->SetLineColor(1);
   entry2_ewq->SetLineWidth(1);
   entry2_ewq->SetMarkerStyle(20);
   entry2_ewq->SetMarkerSize(1);
 
   // B and h
   TLegend *leg2_q = new TLegend(0.47,0.15,0.89,0.3,NULL,"brNDC");
   leg2_q->SetBorderSize(0);
   leg2_q->SetTextFont(62);
   leg2_q->SetTextSize(0.027);
   leg2_q->SetLineColor(1);
   leg2_q->SetLineStyle(1);
   leg2_q->SetLineWidth(1);
   leg2_q->SetFillColor(19);
   leg2_q->SetFillStyle(0);

   TLegendEntry *entry2_q=leg2_q->AddEntry("raah","Charged particles^{^{ }} (0-5%)     |#eta|<1","lpf");
   entry2_q->SetFillColor(TColor::GetColor("#33ccff"));
   entry2_q->SetFillStyle(1001);
   entry2_q->SetLineColor(1);
   entry2_q->SetLineStyle(1);
   entry2_q->SetLineWidth(1);
   entry2_q->SetMarkerStyle(24);
   entry2_q->SetMarkerSize(1);

   entry2_q=leg2_q->AddEntry("raab","b-quarks (0-100%)^{ }                |y|<2.4","lpf");
   entry2_q->SetFillColor(TColor::GetColor("#ee7711"));
   entry2_q->SetFillStyle(1001);
   entry2_q->SetLineColor(1);
   entry2_q->SetLineStyle(1);
   entry2_q->SetLineWidth(1);
   entry2_q->SetMarkerStyle(25);
   entry2_q->SetMarkerSize(1);
   entry2_q=leg2_q->AddEntry("raabjpsi", "(via secondary J/#psi) ","");

   if(bDo2PadZoo)
     {
       TCanvas *pc2 = new TCanvas("pc2","pc2",1200,600);
       TPad *p_0 = new TPad("p_0","p_0",0,0,0.50,1);
       p_0->Draw();
       p_0->cd();
       p_0->SetRightMargin(0);
       
       phAxis_single->Draw("");
       box_lin->Draw();
       
       pgRaaSyst_photon->Draw("2");
       pgRaaSyst_npjpsi->Draw("2");
       pgRaaSyst_h->Draw("2");
       pgRaaSyst_z->Draw("2");
       pgRaaSyst_w->Draw("2");
       
       pgRaa_photon->Draw("P z");
       pgRaa_z->Draw("P z");
       pgRaa_w->Draw("P z");
       pgRaa_h->Draw("P z");
       pgRaa_npjpsi->Draw("P z");
       
       TLine *line_single = new TLine(0.9,1,100,1);
       line_single->SetLineStyle(7);
       line_single->Draw();
       
       // legends 
       leg2_taa->Draw();
       leg2_q->Draw();
       leg2_ewq->Draw();
       
       pc2->cd();
       
       // jet pad
       TPad *p_1 = new TPad("p_1","p_1",0.50,0,1,1);
       p_1->Draw();
       p_1->cd();
       p_1->SetLeftMargin(0.);
       p_1->SetTickx(1);
       p_1->SetTicky(1);
       
       phAxis_jet->Draw("");
       
       pgRaaSyst_bjet->Draw("2");
       pgRaaSyst_jet->Draw("2");
       
       pgRaa_bjet->Draw("P z");
       pgRaa_jet->Draw("P z");
       
       
       TLine *line_jet = new TLine(100,1,300,1);
       line_jet->SetLineStyle(7);
       line_jet->Draw();
       
       leg2->Draw();
       leg2_jet->Draw();
       
       
       p_0->cd();
       gPad->RedrawAxis();
       
       pc2->SaveAs("raaZoo_cms_2pads.pdf");
       pc2->SaveAs("raaZoo_cms_2pads.png");

     }
}
