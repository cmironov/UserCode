#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"

#include "TFile.h"
#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"

#endif
// all minbias numbers
void thermo()
{
  bool bSavePlots = true;
  bool bDoRadius = false;
  bool bDoRatio = false;
  bool bDoNpart  = false;

  bool bDoThermo = true;
  
  // Jpsi, psi(2S), Y(1S), Y(2S), Y(3S)(upper limit)
  double xBindingEn[] = {0.64, 0.05, 1.10, 0.54, 0.20};
  double xRadius[]    = {0.25, 0.45, 0.28, 0.56, 0.78};
  double x_err[]      = {0.0, 0.0, 0.0, 0.0, 0.0};
  double x_sys[]      = {0.025, 0.025, 0.025, 0.025, 0.025};

  double RAA[]     = {0.34, 0.11, 0.56, 0.12, 0.03};
  double RAA_err[] = {0.02, 0.03, 0.08, 0.04, 0.04};
  double RAA_sys[] = {0.04, 0.18, 0.07, 0.02, 0.01};

  double xAx[5];
  double arrowEnd = 0.03;
  if(bDoRadius)
    {
      for(int i=0; i<5; i++)
	xAx[i] = xRadius[i];
    }
  else
    {
      if(bDoRatio)
	{
	  for(int i=0; i<5; i++)
	    xAx[i] = xBindingEn[i]/xRadius[i];
	  arrowEnd = 0.05;
	}
      else
	{
	  for(int i=0; i<5; i++)
	    xAx[i] = xBindingEn[i];
	}
    }
  if(bDoThermo)
    {
      for(int i=0; i<5; i++)
	xAx[i] = 0.5;
      arrowEnd = 0.1;
    }

  // 3S upper limit
  double raaMB_ups3s_limit[1] = {0.1};
  TLine *gUps3_limit = new TLine(xAx[4]-arrowEnd,raaMB_ups3s_limit[0],xAx[4]+arrowEnd,raaMB_ups3s_limit[0]);
  gUps3_limit->SetLineWidth(3);
  gUps3_limit->SetLineColor(kViolet-2);
  TArrow *gUps3_limitArrow = new TArrow( xAx[4],raaMB_ups3s_limit[0], xAx[4],0,0.015,"|>");
  gUps3_limitArrow->SetLineWidth(3);
  gUps3_limitArrow->SetLineColor(kViolet-2);
  
  TGraphErrors *gJpsi = new TGraphErrors(1, &xAx[0], &RAA[0], &x_err[0], &RAA_err[0]);
  TGraphErrors *gPsiP = new TGraphErrors(1, &xAx[1], &RAA[1], &x_err[1], &RAA_err[1]);
  TGraphErrors *gUps1 = new TGraphErrors(1, &xAx[2], &RAA[2], &x_err[2], &RAA_err[2]);
  TGraphErrors *gUps2 = new TGraphErrors(1, &xAx[3], &RAA[3], &x_err[3], &RAA_err[3]);
  TGraphErrors *gUps3 = new TGraphErrors(1, &xAx[4], &RAA[4], &x_err[4], &RAA_err[4]);

  gJpsi->SetName("gJpsi");
  gPsiP->SetName("gPsiP");
  gUps1->SetName("gUps1");
  gUps2->SetName("gUps2");
  gUps3->SetName("gUps3");

  TGraphErrors *gJpsi_sys = new TGraphErrors(1, &xAx[0], &RAA[0], &x_sys[0], &RAA_sys[0]);
  TGraphErrors *gPsiP_sys = new TGraphErrors(1, &xAx[1], &RAA[1], &x_sys[1], &RAA_sys[1]);
  TGraphErrors *gUps1_sys = new TGraphErrors(1, &xAx[2], &RAA[2], &x_sys[2], &RAA_sys[2]);
  TGraphErrors *gUps2_sys = new TGraphErrors(1, &xAx[3], &RAA[3], &x_sys[3], &RAA_sys[3]);
  TGraphErrors *gUps3_sys = new TGraphErrors(1, &xAx[4], &RAA[4], &x_sys[4], &RAA_sys[4]);

  gJpsi_sys->SetName("gJpsi_sys");
  gPsiP_sys->SetName("gPsiP_sys");
  gUps1_sys->SetName("gUps1_sys");
  gUps2_sys->SetName("gUps2_sys");
  gUps3_sys->SetName("gUps3_sys");


  TGraphErrors *gJpsi_p = new TGraphErrors(1, &xAx[0], &RAA[0], &x_err[0], &RAA_err[0]);
  TGraphErrors *gPsiP_p = new TGraphErrors(1, &xAx[1], &RAA[1], &x_err[1], &RAA_err[1]);
  TGraphErrors *gUps1_p = new TGraphErrors(1, &xAx[2], &RAA[2], &x_err[2], &RAA_err[2]);
  TGraphErrors *gUps2_p = new TGraphErrors(1, &xAx[3], &RAA[3], &x_err[3], &RAA_err[3]);
  TGraphErrors *gUps3_p = new TGraphErrors(1, &xAx[4], &RAA[4], &x_err[4], &RAA_err[4]);

  gJpsi_p->SetName("gJpsi_p");
  gPsiP_p->SetName("gPsiP_p");
  gUps1_p->SetName("gUps1_p");
  gUps2_p->SetName("gUps2_p");
  gUps3_p->SetName("gUps3_p");


 

  gJpsi_sys->SetFillColor(kRed-9);
  gPsiP_sys->SetFillColor(kOrange-9);
  gUps1_sys->SetFillColor(kGreen-9);
  gUps2_sys->SetFillColor(kBlue-9);
  gUps3_sys->SetFillColor(kGray);


  gJpsi->SetMarkerStyle(21);
  gJpsi_p->SetMarkerStyle(25);

  gPsiP->SetMarkerStyle(20);
  gPsiP_p->SetMarkerStyle(24);


  gUps1->SetMarkerStyle(29);
  gUps1_p->SetMarkerStyle(30);

  gUps2->SetMarkerStyle(34);
  gUps2_p->SetMarkerStyle(28);

  gUps3->SetMarkerStyle(23);
  gUps3->SetMarkerColor(kViolet);

  gJpsi->SetMarkerColor(kRed+2);
  gPsiP->SetMarkerColor(kOrange+2);
  gUps1->SetMarkerColor(kGreen+2);
  gUps2->SetMarkerColor(kBlue+2);

  gJpsi->SetMarkerSize(1.2);
  gPsiP->SetMarkerSize(1.2);
  gUps1->SetMarkerSize(1.7);
  gUps2->SetMarkerSize(1.7);
  gUps3->SetMarkerSize(1.2);

  gJpsi_p->SetMarkerSize(1.2);
  gPsiP_p->SetMarkerSize(1.2);
  gUps1_p->SetMarkerSize(1.7);
  gUps2_p->SetMarkerSize(1.7);

   // drawing
  TF1 *f1 = new TF1("f1","1",0.1,0.9);
  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("r_{0} [fm]");
  f1->GetYaxis()->SetTitle("R_{AA}");
  f1->GetYaxis()->SetRangeUser(0.0,1.4);
  f1->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f2 = new TF1("f2","1",0,1.2);
  f2->SetLineWidth(1);
  f2->GetXaxis()->SetTitle("#Delta E [GeV]");
  f2->GetYaxis()->SetTitle("R_{AA}");
  f2->GetYaxis()->SetRangeUser(0.0,1.4);
  f2->GetXaxis()->CenterTitle(kTRUE);

  TF1 *f3 = new TF1("f3","1",0,4.5);
  f3->SetLineWidth(1);
  f3->GetXaxis()->SetTitle("#Delta E / r_{0} [GeV/fm]");
  f3->GetYaxis()->SetTitle("R_{AA}");
  f3->GetYaxis()->SetRangeUser(0.0,1.4);
  f3->GetXaxis()->CenterTitle(kTRUE);

  // vs binding energy
  TCanvas *c1 = new TCanvas("c1","c1");
  if(bDoRadius)
    f1->Draw();
  else
    {
      if(bDoRatio)
	f3->Draw();
      else
	f2->Draw();
    }

  gJpsi_sys->Draw("2");
  gPsiP_sys->Draw("2");
  gUps1_sys->Draw("2");
  gUps2_sys->Draw("2");
  // gUps3_sys->Draw("2");

  gJpsi->Draw("P");
  gPsiP->Draw("P");
  gUps1->Draw("P");
  gUps2->Draw("P");
  //  gUps3->Draw("P");

  gJpsi_p->Draw("P");
  gPsiP_p->Draw("P");
  gUps1_p->Draw("P");
  gUps2_p->Draw("P");
  // gUps3_p->Draw("P");

  gUps3_limit->Draw("D");
  gUps3_limitArrow->Draw();

  TLegend *leg1;
  leg1 = new TLegend(0.15,0.45,0.7,0.72);
  leg1->AddEntry(gPsiP,"Inclusive #psi(2S) (6.5 < p_{T} < 30 GeV/c, |y| < 1.6)","P");
  leg1->AddEntry(gUps3,"#varUpsilon(3S) (|y| < 2.4), 95% upper limit","P");
  leg1->AddEntry(gUps2,"#varUpsilon(2S) (|y| < 2.4)","P");
  leg1->AddEntry(gJpsi,"prompt J/#psi (6.5 < p_{T} < 30 GeV/c, |y| < 2.4)","P");
  leg1->AddEntry(gUps1,"#varUpsilon(1S) (|y| < 2.4)","P");

  if(bDoRadius)
    {
      leg1 = new TLegend(0.15,0.53,0.7,0.72);
      leg1->AddEntry(gJpsi,"prompt J/#psi (6.5 < p_{T} < 30 GeV/c, |y| < 2.4)","P");
      leg1->AddEntry(gUps1,"#varUpsilon(1S) (|y| < 2.4)","P");
      leg1->AddEntry(gPsiP,"Inclusive #psi(2S) (6.5 < p_{T} < 30 GeV/c, |y| < 1.6)","P");
      leg1->AddEntry(gUps2,"#varUpsilon(2S) (|y| < 2.4)","P");
      leg1->AddEntry(gUps3,"#varUpsilon(3S) (|y| < 2.4), 95% upper limit","P"); 
    }
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.15);
  leg1->SetTextSize(0.035);
  leg1->Draw();

  TLatex *pre = new TLatex(0.1,1.25,"CMS Preliminary");
  pre->SetTextFont(42);
  pre->SetTextSize(0.05);

 
  if(bDoRadius)
    {
      pre->DrawLatex(0.15,1.25,"CMS Preliminary");
      pre->DrawLatex(0.15,1.15,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
      pre->DrawLatex(0.7,1.25,"0-100%");
      
    }
  else
    {
      if(bDoRatio)
	{
	  pre->DrawLatex(0.4,1.25,"CMS Preliminary");
	  pre->DrawLatex(0.4,1.15,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
	  pre->DrawLatex(3.5,1.25,"0-100%");
	}
      else
	{
	  pre->Draw();
	  pre->DrawLatex(0.1,1.15,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
	  pre->DrawLatex(0.9,1.25,"0-100%");
	}
    }

  gPad->RedrawAxis();

  if(bSavePlots && !bDoThermo)
    {
      if(bDoRadius)
	{
	  c1->SaveAs("raaVSr0_minbias_allOnia.pdf");
	  c1->SaveAs("raaVSr0_minbias_allOnia.png");
	}
      else
	{
	  if(bDoRatio)
	    {
	      c1->SaveAs("raaVSdEr0_minbias_allOnia.pdf");
	      c1->SaveAs("raaVSdEr0_minbias_allOnia.png");
	    }
	  else
	    {
	      c1->SaveAs("raaVSdE_minbias_allOnia.pdf");
	      c1->SaveAs("raaVSdE_minbias_allOnia.png");
	    }
	}
    }



  // make a thermometer, if you can

  TH1 *phAxis = new TH1D("phAxis","phAxis",1,0.,1.5);
  //  phAxis->GetYaxis()->CenterTitle(true);
  phAxis->SetMinimum(0.);
  phAxis->SetMaximum(1.4);
  phAxis->GetXaxis()->SetNdivisions(0);
  //  gPad->SetLogx();
 
 
  TLine *ptBin1 = new TLine(0,1,1.5,1);
 
  ptBin1->SetLineWidth(1);
  ptBin1->SetLineStyle(2);
  // vs binding energy
  TCanvas *pc2 = new TCanvas("pc2","pc2",150,600);
  gPad->SetLeftMargin(0.01);
  gPad->SetRightMargin(0.01);
 
  
  phAxis->Draw("");
  gPad->SetTickx(0);
  // gPad->SetTicky(0);
  
  gJpsi_sys->Draw("2");
  gPsiP_sys->Draw("2");
  gUps1_sys->Draw("2");
  gUps2_sys->Draw("2");
  // gUps3_sys->Draw("2");

  gJpsi->Draw("P");
 
  gUps1->Draw("P");
  gUps2->Draw("P");
  gPsiP->Draw("P");
  //  gUps3->Draw("P");

 
  gPsiP_p->Draw("P");
  gUps1_p->Draw("P");
  gUps2_p->Draw("P");
  // gUps3_p->Draw("P");
 gJpsi_p->Draw("P");

  gUps3_limit->Draw("D");
  gUps3_limitArrow->Draw();

  pre->SetTextSize(0.1);
  pre->DrawLatex(0.15,1.3,"CMS Preliminary");
  pre->SetTextSize(0.095);
  pre->DrawLatex(0.15,1.2,"PbPb #sqrt{s_{NN}} = 2.76 TeV");
  pre->DrawLatex(0.15,1.1,"0-100%");

  pre->SetTextSize(0.1);
  pre->DrawLatex(0.7,0.55,"#varUpsilon(1S)");
  pre->DrawLatex(0.7,0.35,"J/#psi");
  pre->DrawLatex(0.7,0.1,"#psi(2S),#varUpsilon(2S)");
  pre->DrawLatex(0.7,0.06,"#varUpsilon(3S)");

  if(bSavePlots && bDoThermo)
    {
      pc2->SaveAs("raa_minbias_allOnia.pdf");
      pc2->SaveAs("raa_minbias_allOnia.png");
    }
  return;
}
