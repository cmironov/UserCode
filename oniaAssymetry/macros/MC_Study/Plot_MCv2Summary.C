#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
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

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
    const Int_t rows, const Float_t leftOffset=0.,
    const Float_t bottomOffset=0., 
    const Float_t leftMargin=0.2, 
    const Float_t bottomMargin=0.2,
    const Float_t edge=0.05);
void drawDum(float min, float max, double drawXLabel);
void CalEffErr(TGraph *a, double *b);
TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2);
TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2);
void formatTGraph(TGraph* a, int b, int c, int d);
void formatTCanv(TCanvas* a);
void formatTLeg(TLegend* a);
void formatTH1F(TH1* a, int b, int c, int d);

TF1 *fGaus = new TF1("fGaus","0.3*TMath::Gaus(x,15,4,0)",0,30);

void Plot_MCv2Summary(){
  //gROOT->ProcessLine(".x ../rootlogon.C");

  gStyle->SetOptFit(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.10);

  TFile *f1 = new TFile("MCCT_PrJpsi_v2_Cent_EffCor.root");
  TFile *f2 = new TFile("MCCT_PrJpsi_v2_Pts_EffCor.root");
  TFile *f3 = new TFile("MCCT_PrJpsi_v2_Rap_EffCor.root");

  TFile *f4 = new TFile("MCCT_PrJpsi_v2_Cent_NoEffCor.root");
  TFile *f5 = new TFile("MCCT_PrJpsi_v2_Pts_NoEffCor.root");
  TFile *f6 = new TFile("MCCT_PrJpsi_v2_Rap_NoEffCor.root");

  TGraphErrors *gGCent[3], *gRCent[3], *gRCCent[3], gRCentNoEff[3];
  TGraphErrors *gGPts[3], *gRPts[3], *gRCPts[3], gRPtsNoEff[3];
  TGraphErrors *gGRapH[3], *gRRapH[3], *gRCRapH[3], gRRapHNoEff[3];
  TGraphErrors *gGRapL[3], *gRRapL[3], *gRCRapL[3], gRRapLNoEff[3];

  for(int iCat = 0; iCat < 3; iCat++){
    bool bDefault = true;
    bool bCowboy = false;
    bool bSailor = false;

    if(iCat == 0) {bDefault = true;bCowboy = false;bSailor = false;}
    if(iCat == 1) {bDefault = false;bCowboy = true;bSailor = false;}
    if(iCat == 2) {bDefault = false;bCowboy = false;bSailor = true;}

    char cCd[512];
    if(bDefault) sprintf(cCd,"default");
    if(bCowboy) sprintf(cCd,"cowboy");
    if(bSailor) sprintf(cCd,"sailor");

    // gCent_rec_default_ResCorr
    char tmp_genCent[512], tmp_recCent[512], tmp_recCorCent[512];
    sprintf(tmp_genCent,"gCent_gen_%s",cCd);
    sprintf(tmp_recCent,"gCent_rec_%s",cCd);
    sprintf(tmp_recCorCent,"gCent_rec_%s",cCd);
    gGCent[iCat] = (TGraphErrors*) f1->Get(tmp_genCent);
    gRCent[iCat] = (TGraphErrors*) f1->Get(tmp_recCent);
    gRCCent[iCat] = (TGraphErrors*) f4->Get(tmp_recCorCent);

    char tmp_genPts[512], tmp_recPts[512], tmp_recCorPts[512];
    sprintf(tmp_genPts,"gPts_gen_%s",cCd);
    sprintf(tmp_recPts,"gPts_rec_%s",cCd);
    sprintf(tmp_recCorPts,"gPts_rec_%s",cCd);
    gGPts[iCat] = (TGraphErrors*) f2->Get(tmp_genPts);
    gRPts[iCat] = (TGraphErrors*) f2->Get(tmp_recPts);
    gRCPts[iCat] = (TGraphErrors*) f5->Get(tmp_recCorPts);

    // gRapsH_rec_default
    char tmp_genRapH[512], tmp_recRapH[512], tmp_recCorRapH[512];
    sprintf(tmp_genRapH,"gRapsH_gen_%s",cCd);
    sprintf(tmp_recRapH,"gRapsH_rec_%s",cCd);
    sprintf(tmp_recCorRapH,"gRapsH_rec_%s",cCd);
    gGRapH[iCat] = (TGraphErrors*) f3->Get(tmp_genRapH);
    gRRapH[iCat] = (TGraphErrors*) f3->Get(tmp_recRapH);
    gRCRapH[iCat] = (TGraphErrors*) f6->Get(tmp_recCorRapH);

    char tmp_genRapL[512], tmp_recRapL[512], tmp_recCorRapL[512];
    sprintf(tmp_genRapL,"gRapsL_gen_%s",cCd);
    sprintf(tmp_recRapL,"gRapsL_rec_%s",cCd);
    sprintf(tmp_recCorRapL,"gRapsL_rec_%s",cCd);
    gGRapL[iCat] = (TGraphErrors*) f3->Get(tmp_genRapL);
    gRRapL[iCat] = (TGraphErrors*) f3->Get(tmp_recRapL);
    gRCRapL[iCat] = (TGraphErrors*) f6->Get(tmp_recCorRapL);

    gRCent[iCat]->SetMarkerStyle(20);
    gRCent[iCat]->SetMarkerColor(kBlue+2);
    gRCent[iCat]->SetLineColor(kBlue+2);
    gRCent[iCat]->SetMarkerSize(1.8);
    gRPts[iCat]->SetMarkerStyle(20);
    gRPts[iCat]->SetMarkerColor(kBlue+2);
    gRPts[iCat]->SetLineColor(kBlue+2);
    gRPts[iCat]->SetMarkerSize(1.8);
    gRRapH[iCat]->SetMarkerStyle(20);
    gRRapH[iCat]->SetMarkerColor(kBlue+2);
    gRRapH[iCat]->SetLineColor(kBlue+2);
    gRRapH[iCat]->SetMarkerSize(1.8);
    gRRapL[iCat]->SetMarkerStyle(29);
    gRRapL[iCat]->SetMarkerColor(kBlue+2);
    gRRapL[iCat]->SetLineColor(kBlue+2);
    gRRapL[iCat]->SetMarkerSize(1.8);

    gRCCent[iCat]->SetMarkerStyle(25);
    gRCCent[iCat]->SetMarkerColor(kRed+2);
    gRCCent[iCat]->SetLineColor(kRed+2);
    gRCCent[iCat]->SetMarkerSize(1.8);
    gRCPts[iCat]->SetMarkerStyle(25);
    gRCPts[iCat]->SetMarkerColor(kRed+2);
    gRCPts[iCat]->SetLineColor(kRed+2);
    gRCPts[iCat]->SetMarkerSize(1.8);
    gRCRapH[iCat]->SetMarkerStyle(25);
    gRCRapH[iCat]->SetMarkerColor(kRed+2);
    gRCRapH[iCat]->SetLineColor(kRed+2);
    gRCRapH[iCat]->SetMarkerSize(1.8);
    gRCRapL[iCat]->SetMarkerStyle(28);
    gRCRapL[iCat]->SetMarkerColor(kRed+2);
    gRCRapL[iCat]->SetLineColor(kRed+2);
    gRCRapL[iCat]->SetMarkerSize(1.8);

    gGCent[iCat]->SetMarkerStyle(24);
    gGCent[iCat]->SetMarkerColor(kMagenta+2);
    gGCent[iCat]->SetLineColor(kMagenta+2);
    gGCent[iCat]->SetMarkerSize(1.8);
    gGPts[iCat]->SetMarkerStyle(24);
    gGPts[iCat]->SetMarkerColor(kMagenta+2);
    gGPts[iCat]->SetLineColor(kMagenta+2);
    gGPts[iCat]->SetMarkerSize(1.8);
    gGRapH[iCat]->SetMarkerStyle(24);
    gGRapH[iCat]->SetMarkerColor(kMagenta+2);
    gGRapH[iCat]->SetLineColor(kMagenta+2);
    gGRapH[iCat]->SetMarkerSize(1.8);
    gGRapL[iCat]->SetMarkerStyle(30);
    gGRapL[iCat]->SetMarkerColor(kRed+2);
    gGRapL[iCat]->SetLineColor(kRed+2);
    gGRapL[iCat]->SetMarkerSize(1.8);
  }

  TCanvas *c1 = new TCanvas("c1","",1000,400);

  makeMultiPanelCanvas(c1,3,1,0.0,0.0,0.15,0.15,0.05); // 0.12, 0.15, 0.02

  TH1F *hPad1 = new TH1F("hPad1","hPad1;N_{part};v_{2}",10,0.0,400.0);
  TH1F *hPad2 = new TH1F("hPad2","hPad2;p_{T} (GeV/c);v_{2}",20,0.0,40.0);
  TH1F *hPad3 = new TH1F("hPad3","hPad3;|y|;v_{2}",10,0.0,2.4);

  hPad1->GetXaxis()->SetLabelSize(20);
  hPad1->GetXaxis()->SetLabelFont(43);
  hPad1->GetXaxis()->SetTitleSize(24);
  hPad1->GetXaxis()->SetTitleFont(43);
  hPad1->GetXaxis()->SetTitleOffset(1.0);
  hPad1->GetXaxis()->CenterTitle();

  hPad1->GetYaxis()->SetLabelSize(17);
  hPad1->GetYaxis()->SetLabelFont(43);
  hPad1->GetYaxis()->SetTitleSize(30);
  hPad1->GetYaxis()->SetTitleFont(43);
  hPad1->GetYaxis()->SetTitleOffset(0.8);
  hPad1->GetYaxis()->CenterTitle();

  hPad2->GetXaxis()->SetLabelSize(20);
  hPad2->GetXaxis()->SetLabelFont(43);
  hPad2->GetXaxis()->SetTitleSize(24);
  hPad2->GetXaxis()->SetTitleFont(43);
  hPad2->GetXaxis()->SetTitleOffset(1.0);
  hPad2->GetXaxis()->CenterTitle();

  hPad2->GetYaxis()->SetLabelSize(17);
  hPad2->GetYaxis()->SetLabelFont(43);
  hPad2->GetYaxis()->SetTitleSize(30);
  hPad2->GetYaxis()->SetTitleFont(43);
  hPad2->GetYaxis()->SetTitleOffset(0.8);
  hPad2->GetYaxis()->CenterTitle();

  hPad3->GetXaxis()->SetLabelSize(20);
  hPad3->GetXaxis()->SetLabelFont(43);
  hPad3->GetXaxis()->SetTitleSize(24);
  hPad3->GetXaxis()->SetTitleFont(43);
  hPad3->GetXaxis()->SetTitleOffset(1.0);
  hPad3->GetXaxis()->CenterTitle();

  hPad3->GetYaxis()->SetLabelSize(17);
  hPad3->GetYaxis()->SetLabelFont(43);
  hPad3->GetYaxis()->SetTitleSize(30);
  hPad3->GetYaxis()->SetTitleFont(43);
  hPad3->GetYaxis()->SetTitleOffset(0.8);
  hPad3->GetYaxis()->CenterTitle();


  TLatex *lt1 = new TLatex();
  lt1->SetNDC();

  hPad1->SetMaximum(0.2);
  hPad2->SetMaximum(0.2);
  hPad3->SetMaximum(0.2);
  hPad1->SetMinimum(-0.10);
  hPad2->SetMinimum(-0.10);
  hPad3->SetMinimum(-0.10);

  // iSpec : 0 = Cent, 1 = Pt, 2 = Rap
  for(int iSpec = 0; iSpec < 3; iSpec++){
      TLegend *leg1 = new TLegend(0.05,0.17,0.37,0.37);
      TLegend *leg2 = new TLegend(0.45,0.60,0.72,0.80);
      TLegend *leg3 = new TLegend(0.45,0.60,0.72,0.80);
      TLegend *leg4 = new TLegend(0.45,0.60,0.72,0.80);
      formatTLeg(leg1);
      formatTLeg(leg2);
      formatTLeg(leg3);
      formatTLeg(leg4);

      for(int icv = 0; icv < 3; icv++){
        c1->cd(icv+1);
        if(iSpec == 0) hPad1->Draw();
        if(iSpec == 1) hPad2->Draw();
        if(iSpec == 2) hPad3->Draw();

        if(iSpec == 0 && icv == 2){
          leg2->AddEntry(gGCent[icv],"MC Truth","PL");
          leg2->AddEntry(gRCCent[icv],"Reco No Corrected","PL");
          leg2->AddEntry(gRCent[icv],"Reco Corrected","PL");
        }
        if(iSpec == 1 && icv == 2){
          leg3->AddEntry(gGPts[icv],"MC Truth","PL");
          leg3->AddEntry(gRCPts[icv],"Reco No Corrected","PL");
          leg3->AddEntry(gRPts[icv],"Reco Corrected","PL");
        }
        if(iSpec == 2 && icv == 2){
          leg4->SetHeader("6.5 < p_{T} < 40 GeV/c");
          leg4->AddEntry(gGRapH[icv],"MC Truth","PL");
          leg4->AddEntry(gRCRapH[icv],"Reco No Corrected","PL");
          leg4->AddEntry(gRRapH[icv],"Reco Corrected","PL");
          leg1->SetHeader("3 < p_{T} < 6.5 GeV/c");
          leg1->AddEntry(gGRapL[icv],"MC Truth","PL");
          leg1->AddEntry(gRCRapL[icv],"Reco No Corrected","PL");
          leg1->AddEntry(gRRapL[icv],"Reco Corrected","PL");
          //leg1->AddEntry(gGRapH[icv],"6.5 < p_{T} < 40 GeV/c","pl");
          //leg1->AddEntry(gGRapL[icv],"3 < p_{T} < 6.5 GeV/c","pl");
          leg1->Draw();
        }


        if(iSpec == 0){
          gGCent[icv]->Draw("[P]");
          gRCent[icv]->Draw("[P]");
          gRCCent[icv]->Draw("[P]");
        }
        if(iSpec == 0 && icv == 2) leg2->Draw();
        if(iSpec == 1 && icv == 2) leg3->Draw();
        if(iSpec == 2 && icv == 2) leg4->Draw();

        if(icv == 0){
          lt1->SetTextSize(0.05);
          lt1->DrawLatex(0.2,0.87,"CMS Simulation");
          lt1->DrawLatex(0.2,0.80,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
          lt1->SetTextSize(0.05);
          if((iSpec == 0)) {
            lt1->DrawLatex(0.2,0.73,"6.5 < p_{T} < 40 GeV/c");
          }
          if(iSpec == 2) {
            lt1->DrawLatex(0.2,0.73,"|y| < 2.4");
          }
          if(iSpec == 1) {
            lt1->DrawLatex(0.2,0.73,"|y| < 2.4");
            lt1->DrawLatex(0.40,0.18,"PYTHIA+EvtGen+HYDJET");
          }
          lt1->SetTextSize(0.045);
          if(!(iSpec == 1))lt1->DrawLatex(0.2,0.18,"PYTHIA+EvtGen+HYDJET");
        }
        if(icv == 2){
          lt1->SetTextSize(0.055);
          lt1->DrawLatex(0.05,0.87,"Prompt J/#psi");
        }

        if(icv == 0) {
          lt1->SetTextSize(0.05);
          lt1->DrawLatex(0.80,0.87,"Default");
        }
        if(icv == 1) {
          lt1->SetTextSize(0.055);
          if(iSpec == 0) lt1->DrawLatex(0.05,0.87,"|y| < 2.4");
          if(!(iSpec == 0)) lt1->DrawLatex(0.05,0.87,"Cent. 10 - 60 %");
          lt1->DrawLatex(0.78,0.87,"Cowboy");
        }
        if(icv == 2) {
          lt1->SetTextSize(0.055);
          lt1->DrawLatex(0.78,0.87,"Sailor");
        }

        if(iSpec == 1){
          gGPts[icv]->Draw("[P]");
          gRPts[icv]->Draw("[P]");
          gRCPts[icv]->Draw("[P]");
        }

        if(iSpec == 2){
          gGRapH[icv]->Draw("[P]");
          gRRapH[icv]->Draw("[P]");
          gRCRapH[icv]->Draw("[P]");
          gGRapL[iSpec]->Draw("[P]");
          gRRapL[iSpec]->Draw("[P]");
          gRCRapL[iSpec]->Draw("[P]");
          char legs[512];
          //leg1->SetHeader("0.0 < | y | < 1.2");
        }
      }

      //leg2->Delete();
      //leg3->Delete();
      //leg4->Delete();

      char tmp_png[512], tmp_pdf[512];
      char *cats[3] = {"cent","pt","rap"};
      sprintf(tmp_png,"MCv2_MCCT_summary_%s_ZeroV2.png",cats[iSpec]);
      sprintf(tmp_pdf,"MCv2_MCCT_summary_%s_ZeroV2.pdf",cats[iSpec]);
      c1->SaveAs(tmp_png);
      c1->SaveAs(tmp_pdf);
  }
}
//(TH1, color, style, pt, eta, rapidity)
void formatTH1F(TH1* a, int b, int c, int d){
  a->SetLineWidth(2);
  a->SetLineStyle(c);
  a->SetMarkerSize(2);
  a->SetLineColor(b);
  a->SetMarkerColor(b);
  a->GetYaxis()->SetTitle("Single #mu Efficiency");
  if(d == 1){	
    a->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/c)"); 
    //a->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  }else if(d == 2){
    a->GetXaxis()->SetTitle("#eta^{#mu}"); 
  }else if(d == 3){
    a->GetXaxis()->SetTitle("rapidity^{#mu}"); 
  }else if(d == 4){
    a->GetXaxis()->SetTitle("Centrality");
  }else if(d == 5){
    a->GetXaxis()->SetTitle("Centrality (%)");
  }
  a->GetXaxis()->CenterTitle(true);
  //a->GetXaxis()->SetLabelSize(0.05);
  //a->GetXaxis()->SetTitleSize(0.05);
  //a->GetXaxis()->SetTitleOffset(1.1);
  //a->GetYaxis()->SetLabelSize(0.05);
  //a->GetYaxis()->SetTitleSize(0.05);
  //a->GetYaxis()->SetTitleOffset(1.2);

}     

void formatTLeg(TLegend* a){

  a->SetFillStyle(0);  
  a->SetFillColor(0); 
  a->SetBorderSize(0);
  a->SetTextSize(0.04);
  a->SetTextFont(42);
  //a->SetTextSize(20);
}

void formatTCanv(TCanvas* a){
  a->SetBorderSize(2);
  a->SetFrameFillColor(0);
  a->cd();
  a->SetGrid(1);
  a->SetTickx();
  a->SetTicky();
}

void formatTGraph(TGraph* a, int b, int c, int d)
{

  a->SetMarkerStyle(c);
  a->SetMarkerColor(b);
  a->SetMarkerSize(1.0);
  a->SetLineColor(b);
  a->SetLineWidth(1);
  a->GetXaxis()->SetLabelSize(0.05);
  a->GetXaxis()->SetTitleSize(0.06);
  a->GetXaxis()->SetTitleOffset(1.0);
  a->GetYaxis()->SetTitle("Efficiency");
  //a->GetXaxis()->CenterTitle();
  if(d == 1){	
    a->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  }else if(d == 2){
    a->GetXaxis()->SetTitle("eta"); 
  }else if(d == 3){
    a->GetXaxis()->SetTitle("rapidity"); 
  }
  a->GetYaxis()->SetRangeUser(0,1);
  a->GetXaxis()->SetRangeUser(0,10);
  a->GetYaxis()->SetLabelSize(0.05);
  a->GetYaxis()->SetTitleSize(0.05);
  a->GetYaxis()->SetTitleOffset(1.3);


}

TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2)
{
  TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
  gEfficiency->BayesDivide(h1,h2);
  return gEfficiency;
}

TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2)
{

  h1->Sumw2();
  h2->Sumw2();

  TGraphAsymmErrors *result = calcEff(h1,h2);
  return result;
}
void CalEffErr(TGraph *a, double *b){
  const int nbins = 100;
  double x[nbins], y[nbins];
  double sum = 0, errHighSum = 0, errLowSum = 0, sqSumHigh = 0, sqSumLow = 0;
  //double b[3] = 0;

  int nBins = a->GetN();
  for(int i=0;i<a->GetN();i++){
    a->GetPoint(i,x[i],y[i]);
    //cout<<"Eff x = "<<x[i]<<" y = "<<y[i]<<endl;
    double eHigh = a->GetErrorYhigh(i);
    double eLow = a->GetErrorYlow(i);
    //cout<<"Err high = "<<eHigh<<" low = "<<eLow<<endl;
    sum += y[i];
    errHighSum += eHigh;
    sqSumHigh += eHigh*eHigh;
    errLowSum += eLow;
    sqSumLow += eLow*eLow;
  }
  b[0] = sum/nBins;
  b[1] = sqrt(sqSumHigh)/nBins;
  b[2] = sqrt(sqSumLow)/nBins;
  //cout<<"b1 : "<<b[0]<<", b2 : "<<b[1]<<", b3 : "<<b[2]<<endl;

  cout<<b[0]<<"^{"<<b[1]<<"}_{"<<b[2]<<"}"<<endl;
  //return b[3];
}

void makeMultiPanelCanvas(TCanvas*& canv,
    const Int_t columns,
    const Int_t rows,
    const Float_t leftOffset,
    const Float_t bottomOffset,
    const Float_t leftMargin,
    const Float_t bottomMargin,
    const Float_t edge) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth = 
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
        (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
        (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }

  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
          Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);

      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}

void drawDum(float min, float max, double drawXLabel){

  TH1D *hdum = new TH1D("hdum","",20,0,1);
  hdum->SetMaximum(max);

  hdum->SetStats(0);

  if(drawXLabel) hdum->SetXTitle("A_{J} = (E_{T}^{j1}-E_{T}^{j2})/(E_{T}^{j1}+E_{T}^{j2})");
  hdum->GetXaxis()->SetLabelSize(20);
  hdum->GetXaxis()->SetLabelFont(43);
  hdum->GetXaxis()->SetTitleSize(22);
  hdum->GetXaxis()->SetTitleFont(43);
  hdum->GetXaxis()->SetTitleOffset(1.5);
  hdum->GetXaxis()->CenterTitle();

  hdum->GetXaxis()->SetNdivisions(905,true);

  hdum->SetYTitle("Ratio");

  hdum->GetYaxis()->SetLabelSize(20);
  hdum->GetYaxis()->SetLabelFont(43);
  hdum->GetYaxis()->SetTitleSize(20);
  hdum->GetYaxis()->SetTitleFont(43);
  hdum->GetYaxis()->SetTitleOffset(2.5);
  hdum->GetYaxis()->CenterTitle();

  hdum->SetAxisRange(0,0.2,"Y");

  hdum->Draw("");

}

