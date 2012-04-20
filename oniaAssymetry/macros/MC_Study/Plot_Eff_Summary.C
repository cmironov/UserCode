#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TNtuple.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

void Plot_Eff_Summary(){
  //gROOT->ProcessLine(".x ../rootlogon.C");

  for(int iCat = 0; iCat < 3; iCat++){
    bool bDefault = true; // true : default, false : sailor or cowboy
    bool bCowboy  = false; // true : cowboy only, false : salior
    bool bSailor  = false; // true : cowboy only, false : salior

    if(iCat == 0) {bDefault = true;bCowboy = false;bSailor = false;}
    if(iCat == 1) {bDefault = false;bCowboy = true;bSailor = false;}
    if(iCat == 2) {bDefault = false;bCowboy = false;bSailor = true;}

    char cCd[512], pfx[512];
    if(bDefault) {sprintf(cCd, "default"); sprintf(pfx, "Default");}
    if(bCowboy)  {sprintf(cCd, "cowboy"); sprintf(pfx, "Cowboy");}
    if(bSailor)  {sprintf(cCd, "sailor"); sprintf(pfx, "Sailor");}


    gStyle->SetOptFit(0);
    //gStyle->SetTitleFont(62,"xyz");
    //gStyle->SetLabelFont(62,"xyz");
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    // PrJpsi_cowboy_etHFm.root
    char tmp_input[10][512];
    sprintf(tmp_input[0],"PrJpsi_%s.root",cCd);
    TFile *f1 = new TFile(tmp_input[0]);

    const int nCentBins = 6;
    const int nPtBins = 6;
    const int nRapBins = 5;

    double ct_bound[nCentBins] = {0, 4, 8, 12, 24, 40};
    double xct_bound[nCentBins] = {0.0};
    double pt_bound[nPtBins+1] = {3.0, 4.5, 6.5, 7.25, 8.0, 10.0, 40.0};
    double xpt_bound[nPtBins] = {0.0};
    double rap_bound[nRapBins+1] = {-2.4, -1.6, -1.2, 1.2, 1.6, 2.4}; 
    double xrap_bound[nRapBins] = {0.0};


    TH1F *h1D[nCentBins][nPtBins];

    TCanvas *c1 = new TCanvas("c1","",0,0,500,500);
    // eff2D_Pt_0_cowboy
    // eff1D_Pt_0_Rap_0_cowboy
    char tmpIpt[512], tmpIpt_png[512], tmpIrap[512], tmpIrap_png[512];
    for(int ict = 0; ict < nCentBins; ict++){
      for(int ipt = 0; ipt < nPtBins; ipt++){
        if(ict == (nCentBins-1)){
          sprintf(tmpIpt,"eff1D_Cent_%1.f_%1.f_Pt_%0.1f_%0.1f_%s",ct_bound[ict-4],ct_bound[ict-1],pt_bound[ipt],pt_bound[ipt+1],cCd);
        }else{
          sprintf(tmpIpt,"eff1D_Cent_%1.f_%1.f_Pt_%0.1f_%0.1f_%s",ct_bound[ict],ct_bound[ict+1],pt_bound[ipt],pt_bound[ipt+1],cCd);
        }
        h1D[ict][ipt] = (TH1F*)f1->Get(tmpIpt);
        h1D[ict][ipt]->SetMarkerSize(1.8);
        if(ipt == 0) {
          h1D[ict][ipt]->SetMarkerStyle(20);
          h1D[ict][ipt]->SetMarkerColor(kBlue+2);
          h1D[ict][ipt]->SetLineColor(kBlue+2);
        }
        if(ipt == 1) {
          h1D[ict][ipt]->SetMarkerStyle(24);
          h1D[ict][ipt]->SetMarkerColor(kBlue+2);
          h1D[ict][ipt]->SetLineColor(kBlue+2);
        }
        if(ipt == 2) {
          h1D[ict][ipt]->SetMarkerStyle(21);
          h1D[ict][ipt]->SetMarkerColor(kRed+2);
          h1D[ict][ipt]->SetLineColor(kRed+2);
        }
        if(ipt == 3) {
          h1D[ict][ipt]->SetMarkerStyle(25);
          h1D[ict][ipt]->SetMarkerColor(kRed+2);
          h1D[ict][ipt]->SetLineColor(kRed+2);
        }
        if(ipt == 4) {
          h1D[ict][ipt]->SetMarkerStyle(30);
          h1D[ict][ipt]->SetMarkerColor(kGreen+2);
          h1D[ict][ipt]->SetLineColor(kGreen+2);
        }
        if(ipt == 5) {
          h1D[ict][ipt]->SetMarkerStyle(29);
          h1D[ict][ipt]->SetMarkerColor(kGreen+2);
          h1D[ict][ipt]->SetLineColor(kGreen+2);
        }
      }
    }

    //makeMultiPanelCanvas(c1,2,1,0.0,0.0,0.15,0.15,0.05); // 0.12, 0.15, 0.02


    TH1F *hPad1 = new TH1F("hPad1","hPad1;rapidity;Efficiency",10,-2.4,2.4);
    hPad1->GetXaxis()->SetLabelSize(20);
    hPad1->GetXaxis()->SetLabelFont(43);
    hPad1->GetXaxis()->SetTitleSize(24);
    hPad1->GetXaxis()->SetTitleFont(43);
    hPad1->GetXaxis()->SetTitleOffset(1.3);
    hPad1->GetXaxis()->CenterTitle();

    hPad1->GetYaxis()->SetLabelSize(20);
    hPad1->GetYaxis()->SetLabelFont(43);
    hPad1->GetYaxis()->SetTitleSize(30);
    hPad1->GetYaxis()->SetTitleFont(43);
    hPad1->GetYaxis()->SetTitleOffset(1.0);
    hPad1->GetYaxis()->CenterTitle();

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();

    hPad1->SetMaximum(1.2);
    hPad1->SetMinimum(0.0);

    //c1->Divide(3,1);

    //TLegend *leg1 = new TLegend(0.19,0.50,0.41,0.70);
    //formatTLeg(leg1);
    //leg1->SetHeader("0.0 < | y | < 1.2");

    char tmp_leg[512], tmp_1D_png[512];
    char tmp_pt[512], tmp_cCd[512];
    for(int ict = 0; ict < nCentBins; ict++){
      TLegend *leg1 = new TLegend(0.55,0.55,0.90,0.80);
      formatTLeg(leg1);
      hPad1->Draw();
      for(int ipt = 0; ipt < nPtBins; ipt++){
        h1D[ict][ipt]->Draw("elp same");
        sprintf(tmp_leg,"%0.1f < p_{T} < %0.1f",pt_bound[ipt],pt_bound[ipt+1]);
        leg1->AddEntry(h1D[ict][ipt],tmp_leg,"pl");
      }
      if(ict == (nCentBins-1)){
        sprintf(tmp_1D_png,"plot_eff_1D_summary_Cent_%1.f_%1.f_%s.png",ct_bound[ict-4],ct_bound[ict-1],cCd);
      }else{
        sprintf(tmp_1D_png,"plot_eff_1D_summary_Cent_%1.f_%1.f_%s.png",ct_bound[ict],ct_bound[ict+1],cCd);
      }
      leg1->SetTextSize(0.035);
      leg1->Draw();

      lt1->SetTextSize(0.045);
      lt1->DrawLatex(0.2,0.88,"CMS Simulation");
      lt1->DrawLatex(0.2,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
      //lt1->DrawLatex(0.2,0.76,"EP: etHFp+etHFm");
      lt1->SetTextSize(0.04);
      if(ict == (nCentBins-1)){
        sprintf(tmp_pt,"Cent. %1.f - %1.f %%",ct_bound[ict-4]*2.5, ct_bound[ict-1]*2.5);
      }else{
        sprintf(tmp_pt,"Cent. %1.f - %1.f %%",ct_bound[ict]*2.5, ct_bound[ict+1]*2.5);
      }
      lt1->DrawLatex(0.2,0.76,tmp_pt);
      //if(bEtHFm) lt1->DrawLatex(0.2,0.70,"#eta^{J/#psi} > 0");
      //if(!bEtHFm) lt1->DrawLatex(0.2,0.70,"#eta^{J/#psi} < 0");
      lt1->DrawLatex(0.68,0.88,"Prompt J/#psi");
      cout<<"cCd1 : "<<pfx<<endl;
      sprintf(tmp_cCd,"%s",pfx);
      lt1->DrawLatex(0.68,0.82,tmp_cCd);
      lt1->SetTextSize(0.043);
      lt1->DrawLatex(0.2,0.18,"PYTHIA+EvtGen+HYDJET");

      c1->SaveAs(tmp_1D_png);
      leg1->Delete();
      //leg1->Close();
    }
    /*

       lt1->SetTextSize(0.045);
       lt1->DrawLatex(0.2,0.88,"CMS Simulation");
       lt1->DrawLatex(0.2,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
       lt1->DrawLatex(0.2,0.76,"EP: etHFp+etHFm");
       lt1->SetTextSize(0.04);
       lt1->DrawLatex(0.2,0.66,"6.5 < p_{T} < 40 GeV/c");
       lt1->DrawLatex(0.2,0.60,"|y| < 2.4");
       lt1->DrawLatex(0.68,0.88,"Prompt J/#psi");
       lt1->SetTextSize(0.043);
       lt1->DrawLatex(0.2,0.18,"PYTHIA+EvtGen+HYDJET");
    //lt1->SetTextSize(0.04);
    //lt1->DrawLatex(0.54,0.72,"Prompt J/#psi");
    //lt1->DrawLatex(0.7,0.30,"Prompt J/#psi");
    //leg1->Draw();

    TLegend *leg1 = new TLegend(0.19,0.24,0.41,0.34);
    formatTLeg(leg1);
    //leg1->SetHeader("0.0 < | y | < 1.2");
    leg1->AddEntry(gRaps0,"6.5 < p_{T} < 40","pl");
    leg1->AddEntry(gRaps1,"3 < p_{T} < 6.5","pl");

    leg1->Draw();
    c1->SaveAs("MCv2_Raps_0323_default.png");
    c1->SaveAs("MCv2_Raps_0323_default.pdf");
    */

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

