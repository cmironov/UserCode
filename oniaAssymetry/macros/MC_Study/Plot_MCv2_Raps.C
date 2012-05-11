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
void GetEff(TGraph *a, double *b);
void formatTH1F(TH1* a, int b, int c, int d);

// Function determination to fit v2
TF1 *Fit1 = new TF1("Fit1","[0]",0,TMath::Pi()/2);

void GetV2(TGraph *a, double *b, int c);

// Function determination to fit v2
TF1 *v2Fit1 = new TF1("v2Fit1","(2/TMath::Pi())*(1+2*[0]*TMath::Cos(2.0*x))",0,TMath::PiOver2());

void Plot_MCv2_Raps(){
  //gROOT->ProcessLine(".x ../rootlogon.C");

  gStyle->SetOptFit(0);
  //gStyle->SetTitleFont(62,"xyz");
  //gStyle->SetLabelFont(62,"xyz");
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.065);
  gStyle->SetPadLeftMargin(0.16);

  TFile *output = new TFile("MCCT_PrJpsi_v2_Rap_EffCor.root","RECREATE");

  double rBins[] = {0.6, 1.4, 2.0};
  double rBinsErr[] = {0.0, 0.0, 0.0};
  double rBins1[] = {2.0};
  double rBinsErr1[] = {0.0};

  double res = 0.721352;
  double resErr = 0.00192809;

  for(int iCat = 0; iCat < 2; iCat++){
    for(int iSpec = 0; iSpec < 3; iSpec++){
      bool bDefault = true;
      bool bCowboy = false;
      bool bSailor = false;
      bool bReco = true;
      if(iCat == 0) bReco = true;
      if(iCat == 1) bReco = false;

      if(iSpec == 0) {bDefault = true;bCowboy = false;bSailor = false;}
      if(iSpec == 1) {bDefault = false;bCowboy = true;bSailor = false;}
      if(iSpec == 2) {bDefault = false;bCowboy = false;bSailor = true;}

      char cCd[512], cRec[512];
      if(bDefault) sprintf(cCd, "default");
      if(bCowboy)  sprintf(cCd, "cowboy");
      if(bSailor)  sprintf(cCd, "sailor");
      if(bReco) sprintf(cRec, "rec");
      if(!bReco) sprintf(cRec, "gen");

      char tmp_input[10][512];
      sprintf(tmp_input[0],"MCCT_PrJpsi_Raps_0012_dPhi_%s.root",cCd);
      sprintf(tmp_input[1],"MCCT_PrJpsi_Raps_1216_dPhi_%s.root",cCd);
      sprintf(tmp_input[2],"MCCT_PrJpsi_Raps_1624H_dPhi_%s.root",cCd);
      sprintf(tmp_input[3],"MCCT_PrJpsi_Raps_1624L_dPhi_%s.root",cCd);
      TFile *f1 = new TFile(tmp_input[0]);
      TFile *f2 = new TFile(tmp_input[1]);
      TFile *f3 = new TFile(tmp_input[2]);
      TFile *f4 = new TFile(tmp_input[3]);


      TGraphErrors *g[10];
      if(bReco){
        g[0] = (TGraphErrors*) f1->Get("nReco_Jpsi_Raps_0012");
        g[1] = (TGraphErrors*) f2->Get("nReco_Jpsi_Raps_1216");
        g[2] = (TGraphErrors*) f3->Get("nReco_Jpsi_Raps_1624H");
        g[3] = (TGraphErrors*) f4->Get("nReco_Jpsi_Raps_1624L");
      }else{
        g[0] = (TGraphErrors*) f1->Get("nGen_Jpsi_Raps_0012");
        g[1] = (TGraphErrors*) f2->Get("nGen_Jpsi_Raps_1216");
        g[2] = (TGraphErrors*) f3->Get("nGen_Jpsi_Raps_1624H");
        g[3] = (TGraphErrors*) f4->Get("nGen_Jpsi_Raps_1624L");
      }


      TH1F *hPad1 = new TH1F("hPad1","hPad1;| #phi^{J/#psi} -  #Psi_{RP} |(rad);#frac{1}{N_{total J/#psi}} #frac{dN}{d#phi}(rad^{-1})",8,0,TMath::Pi()/2);

      hPad1->GetXaxis()->SetLabelSize(23);
      hPad1->GetXaxis()->SetLabelFont(43);
      hPad1->GetXaxis()->SetTitleSize(20);
      hPad1->GetXaxis()->SetTitleFont(43);
      hPad1->GetXaxis()->SetTitleOffset(1.2);
      hPad1->GetXaxis()->CenterTitle();

      hPad1->GetYaxis()->SetLabelSize(23);
      hPad1->GetYaxis()->SetLabelFont(43);
      hPad1->GetYaxis()->SetTitleSize(20);
      hPad1->GetYaxis()->SetTitleFont(43);
      hPad1->GetYaxis()->SetTitleOffset(1.7);
      hPad1->GetYaxis()->CenterTitle();

      TH1F *hPad2 = (TH1F*)hPad1->Clone();
      TH1F *hPad3 = (TH1F*)hPad1->Clone();

      hPad1->SetMaximum(2.0);
      hPad1->SetMinimum(0.0);

      TLatex *lt1 = new TLatex();
      lt1->SetNDC();

      int cts[4] = {0, 10, 30, 60};

      //TCanvas *c1;
      TCanvas *c1 = new TCanvas("c1","",0,0,500,500);
      //c1->Divide(3,2);

      TLegend *leg1 = new TLegend(0.17,0.53,0.41,0.71);
      formatTLeg(leg1);
      leg1->SetTextSize(0.03);
      char tmp1[512], tmp2[512];
      char tmp3[512], tmp4[512];

      int l1 = 0, l2 = 0, l3 = 0, l4 = 0;
      double b[4] = {0}; double vTwo[10] = {0}; double vTwoErr[10] = {0};
      double chis[4] = {0}; double ndf[4] = {0};
      for(int k = 0; k < 1; k++){
        l1 = 4*k; l2 = 4*k + 1; l3 = 4*k + 2; l4 = 4*k + 3;

        g[l1]->SetMarkerStyle(20);
        g[l2]->SetMarkerStyle(21);
        g[l3]->SetMarkerStyle(34);
        g[l4]->SetMarkerStyle(29);
        g[l1]->SetMarkerSize(1.3);
        g[l2]->SetMarkerSize(1.3);
        g[l3]->SetMarkerSize(1.6);
        g[l4]->SetMarkerSize(1.8);
        g[l1]->SetMarkerColor(kBlue+2);
        g[l2]->SetMarkerColor(kRed+2);
        g[l3]->SetMarkerColor(kGreen+2);
        g[l4]->SetMarkerColor(kMagenta+2);

        GetV2(g[l1],b,862); vTwo[l1]=b[0]; vTwoErr[l1]=b[1]; chis[l1]=b[2]; ndf[l1]=b[3];
        GetV2(g[l2],b,902); vTwo[l2]=b[0]; vTwoErr[l2]=b[1]; chis[l2]=b[2]; ndf[l2]=b[3];
        GetV2(g[l3],b,844); vTwo[l3]=b[0]; vTwoErr[l3]=b[1]; chis[l3]=b[2]; ndf[l3]=b[3];
        GetV2(g[l4],b,871); vTwo[l4]=b[0]; vTwoErr[l4]=b[1]; chis[l4]=b[2]; ndf[l4]=b[3];


        hPad1->Draw();
        g[l1]->Draw("pz");
        g[l2]->Draw("pz");
        g[l3]->Draw("pz");
        g[l4]->Draw("pz");

        lt1->SetTextFont(42);
        lt1->SetTextSize(0.04);
        lt1->DrawLatex(0.18, 0.81, "Cent. 10 - 60 %");
        lt1->DrawLatex(0.18, 0.75, "6.5 < p_{T} < 40.0 GeV/c");
      }

      //v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)
      sprintf(tmp2,"0 < |y| < 1.2 v_{2}: %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)", vTwo[0], vTwoErr[0], chis[0], ndf[0]);
      leg1->AddEntry(g[0],tmp2);
      sprintf(tmp2,"1.2 < |y| < 1.6 v_{2}: %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)", vTwo[1], vTwoErr[1], chis[1], ndf[1]);
      leg1->AddEntry(g[1],tmp2);
      sprintf(tmp2,"1.6 < |y| < 2.4 v_{2}: %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)", vTwo[2], vTwoErr[2], chis[2], ndf[2]);
      leg1->AddEntry(g[2],tmp2);
      sprintf(tmp2,"1.6 < |y| < 2.4 v_{2}: %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)", vTwo[3], vTwoErr[3], chis[3], ndf[3]);
      leg1->AddEntry(g[3],tmp2);
      leg1->Draw();
      lt1->SetTextSize(0.04);
      lt1->DrawLatex(0.645,0.88,"CMS Simulation");
      char tmp_cd[512];
      if(bReco) sprintf(tmp_cd,"Reco %s",cCd);
      if(!bReco) sprintf(tmp_cd,"Gen %s",cCd);
      lt1->DrawLatex(0.70,0.78,tmp_cd);


      lt1->SetTextSize(0.05);
      lt1->DrawLatex(0.18,0.88,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");

      lt1->SetTextSize(0.045);
      lt1->DrawLatex(0.40,0.20,"PYTHIA+EvtGen+HYDJET");

      char tmp_png[512],tmp_pdf[512];
      sprintf(tmp_png,"MCv2_%s_Raps_0405_%s.png", cRec, cCd);
      sprintf(tmp_pdf,"MCv2_%s_Raps_0405_%s.pdf", cRec, cCd);
      c1->SaveAs(tmp_png);
      c1->SaveAs(tmp_pdf);

      for(int i = 0; i < 4; i++){
        cout<<"v_{2} : "<<vTwo[i]<<" #pm "<<vTwoErr[i]<<", #chi^{2} / ndf = "<<chis[i]<<"/"<<ndf[i]<<endl;
      }

      double vTwoH[3] = {vTwo[0], vTwo[1], vTwo[2]};
      double vTwoL[1] = {vTwo[3]};
      double vTwoHErr[3] = {vTwoErr[0], vTwoErr[1], vTwoErr[2]};
      double vTwoLErr[1] = {vTwoErr[3]};
      char tmpRaps[512];
      sprintf(tmpRaps,"gRapsH_%s_%s",cRec,cCd);
      TGraphErrors *gRapsH = new TGraphErrors(3, rBins, vTwoH, rBinsErr, vTwoHErr);
      gRapsH->SetName(tmpRaps);
      output->cd();
      gRapsH->Write();

      sprintf(tmpRaps,"gRapsL_%s_%s",cRec,cCd);
      //TGraphErrors *gRapsL = new TGraphErrors(1, rBins1, vTwo[4], rBinsErr1, vTwoErr[4]);
      cout<<"vTwoL : "<<vTwoL[0]<<endl;
      TGraphErrors *gRapsL = new TGraphErrors(1, rBins1, vTwoL, rBinsErr1, vTwoLErr);
      gRapsL->SetName(tmpRaps);
      output->cd();
      gRapsL->Write();

      double vTwoCorr[4] = {0.0};
      double vTwoCorrErr[4] = {0.0};
      for(int i = 0; i < 4; i++){
        vTwoCorr[i] = (double) vTwo[i]/res; 
        vTwoCorrErr[i] = (double) (sqrt((resErr/res)*(resErr/res)+(vTwoErr[i]/vTwo[i])*(vTwoErr[i]/vTwo[i])));
        vTwoCorrErr[i] = vTwoCorr[i] * vTwoCorrErr[i];
      }
      double vTwoHCorr[3] = {vTwoCorr[0], vTwoCorr[1], vTwoCorr[2]};
      double vTwoLCorr[1] = {vTwoCorr[4]};
      double vTwoHCorrErr[3] = {vTwoCorrErr[0], vTwoCorrErr[1], vTwoCorrErr[2]};
      double vTwoLCorrErr[1] = {vTwoCorrErr[4]};

      sprintf(tmpRaps,"gRapsH_%s_%s_ResCorr",cRec,cCd);
      TGraphErrors *gRapsHCor = new TGraphErrors(3, rBins, vTwoHCorr, rBinsErr, vTwoHCorrErr);
      gRapsHCor->SetName(tmpRaps);
      output->cd();
      gRapsHCor->Write();

      sprintf(tmpRaps,"gRapsL_%s_%s_ResCorr",cRec,cCd);
      TGraphErrors *gRapsLCor = new TGraphErrors(1, rBins1, vTwoLCorr, rBinsErr1, vTwoLCorrErr);
      gRapsLCor->SetName(tmpRaps);
      output->cd();
      gRapsLCor->Write();
    }
  }
  output->Write();
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
  //a->SetTextFont(44);
  //a->SetTextFont(63);
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

void GetEff(TGraph *a, double *b){

  Fit1->SetLineStyle(2);
  Fit1->SetLineWidth(2);
  Fit1->SetLineColor(kOrange+2);

  /*
     a->SetMarkerStyle(20);
     a->SetMarkerColor(kBlack);
     a->SetMarkerSize(1.2);
     */
  a->Fit(Fit1,"ER");
  b[0] = Fit1->GetParameter(0);
  b[1] = Fit1->GetParError(0);
  b[2] = Fit1->GetChisquare();
  b[3] = Fit1->GetNDF();
  //cout<<"Fit : "<<b[0]<<" "<<b[1]<<endl;

}

void GetV2(TGraph *a, double *b, int c){

  v2Fit1->SetLineWidth(1);
  v2Fit1->SetLineColor(c);
  v2Fit1->SetLineStyle(2);

  a->Fit(v2Fit1,"R");
  //a->Fit(v2Fit1,"qr");
  b[0] = v2Fit1->GetParameter(0);
  b[1] = v2Fit1->GetParError(0);
  b[2] = v2Fit1->GetChisquare();
  b[3] = v2Fit1->GetNDF();

}
