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
    //gStyle->SetTitleFont(62,"xyz");
    //gStyle->SetLabelFont(62,"xyz");
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.10);

    double cBins[] = {355.4, 261.4178, 187.1470, 89.9};
    double cBinsErr[] = {0.0, 0.0, 0.0, 0.0};
    double ptBins[] = {7.3, 9.0, 13.4};
    double ptBinsErr[] = {0.0, 0.0, 0.0};
    double rBins[] = {0.6, 1.4, 2.0};
    double rBinsErr[] = {0.0, 0.0, 0.0};
    double rBins1[] = {2.0};
    double rBinsErr1[] = {0.0};

    /*
    double rCt0[] = {0.0178, -0.0098, -0.0041, 0.0291};
    double rCtErr0[] = {0.0189, 0.0146, 0.0135, 0.0104};
    double rPt0[] = {-0.0125, 0.0011, 0.0076};
    double rPtErr0[] = {0.0227, 0.0197, 0.0086};

    double rRap0[] = {0.0191, 0.0027, -0.0160};
    double rRapErr0[] = {0.0142, 0.0213, 0.0193};
    double rRap1[] = {0.0676};
    double rRapErr1[] = {0.0353};
    */
    /*
    double gCt0[] = {0.1390, 0.1365, 0.1318, 0.1348};
    double gCtErr0[] = {0.0152, 0.0164, 0.0179, 0.0130};
    double gPt0[] = {0.0213, 0.1114, 0.2674};
    double gPtErr0[] = {0.0127, 0.0143, 0.0224};

    double gRap0[] = {0.1396, 0.1013, 0.1237};
    double gRapErr0[] = {0.0144, 0.0155, 0.0153};
    double gRap1[] = {0.0570};
    double gRapErr1[] = {0.0155};
    */

    // After Eff Correction and No Res corr
    double rCt0[] = {0.0951,0.1043,0.0785,0.1262};
    double rCtErr0[] = {0.0198,0.0210,0.0196,0.0183};
    double rPt0[] = {0.0398,0.094,0.2155};
    double rPtErr0[] = {0.0197,0.0191,0.0254};

    double rRap0[] = {0.1271,0.1045,0.0788};
    double rRapErr0[] = {0.0196,0.0225,0.0200};
    double rRap1[] = {0.0884};
    double rRapErr1[] = {0.0344};

    // After Eff Correction and resolution 
    double rACt0[] = {0.1572,0.1330,0.0998,0.1966};
    double rACtErr0[] = {0.0231,0.0189,0.0176,0.0202};
    double rAPt0[] = {0.0552,0.1303,0.2987};
    double rAPtErr0[] = {0.0193,0.0187,0.0249};

    double rARap0[] = {0.1762,0.1449,0.1092};
    double rARapErr0[] = {0.0192,0.0221,0.0196};
    double rARap1[] = {0.0884};
    double rARapErr1[] = {0.0344};

    // After Eff and No resolution correction with Cent Free
    double rBCt0[] = {0.0951,0.1043,0.0785,0.1262};
    double rBCtErr0[] = {0.0231,0.0189,0.0176,0.0202};
    double rBPt0[] = {0.035,0.0944,0.2152};
    double rBPtErr0[] = {0.0196,0.0197,0.0254};

    double rBRap0[] = {0.1275,0.1022,0.0775};
    double rBRapErr0[] = {0.0193,0.0226,0.0204};
    double rBRap1[] = {0.0891};
    double rBRapErr1[] = {0.0382};

    /*
    // After Eff and Resolution correction with Cent Free
    double rBCt0[] = {0.1572,0.1330,0.0998,0.1966};
    double rBCtErr0[] = {0.0231,0.0189,0.0176,0.0202};
    double rBPt0[] = {0.0485,0.1309,0.2983};
    double rBPtErr0[] = {0.0192,0.0193,0.0249};

    double rBRap0[] = {0.1768,0.1417,0.1074};
    double rBRapErr0[] = {0.0189,0.0222,0.0200};
    double rBRap1[] = {0.0891};
    double rBRapErr1[] = {0.0382};
    */


    // Gen
    double gCt0[] = {0.1148, 0.1006, 0.0830, 0.1056};
    double gCtErr0[] = {0.0156, 0.0146, 0.0131, 0.0133};
    double gPt0[] = {0.0313, 0.0871, 0.2097};
    double gPtErr0[] = {0.0104, 0.0132, 0.0243};

    double gRap0[] = {0.1090, 0.0782, 0.0966};
    double gRapErr0[] = {0.0143, 0.0137, 0.0143};
    double gRap1[] = {0.0438};
    double gRapErr1[] = {0.0125};

    TGraphErrors *gCent0 = new TGraphErrors(4, cBins, gCt0, cBinsErr, gCtErr0);
    TGraphErrors *gPts0 = new TGraphErrors(3, ptBins, gPt0, ptBinsErr, gPtErr0);
    TGraphErrors *gRaps0 = new TGraphErrors(3, rBins, gRap0, rBinsErr, gRapErr0);
    TGraphErrors *gRaps1 = new TGraphErrors(1, rBins1, gRap1, rBinsErr1, gRapErr1);

    TGraphErrors *rCent0 = new TGraphErrors(4, cBins, rCt0, cBinsErr, rCtErr0);
    TGraphErrors *rPts0 = new TGraphErrors(3, ptBins, rPt0, ptBinsErr, rPtErr0);
    TGraphErrors *rRaps0 = new TGraphErrors(3, rBins, rRap0, rBinsErr, rRapErr0);
    TGraphErrors *rRaps1 = new TGraphErrors(1, rBins1, rRap1, rBinsErr1, rRapErr1);

    TGraphErrors *rACent0 = new TGraphErrors(4, cBins, rACt0, cBinsErr, rACtErr0);
    TGraphErrors *rAPts0 = new TGraphErrors(3, ptBins, rAPt0, ptBinsErr, rAPtErr0);
    TGraphErrors *rARaps0 = new TGraphErrors(3, rBins, rARap0, rBinsErr, rARapErr0);
    TGraphErrors *rARaps1 = new TGraphErrors(1, rBins1, rARap1, rBinsErr1, rARapErr1);

    TGraphErrors *rBCent0 = new TGraphErrors(4, cBins, rBCt0, cBinsErr, rBCtErr0);
    TGraphErrors *rBPts0 = new TGraphErrors(3, ptBins, rBPt0, ptBinsErr, rBPtErr0);
    TGraphErrors *rBRaps0 = new TGraphErrors(3, rBins, rBRap0, rBinsErr, rBRapErr0);
    TGraphErrors *rBRaps1 = new TGraphErrors(1, rBins1, rBRap1, rBinsErr1, rBRapErr1);

    rBCent0->SetMarkerStyle(30);
    rBCent0->SetMarkerColor(kBlue+2);
    rBCent0->SetLineColor(kBlue+2);
    rBCent0->SetMarkerSize(1.8);
    rBPts0->SetMarkerStyle(30);
    rBPts0->SetMarkerColor(kBlue+2);
    rBPts0->SetLineColor(kBlue+2);
    rBPts0->SetMarkerSize(1.8);
    rBRaps0->SetMarkerStyle(30);
    rBRaps0->SetMarkerColor(kBlue+2);
    rBRaps0->SetLineColor(kBlue+2);
    rBRaps0->SetMarkerSize(1.8);
    rBRaps1->SetMarkerStyle(30);
    rBRaps1->SetMarkerColor(kBlue+2);
    rBRaps1->SetLineColor(kBlue+2);
    rBRaps1->SetMarkerSize(1.8);

    rACent0->SetMarkerStyle(25);
    rACent0->SetMarkerColor(kRed+2);
    rACent0->SetLineColor(kRed+2);
    rACent0->SetMarkerSize(1.8);
    rAPts0->SetMarkerStyle(25);
    rAPts0->SetMarkerColor(kRed+2);
    rAPts0->SetLineColor(kRed+2);
    rAPts0->SetMarkerSize(1.8);
    rARaps0->SetMarkerStyle(25);
    rARaps0->SetMarkerColor(kRed+2);
    rARaps0->SetLineColor(kRed+2);
    rARaps0->SetMarkerSize(1.8);
    rARaps1->SetMarkerStyle(25);
    rARaps1->SetMarkerColor(kRed+2);
    rARaps1->SetLineColor(kRed+2);
    rARaps1->SetMarkerSize(1.8);


    rCent0->SetMarkerStyle(20);
    rCent0->SetMarkerColor(kBlue+2);
    rCent0->SetLineColor(kBlue+2);
    rCent0->SetMarkerSize(1.8);
    rPts0->SetMarkerStyle(20);
    rPts0->SetMarkerColor(kBlue+2);
    rPts0->SetLineColor(kBlue+2);
    rPts0->SetMarkerSize(1.8);
    rRaps0->SetMarkerStyle(20);
    rRaps0->SetMarkerColor(kBlue+2);
    rRaps0->SetLineColor(kBlue+2);
    rRaps0->SetMarkerSize(1.8);
    rRaps1->SetMarkerStyle(21);
    rRaps1->SetMarkerColor(kBlue+2);
    rRaps1->SetLineColor(kBlue+2);
    rRaps1->SetMarkerSize(1.8);

    gCent0->SetMarkerStyle(24);
    gCent0->SetMarkerColor(kMagenta+2);
    gCent0->SetLineColor(kMagenta+2);
    gCent0->SetMarkerSize(1.8);
    gPts0->SetMarkerStyle(24);
    gPts0->SetMarkerColor(kMagenta+2);
    gPts0->SetLineColor(kMagenta+2);
    gPts0->SetMarkerSize(1.8);
    gRaps0->SetMarkerStyle(24);
    gRaps0->SetMarkerColor(kMagenta+2);
    gRaps0->SetLineColor(kMagenta+2);
    gRaps0->SetMarkerSize(1.8);
    gRaps1->SetMarkerStyle(25);
    gRaps1->SetMarkerColor(kRed+2);
    gRaps1->SetLineColor(kRed+2);
    gRaps1->SetMarkerSize(1.8);

    TCanvas *c1 = new TCanvas("c1","",500,500);

    //makeMultiPanelCanvas(c1,2,1,0.0,0.0,0.15,0.15,0.05); // 0.12, 0.15, 0.02


    TH1F *hPad1 = new TH1F("hPad1","hPad1;N_{part};v_{2}",10,0.0,400.0);
    TH1F *hPad2 = new TH1F("hPad2","hPad2;p_{T} (GeV/c);v_{2}",20,0.0,40.0);
    TH1F *hPad3 = new TH1F("hPad3","hPad3;|y|;v_{2}",10,0.0,2.4);

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

    hPad2->GetXaxis()->SetLabelSize(20);
    hPad2->GetXaxis()->SetLabelFont(43);
    hPad2->GetXaxis()->SetTitleSize(24);
    hPad2->GetXaxis()->SetTitleFont(43);
    hPad2->GetXaxis()->SetTitleOffset(1.3);
    hPad2->GetXaxis()->CenterTitle();

    hPad2->GetYaxis()->SetLabelSize(20);
    hPad2->GetYaxis()->SetLabelFont(43);
    hPad2->GetYaxis()->SetTitleSize(30);
    hPad2->GetYaxis()->SetTitleFont(43);
    hPad2->GetYaxis()->SetTitleOffset(1.0);
    hPad2->GetYaxis()->CenterTitle();

    hPad3->GetXaxis()->SetLabelSize(20);
    hPad3->GetXaxis()->SetLabelFont(43);
    hPad3->GetXaxis()->SetTitleSize(24);
    hPad3->GetXaxis()->SetTitleFont(43);
    hPad3->GetXaxis()->SetTitleOffset(1.3);
    hPad3->GetXaxis()->CenterTitle();

    hPad3->GetYaxis()->SetLabelSize(20);
    hPad3->GetYaxis()->SetLabelFont(43);
    hPad3->GetYaxis()->SetTitleSize(30);
    hPad3->GetYaxis()->SetTitleFont(43);
    hPad3->GetYaxis()->SetTitleOffset(1.0);
    hPad3->GetYaxis()->CenterTitle();

    TLegend *leg2 = new TLegend(0.5,0.50,0.72,0.70);
    TLegend *leg3 = new TLegend(0.5,0.45,0.72,0.70);
    TLegend *leg4 = new TLegend(0.5,0.50,0.72,0.70);
    formatTLeg(leg2);
    formatTLeg(leg3);
    formatTLeg(leg4);

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();

    hPad1->SetMaximum(0.5);
    hPad2->SetMaximum(0.5);
    hPad3->SetMaximum(0.5);
    hPad1->SetMinimum(0);
    hPad2->SetMinimum(0);
    hPad3->SetMinimum(0);
    
    leg2->AddEntry(gCent0,"Gen","PL");
    leg2->AddEntry(rCent0,"Reco Eff Corr","PL");
    leg2->AddEntry(rACent0,"Reco Eff+Res Corr","PL");
    //leg2->AddEntry(rBCent0,"Reco Eff Corr(Cent Free)","PL");
    //c1->Divide(3,1);

    // Centrality
    hPad1->Draw();
    rCent0->Draw("[P]");
    //rBCent0->Draw("[P]");
    rACent0->Draw("[P]");
    gCent0->Draw("[P]");
    leg2->Draw();

    lt1->SetTextSize(0.045);
    lt1->DrawLatex(0.2,0.88,"CMS Simulation");
    lt1->DrawLatex(0.2,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    //lt1->DrawLatex(0.2,0.76,"EP: etHFp+etHFm");
    lt1->SetTextSize(0.04);
    lt1->DrawLatex(0.2,0.70,"6.5 < p_{T} < 40 GeV/c");
    lt1->DrawLatex(0.2,0.64,"|y| < 2.4");
    lt1->DrawLatex(0.68,0.88,"Prompt J/#psi");
    lt1->DrawLatex(0.68,0.82,"Default");
    lt1->SetTextSize(0.043);
    lt1->DrawLatex(0.2,0.18,"PYTHIA+EvtGen+HYDJET");
    //lt1->SetTextSize(0.04);
    //lt1->DrawLatex(0.54,0.72,"Prompt J/#psi");
    //lt1->DrawLatex(0.7,0.30,"Prompt J/#psi");
    //leg1->Draw();

    c1->SaveAs("MCv2_Cents_0420_default.png");
    c1->SaveAs("MCv2_Cents_0420_default.pdf");
 
    leg3->AddEntry(gPts0,"Gen","PL");
    leg3->AddEntry(rPts0,"Reco Eff Corr","PL");
    leg3->AddEntry(rAPts0,"Reco Eff+Res Corr","PL");
    //leg3->AddEntry(rBPts0,"Reco Eff Corr(Cent Free)","PL");
    leg3->AddEntry(fGaus,"Gaussian","L");

    fGaus->SetLineColor(871);
    fGaus->SetLineWidth(2);
    fGaus->SetLineStyle(1);
    hPad2->Draw();
    rPts0->Draw("[P]");
    rAPts0->Draw("[P]");
    //rBPts0->Draw("[P]");
    gPts0->Draw("[P]");
    fGaus->Draw("same");

    leg3->Draw();
    
    lt1->SetTextSize(0.045);
    lt1->DrawLatex(0.2,0.88,"CMS Simulation");
    lt1->DrawLatex(0.2,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    //lt1->DrawLatex(0.2,0.76,"EP: etHFp+etHFm");
    lt1->SetTextSize(0.04);
    lt1->DrawLatex(0.2,0.70,"Cent. 10 - 60 %");
    lt1->DrawLatex(0.2,0.64,"|y| < 2.4");
    lt1->DrawLatex(0.68,0.88,"Prompt J/#psi");
    lt1->DrawLatex(0.68,0.82,"Default");
    lt1->SetTextSize(0.043);
    lt1->DrawLatex(0.4,0.18,"PYTHIA+EvtGen+HYDJET");


    c1->SaveAs("MCv2_Pts_0420_default.png");
    c1->SaveAs("MCv2_Pts_0420_default.pdf");

    hPad3->Draw();
    rRaps0->Draw("[P]");
    rARaps0->Draw("[P]");
    //rBRaps0->Draw("[P]");
    //rRaps1->Draw("[P]");
    gRaps0->Draw("[P]");
    //gRaps1->Draw("[P]");
    
    leg4->AddEntry(gRaps0,"Gen","PL");
    leg4->AddEntry(rRaps0,"Reco Eff Corr","PL");
    leg4->AddEntry(rARaps0,"Reco Eff+Res Corr","PL");
    //leg4->AddEntry(rBRaps0,"Reco Eff Corr(Cent Free)","PL");

    lt1->SetTextSize(0.045);
    lt1->DrawLatex(0.2,0.88,"CMS Simulation");
    lt1->DrawLatex(0.2,0.82,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
    //lt1->DrawLatex(0.2,0.76,"EP: etHFp+etHFm");
    lt1->SetTextSize(0.04);
    lt1->DrawLatex(0.2,0.70,"Cent. 10 - 60 %");
    //lt1->DrawLatex(0.2,0.60,"6.5 < p_{T} < 40 GeV/c");
    lt1->DrawLatex(0.68,0.88,"Prompt J/#psi");
    lt1->DrawLatex(0.68,0.82,"Default");
    lt1->SetTextSize(0.043);
    lt1->DrawLatex(0.2,0.18,"PYTHIA+EvtGen+HYDJET");

    char legs[512];
    TLegend *leg1 = new TLegend(0.19,0.24,0.41,0.34);
    formatTLeg(leg1);
    //leg1->SetHeader("0.0 < | y | < 1.2");
    leg1->AddEntry(rRaps0,"6.5 < p_{T} < 40","pl");
    leg1->AddEntry(rRaps1,"3 < p_{T} < 6.5","pl");

    //leg1->Draw();
    leg4->Draw();

    c1->SaveAs("MCv2_Raps_0420_default.png");
    c1->SaveAs("MCv2_Raps_0420_default.pdf");

   // 
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

