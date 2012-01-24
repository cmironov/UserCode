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

void makeMultiPanelCanvas(TCanvas*& canv, 
			  const Int_t columns, 
			  const Int_t rows, 
			  const Float_t leftOffset=0.,
			  const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
			  const Float_t edge=0.05);
void drawDum(float min, float max, double drawXLabel);
void CalEffErr(TGraphErrors *a, double *b);
TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2);
TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2);
void formatTGraph(TGraphErrors* a, int b, int c, int d);
void formatTCanv(TCanvas* a);
void formatTLeg(TLegend* a);
void formatTH1F(TH1* a, int b, int c, int d);
void ptInt_V2Ana();
void GetV2(TGraphErrors *a, double *b);


//__________________________________________________________________________
void ptInt_V2Ana()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gStyle->SetOptFit(0);
 
 //  TFile *f1 = new TFile("inclusive150_jpsi/InclusiveYield/ep21/xopt21_mbFixedIncl/summary/saved_histo.root");
//   const char* nameoutfile[2] = {"",
// 				"ep21_xopt21_mbFixedIncl"
//   };

//   const char* eventPlane[2] = {"",
// 			        "EP: etHF"
//   };

//  TFile *f1 = new TFile("inclusive150_jpsi/InclusiveYield/ep23/xopt21_allMuTrig/summary/saved_histo.root");
//   const char* nameoutfile[2] = {"",
// 				"ep23_xopt21_allMuTrig"
//   };

//   const char* eventPlane[2] = {"",
// 			       "EP: etHFm"
//   };


 TFile *f1 = new TFile("inclusive150_jpsi/InclusiveYield/ep22/xopt21_bit2/summary/saved_histo.root");
  const char* nameoutfile[2] = {"",
				"ep22_xopt21_bit2"
  };

  const char* eventPlane[2] = {"",
			       "EP: etHFp"
  };
  
  
  // pt integrated
  const int ncentbins = 5;
  const int nrapbins  = 2;
  const int nptbins   = 2; 
  int    cts[ncentbins]  = {0,  5,   10,   30,   60};
  double raps[nrapbins]  = {0.0, 2.4};
  double pts[nptbins]    = {6.5,40.0};

  const char* signal[4]      = {"","NSig","NPr","NNp"};
  const char* legend[4]      = {"","J/#psi","Prompt J/#psi","Non-prompt J/#psi"};
  int choseSignal            = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal   = signal[choseSignal];

  // options
  bool bDoDebug   = false;
  bool bSavePlots = true; 
  bool doCentIntegrated = true;
  
  int centrality_start = 0;
  int centrality_end   = 4; 

  int y_start = 0;
  int y_end   = 1;

  int pt_start = 0;
  int pt_end   = 1;
  
  int nPads = centrality_end - centrality_start;  


  TGraphErrors *g[10][10][10];
  double pr[10][10][10], ep[10][10][10], chi[10][10][10], ndf[10][10][10];
  double vals[4];
   
  char gTmp[512];
  char tmp0[512],tmp2[512], tmp3[512], tmp4[512], tmp5[512], tmp6[512];
  double vraps1, vraps2, vpts1, vpts2;
  int vcts1, vcts2;
  
  for(int mcent = centrality_start; mcent < centrality_end; mcent++)
    {
      vcts1 = cts[mcent]; 
      vcts2 = cts[mcent+1];
      cout<<"Centrlality bin: "<<vcts1<<"-"<<vcts2<<endl;
      
      for(int iy = y_start; iy < y_end; iy++)
	{
	  vraps1 = raps[iy];
	  vraps2 = raps[iy+1];
	  cout<<"Rapidity bin: "<<vraps1<<"-"<<vraps2<<endl;
	  // if(iy==0) {pt_start = 1;}
	  // else pt_start = 0;}
	  for(int jpt = pt_start; jpt < pt_end; jpt++)
	    {
	      cout<<" producing the TGraphs : "<<mcent<<" "<<iy<<" "<<jpt<<" "<<endl;
	      
	      sprintf(gTmp,"rap%.1f-%.1f_cent%d-%d_pT%.1f-%.1f_%s",vraps1,vraps2,vcts1,vcts2,pts[jpt],pts[jpt+1],chosenSignal);
	      cout<<"TGraph name : "<<gTmp<<endl;
	      TGraphErrors *pgTemp =  (TGraphErrors*)f1->Get(gTmp);
	      g[mcent][iy][jpt]  = pgTemp;
	      if(!g[mcent][iy][jpt]) {cout<<"Warning: No graph found !!!!"<<endl;continue;}
	      cout<<g[mcent][iy][jpt]<<endl;
	      double b[4];
	      GetV2(g[mcent][iy][jpt], b);
	      
	      pr[mcent][iy][jpt]  = b[0];
	      ep[mcent][iy][jpt]  = b[1];
	      chi[mcent][iy][jpt] = b[2];
	      ndf[mcent][iy][jpt] = b[3];
	      
	    }//pt bin loop
	}//rapidity loop
    }//centrlaity loop
  
  // #### Drawing: 
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  
  // frawing: 
  TCanvas *pc1 = new TCanvas("pc1",Form("pcV2_intPt"),0,0,1200,500);
  TH1F    *pp  = new TH1F("pp",
			  Form(";|#phi_{J/#psi} - #Psi_{RP}|;#frac{1}{N_{total J/#psi}} #frac{dN}{d#phi}"),
			  4,0,TMath::PiOver2());
  
  makeMultiPanelCanvas(pc1,nPads,1,0.0,0.0,0.2,0.15,0.02);
  int ind = 0;
  for(int mcent = centrality_start; mcent < centrality_end; mcent++)
    {
      vcts1  = cts[mcent]; 
      vcts2  = cts[mcent+1];
      cout<<"Canvas Centrality: "<<vcts1<<"-"<<vcts2<<endl;
           
     
      pc1->cd(ind+1);
      for(int ky = y_start; ky < y_end; ky++)
	{
	  vraps1 = raps[ky];
	  vraps2 = raps[ky+1];
	  cout << "Canvas Rapidity = "<< vraps1 << "\t"<<vraps2<<endl;
	
	  for(int lpt = pt_start; lpt < pt_end; lpt++)
	    {
	      cout<<"ind : "<<ind<<endl;
	      vpts1 = pts[lpt]; 
	      vpts2 = pts[lpt+1];
	   	    
	      pp->SetMaximum(1.2);
	      pp->SetMinimum(0.2);
	      pp->Draw();
	   
	      lt1->SetTextSize(0.05);
	      if(ind == 0) 
		{
		  lt1->DrawLatex(0.75,0.8,Form("|y| < %.1f",vraps2));       // rapidity
		  lt1->DrawLatex(0.25,0.2,Form("%.1f < p_{T} <%.1f", vpts1, vpts2)); 
		}
	      
	      if(!g[mcent][ky][lpt]) 
		{
		  cout<<"No graph! continued !!!!"<<endl;
		  continue;
		}
	      cout<<"#### Drawing: cent:"<<vcts1<<"-"<<vcts2<<"\t rapidity"<<vraps1<<"-"<<vraps2<<"\t pt"<<vpts1<<"-"<<vpts2<<endl;
	      g[mcent][ky][lpt]->Draw("pz");
	      if(ind==0 )
 		{
		  // drapwing v2 value and chi2
		  lt1->SetTextSize(0.04);
		  lt1->DrawLatex(0.23,0.3,Form("v_{2}: %.4f #pm %.4f (#chi^{2} / ndof: %.3f / %.0f)", pr[mcent][ky][lpt], ep[mcent][ky][lpt], chi[mcent][ky][lpt], ndf[mcent][ky][lpt]));
		  lt1->DrawLatex(0.75,0.9,Form("%d - %d %%",vcts1, vcts2)); // centrality 
		}
 	      else
 		{
 		  lt1->DrawLatex(0.05,0.3,Form("v_{2}: %.4f #pm %.4f (#chi^{2} / ndof: %.3f / %.0f)", pr[mcent][ky][lpt], ep[mcent][ky][lpt], chi[mcent][ky][lpt], ndf[mcent][ky][lpt]));
		  lt1->DrawLatex(0.75,0.9,Form("%d - %d %%",vcts1, vcts2));
 		}
	    }//pt bin loop
	}//rapidity bin loop
        ind++;
    }
 //_______ stuff to write
  pc1->cd(1);
  TLatex *tex1 = new TLatex(0.23,0.95,"CMS Preliminary");
  tex1->SetNDC();
  tex1->SetTextAlign(13);
  tex1->SetTextFont(43);
  tex1->SetTextSize(23);
  tex1->SetLineWidth(1);
  tex1->Draw();

  TLatex *tex2 = new TLatex(0.23,0.89,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(43);
  tex2->SetTextSize(18);
  tex2->SetLineWidth(2);
  tex2->Draw();

  TLatex *tex3 = new TLatex(0.23,0.83,"L_{int} = 150 #mub^{-1}");
  tex3->SetNDC();
  tex3->SetTextAlign(13);
  tex3->SetTextFont(43);
  tex3->SetTextSize(18);
  tex3->SetLineWidth(2);
  tex3->Draw();

  lt1->DrawLatex(0.25,0.7,Form("%s",legend[choseSignal]));  // what signal is
 
  //_______
  if(bSavePlots)
    {
      pc1->SaveAs(Form("%s_%s_Jpsi_%d_%d_150mub.png",chosenSignal,nameoutfile[1],pt_start,pt_end));
      pc1->SaveAs(Form("%s_%s_Jpsi_%d_%d_150mub.pdf",chosenSignal,nameoutfile[1],pt_start,pt_end));
      
      pc1->Clear();
    }
    
  // ####### SUMMARY PLOT!!!
  TCanvas *c2 = new TCanvas("c2","c2");
  //  makeMultiPanelCanvas(c2,1,1,0.0,0.0,0.2,0.15,0.02);
  char tmp7[512];
  TGraphErrors *gPtBarr[3];

  TH1F *hPad2 = new TH1F("hPad2",";N_{coll};Uncorrected v2;",1500,100,1700);
  hPad2->SetMaximum(0.12);
  hPad2->SetMinimum(-0.02);
  hPad2->Draw();
  TLegend *leg1 = new TLegend(0.25,0.75,0.5,0.95);
  leg1->SetFillColor(0);
  c2->cd(1);

  TLatex *tex4 = new TLatex(0.67,0.95,"CMS Preliminary");
  tex4->SetNDC();
  tex4->SetTextAlign(13);
  tex4->SetTextFont(43);
  tex4->SetTextSize(23);
  tex4->SetLineWidth(1);
  tex4->Draw();

  TLatex *tex5 = new TLatex(0.67,0.89,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  tex5->SetNDC();
  tex5->SetTextAlign(13);
  tex5->SetTextFont(43);
  tex5->SetTextSize(18);
  tex5->SetLineWidth(2);
  tex5->Draw();

  TLatex *tex6 = new TLatex(0.67,0.83,"L_{int} = 150 #mub^{-1}");
  tex6->SetNDC();
  tex6->SetTextAlign(13);
  tex6->SetTextFont(43);
  tex6->SetTextSize(18);
  tex6->SetLineWidth(2);
  tex6->Draw();
  lt1->DrawLatex(0.8,0.6,Form("%s",legend[choseSignal]));  // what signal is

 
  double ptBins[1]      = {16.75};
  double ptErrs[2]      = {0.0};
  //  double ncoll[3]       = {1483.89,744.34, 196.15};
  double ncoll[4]       = {1600.00, 1300.00,744.34, 196.15};
  for(int ic = centrality_end-1; ic >= centrality_start; ic--)
    {
     
      // [centrality][rapidity][pt]
      double v2PtBarr[1]    = {pr[ic][0][0]};
      double v2PtBarrErr[1] = {ep[ic][0][0]};
      double cBin[1]        = {ncoll[ic]};
      double cBinErr[1]     = {0.0};

      gPtBarr[ic] = new TGraphErrors(1, cBin, v2PtBarr, cBinErr, v2PtBarrErr);  
      gPtBarr[ic]->SetMarkerStyle(20+ic);
      gPtBarr[ic]->SetMarkerSize(1.2);
      gPtBarr[ic]->SetMarkerColor(3*ic);
      if(ic==0) gPtBarr[ic]->SetMarkerColor(kBlue+2);
       
      gPtBarr[ic]->Draw("pz");
      vcts1 = cts[ic]; vcts2 = cts[ic+1];
      leg1->AddEntry(gPtBarr[ic],Form("%d-%d %%",cts[ic],cts[ic+1]),"P");
    
      cout<<" producing the v2 plots in the Canvas : "<<ic<<" "<<endl;
      c2->Update();
    }
  leg1->Draw("same");
  lt1->DrawLatex(0.35,0.4,Form("|y| < %.1f",vraps2));       // rapidity
  lt1->DrawLatex(0.35,0.3,Form("%.1f < p_{T} <%.1f", vpts1, vpts2)); 
  lt1->DrawLatex(0.35,0.2,Form("%s",eventPlane[1]));
  c2->SaveAs(Form("%s_jpsi_%s_v2_pt.png",chosenSignal,nameoutfile[1]));
  c2->SaveAs(Form("%s_jpsi_%s_v2_pt.pdf",chosenSignal,nameoutfile[1]));
}


//____________________________________________________________________________
void formatTH1F(TH1* a, int b, int c, int d)
{
  //(TH1, color, style, pt, eta, rapidity)

  a->SetLineWidth(2);
  a->SetLineStyle(c);
  a->SetMarkerSize(2);
  a->SetLineColor(b);
  a->SetMarkerColor(b);
  a->GetYaxis()->SetTitle("Single #mu Efficiency");
  if(d == 1)
    {	
      a->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/c)"); 
      //a->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
    }
  else if(d == 2)
    {
      a->GetXaxis()->SetTitle("#eta^{#mu}"); 
    }
  else if(d == 3)
    {
      a->GetXaxis()->SetTitle("rapidity^{#mu}"); 
    }
  else if(d == 4)
    {
      a->GetXaxis()->SetTitle("Centrality");
    }
  else if(d == 5)
    {
      a->GetXaxis()->SetTitle("Centrality (%)");
    }
  a->GetXaxis()->CenterTitle(true);
   
}     


//__________________________________________________________________________
void formatTLeg(TLegend* a)
{
  a->SetFillColor(0); 
  a->SetBorderSize(0);
  a->SetTextSize(0.03);
  a->SetTextFont(63);
}


//__________________________________________________________________________
void formatTCanv(TCanvas* a)
{
  a->SetBorderSize(2);
  a->SetFrameFillColor(0);
  a->cd();
  a->SetGrid(1);
  a->SetTickx();
  a->SetTicky();
}


//__________________________________________________________________________
void formatTGraph(TGraphErrors* a, int b, int c, int d)
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
  if(d == 1)
    {	
      a->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
    }
  else if(d == 2)
    {
      a->GetXaxis()->SetTitle("eta"); 
    }
  else if(d == 3)
    {
      a->GetXaxis()->SetTitle("rapidity"); 
    }
  a->GetYaxis()->SetRangeUser(0,1);
  a->GetXaxis()->SetRangeUser(0,10);
  a->GetYaxis()->SetLabelSize(0.05);
  a->GetYaxis()->SetTitleSize(0.05);
  a->GetYaxis()->SetTitleOffset(1.3);
  
  
}


//__________________________________________________________________________
TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2)
{
  TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
  gEfficiency->BayesDivide(h1,h2);
  return gEfficiency;
}


//__________________________________________________________________________
TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2)
{
  
  h1->Sumw2();
  h2->Sumw2();
  
  TGraphAsymmErrors *result = calcEff(h1,h2);
  return result;
}


//__________________________________________________________________________
void CalEffErr(TGraphErrors *a, double *b)
{
  const int nbins = 100;
  double x[nbins], y[nbins];
  double sum = 0, errHighSum = 0, errLowSum = 0, sqSumHigh = 0, sqSumLow = 0;
  //double b[3] = 0;
  
  int nBins = a->GetN();
  for(int i=0;i<a->GetN();i++)
    {
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


//__________________________________________________________________________
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) 
{
   if (canv==0) 
     {
       Error("makeMultiPanelCanvas","Got null canvas.");
       return;
     }
   canv->Clear();
   
   TPad* pad[columns][rows];
   
   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth   = (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
					  (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight  = (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
					    (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0]            = leftOffset;
   Xup[0]             = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1]     = 1;
   Xlow[columns-1]    = 1.0-PadWidth/(1.0-edge);
   
   Yup[0]             = 1;
   Ylow[0]            = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1]       = bottomOffset;
   Yup[rows-1]        = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) 
     {
       Xlow[i] = Xup[0] + (i-1)*PadWidth;
       Xup[i]  = Xup[0] + (i)*PadWidth;
     }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) 
     {
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i]  = Yup[rows-1] + (ct+1)*PadHeight;
       ct++;
     }
   
   TString padName;
   for(Int_t i=0;i<columns;i++) 
     {
       for(Int_t j=0;j<rows;j++) 
	 {
	   canv->cd();
	   padName   = Form("p_%d_%d",i,j);
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


//__________________________________________________________________________
void drawDum(float min, float max, double drawXLabel)
{
  
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


//__________________________________________________________________________
void GetV2(TGraphErrors *a, double *b)
{
  // fit to extract v2 values
// Function determination to fit v2
//TF1 *v2Fit1 = new TF1("v2Fit1","[0]*(1+2*[1]*TMath::Cos(2.0*x))",-TMath::PiOver2(),TMath::PiOver2()); 
  
  TF1 *v2Fit1 = new TF1("v2Fit1","(2/TMath::Pi())*(1+2*[0]*TMath::Cos(2.0*x))",0,TMath::PiOver2()); 
 
  v2Fit1->SetLineStyle(1);
  v2Fit1->SetLineColor(kOrange+2);

  v2Fit1->SetLineStyle(1);
  v2Fit1->SetLineColor(kOrange+2);
  
  a->SetMarkerStyle(20);
  a->SetMarkerColor(kBlack); 
  a->SetMarkerSize(1.2);
  a->Fit(v2Fit1,"qrm");
  
  b[0] = v2Fit1->GetParameter(0);
  b[1] = v2Fit1->GetParError(0);
  b[2] = v2Fit1->GetChisquare();
  b[3] = v2Fit1->GetNDF();
  
}
//__________________________________________________________________________
