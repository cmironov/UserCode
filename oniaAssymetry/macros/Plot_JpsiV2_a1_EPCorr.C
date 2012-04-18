#if !defined(__CINT__) || defined(__MAKECINT__)

#include <fstream>
#include <Riostream.h>
#include <iostream>
#include <map>
#include <utility>

#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TApplication.h>
#include <TInterpreter.h>

#endif


void makeMultiPanelCanvas(TCanvas*& canv, 
			  const Int_t columns, 
			  const Int_t rows, 
			  const Float_t leftOffset=0.,
			  const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
			  const Float_t edge=0.05);
void GetV2(TGraphErrors *a, double *b);
void drawDum(float min, float max, double drawXLabel);
void CalEffErr(TGraphErrors *a, double *b);
TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2);
TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2);
void formatTGraph(TGraphErrors* a, int b, int c, int d);
void formatTCanv(TCanvas* a);
void formatTLeg(TLegend* a);
void formatTH1F(TH1* a, int b, int c, int d);
void TGetPoints(TGraphErrors *a, double *b, double *c);
void getEPCorrection(int epType, int centLow, int centHigh, double *corrVal, double *corrErr) ;

//__________________________________________________________________________
void Plot_JpsiV2_a1_EPCorr(const char* inDirName = "./")
{
  gROOT->Macro("./rootlogon.C");
  //gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C");
  gStyle->SetOptFit(0);
  
 
//  const int nPrefix = 10;
//  const char *prefixarr[nPrefix] = {"nominal", "polFunct", "constrained", "signalCB3WN", "cowboy", "sailor", "bit1", "noFlat", "zVtxLT10", "autoCorr"};
  const int nPrefix = 3;
  const char *prefixarr[nPrefix] = {"3DEff_nominal", "3DEff_cowboy", "3DEff_sailor"};

  const char* signal[4]      = {"NSig","NBkg","NPr","NNp"};
  const char* legend[4]      = {"Inclusive J/#psi","Background","Prompt J/#psi","Non-prompt J/#psi"};

  ofstream output("./a1_v2_Result.txt");
  if(!output.is_open()) { cout << "cannot open a1_v2_Result.txt. Exit\n"; return ;}
 
  const int ncentbins = 4; const int cts[ncentbins+1]    = {0, 10, 20, 30, 60};
  const int nrapbins  = 1; const double raps[nrapbins+1] = {0.0, 2.4};
  const int nptbins   = 1; const double pts[nptbins+1]   = {6.5, 40.0};

  double ncoll[4]       = {355.4, 261.4178, 187.1470, 89.9};

   // options
  bool bSavePlots = true; 

  int prefix_start     = 0; // which setting for v2
  int prefix_end       = nPrefix;
  int signal_start     = 0;// sgn, bkg, pr, npr
  int signal_end       = 1;
  int centrality_start = 0;
  int centrality_end   = 4; 
  int y_start  = 0;
  int y_end    = 1;
  int pt_start = 0;
  int pt_end   = 1;

  double ptBins[1]      = {16.75};
  double ptErrs[2]      = {0.0};
  int nPads = centrality_end - centrality_start;
  
  // 1st column: Different fit method or datasets (prefixarr contains all set)
  // 2nd column: [0] etHFm, [1] etHFp, [2] etHF
  // 3rd column: [0] inclusive yields, [1] bkg, [2] Prompt, [3]Non-prompt
  double v2[nPrefix][3][4][ncentbins][nrapbins][nptbins] = {{{{{{0.0}}}}}};
  double v2Err[nPrefix][3][4][ncentbins][nrapbins][nptbins] = {{{{{{0.0}}}}}};
  double chi[nPrefix][3][4][ncentbins][nrapbins][nptbins] = {{{{{{0.0}}}}}};
  double ndf[nPrefix][3][4][ncentbins][nrapbins][nptbins] = {{{{{{0.0}}}}}};

  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
      for(int iCat = 0; iCat < 2; iCat++)
	{ // etHFp and etHFm

	  TFile *f1;
	  char nameoutfile[512];
	  char eventPlane[512];
	  if(iCat == 0) {
	    f1 = new TFile(Form("%s/ep23_%s/summary/saved_histo.root",inDirName,prefixarr[prefix]));
	    sprintf(nameoutfile, "etHFm");
	    sprintf(eventPlane, "#eta_{J/#psi} > 0");
	  }
	  if(iCat == 1) {
	    f1 = new TFile(Form("%s/ep22_%s/summary/saved_histo.root",inDirName,prefixarr[prefix]));
	    sprintf(nameoutfile, "etHFp");
	    sprintf(eventPlane, "#eta_{J/#psi} < 0");
	  }
	  if(iCat == 2) 
	    {
	      f1 = new TFile("etHF/summary/saved_histo.root");
	      sprintf(nameoutfile, "etHF");
	      sprintf(eventPlane, "EP: etHF");
	    }
	  if(!f1) continue;

	  // Draw and get v2 plots signal/bkg/Prompt/Non-prompt
      for(int choseSignal = signal_start; choseSignal<signal_end; choseSignal++)
	{
	  const char* chosenSignal = signal[choseSignal];
	  
	  TGraphErrors *g[ncentbins][nrapbins][nptbins];
	  char gTmp[512];
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
		  
		  for(int jpt = pt_start; jpt < pt_end; jpt++)
		    {
		      cout<<" producing the TGraphs : "<<mcent<<" "<<iy<<" "<<jpt<<" "<<endl;
		      
		      sprintf(gTmp,"rap%.1f-%.1f_cent%d-%d_pT%.1f-%.1f_%s",vraps1,vraps2,vcts1,vcts2,pts[jpt],pts[jpt+1],chosenSignal);
		      cout<<"TGraph name : "<<gTmp<<endl;
		      TGraphErrors *pgTemp =  (TGraphErrors*)f1->Get(gTmp);
		      g[mcent][iy][jpt]  = pgTemp;
		      if(!g[mcent][iy][jpt]) {cout<<"Warning: No graph found !!!!"<<endl;continue;}
		      cout<<g[mcent][iy][jpt]<<endl;
		      
		      double c[4] = {0.0, 0.0, 0.0, 0.0};
		      
		      GetV2(g[mcent][iy][jpt], c);
		      
		      v2[prefix][iCat][choseSignal][mcent][iy][jpt]  = c[0];
		      v2Err[prefix][iCat][choseSignal][mcent][iy][jpt]  = c[1];
		      chi[prefix][iCat][choseSignal][mcent][iy][jpt] = c[2];
		      ndf[prefix][iCat][choseSignal][mcent][iy][jpt] = c[3];
		      
		    }//pt bin loop
		}//rapidity loop
	    }//centrlaity loop
	  
	  // #### Drawing: 
	  TLatex *lt1 = new TLatex(); lt1->SetNDC();
	  
	  // Drawing: 
	  TCanvas *pc1 = new TCanvas("pc1",Form("pcV2_intPt"),0,0,1200,350);
	  TH1F    *pp  = new TH1F("pp", Form(";| #phi_{J/#psi} - #Psi_{EP} | (rad);#frac{1}{N_{total J/#psi}} #frac{dN}{d#phi} (rad^{-1})"),4,0,TMath::PiOver2());
	  
	  pp->GetXaxis()->SetLabelSize(20);
	  pp->GetXaxis()->SetLabelFont(43);
	  pp->GetXaxis()->SetTitleSize(18);
	  pp->GetXaxis()->SetTitleFont(43);
	  pp->GetXaxis()->SetTitleOffset(1.2);
	  pp->GetXaxis()->CenterTitle();
	  
	  pp->GetYaxis()->SetLabelSize(20);
	  pp->GetYaxis()->SetLabelFont(43);
	  pp->GetYaxis()->SetTitleSize(19);
	  pp->GetYaxis()->SetTitleFont(43);
	  pp->GetYaxis()->SetTitleOffset(1.43);
	  pp->GetYaxis()->CenterTitle();
	  
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
		      if(ind == 2) 
			{
			  //lt1->DrawLatex(0.25,0.2,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
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
			  lt1->SetTextSize(0.045);
			  lt1->DrawLatex(0.24,0.23,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[prefix][iCat][choseSignal][mcent][ky][lpt],
							v2Err[prefix][iCat][choseSignal][mcent][ky][lpt],
							chi[prefix][iCat][choseSignal][mcent][ky][lpt],
							ndf[prefix][iCat][choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.055);
			  //lt1->DrawLatex(0.78,0.89,Form("%d - %d %%",vcts1, vcts2)); // centrality 
			  lt1->DrawLatex(0.24,0.3,Form("Cent. %d - %d %%",vcts1, vcts2)); // centrality 
			}
		      else
			{
			  lt1->SetTextSize(0.05);
			  lt1->DrawLatex(0.05,0.23,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[prefix][iCat][choseSignal][mcent][ky][lpt],
							v2Err[prefix][iCat][choseSignal][mcent][ky][lpt],
							chi[prefix][iCat][choseSignal][mcent][ky][lpt],
							ndf[prefix][iCat][choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.06);
			  if(ind == 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  if(ind > 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  //if(ind == 1) lt1->DrawLatex(0.73,0.89,Form("%d - %d %%",vcts1, vcts2));
			  //if(ind > 1) lt1->DrawLatex(0.71,0.89,Form("%d - %d %%",vcts1, vcts2));
			  if(ind == 3){
			    lt1->SetTextSize(0.06);
			    if (!strcmp(nameoutfile,"etHFm")) 
			      lt1->DrawLatex(0.05,0.75,Form("%.1f < y < %.1f",vraps1,vraps2));       // rapidity
			    else if (!strcmp(nameoutfile,"etHFp")) {
			      if (vraps1 == 0)
				lt1->DrawLatex(0.05,0.75,Form("-%.1f < y < %.1f",vraps2,vraps1));       // rapidity
			      else
				lt1->DrawLatex(0.05,0.75,Form("-%.1f < y < -%.1f",vraps2,vraps1));       // rapidity
			    } else
			      lt1->DrawLatex(0.05,0.75,Form("|y| < %.1f",vraps2));       // rapidity
			  }
			  if(ind == 3) lt1->DrawLatex(0.05,0.68,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			}
		    }//pt bin loop
		}//rapidity bin loop
	      ind++;
	    }
        //_______ stuff to write
	  pc1->cd(1);
	  TLatex *tex1 = new TLatex(0.23,0.93,"CMS Preliminary");
	  tex1->SetNDC();
	  tex1->SetTextAlign(13);
	  tex1->SetTextFont(43);
	  tex1->SetTextSize(23);
	  tex1->SetLineWidth(1);
	  tex1->Draw();
	  
	  TLatex *tex2 = new TLatex(0.23,0.84,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
	  tex2->SetNDC();
	  tex2->SetTextAlign(13);
	  tex2->SetTextFont(43);
	  tex2->SetTextSize(18);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	  
	  pc1->cd(1);
	  TLatex *tex3 = new TLatex(0.23,0.75,"L_{int} = 150 #mub^{-1}");
	  tex3->SetNDC();
	  tex3->SetTextAlign(13);
	  tex3->SetTextFont(43);
	  tex3->SetTextSize(18);
	  tex3->SetLineWidth(2);
	  tex3->Draw();
	  
	  pc1->cd(4);
	  lt1->SetTextSize(0.06);
	  lt1->DrawLatex(0.05,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  
	  lt1->SetTextSize(0.06);
	  //        lt1->DrawLatex(0.05,0.82,Form("%s",eventPlane));
	  //_______
	  if(bSavePlots)
	    {
	      if(iCat == 0){
		gSystem->mkdir(Form("./plots/etHFm_%s",prefixarr[prefix]),kTRUE);
		pc1->SaveAs(Form("./plots/etHFm_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.png",prefixarr[prefix],chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
		pc1->SaveAs(Form("./plots/etHFm_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.pdf",prefixarr[prefix],chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      }
	      if(iCat == 1){
		gSystem->mkdir(Form("./plots/etHFp_%s",prefixarr[prefix]),kTRUE);
		pc1->SaveAs(Form("./plots/etHFp_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.png",prefixarr[prefix],chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
		pc1->SaveAs(Form("./plots/etHFp_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.pdf",prefixarr[prefix],chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      }
	      if(iCat == 2){
		gSystem->mkdir("./plots/etHF",kTRUE);
		pc1->SaveAs(Form("./plots/etHF/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.png",chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
		pc1->SaveAs(Form("./plots/etHF/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.pdf",chosenSignal,nameoutfile,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      }
	      
	      pc1->Clear();
	    }
	  
	  // ####### SUMMARY PLOT!!! (UnCorrected)
	  TCanvas *c2 = new TCanvas("c2","c2");
	  TGraphErrors *gPtBarr[3];
	  
	  TH1F *hPad2 = new TH1F("hPad2",";N_{part};Uncorrected v_{2};",400,0,400);
	  hPad2->GetXaxis()->SetLabelSize(20);
	  hPad2->GetXaxis()->SetLabelFont(43);
	  hPad2->GetXaxis()->SetTitleSize(30);
	  hPad2->GetXaxis()->SetTitleFont(43);
	  hPad2->GetXaxis()->SetTitleOffset(1.1);
	  hPad2->GetXaxis()->CenterTitle();
	  
	  hPad2->GetYaxis()->SetLabelSize(20);
	  hPad2->GetYaxis()->SetLabelFont(43);
	  hPad2->GetYaxis()->SetTitleSize(32);
	  hPad2->GetYaxis()->SetTitleFont(43);
	  hPad2->GetYaxis()->SetTitleOffset(1.1);
	  hPad2->GetYaxis()->CenterTitle();
	  
	  hPad2->SetMaximum(0.25);
	  hPad2->SetMinimum(-0.1);
	  
	  hPad2->Draw();
	  TLegend *leg1 = new TLegend(0.1845638,0.689021,0.4412752,0.9230769);
	  leg1->SetHeader("Centrality");
	  leg1->SetFillColor(0);
	  leg1->SetBorderSize(0);
	  leg1->SetTextSize(0.03);
	  
	  c2->cd(1);
	  
	  TLatex *tex4 = new TLatex(0.53,0.92,"CMS Preliminary");
	  tex4->SetNDC();
	  tex4->SetTextAlign(13);
	  tex4->SetTextFont(43);
	  tex4->SetTextSize(25);
	  tex4->SetLineWidth(1);
	  tex4->Draw();
	  
	  TLatex *tex5 = new TLatex(0.53,0.86,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
	  tex5->SetNDC();
	  tex5->SetTextAlign(13);
	  tex5->SetTextFont(43);
	  tex5->SetTextSize(25);
	  tex5->SetLineWidth(2);
	  tex5->Draw();
	  
	  TLatex *tex6 = new TLatex(0.53,0.80,"L_{int} = 150 #mub^{-1}");
	  tex6->SetNDC();
	  tex6->SetTextAlign(13);
	  tex6->SetTextFont(43);
	  tex6->SetTextSize(25);
	  tex6->SetLineWidth(2);
	  tex6->Draw();
	  
	  for(int ic = centrality_end-1; ic >= centrality_start; ic--)
	    {
	      double v2PtBarr[1]    = {v2[prefix][iCat][choseSignal][ic][0][0]};
	      double v2PtBarrErr[1] = {v2Err[prefix][iCat][choseSignal][ic][0][0]};
	      double cBin[1]        = {ncoll[ic]};
	      double cBinErr[1]     = {0.0};
	      int stl = 20;
	      double sz = 0.0;
	      int scol = 0;
	      stl = 20; sz = 1.8, scol = 602;
	      gPtBarr[ic] = new TGraphErrors(1, cBin, v2PtBarr, cBinErr, v2PtBarrErr);  
	      gPtBarr[ic]->SetMarkerStyle(stl);
	      gPtBarr[ic]->SetMarkerSize(sz);
	      gPtBarr[ic]->SetMarkerColor(scol);
	      
	      gPtBarr[ic]->Draw("pz");
	      leg1->AddEntry(gPtBarr[ic],Form("%d-%d %%",cts[ic],cts[ic+1]),"P");
	      
	      cout<<" producing the v2 plots in the Canvas : "<<ic<<" "<<endl;
	      c2->Update();
	    } 
	  
	  cout<<""<<endl;
	  cout<<prefixarr[prefix]<< endl;
	  cout<<"%%%%% Category : "<<signal[choseSignal]<<", "<<eventPlane<<" %%%%%"<<endl;
	  cout<<"|  Cent.   |  v2  |  error  |"<<endl;
	  cout<<"|  0-10%   |  "<<v2[prefix][iCat][choseSignal][0][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][0][0][0]<<"  |  "<<endl;
	  cout<<"|  10-20%  |  "<<v2[prefix][iCat][choseSignal][1][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][1][0][0]<<"  |  "<<endl;
	  cout<<"|  20-30%  |  "<<v2[prefix][iCat][choseSignal][2][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][2][0][0]<<"  |  "<<endl;
	  cout<<"|  30-60%  |  "<<v2[prefix][iCat][choseSignal][3][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][3][0][0]<<"  |  "<<endl;
	  cout<<""<<endl;
	  
	  output<<"%%%%% Fit : " << eventPlane << " " << prefixarr[prefix] << " Category : "<<signal[choseSignal]<<" %%%%%"<<endl;
	  output<<"|  Cent.   |  v2  |  error  |"<<endl;
	  output<<"|  0-10%   |  "<<v2[prefix][iCat][choseSignal][0][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][0][0][0]<<"  |  "<<endl;
	  output<<"|  10-20%  |  "<<v2[prefix][iCat][choseSignal][1][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][1][0][0]<<"  |  "<<endl;
	  output<<"|  20-30%  |  "<<v2[prefix][iCat][choseSignal][2][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][2][0][0]<<"  |  "<<endl;
	  output<<"|  30-60%  |  "<<v2[prefix][iCat][choseSignal][3][0][0]<<"  |  "<<v2Err[prefix][iCat][choseSignal][3][0][0]<<"  |  "<<endl;
	  output<<endl;
	  
	  //leg1->Draw("same");
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  lt1->SetTextSize(0.038);
	  if (!strcmp(nameoutfile,"etHFm")) {
	    lt1->DrawLatex(0.18,0.83,Form("%.1f < y < %.1f",vraps1,vraps2));       // rapidity
	  } else if (!strcmp(nameoutfile,"etHFp")) {
	    if (vraps1 == 0)
	      lt1->DrawLatex(0.18,0.83,Form("-%.1f < y < %.1f",vraps2,vraps1));       // rapidity
	    else
	      lt1->DrawLatex(0.18,0.83,Form("-%.1f < y < -%.1f",vraps2,vraps1));       // rapidity
	  } else
            lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",vraps2));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
	  //        lt1->DrawLatex(0.18,0.71,Form("%s",eventPlane));
	  if(iCat == 0){
	    c2->SaveAs(Form("./plots/etHFm_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.png",prefixarr[prefix],chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	    c2->SaveAs(Form("./plots/etHFm_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.pdf",prefixarr[prefix],chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	  }
	  if(iCat == 1){
	    c2->SaveAs(Form("./plots/etHFp_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.png",prefixarr[prefix],chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	    c2->SaveAs(Form("./plots/etHFp_%s/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.pdf",prefixarr[prefix],chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	  }
	  if(iCat == 2){
	    c2->SaveAs(Form("./plots/etHF/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.png",chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	    c2->SaveAs(Form("./plots/etHF/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.pdf",chosenSignal,nameoutfile,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	  }
	} // end of loop for signal/bkg/prompt/non-prompt loop
	} //end of loop for all categories
    }// end of loop for all prefixes
  
  //__________________________________________________________________________________________
  // Get Event Plane correction number and apply it to uncorrected v2
  double corrEPV2[nPrefix][3][4][ncentbins][nrapbins][nptbins] = {{{{{{0.0}}}}}};
  // Get final etHFp + etHFm combined v2 from EP corrected v2
  double finalV2[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  double finalV2Err[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  // Get RMS of all v2
  double RMSV2[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};


  TFile *rootoutput = new TFile("a1_corrV2.root","recreate");
  if (!rootoutput->IsOpen()) { cout << "cannot open result root file. exit.\n"; return;}

  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
      for (int choseSignal= signal_start; choseSignal < signal_end; choseSignal++) 
	{
	  for (int cent =  centrality_start; cent <  centrality_end; cent++) 
	    {
	      // Event plane correction factor
	      double corrVal_etHFp = 0, corrErr_etHFp = 0;
	      double corrVal_etHFm = 0, corrErr_etHFm = 0;
	      
	      getEPCorrection(0,cts[cent],cts[cent+1],&corrVal_etHFm,&corrErr_etHFm);
	      getEPCorrection(1,cts[cent],cts[cent+1],&corrVal_etHFp,&corrErr_etHFp);
	      
	      corrEPV2[prefix][0][choseSignal][cent][0][0] = v2[prefix][0][choseSignal][cent][0][0]/corrVal_etHFm;
	      corrEPV2[prefix][1][choseSignal][cent][0][0] = v2[prefix][1][choseSignal][cent][0][0]/corrVal_etHFp;
	      
	      double v2_etHFm = v2[prefix][0][choseSignal][cent][0][0];
	      double v2_etHFp = v2[prefix][1][choseSignal][cent][0][0];
	      double v2Err_etHFm = v2Err[prefix][0][choseSignal][cent][0][0];
	      double v2Err_etHFp = v2Err[prefix][1][choseSignal][cent][0][0];
	      double corrEPV2_etHFm = corrEPV2[prefix][0][choseSignal][cent][0][0];
	      double corrEPV2_etHFp = corrEPV2[prefix][1][choseSignal][cent][0][0];
	      
	      finalV2[prefix][choseSignal][cent][0][0] = ( corrEPV2_etHFm + corrEPV2_etHFp )/2;
	      finalV2Err[prefix][choseSignal][cent][0][0] = 0.5*( sqrt(
								       pow(corrEPV2_etHFm,2)*( pow(v2Err_etHFm/v2_etHFm,2) + pow(corrErr_etHFm/corrVal_etHFm,2) ) +
								       pow(corrEPV2_etHFp,2)*( pow(v2Err_etHFp/v2_etHFp,2) + pow(corrErr_etHFp/corrVal_etHFp,2) ) 
								       ) );
	      
	      output<<"%%%%% Corrected v2 : " << prefixarr[prefix] << " Category : "<<signal[choseSignal]<<" %%%%%"<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Correction etHFm: " << corrVal_etHFm <<"  |  Correction etHFp:"<< corrVal_etHFp <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Uncorrected etHFm: " << v2_etHFm <<"  |  Uncorrected etHFp:"<< v2_etHFp <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Correction err etHFm: " << corrErr_etHFm <<"  |  Correction err etHFp:"<< corrErr_etHFp <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Uncorrected err etHFm: " << v2Err_etHFm <<"  |  Uncorrected err etHFp:"<< v2Err_etHFp <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Corrected etHFm: " << corrEPV2_etHFm <<"  |  Corrected etHFp:"<< corrEPV2_etHFp <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  " << finalV2[prefix][choseSignal][cent][0][0]<<"  |  "<<finalV2Err[prefix][choseSignal][cent][0][0]<<"  |  "<<endl;
	      output<<endl;
	      
	    } // end of centrality bins
	  
	  rootoutput->cd();
	  
	  char histname[200];
	  sprintf(histname,"%s_%s_rap%.1f-%.1f_pT%.1f-%.1f",prefixarr[prefix],signal[choseSignal],raps[0],raps[1],pts[0],pts[1]);
	  double cts_bound[ncentbins] = {0.0};
	  double finalV2Cent[ncentbins] ={0.0} , finalV2ErrCent[ncentbins] = {0.0};
	  for (int cent = centrality_end-1; cent >= centrality_start; cent--) 
	    {
	      cts_bound[cent] = ncoll[cent];
	      finalV2Cent[cent] = finalV2[prefix][choseSignal][cent][0][0];
	      finalV2ErrCent[cent] = finalV2Err[prefix][choseSignal][cent][0][0];
	      cout << finalV2Cent[cent] << " " << finalV2ErrCent[cent] << endl;
	    }
	  TGraphErrors hFinal(ncentbins,cts_bound,finalV2Cent,0,finalV2ErrCent);
	  hFinal.SetName(histname);
	  hFinal.Write();
	  
	  // ####### SUMMARY PLOT!!! (Corrected)
	  TCanvas *c22 = new TCanvas("c22","c22");
	  char tmp7[512];
	  TGraphErrors *gPtBarrCorr[3];
	  
	  TH1F *hPad22 = new TH1F("hPad22",";N_{part};Corrected v_{2};",400,0,400);
	  hPad22->GetXaxis()->SetLabelSize(20);
	  hPad22->GetXaxis()->SetLabelFont(43);
	  hPad22->GetXaxis()->SetTitleSize(30);
	  hPad22->GetXaxis()->SetTitleFont(43);
	  hPad22->GetXaxis()->SetTitleOffset(1.1);
	  hPad22->GetXaxis()->CenterTitle();
	  
	  hPad22->GetYaxis()->SetLabelSize(20);
	  hPad22->GetYaxis()->SetLabelFont(43);
	  hPad22->GetYaxis()->SetTitleSize(32);
	  hPad22->GetYaxis()->SetTitleFont(43);
	  hPad22->GetYaxis()->SetTitleOffset(1.1);
	  hPad22->GetYaxis()->CenterTitle();
	  
	  hPad22->SetMaximum(0.25);
	  hPad22->SetMinimum(-0.1);
	  
	  hPad22->Draw();
	  
	  c22->cd(1);
	  
	  TLatex *tex4 = new TLatex(0.53,0.92,"CMS Preliminary");
	  tex4->SetNDC();
	  tex4->SetTextAlign(13);
	  tex4->SetTextFont(43);
	  tex4->SetTextSize(25);
	  tex4->SetLineWidth(1);
	  tex4->Draw();
	  
	  TLatex *tex5 = new TLatex(0.53,0.86,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
	  tex5->SetNDC();
	  tex5->SetTextAlign(13);
	  tex5->SetTextFont(43);
	  tex5->SetTextSize(25);
	  tex5->SetLineWidth(2);
	  tex5->Draw();
	  
	  TLatex *tex6 = new TLatex(0.53,0.80,"L_{int} = 150 #mub^{-1}");
	  tex6->SetNDC();
	  tex6->SetTextAlign(13);
	  tex6->SetTextFont(43);
	  tex6->SetTextSize(25);
	  tex6->SetLineWidth(2);
	  tex6->Draw();
	  
	  for(int ic = centrality_end-1; ic >= centrality_start; ic--)
	    {
	      double v2PtBarr[1]    = {finalV2[prefix][choseSignal][ic][0][0]};
	      double v2PtBarrErr[1] = {finalV2Err[prefix][choseSignal][ic][0][0]};
	      double cBin[1]        = {ncoll[ic]};
	      double cBinErr[1]     = {0.0};
	      int stl = 20;
	      double sz = 0.0;
	      int scol = 0;
	      stl = 20; sz = 1.8, scol = 602;
	      gPtBarrCorr[ic] = new TGraphErrors(1, cBin, v2PtBarr, cBinErr, v2PtBarrErr);  
	      gPtBarrCorr[ic]->SetMarkerStyle(stl);
	      gPtBarrCorr[ic]->SetMarkerSize(sz);
	      gPtBarrCorr[ic]->SetMarkerColor(scol);
	      
	      gPtBarrCorr[ic]->Draw("pz");
	      c22->Update();
	    }
	  
	  TLatex *lt1 = new TLatex(); lt1->SetNDC();
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  lt1->SetTextSize(0.038);
	  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",raps[1]));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.1f GeV/c", pts[0], pts[1])); 
	  
	  gSystem->mkdir("./plots/etHFp_etHFm_combined",kTRUE);
    c22->SaveAs(Form("./plots/etHFp_etHFm_combined/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.png",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
    c22->SaveAs(Form("./plots/etHFp_etHFm_combined/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.pdf",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	  
	  delete hPad22;
	  delete c22;
	  
	} // end of choseSignal (NSig, NBkg, NPr, NNp)
    } // end of prefix (Different fit methods, datasets and etc)
  
  
  /*  for (int choseSignal= 0; choseSignal < 2; choseSignal++) {
    for (int cent = 0; cent < ncentbins; cent++) {
      for(int prefix=0; prefix<nPrefix; prefix++) {
        RMSV2[choseSignal][cent][0][0] += RMSV2[choseSignal][cent][0][0] + pow(finalV2[prefix][choseSignal][cent][0][0],2);
      }
      RMSV2[choseSignal][cent][0][0] = sqrt ( 1.0/(double)nPrefix * RMSV2[choseSig"  |  "nal][cent][0][0] ); // Get RMS
    }
    cout << endl;
    cout << "|  " << signal[choseSignal] << "  |\n";
    cout << "|  Cent.   |  RMS of v2  |  v2 error  |\n";
    cout << "|  0-10%   |  "<<RMSV2[choseSignal][0][0][0]<<"  |  "<<RMSV2[choseSignal][0][0][0]<<"  |\n";
    cout << "|  10-20%  |  "<<RMSV2[choseSignal][1][0][0]<<"  |  "<<RMSV2[choseSignal][1][0][0]<<"  |\n";
    cout << "|  20-30%  |  "<<RMSV2[choseSignal][2][0][0]<<"  |  "<<RMSV2[choseSignal][2][0][0]<<"  |\n";
    cout << "|  30-60%  |  "<<RMSV2[choseSignal][3][0][0]<<"  |  "<<RMSV2[choseSignal][3][0][0]<<"  |\n";
  } */

  rootoutput->Close();
  output.close();
  gApplication->Terminate();
 
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


void GetV2(TGraphErrors *a, double *b)
{
  // fit to extract v2 values
  // Function determination to fit v2
  //TF1 *v2Fit1 = new TF1("v2Fit1","[0]*(1+2*[1]*TMath::Cos(2.0*x))",-TMath::PiOver2(),TMath::PiOver2()); 

  TF1 *v2Fit1 = new TF1("v2Fit1","(2/TMath::Pi())*(1+2*[0]*TMath::Cos(2.0*x))",0,TMath::PiOver2()); 

  v2Fit1->SetLineStyle(2);
  v2Fit1->SetLineWidth(1);
  v2Fit1->SetLineColor(kViolet+7); // kRed-3

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


void TGetPoints(TGraphErrors *a, double *b, double *c)
{
  int na = a->GetN();
  for(int i = 0; i < na; i++){
    a->GetPoint(i, b[i], c[i]);
  }

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
void getEPCorrection(int epType, int centLow, int centHigh, double *corrVal, double *corrErr) {
  const int nBins =20;
  string centBins[nBins] = {
    "0-5",
    "5-10",
    "10-15",
    "15-20",
    "20-25",
    "25-30",
    "30-35",
    "35-40",
    "40-50",
    "50-60",
    "60-70",
    "70-80",
    "80-90",
    "90-100",
    "10-60",
    "10-30",
    "30-60",
    "0-10",
    "10-20",
    "20-30"
  };

  double etHF[nBins*2] = {
0.660599,0.003897,
0.823229,0.002848,
0.883223,0.002547,
0.905866,0.002477,
0.911496,0.002443,
0.910912,0.00243,
0.901351,0.002494,
0.8816,0.002547,
0.842535,0.001955,
0.744624,0.002426,
0.574144,0.003535,
0.367572,0.006029,
0.198353,0.010795,
0.142379,0.032124,
0.872701,0.000856,
0.902874,0.001017,
0.842528,0.001339,
0.741914,0.002413,
0.8945445,0.001776,
0.911204,0.001723
  };

  double etHFp[nBins*2] = {
0.533455,0.002401,
0.723272,0.002093,
0.798937,0.002021,
0.832178,0.002045,
0.842435,0.002047,
0.837015,0.00202,
0.821521,0.002034,
0.798025,0.002033,
0.737325,0.001448,
0.609293,0.001542,
0.432927,0.001939,
0.251981,0.003116,
0.100606,0.005831,
0.05576,0.012724,
0.784591,0.000676,
0.827641,0.001017,
0.741541,0.000879,
0.628364,0.001593,
0.8155575,0.001438,
0.839725,0.001438,
  };

  double etHFm[nBins*2] = {
0.53945,0.002428,
0.722841,0.002092,
0.801483,0.002027,
0.831597,0.002044,
0.839228,0.002039,
0.835527,0.002016,
0.821654,0.002035,
0.796565,0.002029,
0.736792,0.001446,
0.609372,0.001542,
0.432557,0.001937,
0.24412,0.003019,
0.101207,0.005866,
0.046039,0.010545,
0.784027,0.000676,
0.826959,0.001016,
0.741096,0.000878,
0.631146,0.001602,
0.81654,0.001439,
0.8373775,0.001434
  };

  map<string, pair<double,double> > Corr_etHF;
  map<string, pair<double,double> > Corr_etHFp;
  map<string, pair<double,double> > Corr_etHFm;
  map<string, pair<double,double> >::iterator it;

  for (int i=0; i<nBins*2; i=i+2) {
    Corr_etHF[centBins[i/2]] = make_pair(etHF[i],etHF[i+1]);
    Corr_etHFp[centBins[i/2]] = make_pair(etHFp[i],etHFp[i+1]);
    Corr_etHFm[centBins[i/2]] = make_pair(etHFm[i],etHFm[i+1]);
  }

  char binName[20] = {0};
  sprintf(binName,"%d-%d",centLow,centHigh);
  if (epType == 0) {  //etHFm
    it = Corr_etHFm.find(binName);
  } else if (epType == 1) { //etHFp
    it = Corr_etHFp.find(binName);
  } else if (epType == 2) { //etHF
    it = Corr_etHF.find(binName);
  }

/*  cout << "Correction factor bin: " << it->first << endl;
  cout << "Correction factor value: " << it->second.first << endl;
  cout << "Correction factor value: " << it->second.second << endl;*/

  *corrVal = it->second.first;
  *corrErr = it->second.second;

}
