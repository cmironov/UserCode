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
void Plot_JpsiV2_a1_EPCorr(const char* inDirName = "./") //pnp //6bin
{
  // gROOT->Macro("./rootlogon.C");
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C");
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(5);

  
  bool bSavePlots      = true; 
  bool bIncludeSystem  = true;
  bool doDebug         = false;
  bool bExtraStudy     = false;

  bool bDoPNp = false;

  bool bDo6Bin = false;

  int signal_start     = 0;// sgn, bkg, pr, npr
  int signal_end       = 1;


  int nFits  = 3; // number of fits systm; need a counting for calcualting the rms
  int nDecay = 2; // cowboys-sailor; need a counter that it's not hard-coded
  //!!!!!!!!!!!! ############## DO NOT CHANGE THE ORDER IN THE  prefixarr[] aARRAY!!!!!!!! DO NOT CHANGE THE ORDER!!!!!!
  //**************************** modified errors
  // %%%%%%%% inclusive
   const int nPrefix = 9;
  // first in this case is always the nominal case 
   const char *prefixarr[nPrefix] = {"default_bit1","default_sailor","default_cowboy","default_constrained","default_polFunct","default_signalCB3WN","noFlat_bit1","zVtxLT10_bit1","default_bit1_weight"};
  
 // %%%%%%%% prompt - NPr
 //  const int nPrefix = 11;
//   const char *prefixarr[nPrefix] = {"default_bit1","default_sailor","default_cowboy","default_constrained","default_polFunct","default_signalCB3WN","default_bit1_1GaussResol","default_bit1_ResolFixToPRMC","noFlat_bit1","zVtxLT10_bit1","default_bit1_weight"};
 //**************************** extra Studies
  // bExtraStudy
 //  const int nPrefix = 19;
//   const char *prefixarr[nPrefix] = {"default_bit1","default_sailor","default_cowboy","default_bit1_weight","default_cowboy_weight","default_sailor_weight","nMuValHits12_bit1","nMuValHits12_cowboy","nMuValHits12_sailor","singleMuLTeta1.2_bit1","singleMuLTeta1.2_cowboy","singleMuLTeta1.2_sailor","zVtxLT10_bit1","zVtxLT10_cowboy","zVtxLT10_sailor","autoCorr_bit1","segMatchGT1_sailor","segMatchGT1_cowboy","segMatchGT1_bit1"};

  // bExtraStudy && bDoPNp
  // const int nPrefix = 13;
  //  const char *prefixarr[nPrefix] = {"default_bit1","default_sailor","default_cowboy","default_bit1_weight","default_cowboy_weight","default_sailor_weight","singleMuLTeta1.2_bit1","singleMuLTeta1.2_cowboy","singleMuLTeta1.2_sailor","zVtxLT10_bit1","zVtxLT10_cowboy","zVtxLT10_sailor","autoCorr_bit1","segMatchGT1_sailor","segMatchGT1_cowboy","segMatchGT1_bit1"};

   // options
  int prefix_start     = 0; // which setting for v2
  int prefix_end       = nPrefix;
  int centrality_start = 0;
  int centrality_end   = 4; 
  int y_start  = 0;
  int y_end    = 1;
  int pt_start = 0;
  int pt_end   = 1;

  int nPads = centrality_end - centrality_start;

  const char* signal[4]      = {"NSig","NBkg","NPr","NNp"};
  const char* legend[4]      = {"Inclusive J/#psi","Background","Prompt J/#psi","Non-prompt J/#psi"};

  ofstream output;
  if(bDoPNp)
    {
      if(bExtraStudy)
	{
	  gSystem->mkdir("./extrastud",kTRUE);
	  output.open("./extrastud/a1_v2_Result_pnp.txt");
	}
      else
	output.open("./a1_v2_Result_pnp.txt");
    }
  else
    {
      if(bExtraStudy)
	{
	  gSystem->mkdir("./extrastud",kTRUE);
	  output.open("./extrastud/a1_v2_Result.txt");
	}
      else
	output.open("./a1_v2_Result.txt");
    }
  if(bDo6Bin) output.open("./a1_v2_Result_6bin.txt");
   
  if(!output.is_open()) { cout << "cannot open a1_v2_Result.txt. Exit\n"; return ;}
 
  const int ncentbins = 4; const int cts[ncentbins+1]    = {0, 10, 20, 30, 60};
  const int nrapbins  = 1; const double raps[nrapbins+1] = {0.0, 2.4};
  const int nptbins   = 1; const double pts[nptbins+1]   = {6.5, 40.0};

  double ncoll[4]     = {355.4, 261.4178, 187.1470, 89.9};

  

  // 1st column: Different fit method or datasets (prefixarr contains all set)
  // 2nd column: [0] etHFm, [1] etHFp, [2] etHF
  // 3rd column: [0] inclusive yields, [1] bkg, [2] Prompt, [3]Non-prompt

 // Get Event Plane correction number and apply it to uncorrected v2 --- stat uncertainties only
  double corrEPV2[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};

 
  // stat only
  double v2[nPrefix][4][ncentbins][nrapbins][nptbins]    = {{{{{0.0}}}}};
  double v2Err[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  double chi[nPrefix][4][ncentbins][nrapbins][nptbins]      = {{{{{0.0}}}}};
  double ndf[nPrefix][4][ncentbins][nrapbins][nptbins]      = {{{{{0.0}}}}};

  // stat + res corrections
  double finalV2[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  double finalV2Err[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};

  // stat+syst+res included
  double v2_final[4][ncentbins][nrapbins][nptbins]    = {{{{0.0}}}};
  double v2Err_final[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
  double chi_final[4][ncentbins][nrapbins][nptbins]   = {{{{0.0}}}};
  double ndf_final[4][ncentbins][nrapbins][nptbins]   = {{{{0.0}}}};

  // stat+syst
  double v2_final_sansRes[4][ncentbins][nrapbins][nptbins]    = {{{{0.0}}}};
  double v2Err_final_sansRes[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};

  // stat only
  double v2Err_statOnly[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};

  // syst only
  double v2Err_systOnly[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
 
  // each yieldType, sgn type, centrality, rapidity, pt bin
  TGraphErrors *g[nPrefix][4][ncentbins][nrapbins][nptbins];
  TGraphErrors *pgFinal[4][ncentbins][nrapbins][nptbins];

  // some drawing stuff
  TH1F *pp  = new TH1F("pp", Form(";| #phi_{J/#psi} - #Psi_{EP} | (rad);#frac{1}{N_{total J/#psi}} #frac{dN}{d#phi} (rad^{-1})"),4,0,TMath::PiOver2());
	  
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

  TLatex *lt1 = new TLatex(); lt1->SetNDC();
	  
  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
      char eventPlane[512];
      TFile *f1 = new TFile(Form("%s/%s/summary/saved_histo.root",inDirName,prefixarr[prefix]));
     
      for(int choseSignal = signal_start; choseSignal<signal_end; choseSignal++)
	{
	  const char* chosenSignal = signal[choseSignal];
	  
	  
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
		      g[prefix][choseSignal][mcent][iy][jpt]  = pgTemp;
		      if(!g[prefix][choseSignal][mcent][iy][jpt]) {cout<<"Warning: No graph found !!!!"<<endl;continue;}
		      cout<<g[prefix][choseSignal][mcent][iy][jpt]<<endl;
		      
		      double c[4] = {0.0, 0.0, 0.0, 0.0};
		      
		      GetV2(g[prefix][choseSignal][mcent][iy][jpt], c);
		      
		      // this gives the v2 with statistical uncertainties only
		      v2[prefix][choseSignal][mcent][iy][jpt]     = c[0];
		      v2Err[prefix][choseSignal][mcent][iy][jpt]  = c[1];
		      chi[prefix][choseSignal][mcent][iy][jpt]    = c[2];
		      ndf[prefix][choseSignal][mcent][iy][jpt]    = c[3];

		      if(doDebug && prefix==0) cout<<"########### Stat only, cent bin "<<mcent<<"\t v2= "<< v2[prefix][choseSignal][mcent][iy][jpt]<<"\t +- "<<v2Err[prefix][choseSignal][mcent][iy][jpt] <<endl;
		      
		    }//pt bin loop
		}//rapidity loop
	    }//centrlaity loop
	  
	  // #### Drawing: 
	 
	  // Drawing jsut the v2 that has stat uncertainties with it
	  TCanvas *pc1 = new TCanvas("pc1",Form("pcV2_intPt"),0,0,1200,350);
	  
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
		      
		      pp->SetMaximum(1.);
		      pp->SetMinimum(0.3);
		      pp->Draw();
		      
		      lt1->SetTextSize(0.05);
		      if(ind == 2) 
			{
			  //lt1->DrawLatex(0.25,0.2,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			}
		      
		      if(!g[prefix][choseSignal][mcent][ky][lpt]) 
			{
			  cout<<"No graph! continued !!!!"<<endl;
			  continue;
			}
		      cout<<"#### Drawing: cent:"<<vcts1<<"-"<<vcts2<<"\t rapidity"<<vraps1<<"-"<<vraps2<<"\t pt"<<vpts1<<"-"<<vpts2<<endl;
		      g[prefix][choseSignal][mcent][ky][lpt]->Draw("pz");
		      if(ind==0 )
			{
			  // drapwing v2 value and chi2
			  lt1->SetTextSize(0.045);
			  lt1->DrawLatex(0.24,0.23,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[prefix][choseSignal][mcent][ky][lpt],
							v2Err[prefix][choseSignal][mcent][ky][lpt],
							chi[prefix][choseSignal][mcent][ky][lpt],
							ndf[prefix][choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.055);
			  //lt1->DrawLatex(0.78,0.89,Form("%d - %d %%",vcts1, vcts2)); // centrality 
			  lt1->DrawLatex(0.24,0.3,Form("Cent. %d - %d %%",vcts1, vcts2)); // centrality 
			}
		      else
			{
			  lt1->SetTextSize(0.05);
			  lt1->DrawLatex(0.05,0.23,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[prefix][choseSignal][mcent][ky][lpt],
							v2Err[prefix][choseSignal][mcent][ky][lpt],
							chi[prefix][choseSignal][mcent][ky][lpt],
							ndf[prefix][choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.06);
			  if(ind == 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  if(ind > 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  
			  if(ind == 3)
			    {
			      lt1->SetTextSize(0.06);
			      lt1->DrawLatex(0.05,0.75,Form("|y| < %.1f",vraps2));       // rapidity
			      lt1->DrawLatex(0.05,0.68,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			    }
			}
		    }//pt bin loop
		}//rapidity bin loop
	      ind++;
	    }//centrlaity
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

	  pc1->Update();
	  //_______
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/centdependence",kTRUE);
	      pc1->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.png",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      pc1->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_elements.pdf",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	     
	      pc1->Clear();
	    }
	  
	  // ####### SUMMARY PLOT for stat only!!! (UnCorrected for EP resolution)
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
	      double v2PtBarr[1]    = {v2[prefix][choseSignal][ic][0][0]};
	      double v2PtBarrErr[1] = {v2Err[prefix][choseSignal][ic][0][0]};
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
	  cout<<"|  Cent.   |  v2_raw_statError  |  error  |"<<endl;
	  cout<<"|  0-10%   |  "<<v2[prefix][choseSignal][0][0][0]<<"  |  "<<v2Err[prefix][choseSignal][0][0][0]<<"  |  "<<endl;
	  cout<<"|  10-20%  |  "<<v2[prefix][choseSignal][1][0][0]<<"  |  "<<v2Err[prefix][choseSignal][1][0][0]<<"  |  "<<endl;
	  cout<<"|  20-30%  |  "<<v2[prefix][choseSignal][2][0][0]<<"  |  "<<v2Err[prefix][choseSignal][2][0][0]<<"  |  "<<endl;
	  cout<<"|  30-60%  |  "<<v2[prefix][choseSignal][3][0][0]<<"  |  "<<v2Err[prefix][choseSignal][3][0][0]<<"  |  "<<endl;
	  cout<<""<<endl;
	 	  
	  //leg1->Draw("same");
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  lt1->SetTextSize(0.038);

	  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",vraps2));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
	  	 
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/centdependence",kTRUE);
	      c2->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.png",prefixarr[prefix],chosenSignal,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	      c2->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Uncorr.pdf",prefixarr[prefix],chosenSignal,vraps1,vraps2,pts[pt_start],pts[pt_end]));
	    }
	 
	}//choseSignal end of loop for signal/bkg/prompt/non-prompt 
    }//prefix loop (different input files)
  
  //__________________________________________________________________________________________
 
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
  
  TFile *rootoutput;
  if(bDoPNp)
    {
      if(!bExtraStudy) rootoutput = new TFile("a1_corrV2_pnp.root","recreate"); 
      else rootoutput = new TFile("./extrastud/a1_corrV2_pnp.root","recreate");
    }
  else
    {
      if(!bExtraStudy) rootoutput = new TFile("a1_corrV2.root","recreate"); 
      else rootoutput = new TFile("./extrastud/a1_corrV2.root","recreate");
    }
  if(bDo6Bin) rootoutput = new TFile("a1_corrV2_6bin.root","recreate"); 

  if (!rootoutput->IsOpen()) { cout << "cannot open result root file. exit.\n"; return;}

  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
      for (int choseSignal= signal_start; choseSignal < signal_end; choseSignal++) 
	{
	  for (int cent =  centrality_start; cent <  centrality_end; cent++) 
	    {
	      // Event plane correction factor
	      double corrVal_etHFpm = 0, corrErr_etHFpm = 0;
	     	      
	      getEPCorrection(0,cts[cent],cts[cent+1],&corrVal_etHFpm,&corrErr_etHFpm);
	      corrEPV2[prefix][choseSignal][cent][0][0] = v2[prefix][choseSignal][cent][0][0]/corrVal_etHFpm;

	      double v2_etHFpm       = v2[prefix][choseSignal][cent][0][0];
	      double v2Err_etHFpm    = v2Err[prefix][choseSignal][cent][0][0];
	      double corrEPV2_etHFpm = corrEPV2[prefix][choseSignal][cent][0][0];
	      
	      finalV2[prefix][choseSignal][cent][0][0]    = corrEPV2_etHFpm;
	      finalV2Err[prefix][choseSignal][cent][0][0] = TMath::Abs(corrEPV2_etHFpm) * sqrt(pow(v2Err_etHFpm/v2_etHFpm,2) + pow(corrErr_etHFpm/corrVal_etHFpm,2));
	      
	      output<<"%%%%% Res Corrected v2_statError : " << prefixarr[prefix] << " Category : "<<signal[choseSignal]<<" %%%%%"<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Correction EP" << corrVal_etHFpm <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Correction err EP: " << corrErr_etHFpm  <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Uncorrected v2: " << v2_etHFpm <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Uncorrected v2 err: " << v2Err_etHFpm  <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  Corrected v2: " << corrEPV2_etHFpm  <<"  |  "<<endl;
	      output<<"|  cent: " << cts[cent] << "-" << cts[cent+1] << "  |  " << finalV2[prefix][choseSignal][cent][0][0]<<"  |  "<<finalV2Err[prefix][choseSignal][cent][0][0]<<"  |  "<<endl;
	      output<<endl;
	      
	    } // end of centrality bins
	  
	  // write out the tgraph errors ---- this includes stat+res corrections!!!!!!!
	  rootoutput->cd();
	  char histname[200];
	  sprintf(histname,"%s_%s_rap%.1f-%.1f_pT%.1f-%.1f",prefixarr[prefix],signal[choseSignal],raps[0],raps[1],pts[0],pts[1]);
	  double cts_bound[ncentbins]   = {0.0};
	  double finalV2Cent[ncentbins] = {0.0} , finalV2ErrCent[ncentbins] = {0.0};
	  for (int cent = centrality_end-1; cent >= centrality_start; cent--) 
	    {
	      cts_bound[cent]      = ncoll[cent];
	      finalV2Cent[cent]    = finalV2[prefix][choseSignal][cent][0][0];
	      finalV2ErrCent[cent] = finalV2Err[prefix][choseSignal][cent][0][0];
	      cout << finalV2Cent[cent] << " " << finalV2ErrCent[cent] << endl;
	    }
	  TGraphErrors hFinal(ncentbins,cts_bound,finalV2Cent,0,finalV2ErrCent);
	  hFinal.SetName(histname);
	  hFinal.Write();
	  // write out the tgraph errors ---- this includes stat+res corrections!!!!!!!


	  // ####### SUMMARY PLOT!!! (Corrected)
	  TCanvas *c22 = new TCanvas("c22","c22");
	  TGraphErrors *gPtBarrCorr[3];

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
		  
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  lt1->SetTextSize(0.038);
	  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",raps[1]));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.1f GeV/c", pts[0], pts[1])); 
	  
	  c22->Update();

	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/centdependence",kTRUE);
	      c22->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.png",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      c22->SaveAs(Form("./plots/centdependence/%s_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.pdf",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	    }
 
	} // end of choseSignal (NSig, NBkg, NPr, NNp)
    } // end of prefix (Different fit methods, datasets and etc)
  
  //_========================================================================
  // include only the first nOnlyForSystm from the list of prefixes!!!!!
  if(bIncludeSystem)
    {
      for (int choseSignal= signal_start; choseSignal < signal_end; choseSignal++) 
	{
	  for (int icent = 0; icent < ncentbins; icent++) 
	    {
	      TGraphErrors *pgTemp_nominal =  (TGraphErrors*)g[0][choseSignal][icent][0][0];// nominal value
	      double yield[4]              = {0};
	      double yield_newErr[4]       = {0};
	      double xcenter[4]            = {0};
	      double xerr[4]               = {0};

	      for (int iphi=0; iphi<4; iphi++)
		{
		  double y_nom,xbin;
		  int bin1        = pgTemp_nominal->GetPoint(iphi,xbin,y_nom);
		  double yerr_nom = pgTemp_nominal->GetErrorY(iphi);
		  double rms_fits = 0;
		  double rms_cs = 0;
		  double rms_rest = 0;
		  //		  int nFits  = 3; // number of fits systm; need a counting for calcualting the rms
		  // int nDecay = 2; // cowboys-sailor; need a counter that it's not hard-coded
		  if(choseSignal >0) nFits = 5;
		  for(int prefix=1; prefix<prefix_end; prefix++) 
		    {
		      TGraphErrors *pgTemp_systm =  (TGraphErrors*)g[prefix][choseSignal][icent][0][0];
		      double y,x2;
		      int bin2 = pgTemp_systm->GetPoint(iphi,x2,y);
		      
		      if(prefix<nDecay+1)
			rms_cs+=pow(y_nom-y,2);
		      
		      if(prefix>=nDecay+1 && prefix<nFits+1)
			rms_fits+=pow(y_nom-y,2);
			
		      if(prefix >= nDecay+nFits+1) rms_rest+=pow(y_nom-y,2);
		      if(doDebug) cout <<"prefix = "<<prefixarr[prefix]<<"\t rms_cs="<<rms_cs<<"\t "<<"\t rms_fits="<<rms_fits<<"\t rms_rest="<<rms_rest<<endl;
		    }
		  rms_cs   = rms_cs/(nDecay);
		  rms_fits = rms_fits/nFits;

		  xcenter[iphi]       = xbin;
		  yield[iphi]         = y_nom;
		  yield_newErr[iphi]  = sqrt(yerr_nom*yerr_nom + rms_cs+rms_fits+rms_rest);
		  if(doDebug)
		    cout<<" phiBin = "<<iphi<<"\t yeild= "<<yield[iphi]<<"\t oldError= "<< yerr_nom <<"\t newError= "<<yield_newErr[iphi]<<endl;
		}// phi bins
	      // define new TGraph
	      TGraphErrors *pgNew = new TGraphErrors(4,xcenter,yield,xerr,yield_newErr);
	      pgFinal[choseSignal][icent][0][0] = pgNew;
	    
	      // this gets the v2 with systm+stat errors
	      double c2[4] = {0.0, 0.0, 0.0, 0.0};
	      GetV2(pgNew, c2);
	      v2_final_sansRes[choseSignal][icent][0][0]     = c2[0];
	      v2Err_final_sansRes[choseSignal][icent][0][0]  = c2[1];
	      chi_final[choseSignal][icent][0][0]            = c2[2];
	      ndf_final[choseSignal][icent][0][0]            = c2[3];

	      double errStat2  = v2Err[0][choseSignal][icent][0][0]*v2Err[0][choseSignal][icent][0][0]; // prefix=0 -- the nominal, first value in array
	      double errFinal2 = v2Err_final_sansRes[choseSignal][icent][0][0] * v2Err_final_sansRes[choseSignal][icent][0][0];
	      // calcualte separatelly the systm only: total_woRes-stat
	      v2Err_systOnly[choseSignal][icent][0][0] = sqrt(errFinal2 - errStat2);

	      // apply resolution corrections on the new v2:
	      double resCorrection = 0, resCorrection_error = 0;
	     	      
	      getEPCorrection(0,cts[icent],cts[icent+1],&resCorrection,&resCorrection_error);
	      double resCorrectedV2       = c2[0]/resCorrection;
	      double resCorrectedV2_error = TMath::Abs(resCorrectedV2) * sqrt( pow(c2[1]/c2[0],2) + pow(resCorrection_error/resCorrection,2));
	      
	      cout<<"Bin "<<icent<<"\t resCorrection "<<resCorrection<<"\t "<<resCorrection_error<<endl;
	      cout<<"Bin "<<icent<<"\t v2 "<<resCorrectedV2<<"\t "<<resCorrectedV2_error<<endl;

	      v2_final[choseSignal][icent][0][0]     = resCorrectedV2;
	      v2Err_final[choseSignal][icent][0][0]  = resCorrectedV2_error;

	    }//centrality
	  cout << endl;
	  output << "%%%%%%%%%  " << signal[choseSignal] << "  |\n";
	
	  output << "|  Cent.   |  v2_statOnly_noRes  |  v2_statOnly_noRes error (first in the prefix list)  |\n";
	  output << "|  0-10%   |  "<<v2[0][choseSignal][0][0][0]<<"  |  "<<v2Err[0][choseSignal][0][0][0]<<"  |\n";
	  output << "|  10-20%  |  "<<v2[0][choseSignal][1][0][0]<<"  |  "<<v2Err[0][choseSignal][1][0][0]<<"  |\n";
	  output << "|  20-30%  |  "<<v2[0][choseSignal][2][0][0]<<"  |  "<<v2Err[0][choseSignal][2][0][0]<<"  |\n";
	  output << "|  30-60%  |  "<<v2[0][choseSignal][3][0][0]<<"  |  "<<v2Err[0][choseSignal][3][0][0]<<"  |\n";

	  output << "|  Cent.   |  v2_systematic only error  |\n";
	  output << "|  0-10%   |  "<<v2Err_systOnly[choseSignal][0][0][0]<<"  |\n";
	  output << "|  10-20%  |  "<<v2Err_systOnly[choseSignal][1][0][0]<<"  |\n";
	  output << "|  20-30%  |  "<<v2Err_systOnly[choseSignal][2][0][0]<<"  |\n";
	  output << "|  30-60%  |  "<<v2Err_systOnly[choseSignal][3][0][0]<<"  |\n";

	  output << "|  Cent.   |  v2_systStat_noRes  |  v2_noRes error  |\n";
	  output << "|  0-10%   |  "<<v2_final_sansRes[choseSignal][0][0][0]<<"  |  "<<v2Err_final_sansRes[choseSignal][0][0][0]<<"  |\n";
	  output << "|  10-20%  |  "<<v2_final_sansRes[choseSignal][1][0][0]<<"  |  "<<v2Err_final_sansRes[choseSignal][1][0][0]<<"  |\n";
	  output << "|  20-30%  |  "<<v2_final_sansRes[choseSignal][2][0][0]<<"  |  "<<v2Err_final_sansRes[choseSignal][2][0][0]<<"  |\n";
	  output << "|  30-60%  |  "<<v2_final_sansRes[choseSignal][3][0][0]<<"  |  "<<v2Err_final_sansRes[choseSignal][3][0][0]<<"  |\n";

	  output << "|  Cent.   |  v2_final  |  v2_final error  |\n";
	  output << "|  0-10%   |  "<<v2_final[choseSignal][0][0][0]<<"  |  "<<v2Err_final[choseSignal][0][0][0]<<"  |\n";
	  output << "|  10-20%  |  "<<v2_final[choseSignal][1][0][0]<<"  |  "<<v2Err_final[choseSignal][1][0][0]<<"  |\n";
	  output << "|  20-30%  |  "<<v2_final[choseSignal][2][0][0]<<"  |  "<<v2Err_final[choseSignal][2][0][0]<<"  |\n";
	  output << "|  30-60%  |  "<<v2_final[choseSignal][3][0][0]<<"  |  "<<v2Err_final[choseSignal][3][0][0]<<"  |\n";
	  output<<endl;


	}//signal type

      // let's do some plotting please
      // ****** yields distributions, with statOnly and stat+systm uncertainties
      // Drawing jsut the v2 that has stat uncertainties with it
      
      
      for(int choseSignal = signal_start; choseSignal<signal_end; choseSignal++)
	{
	  TCanvas *pc7 = new TCanvas("pc7",Form("pc_dPhiDistribStatSyst_%s",signal[choseSignal]),0,0,1200,350);
	  makeMultiPanelCanvas(pc7,nPads,1,0.0,0.0,0.2,0.15,0.02);
	  int ind = 0;
	  for(int mcent = centrality_start; mcent < centrality_end; mcent++)
	    {
	      int vcts1  = cts[mcent]; 
	      int vcts2  = cts[mcent+1];
	      cout<<"Canvas Centrality: "<<vcts1<<"-"<<vcts2<<endl;
	      
	      pc7->cd(ind+1);
	      for(int ky = y_start; ky < y_end; ky++)
		{
		  double vraps1 = raps[ky];
		  double vraps2 = raps[ky+1];
		  cout << "Canvas Rapidity = "<< vraps1 << "\t"<<vraps2<<endl;
		  
		  for(int lpt = pt_start; lpt < pt_end; lpt++)
		    {
		      cout<<"ind : "<<ind<<endl;
		      double vpts1 = pts[lpt]; 
		      double vpts2 = pts[lpt+1];
		      
		      pp->SetMaximum(1.0);
		      pp->SetMinimum(0.4);
		      pp->Draw();
		      
		      lt1->SetTextSize(0.05);
		      if(!g[0][choseSignal][mcent][ky][lpt]) 
			{
			  cout<<"No graph! continued !!!!"<<endl;
			  continue;
			}
		      cout<<"#### Drawing: cent:"<<vcts1<<"-"<<vcts2<<"\t rapidity"<<vraps1<<"-"<<vraps2<<"\t pt"<<vpts1<<"-"<<vpts2<<endl;		     
		      TGraphErrors *pgTemp2 = (TGraphErrors *)pgFinal[choseSignal][mcent][0][0];
		      pgTemp2->SetMarkerColor(2);
		      pgTemp2->SetLineColor(2);
		      pgTemp2->SetMarkerStyle(20);
		      pgTemp2->Draw("[PZ]");

		      TGraphErrors *pgTemp1 = (TGraphErrors *)g[0][choseSignal][mcent][ky][lpt];
		      pgTemp1->SetMarkerStyle(24);
		      pgTemp1->Draw("[P]");
		      if(ind==0 )
			{
			  // drapwing v2 value and chi2
			  lt1->SetTextSize(0.045);
			  lt1->DrawLatex(0.24,0.23,Form("v_{2}^{stat} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[0][choseSignal][mcent][ky][lpt],
							v2Err[0][choseSignal][mcent][ky][lpt],
							chi[0][choseSignal][mcent][ky][lpt],
							ndf[0][choseSignal][mcent][ky][lpt]));
			  lt1->DrawLatex(0.24,0.18,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2_final[choseSignal][mcent][ky][lpt],
							v2Err_final[choseSignal][mcent][ky][lpt],
							chi_final[choseSignal][mcent][ky][lpt],
							ndf_final[choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.055);
			  lt1->DrawLatex(0.24,0.3,Form("Cent. %d - %d %%",vcts1, vcts2)); // centrality 
			}
		      else
			{
			  lt1->SetTextSize(0.05);
			  lt1->DrawLatex(0.05,0.23,Form("v_{2}^{stat} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2[0][choseSignal][mcent][ky][lpt],
							v2Err[0][choseSignal][mcent][ky][lpt],
							chi[0][choseSignal][mcent][ky][lpt],
							ndf[0][choseSignal][mcent][ky][lpt]));
			  lt1->DrawLatex(0.05,0.18,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							v2_final[choseSignal][mcent][ky][lpt],
							v2Err_final[choseSignal][mcent][ky][lpt],
							chi_final[choseSignal][mcent][ky][lpt],
							ndf_final[choseSignal][mcent][ky][lpt]));
			  lt1->SetTextSize(0.06);
			  if(ind == 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  if(ind > 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
			  
			  if(ind == 3)
			    {
			      lt1->DrawLatex(0.05,0.8,Form("|y| < %.1f",vraps2));       // rapidity
			      if(ind == 3) lt1->DrawLatex(0.05,0.72,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			    }
			}
		    }//pt loop
		}//rapidity loop
	       ind++;
	    }//centrality loop

        //_______ stuff to write
	  pc7->cd(1);
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
	  
	  TLatex *tex3 = new TLatex(0.23,0.75,"L_{int} = 150 #mub^{-1}");
	  tex3->SetNDC();
	  tex3->SetTextAlign(13);
	  tex3->SetTextFont(43);
	  tex3->SetTextSize(18);
	  tex3->SetLineWidth(2);
	  tex3->Draw();
	  
	  pc7->cd(4);
	  lt1->SetTextSize(0.06);
	  lt1->DrawLatex(0.05,0.89,Form("%s",legend[choseSignal]));  // what signal is

	  //_______
	
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/final",kTRUE);
	      pc7->SaveAs(Form("./plots/final/cent_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.png",signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      pc7->SaveAs(Form("./plots/final/cent_%s_rap%.1f_%.1f_pT%.1f-%.1f_a1_Corr.pdf",signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	     
	    }
	
      	  // write out the tgraph errors
	  // v2_statOnly, v2_statSyst_sansResolution, v2_statSyst_final, v2_syst_final?
	  rootoutput->cd();
	  double cent_bound[ncentbins]    = {0.0};
	  double cent_err[ncentbins]      = {5., 5., 5.,5.};
	  double v2SystOnlyErr[ncentbins] = {0.0};
	  double v2StatOnlyErr[ncentbins] = {0.0};
	  
	  double v2finalErr[ncentbins]    = {0.0};
	  double v2final[ncentbins]       = {0.0};

	  for (int cent = centrality_end-1; cent >= centrality_start; cent--) 
	    {
	      cent_bound[cent]  = ncoll[cent];
	      v2final[cent]     = v2_final[choseSignal][cent][0][0];    // stat+syst+res
	     
	      
	      v2StatOnlyErr[cent]  = v2Err[0][choseSignal][cent][0][0];          // statistical only, prefix=0, nominal, first entry in the array
	      v2SystOnlyErr[cent]  = v2Err_systOnly[choseSignal][cent][0][0]; // systematic errors only
	      v2finalErr[cent]     = v2Err_final[choseSignal][cent][0][0];    // stat+syst+res

	      cout << "Signal:"<< signal[choseSignal]<<"\t"<<v2final[cent] << " Last chance to screw up!!!!! " << v2finalErr[cent] << endl;
	    }
	
	  TGraphErrors pgWrite(ncentbins,cent_bound,v2final,cent_err,v2finalErr);
	  pgWrite.Write(Form("final_centDependence_%s_rap%.1f-%.1f_pT%.1f-%.1f",signal[choseSignal],raps[0],raps[1],pts[0],pts[1]));

	  TGraphErrors pgWrite_stat(ncentbins,cent_bound,v2final,cent_err,v2StatOnlyErr);
	  pgWrite_stat.Write(Form("final_centDependence_statErr_%s_rap%.1f-%.1f_pT%.1f-%.1f",signal[choseSignal],raps[0],raps[1],pts[0],pts[1]));

	  TGraphErrors pgWrite_syst(ncentbins,cent_bound,v2final,cent_err,v2SystOnlyErr);
	  pgWrite_syst.Write(Form("final_centDependence_systErr_%s_rap%.1f-%.1f_pT%.1f-%.1f",signal[choseSignal],raps[0],raps[1],pts[0],pts[1]));

	  if(doDebug)
	    {
	      double x,y;
	      cout<<pgWrite.GetPoint(1,x,y);
	      cout<<"Full"<<x<<"\t "<<y<<endl;

	      cout<<pgWrite_stat.GetPoint(1,x,y);
	      cout<<"Stat"<<x<<"\t "<<y<<endl;
	      
	      cout<<pgWrite_syst.GetPoint(1,x,y);
	      cout<<"Syst"<<x<<"\t "<<y<<endl;

	    }
	  //----------------- done writing

	}//chose signal
    } // include systm

  rootoutput->Close();
  output.close();

  //  gApplication->Terminate();
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


//__________________________________________________________________________________________
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


//__________________________________________________________________________________________
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
void getEPCorrection(int epType, int centLow, int centHigh, double *corrVal, double *corrErr) 
{
  // retrieve the right resolution correction factor
  const int nBins =5;
  string centBins[nBins] = {"0-10",
			    "10-20",
			    "20-30",
			    "30-60",
			    "10-60"};

  double etHFpm[nBins*2] = {0.622758,	0.00139869,
			    0.816579,	0.00178952,
			    0.841823,	0.00251005,
			    0.781244,	0.00292897,
			    0.816756,	0.00130408
  };

  map<string, pair<double,double> > Corr_etHFpm;
  map<string, pair<double,double> >::iterator it;

  for (int i=0; i<nBins*2; i=i+2) 
    {
      Corr_etHFpm[centBins[i/2]] = make_pair(etHFpm[i],etHFpm[i+1]);
    }

  char binName[20] = {0};
  sprintf(binName,"%d-%d",centLow,centHigh);
  if (epType == 0 ||epType == 1 ) {  //etHFm or etHFp
    it = Corr_etHFpm.find(binName);
  } 

/*  cout << "Correction factor bin: " << it->first << endl;
  cout << "Correction factor value: " << it->second.first << endl;
  cout << "Correction factor value: " << it->second.second << endl;*/

  *corrVal = it->second.first;
  *corrErr = it->second.second;

}

