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
void Plot_JpsiV2_a2_EPCorr(const char* inDirName = "./pnp")
{
  gROOT->Macro("./rootlogon.C");
  //gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C");
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(5);
  bool bSavePlots      = true; 
  bool bIncludeSystem  = true;
  bool doDebug         = false;
  bool bExtraStudy     = false;

  bool bDoPNp = false;  
  bool bDo6Bin = false;
  

  int nFits  = 3; // number of fits systm; need a counting for calcualting the rms
  int nDecay = 2; // cowboys-sailor; need a counter that it's not hard-coded
  //!!!!!!!!!!!! ############## DO NOT CHANGE THE ORDER IN THE  prefixarr[] aARRAY!!!!!!!! DO NOT CHANGE THE ORDER!!!!!!
 //**************************** modified errors
  // %%%%%%%% inclusive
  //  const int nPrefix = 9;
// //   // //  //  // first in this case is always the 
  // const char *prefixarr[nPrefix] = {"default_bit1","noFlat_bit1","zVtxLT10_bit1","default_sailor","default_cowboy","default_bit1_weight","default_constrained","default_polFunct","default_signalCB3WN"};
  
  // %%%%%%%% prompt - NPr
  const int nPrefix = 11;
// // //   // first in this case is always the nominal case 
// // //   // "default_bit1_weight" -- this is 3D correction: pt, y, centrality (averaged over dPhi)
  const char *prefixarr[nPrefix] = {"default_bit1","noFlat_bit1","zVtxLT10_bit1","default_sailor","default_cowboy","default_bit1_weight","default_constrained","default_polFunct","default_signalCB3WN","default_bit1_1GaussResol","default_bit1_ResolFixToPRMC"};

 //**************************** extra Studies
//   const int nPrefix = 19;
//   const char *prefixarr[nPrefix] = {"default_bit1","default_sailor","default_cowboy","default_bit1_weight","default_bit1_4Dweight","default_cowboy_weight","default_cowboy_4Dweight","default_sailor_4Dweight","default_sailor_weight","nMuValHits12_bit1","nMuValHits12_cowboy","nMuValHits12_sailor","singleMuLTeta1.2_bit1","singleMuLTeta1.2_cowboy","singleMuLTeta1.2_sailor","zVtxLT10_bit1","zVtxLT10_cowboy","zVtxLT10_sailor","autoCorr_bit1"};



  const char* signal[4]          = {"NSig","NBkg","NPr","NNp"};
  const char* legend[4]          = {"Inclusive J/#psi","Background","Prompt J/#psi","Non-prompt J/#psi"};
 
  const int ncentbins = 1; int cts[ncentbins+1]    = {10, 60};
  const int nrapbins  = 1; double raps[nrapbins+1] = {0.0, 2.4};
  const int nptbins   = 3; double pts[nptbins+1]   = {6.5, 8.0, 10.0, 40.0};

  double pts_bound[nptbins] = {7.3, 9.0, 13.4}; // <pt> values in the AN

  int prefix_start     = 0; // which setting for v2
  int prefix_end       = nPrefix;
  int signal_start     = 2;// sgn, bkg, pr, npr
  int signal_end       = 4;
  int centrality_start  = 0;
  int centrality_end    = 1; 
  
  int y_start  = 0;
  int y_end    = 1;
  
  int pt_start = 0;
  int pt_end   = 3;
  
  int nPads    = pt_end - pt_start;  

  ofstream output;
 if(bDoPNp)
    {
      if(bExtraStudy)
	{
	  gSystem->mkdir("./extrastud",kTRUE);
	  output.open("./extrastud/a2_v2_Result_pnp.txt");
	}
      else
	output.open("./a2_v2_Result_pnp.txt");
    }
  else
    {
      if(bExtraStudy)
	{
	  gSystem->mkdir("./extrastud",kTRUE);
	  output.open("./extrastud/a2_v2_Result.txt");
	}
      else
	output.open("./a2_v2_Result.txt");
    }
  
  if(!output.is_open()) { cout << "cannot open a3_v2_Result.txt. Exit\n"; return ;}
   
 
    

  // 1st column: Different fit method or datasets (prefixarr contains all set)
  // 2nd column: [0] etHFm, [1] etHFp, [2] etHF
  // 3rd column: [0] inclusive yields, [1] bkg, [2] Prompt, [3]Non-prompt
  double v2[nPrefix][4][ncentbins][nrapbins][nptbins]    = {{{{{0.0}}}}};
  double v2Err[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  double chi[nPrefix][4][ncentbins][nrapbins][nptbins]   = {{{{{0.0}}}}};
  double ndf[nPrefix][4][ncentbins][nrapbins][nptbins]   = {{{{{0.0}}}}};

  // with systm included
  double v2_final[4][ncentbins][nrapbins][nptbins]    = {{{{0.0}}}};
  double v2Err_final[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
  double chi_final[4][ncentbins][nrapbins][nptbins]   = {{{{0.0}}}};
  double ndf_final[4][ncentbins][nrapbins][nptbins]   = {{{{0.0}}}};

  // stat+syst+res
  double v2_final_sansRes[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
  double v2Err_final_sansRes[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
 
  // stat only
  double v2Err_statOnly[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};
  // syst only
  double v2Err_systOnly[4][ncentbins][nrapbins][nptbins] = {{{{0.0}}}};

 
  // each yieldType, centrality, rapidity, pt bin
  TGraphErrors *g[nPrefix][4][ncentbins][nrapbins][nptbins];
  TGraphErrors *pgFinal[4][ncentbins][nrapbins][nptbins];
  
  // some drawing stuff
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  TH1F    *pp  = new TH1F("pp", Form(";| #phi_{J/#psi} - #Psi_{EP} | (rad);#frac{1}{N_{total J/#psi}} #frac{dN}{d#phi} (rad^{-1})"),4,0,TMath::PiOver2());
	      
  pp->GetXaxis()->SetLabelSize(20);
  pp->GetXaxis()->SetLabelFont(43);
  pp->GetXaxis()->SetTitleSize(14);
  pp->GetXaxis()->SetTitleFont(43);
  pp->GetXaxis()->SetTitleOffset(1.5);
  pp->GetXaxis()->CenterTitle();
  
  pp->GetYaxis()->SetLabelSize(20);
  pp->GetYaxis()->SetLabelFont(43);
  pp->GetYaxis()->SetTitleSize(14);
  pp->GetYaxis()->SetTitleFont(43);
  pp->GetYaxis()->SetTitleOffset(1.8);
  pp->GetYaxis()->CenterTitle();
  

  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
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
		      if(!g[prefix][choseSignal][mcent][iy][jpt]) {cout<<"@@@@@@@@@@@@@@@@ Warning: No graph found !!!!"<<endl;continue;}
		      cout<<g[prefix][choseSignal][mcent][iy][jpt]<<endl;
		      
		      double c[4] = {0.0, 0.0, 0.0, 0.0};
		      
		      GetV2(g[prefix][choseSignal][mcent][iy][jpt], c);
		      
		      v2[prefix][choseSignal][mcent][iy][jpt]     = c[0];
		      v2Err[prefix][choseSignal][mcent][iy][jpt]  = c[1];
		      chi[prefix][choseSignal][mcent][iy][jpt]    = c[2];
		      ndf[prefix][choseSignal][mcent][iy][jpt]    = c[3];

		      if(doDebug && prefix==0) cout<<"########### Stat only, cent bin "<<mcent<<"\t v2= "<< v2[prefix][choseSignal][mcent][iy][jpt]<<"\t +- "<<v2Err[prefix][choseSignal][mcent][iy][jpt] <<endl;
		      
		    }//pt bin loop
		}//rapidity loop
	    }//centrlaity loop
	      
	  // ###### Drawing: 
	  TCanvas *pc1 = new TCanvas("pc1",Form("pcV2_intPt"),0,0,850,300);
	  makeMultiPanelCanvas(pc1,nPads,1,0.0,0.0,0.2,0.15,0.02);
	  int ind = 0;
	  for(int mpt = pt_start; mpt < pt_end; mpt++)
	    {
	      cout<<"ind : "<<ind<<endl;
	      vpts1 = pts[mpt]; 
	      vpts2 = pts[mpt+1];
	      
	      pc1->cd(ind+1);
	      for(int ky = y_start; ky < y_end; ky++)
		{
		  vraps1 = raps[ky];
		  vraps2 = raps[ky+1];
		  cout << "Canvas Rapidity = "<< vraps1 << "\t"<<vraps2<<endl;
		  
		  for(int lc = centrality_start; lc < centrality_end; lc++)
		    {
		      vcts1  = cts[lc]; 
		      vcts2  = cts[lc+1];
		      cout<<"Canvas Centrality: "<<vcts1<<"-"<<vcts2<<endl;
		      
		      pp->SetMaximum(1.2);
		      pp->SetMinimum(0.2);
		      pp->Draw();
		      
		      lt1->SetTextSize(0.05);
		      if(ind == 2) 
			    {
			      lt1->DrawLatex(0.62,0.90,Form("Cent. %d - %d %%",vcts1, vcts2)); // centrality
			      lt1->DrawLatex(0.07,0.74,Form("|y| < %.1f",vraps2));       // rapidity
			    }
		      
		      if(!g[prefix][choseSignal][lc][ky][mpt]) 
			{
			  cout<<"No graph! continued !!!!"<<endl;
			  continue;
			}
		      cout<<"#### Drawing: cent:"<<vcts1<<"-"<<vcts2<<"\t rapidity"<<vraps1<<"-"<<vraps2<<"\t pt"<<vpts1<<"-"<<vpts2<<endl;
			  g[prefix][choseSignal][lc][ky][mpt]->Draw("pz");
			  if(ind==0 )
			    {
			      // drapwing v2 value and chi2
			      lt1->SetTextSize(0.04);
			      lt1->DrawLatex(0.25,0.22,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							    v2[prefix][choseSignal][lc][ky][mpt],
							    v2Err[prefix][choseSignal][lc][ky][mpt],
							    chi[prefix][choseSignal][lc][ky][mpt],
							    ndf[prefix][choseSignal][lc][ky][mpt]));
			      lt1->DrawLatex(0.25,0.3,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			    }
			  else
			    {
			      lt1->SetTextSize(0.046);
			      lt1->DrawLatex(0.07,0.22,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
							    v2[prefix][choseSignal][lc][ky][mpt],
							    v2Err[prefix][choseSignal][lc][ky][mpt],
							    chi[prefix][choseSignal][lc][ky][mpt],
							    ndf[prefix][choseSignal][lc][ky][mpt]));
			      lt1->DrawLatex(0.07,0.3,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
			    }
		    }//centrality bin loop
		    }//rapidity bin loop
	      ind++;
	    }//pt
	  //_______ stuff to write
	  pc1->cd(1);
	  TLatex *tex1 = new TLatex(0.25,0.93,"CMS Preliminary");
	  tex1->SetNDC();
	  tex1->SetTextAlign(13);
	  tex1->SetTextFont(43);
	  tex1->SetTextSize(16);
	  tex1->SetLineWidth(1);
	  tex1->Draw();
	  
	  TLatex *tex2 = new TLatex(0.25,0.85,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
	  tex2->SetNDC();
	  tex2->SetTextAlign(13);
	  tex2->SetTextFont(43);
	  tex2->SetTextSize(16);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	      
	  TLatex *tex3 = new TLatex(0.25,0.78,"L_{int} = 150 #mub^{-1}");
	  tex3->SetNDC();
	  tex3->SetTextAlign(13);
	  tex3->SetTextFont(43);
	  tex3->SetTextSize(16);
	  tex3->SetLineWidth(2);
	  tex3->Draw();
	  
	  pc1->cd(3);
	  lt1->SetTextSize(0.055);
	  lt1->DrawLatex(0.07,0.90,Form("%s",legend[choseSignal]));  // what signal is
	  	  
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/ptdependence",kTRUE);
	      pc1->SaveAs(Form("./plots/etHFpm_%s/%s_rap%.1f_%.1f_cent%d-%d_a2_elements.png",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	      pc1->SaveAs(Form("./plots/etHFpm_%s/%s_rap%.1f_%.1f_cent%d-%d_a2_elements.pdf",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	      pc1->Clear();
	    }
	  	      
	  // ####### SUMMARY PLOT!!!
	  TCanvas *c2 = new TCanvas("c2","c2");
	  TH1F *hPad2 = new TH1F("hPad2",";p_{T} (GeV/c);Uncorrected v_{2};",100,0,40);
	  hPad2->GetXaxis()->SetLabelSize(20);
	  hPad2->GetXaxis()->SetLabelFont(43);
	  hPad2->GetXaxis()->SetTitleSize(27);
	  hPad2->GetXaxis()->SetTitleFont(43);
	  hPad2->GetXaxis()->SetTitleOffset(1.2);
	  hPad2->GetXaxis()->CenterTitle();
	  
	  hPad2->GetYaxis()->SetLabelSize(20);
	  hPad2->GetYaxis()->SetLabelFont(43);
	  hPad2->GetYaxis()->SetTitleSize(32);
	  hPad2->GetYaxis()->SetTitleFont(43);
	  hPad2->GetYaxis()->SetTitleOffset(1.1);
	  hPad2->GetYaxis()->CenterTitle();
	  
	  hPad2->SetMaximum(0.25);
	  hPad2->SetMinimum(-0.10);
	  hPad2->Draw();
	  TLegend *leg1 = new TLegend(0.5100671,0.2255245,0.7667785,0.3304196);
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
	  
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  
	  TGraphErrors *gPtBarr[3];
	  double ptErrs[3]      = {0.0};
	
	  for(int ipt = pt_start; ipt < pt_end; ipt++)
	    {
	      // [centrality][rapidity][pt]
	      double v2PtBarr[1]       = {v2[prefix][choseSignal][0][0][ipt]};
	      double v2PtBarrErr[1]    = {v2Err[prefix][choseSignal][0][0][ipt]};
	      double ptbinlocal[1]     = {pts_bound[ipt]};
	      double ptbinlocalerr[1]  = {ptErrs[ipt]};
	      //cout<<"ptbinlocal[nptbins-5] : "<<ptbinlocal[nptbins-5]<<endl;
	      //cout<<"ptBins[ipt] : "<<ptBins[ipt]<<endl;
	      gPtBarr[ipt] = new TGraphErrors(1, ptbinlocal, v2PtBarr, ptbinlocalerr, v2PtBarrErr);  
	      gPtBarr[ipt]->SetMarkerStyle(20);
	      gPtBarr[ipt]->SetMarkerSize(1.8);
	      gPtBarr[ipt]->SetMarkerColor(kBlue+2);
	      
	      gPtBarr[ipt]->Draw("pz");
	      vpts1 = pts[ipt]; vpts2 = pts[ipt+1];
	      leg1->AddEntry(gPtBarr[ipt],Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2),"P");
	      
	      cout<<" producing the v2 plots in the Canvas : "<<ipt<<" "<<endl;
	      c2->Update();
	    }
	  cout<<""<<endl;
	  output<<prefixarr[prefix]<< endl;
	  output<<"%%%%% Category : "<<signal[choseSignal]<<" %%%%%"<<endl;
	  output<<"|  (pT)  |  v2  |  error  |\n";
	  output<<"|  6.5-8.0  |  "<<v2[prefix][choseSignal][0][0][0]<<"  |  "<<v2Err[prefix][choseSignal][0][0][0]<<"  |"<<endl;
	  output<<"|  8.0-10.0  |  "<<v2[prefix][choseSignal][0][0][1]<<"  |  "<<v2Err[prefix][choseSignal][0][0][1]<<"  |"<<endl;
	  output<<"|  10.0-40.0  |"<<v2[prefix][choseSignal][0][0][2]<<"  |  "<<v2Err[prefix][choseSignal][0][0][2]<<"  |"<<endl;
	 
	  output<<endl;
	  
	  //leg1->Draw("same");
	  lt1->SetTextSize(0.038);
	  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",vraps2));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("Cent. 10 - 60 %%")); 
	  
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/ptdependence",kTRUE);
	      c2->SaveAs(Form("./plots/ptdependence/%s_%s_rap%.1f_%.1f_cent%d-%d_a2_Uncorr.png",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	      c2->SaveAs(Form("./plots/ptdependence/%s_%s_rap%.1f_%.1f_cent%d-%d_a2_Uncorr.pdf",prefixarr[prefix],chosenSignal,raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	    }
	  
	} // end of loop for signal/bkg/prompt/non-prompt loop
    }// end of loop for all prefixes

  //__________________________________________________________________________________________
  // Get Event Plane correction number and apply it to uncorrected v2
  double corrEPV2[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  // Get final etHFp + etHFm combined v2 from EP corrected v2
  double finalV2[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  double finalV2Err[nPrefix][4][ncentbins][nrapbins][nptbins] = {{{{{0.0}}}}};
  
  TH1F *hPad22 = new TH1F("hPad2",";p_{T} (GeV/c);Corrected v_{2};",100,0,40);
  hPad22->GetXaxis()->SetLabelSize(20);
  hPad22->GetXaxis()->SetLabelFont(43);
  hPad22->GetXaxis()->SetTitleSize(27);
  hPad22->GetXaxis()->SetTitleFont(43);
  hPad22->GetXaxis()->SetTitleOffset(1.2);
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
      if(!bExtraStudy) rootoutput = new TFile("a2_corrV2_pnp.root","recreate"); 
      else rootoutput = new TFile("./extrastud/a2_corrV2_pnp.root","recreate");
    }
  else
    {
      if(!bExtraStudy) rootoutput = new TFile("a2_corrV2.root","recreate"); 
      else rootoutput = new TFile("./extrastud/a2_corrV2.root","recreate");
    }
  if(bDo6Bin) rootoutput = new TFile("a2_corrV2_6bin.root","recreate"); 

  for(int prefix=prefix_start; prefix<prefix_end; prefix++) 
    {
      for (int choseSignal= signal_start; choseSignal < signal_end; choseSignal++) 
	{
	  for (int ipt = pt_start; ipt < pt_end; ipt++) 
	    {
	      // Event plane correction factor
	      double corrVal_etHFpm = 0, corrErr_etHFpm = 0;
	      
	      int cent =0;  //10-60% centrality
	      getEPCorrection(0,cts[0],cts[1],&corrVal_etHFpm,&corrErr_etHFpm);
	      corrEPV2[prefix][choseSignal][0][0][ipt] = v2[prefix][choseSignal][cent][0][ipt]/corrVal_etHFpm;

	      double v2_etHFpm       = v2[prefix][choseSignal][cent][0][ipt];
	      double v2Err_etHFpm    = v2Err[prefix][choseSignal][cent][0][ipt];
	      double corrEPV2_etHFpm = corrEPV2[prefix][choseSignal][cent][0][ipt];
	      
	      finalV2[prefix][choseSignal][cent][0][ipt]    = corrEPV2_etHFpm;
	      finalV2Err[prefix][choseSignal][cent][0][ipt] = TMath::Abs(corrEPV2_etHFpm) * sqrt(pow(v2Err_etHFpm/v2_etHFpm,2) + pow(corrErr_etHFpm/corrVal_etHFpm,2));

	      output<<"%%%%% Res Corrected v2_statError : " << prefixarr[prefix] << " Category : "<<signal[choseSignal]<<" %%%%%"<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  Correction EP" << corrVal_etHFpm <<"  |  "<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  Correction err EP: " << corrErr_etHFpm  <<"  |  "<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  Uncorrected v2: " << v2_etHFpm <<"  |  "<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  Uncorrected v2 err: " << v2Err_etHFpm  <<"  |  "<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  Corrected v2: " << corrEPV2_etHFpm  <<"  |  "<<endl;
	      output<<"|  pt: " << pts[ipt] << "-" << pts[ipt+1] << "  |  " << finalV2[prefix][choseSignal][0][0][ipt]<<"  |  "<<finalV2Err[prefix][choseSignal][0][0][ipt]<<"  |  "<<endl;
	      output<<endl;
	    } // end of pt bins
	  
	  rootoutput->cd();
	  char histname[200];
	  sprintf(histname,"%s_%s_rap%.1f-%.1f_cent%d-%d",prefixarr[prefix],signal[choseSignal],raps[0],raps[1],cts[0],cts[1]);
	  double finalV2Pt[nptbins] = {0.0} , finalV2ErrPt[nptbins] = {0.0};
	  
	  for (int pt = pt_start; pt < pt_end; pt++) 
	    {
	      int cent =0;  //10-60% centrality
	      finalV2Pt[pt]    = finalV2[prefix][choseSignal][cent][0][pt];
	      finalV2ErrPt[pt] = finalV2Err[prefix][choseSignal][cent][0][pt];
	      cout << finalV2Pt[pt] << " " << finalV2ErrPt[pt] << endl;
	    }
	  TGraphErrors hFinal(nptbins,pts_bound,finalV2Pt,0,finalV2ErrPt);
	  hFinal.SetName(histname);
	  hFinal.Write();
	  
	  //####### SUMMARY PLOT!!! (Corrected)
	  TCanvas *c22 = new TCanvas("c22","c22");
	  hPad22->SetMaximum(0.25);
	  hPad22->SetMinimum(-0.10);
	  hPad22->Draw();
	  TLegend *leg1 = new TLegend(0.5100671,0.2255245,0.7667785,0.3304196);
	  leg1->SetFillColor(0);
	  leg1->SetBorderSize(0);
	  leg1->SetTextSize(0.03);
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
	  
	  lt1->SetTextSize(0.04);
	  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
	  lt1->SetTextSize(0.038);
	  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",raps[1]));       // rapidity
	  lt1->DrawLatex(0.18,0.77,Form("Cent. %d - %d %%",cts[0],cts[1]));
	  
	  TGraphErrors *gPtBarrCorr[nptbins];
 	  double ptErrs[nptbins]= {0.0};

	  for(int ipt = pt_start; ipt < pt_end; ipt++)
	    {
	      // [centrality][rapidity][pt]
	      double v2PtBarr[1]    = {finalV2[prefix][choseSignal][0][0][ipt]};
	      double v2PtBarrErr[1] = {finalV2Err[prefix][choseSignal][0][0][ipt]};
	      double ptbinlocal[1]  = {pts_bound[ipt]};
	      double ptbinlocalerr[1]  = {ptErrs[ipt]};
	      //cout<<"ptbinlocal[nptbins-5] : "<<ptbinlocal[nptbins-5]<<endl;
	      //cout<<"ptBins[ipt] : "<<ptBins[ipt]<<endl;
	      gPtBarrCorr[ipt] = new TGraphErrors(1, ptbinlocal, v2PtBarr, ptbinlocalerr, v2PtBarrErr);  
	      gPtBarrCorr[ipt]->SetMarkerStyle(20);
	      gPtBarrCorr[ipt]->SetMarkerSize(1.8);
	      gPtBarrCorr[ipt]->SetMarkerColor(kBlue+2);
	      
	      gPtBarrCorr[ipt]->Draw("pz");
	      //leg1->AddEntry(gPtBarrCorr[ipt],Form("%.1f < p_{T} < %.1f GeV/c", pts[ipt], pts[ipt+1]),"P");
	      c22->Update();
	    }
	  
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/ptdependence",kTRUE);
	      c22->SaveAs(Form("./plots/ptdependence/%s_%s_rap%.1f_%.1f_cent%d-%d_a2_Corr.png",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	      c22->SaveAs(Form("./plots/ptdependence/%s_%s_rap%.1f_%.1f_cent%d-%d_a2_Corr.pdf",prefixarr[prefix],signal[choseSignal],raps[y_start],raps[y_end],cts[centrality_start],cts[centrality_end]));
	    }

	  	  
	} // end of choseSignal (NSig, NBkg, NPr, NNp)
    } // end of prefix (Different fit methods, datasets and etc)

  //=================================================================================================
  // include only the first nOnlyForSystm from the list of prefixes!!!!!
  if(bIncludeSystem)
    {
      for (int choseSignal= signal_start; choseSignal < signal_end; choseSignal++) 
	{
	  for (int ipt = 0; ipt < nptbins; ipt++) 
	    {
	      TGraphErrors *pgTemp_nominal =  (TGraphErrors*)g[0][choseSignal][0][0][ipt]; // the nominal value
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
		      TGraphErrors *pgTemp_systm =  (TGraphErrors*)g[prefix][choseSignal][0][0][ipt];
		      double y,x2;
		      int bin2 = pgTemp_systm->GetPoint(iphi,x2,y);
		      if(prefix<nDecay+1)
			rms_cs+=pow(y_nom-y,2);
		      
		      if(prefix>=nDecay+1 && prefix<nFits+1)
			rms_fits+=pow(y_nom-y,2);
		      
		      if(prefix >= nDecay+nFits+1) rms_rest+=pow(y_nom-y,2);
		      
		      //    if(doDebug) cout <<"prefix = "<<prefixarr[prefix]<<"\t yield_nom= "<< y_nom<< "\t yield_pre = "<<y <<"\t rms_cs="<<rms_cs<<"\t "<<"\t rms_fits="<<rms_fits<<"\t rms_rest="<<rms_rest<<endl;
		      if(doDebug) cout <<"prefix = "<<prefixarr[prefix]<<"\t rms_cs="<<rms_cs<<"\t "<<"\t rms_fits="<<rms_fits<<"\t rms_rest="<<rms_rest<<endl;
		    }

		  rms_cs   = rms_cs/(nDecay);
		  rms_fits = rms_fits/nFits;
		  
		  xcenter[iphi]      = xbin;
		  yield[iphi]        = y_nom;
		  yield_newErr[iphi] = sqrt(yerr_nom*yerr_nom + rms_cs+rms_fits+rms_rest);
		}// phi bins
	      // define new TGraph
	      TGraphErrors *pgNew = new TGraphErrors(4,xcenter,yield,xerr,yield_newErr);
	      pgFinal[choseSignal][0][0][ipt] = pgNew;
	    
	      // get v2 with the new errors
	      double c2[4] = {0.0, 0.0, 0.0, 0.0};
	      GetV2(pgNew, c2);
	      v2_final_sansRes[choseSignal][0][0][ipt]     = c2[0];
	      v2Err_final_sansRes[choseSignal][0][0][ipt]  = c2[1];
	      chi_final[choseSignal][0][0][ipt]            = c2[2];
	      ndf_final[choseSignal][0][0][ipt]            = c2[3];

	      // get the systematics away from stat
	      double errStat2  = v2Err[0][choseSignal][0][0][ipt]*v2Err[0][choseSignal][0][0][ipt]; // prefix=0 -- the nominal, first value in array
	      double errFinal2 = v2Err_final_sansRes[choseSignal][0][0][ipt] * v2Err_final_sansRes[choseSignal][0][0][ipt];
	      // calcualte separatelly the systm only: total_woRes-stat
	      v2Err_systOnly[choseSignal][0][0][ipt] = sqrt(errFinal2 - errStat2);

	      // apply resolution corrections on the new v2:
	      double resCorrection = 0, resCorrection_error = 0;
	     	      
	      getEPCorrection(0,cts[0],cts[1],&resCorrection,&resCorrection_error);
	      double resCorrectedV2       = c2[0]/resCorrection;
	      double resCorrectedV2_error = TMath::Abs(resCorrectedV2) * sqrt( pow(c2[1]/c2[0],2) + pow(resCorrection_error/resCorrection,2));
	      
	      cout<<"Bin "<<ipt<<"\t resCorrection "<<resCorrection<<"\t "<<resCorrection_error<<endl;
	      cout<<"Bin "<<ipt<<"\t v2 "<<resCorrectedV2<<"\t "<<resCorrectedV2_error<<endl;

	      v2_final[choseSignal][0][0][ipt]     = resCorrectedV2;
	      v2Err_final[choseSignal][0][0][ipt]  = resCorrectedV2_error;

	    }//centrality
	  cout << endl;
	  output << "%%%%%%%%%%%  " << signal[choseSignal] << "  |\n";
	  
	  output << "|  pT.   |  v2_statOnly_noRes  |  v2_statOnly_noRes error (first in the prefix list)  |\n";
	  output << "|   6.5-8.0   |  "<<v2[0][choseSignal][0][0][0]<<"  |  "<<v2Err[0][choseSignal][0][0][0]<<"  |\n";
	  output << "|  8.0-10.0   |  "<<v2[0][choseSignal][0][0][1]<<"  |  "<<v2Err[0][choseSignal][0][0][1]<<"  |\n";
	  output << "|  10.0-40.0  |  "<<v2[0][choseSignal][0][0][2]<<"  |  "<<v2Err[0][choseSignal][0][0][2]<<"  |\n";

	  output << "|  pT.   |  v2_systematic only error  |\n";
	  output << "|  6.5-8.0    |  "<<v2Err_systOnly[choseSignal][0][0][0]<<"  |\n";
	  output << "|  8.0-10.0   |  "<<v2Err_systOnly[choseSignal][0][0][1]<<"  |\n";
	  output << "|  10.0-40.0  |  "<<v2Err_systOnly[choseSignal][0][0][2]<<"  |\n";
	
	  output << "|  pT   |  v2_systStat_noRes  |  v2_noRes error  |\n";
	  output << "|   6.5-8.0   |  "<<v2_final_sansRes[choseSignal][0][0][0]<<"  |  "<<v2Err_final_sansRes[choseSignal][0][0][0]<<"  |\n";
	  output << "|  8.0-10.0   |  "<<v2_final_sansRes[choseSignal][0][0][1]<<"  |  "<<v2Err_final_sansRes[choseSignal][0][0][1]<<"  |\n";
	  output << "|  10.0-40.0  |  "<<v2_final_sansRes[choseSignal][0][0][2]<<"  |  "<<v2Err_final_sansRes[choseSignal][0][0][2]<<"  |\n";
	 

	  output << "|  pT   |  v2_final  |  v2_final error  |\n";
	  output << "|   6.5-8.0   |  "<<v2_final[choseSignal][0][0][0]<<"  |  "<<v2Err_final[choseSignal][0][0][0]<<"  |\n";
	  output << "|  8.0-10.0   |  "<<v2_final[choseSignal][0][0][1]<<"  |  "<<v2Err_final[choseSignal][0][0][1]<<"  |\n";
	  output << "|  10.0-40.0  |  "<<v2_final[choseSignal][0][0][2]<<"  |  "<<v2Err_final[choseSignal][0][0][2]<<"  |\n";
	 
	  output<<endl;


	}//signal type

      // let's do some plotting please
      // ****** yields distributions, with statOnly and stat+systm uncertainties
      // Drawing jsut the v2 that has stat uncertainties with it
      TCanvas *pc7 = new TCanvas("pc7",Form("pc_dPhiDistribStatSyst"),0,0,1200,350);
      makeMultiPanelCanvas(pc7,nPads,1,0.0,0.0,0.2,0.15,0.02);

      int ind = 0;
      for(int choseSignal = signal_start; choseSignal<signal_end; choseSignal++)
	{
	  int vcts1  = cts[0]; 
	  int vcts2  = cts[1];
	  cout<<"Canvas Centrality: "<<vcts1<<"-"<<vcts2<<endl;
	  for(int ky = y_start; ky < y_end; ky++)
	    {
	      double vraps1 = raps[ky];
	      double vraps2 = raps[ky+1];
	      cout << "Canvas Rapidity = "<< vraps1 << "\t"<<vraps2<<endl;
	      
	      for(int lpt = pt_start; lpt < pt_end; lpt++)
		{
		  pc7->cd(ind+1);
		  cout<<"ind : "<<ind<<endl;
		  double vpts1 = pts[lpt]; 
		  double vpts2 = pts[lpt+1];
		  
		  pp->SetMaximum(1.);
		  pp->SetMinimum(0.3);
		  pp->Draw();
		  
		  lt1->SetTextSize(0.05);
		  if(!g[0][choseSignal][0][ky][lpt]) 
		    {
		      cout<<"No graph! continued !!!!"<<endl;
		      continue;
		    }
		  cout<<"#### Drawing: cent:"<<vcts1<<"-"<<vcts2<<"\t rapidity"<<vraps1<<"-"<<vraps2<<"\t pt"<<vpts1<<"-"<<vpts2<<endl;		     
		  TGraphErrors *pgTemp2 = (TGraphErrors *)pgFinal[choseSignal][0][0][lpt];
		  pgTemp2->SetMarkerColor(2);
		  pgTemp2->SetLineColor(2);
		  pgTemp2->SetMarkerStyle(20);
		  pgTemp2->Draw("[PZ]");
		  
		  TGraphErrors *pgTemp1 = (TGraphErrors *)g[0][choseSignal][0][ky][lpt];
		  pgTemp1->SetMarkerStyle(24);
		  pgTemp1->Draw("[P]");
		  if(ind==0 )
		    {
			  // drapwing v2 value and chi2
		      lt1->SetTextSize(0.06);
		      lt1->DrawLatex(0.75,0.89,Form("%s",legend[choseSignal]));  // what signal is

		      lt1->SetTextSize(0.05);
		      lt1->DrawLatex(0.24,0.23,Form("v_{2}^{stat} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
						    v2[0][choseSignal][0][ky][lpt],
						    v2Err[0][choseSignal][0][ky][lpt],
						    chi[0][choseSignal][0][ky][lpt],
						    ndf[0][choseSignal][0][ky][lpt]));
		      lt1->DrawLatex(0.24,0.18,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
						    v2_final[choseSignal][0][ky][lpt],
						    v2Err_final[choseSignal][0][ky][lpt],
						    chi_final[choseSignal][0][ky][lpt],
						    ndf_final[choseSignal][0][ky][lpt]));
		      lt1->SetTextSize(0.06);
		      lt1->DrawLatex(0.24,0.35,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); // centrality 
		     
		    }
		  else
		    {
		      lt1->SetTextSize(0.05);
		      lt1->DrawLatex(0.05,0.23,Form("v_{2}^{stat} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
						    v2[0][choseSignal][0][ky][lpt],
							v2Err[0][choseSignal][0][ky][lpt],
						    chi[0][choseSignal][0][ky][lpt],
						    ndf[0][choseSignal][0][ky][lpt]));
		      lt1->DrawLatex(0.05,0.18,Form("v_{2} = %.4f #pm %.4f (#chi^{2} / ndf = %.3f / %.0f)",
						    v2_final[choseSignal][0][ky][lpt],
						    v2Err_final[choseSignal][0][ky][lpt],
						    chi_final[choseSignal][0][ky][lpt],
						    ndf_final[choseSignal][0][ky][lpt]));
		      lt1->SetTextSize(0.06);
		      lt1->DrawLatex(0.05,0.35,Form("%.1f < p_{T} < %.1f GeV/c", vpts1, vpts2)); 
		      		      
		      if(ind == 2)
			{
			  lt1->DrawLatex(0.6,0.82,Form("|y| < %.1f",vraps2));       // rapidity
			  lt1->DrawLatex(0.6,0.89,Form("Cent. 10 - 60 %%")); 
			}
		    }
		  ind++;
		}//pt loop
	    }//rapidity loop

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
		
	  
	  //_______
	  const char* chosenSignal = signal[choseSignal];
	  if(bSavePlots)
	    {
	      gSystem->mkdir("./plots/final",kTRUE);
	      pc7->SaveAs(Form("./plots/final/pt_%s_rap%.1f_%.1f_pT%.1f-%.1f_a2_Corr.png",signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	      pc7->SaveAs(Form("./plots/final/pt_%s_rap%.1f_%.1f_pT%.1f-%.1f_a2_Corr.pdf",signal[choseSignal],raps[y_start],raps[y_end],pts[pt_start],pts[pt_end]));
	    }

      	  // write out the tgraph errors
	  
	  rootoutput->cd();
	  double pt_err[nptbins]     = {0.5,0.5,0.5};
	  double v2final[nptbins]  = {0.0} , v2finalErr[nptbins] = {0.0};

	  double v2SystOnlyErr[nptbins] = {0.0};
	  double v2StatOnlyErr[nptbins] = {0.0};

	  for (int ipt = pt_start; ipt < pt_end; ipt++) 
	    {
	      v2final[ipt]     = v2_final[choseSignal][0][0][ipt];
	      v2finalErr[ipt]  = v2Err_final[choseSignal][0][0][ipt];

	      v2StatOnlyErr[ipt]  = v2Err[0][choseSignal][0][0][ipt];          // statistical only, prefix=0, nominal, first entry in the array
	      v2SystOnlyErr[ipt]  = v2Err_systOnly[choseSignal][0][0][ipt]; // systematic errors only

	      cout << v2final[ipt] << " Last chance to screw up!!!!! " << v2finalErr[ipt] << endl;
	    }
	  
	  TGraphErrors pgWrite(nptbins,pts_bound,v2final,pt_err,v2finalErr);
	  pgWrite.Write(Form("final_ptDependence_%s_rap%.1f-%.1f_cent%d-%d",signal[choseSignal],raps[0],raps[1],cts[0],cts[1]));

	  TGraphErrors pgWrite_stat(nptbins,pts_bound,v2final,pt_err,v2StatOnlyErr);
	  pgWrite_stat.Write(Form("final_ptDependence_statErr_%s_rap%.1f-%.1f_cent%d-%d",signal[choseSignal],raps[0],raps[1],cts[0],cts[1]));

	  TGraphErrors pgWrite_syst(nptbins,pts_bound,v2final,pt_err,v2SystOnlyErr);
	  pgWrite_syst.Write(Form("final_ptDependence_systErr_%s_rap%.1f-%.1f_cent%d-%d",signal[choseSignal],raps[0],raps[1],cts[0],cts[1]));


	  if(doDebug)
	    {
	      double x,y;
	      cout<<pgWrite.GetPoint(1,x,y);
	      cout<<x<<"\t "<<y<<endl;

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
  gApplication->Terminate();
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

void TGetPoints(TGraphErrors *a, double *b, double *c)
{

  int na = a->GetN();
  for(int i = 0; i < na; i++){
    a->GetPoint(i, b[i], c[i]);
  }

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

//__________________________________________________________________________
//
