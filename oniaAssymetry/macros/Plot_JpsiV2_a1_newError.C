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
void GetV2(TGraphErrors *a, double *b);
void TGetPoints(TGraphErrors *a, double *b, double *c);
void assignArray(double adInput[4],double adOutput[4]);

 double dPhiBinWidth =  TMath::PiOver2()/4;
//__________________________________________________________________________
void Plot_JpsiV2_a1_newError()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gStyle->SetOptFit(0);
 
  // options
  bool bDoDebug   = false;
  bool bSavePlots = true; 

  const int ncentbins    = 4;
  int    cts[ncentbins+1]  = {0, 10, 20, 30, 60};
  int centrality_start   = 0;
  int centrality_end     = 4; 
  int nPads = centrality_end - centrality_start; 
 
  const int nsignals           = 3;
  const char* signals[nsignals]      = {"NSig","NPr","NNp"};
  const char* legend[nsignals] = {"Inclusive J/#psi","Prompt J/#psi","Non-prompt J/#psi"};
  int signal_start = 0;
  int signal_end   = 1; 
  
  const int nrapbins     = 2;
  double raps[nrapbins]  = {0.0, 2.4};
 
  const int nptbins      = 2; 
  double pts[nptbins]    = {6.5, 40.0};
  int pt_start = 0;
  int pt_end   = 1; 
 
  // phi bins
  const int nphibins = 4;
  float dPhiLimits[nphibins+1] = {0.0, TMath::PiOver2()/4, 2*TMath::PiOver2()/4, 3*TMath::PiOver2()/4, 4*TMath::PiOver2()/4};
  double dPhiBinCenters[nphibins];
  double dPhiBinCentersErr[nphibins];
 
  for(int ibin=0; ibin<nphibins; ibin++)
    {
      dPhiBinCenters[ibin]	= dPhiLimits[ibin]+dPhiBinWidth/2;
      dPhiBinCentersErr[ibin]   = 0.0;
      cout<<"bin limits"<< dPhiLimits[ibin] << "\t"<< dPhiLimits[ibin+1]<<"\t"<<dPhiBinCenters[ibin]<<endl;
    }
 

 // inclusive
   
  double gYieldP_0010_nsig[nphibins] = {0.6791,
					 0.6419,
					 0.6078,
					 0.6177}; // yield of etHFp for 0 - 10 %
  double gYieldP_1020_nsig[nphibins] = {0.6788,
					0.6745,
					0.6424,
					0.5508}; // yield of etHFp for 10 - 20 %
  double gYieldP_2030_nsig[nphibins] = {0.7705,
					0.5811,
					0.6172,
					0.5776}; // yield of etHFp for 20 - 30 %
  double gYieldP_3060_nsig[nphibins] = {0.6405,
					0.6994,
					0.5997,
					0.6069}; // yield of etHFp for 30 - 60 %
  
  double gYieldM_0010_nsig[nphibins] = {0.6796,
					 0.6196,
					 0.6358,
					 0.6115}; // yield of etHFm for 0 - 10 %
  double gYieldM_1020_nsig[nphibins] = {0.7550,
					0.6097,
					0.5786,
					0.6032}; // yield of etHFm for 10 - 20 %
  double gYieldM_2030_nsig[nphibins] = {0.7167,
					0.6719,
					0.5886,
					0.5693}; // yield of etHFm for 20 - 30 %
  double gYieldM_3060_nsig[nphibins] = {0.7554,
					0.6385,
					0.5759,
					0.5767}; // yield of etHFm for 30 - 60 %


  double gYieldP_0010_nsig_err[nphibins] = {0.0409,
					     0.0404,
					     0.0389,
					     0.0382}; // error of etHFp for 0 - 10 %
  double gYieldP_1020_nsig_err[nphibins] = {0.0436,
					    0.0427,
					    0.0413,
					    0.0386}; // error of etHFp for 10 - 20 %
  double gYieldP_2030_nsig_err[nphibins] = {0.0523,
					    0.0461,
					    0.0467,
					    0.0449}; // error of etHFp for 20 - 30 %
  double gYieldP_3060_nsig_err[nphibins] = {0.0396,
					    0.0409,
					    0.0378,
					    0.0379}; // error of etHFp for 30 - 60 %
  
  double gYieldM_0010_nsig_err[nphibins] = {0.0419,
					     0.0405,
					     0.0403,
					     0.0399}; // error of etHFm for 0 - 10 %
  double gYieldM_1020_nsig_err[nphibins] = {0.0469,
					    0.0425,
					    0.0414,
					    0.0419}; // error of etHFm for 10 - 20 %
  double gYieldM_2030_nsig_err[nphibins] = {0.0504,
					    0.0485,
					    0.0454,
					    0.0448}; // error of etHFm for 20 - 30 %
  double gYieldM_3060_nsig_err[nphibins] = {0.0423,
					    0.0389,
					    0.0369,
					    0.0371}; // error of etHFm for 30 - 60 %

   // non-prompt
  double gYieldP_0010_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 0 - 10 %
  double gYieldP_1020_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 10 - 20 %
  double gYieldP_2030_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 20 - 30 %
  double gYieldP_3060_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 30 - 60 %
  
  double gYieldM_0010_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 0 - 10 %
  double gYieldM_1020_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 10 - 20 %
  double gYieldM_2030_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 20 - 30 %
  double gYieldM_3060_nnp[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 30 - 60 %

  double gYieldP_0010_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 0 - 10 %
  double gYieldP_1020_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 10 - 20 %
  double gYieldP_2030_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 20 - 30 %
  double gYieldP_3060_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 30 - 60 %
  
  double gYieldM_0010_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 0 - 10 %
  double gYieldM_1020_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 10 - 20 %
  double gYieldM_2030_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 20 - 30 %
  double gYieldM_3060_nnp_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 30 - 60 %

  // prompt
  double gYieldP_0010_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 0 - 10 %
  double gYieldP_1020_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 10 - 20 %
  double gYieldP_2030_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 20 - 30 %
  double gYieldP_3060_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFp for 30 - 60 %
  
  double gYieldM_0010_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 0 - 10 %
  double gYieldM_1020_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 10 - 20 %
  double gYieldM_2030_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 20 - 30 %
  double gYieldM_3060_npr[nphibins] = {0.05, 0.05, 0.05, 0.05}; // yield of etHFm for 30 - 60 %

  double gYieldP_0010_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 0 - 10 %
  double gYieldP_1020_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 10 - 20 %
  double gYieldP_2030_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 20 - 30 %
  double gYieldP_3060_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFp for 30 - 60 %
  
  double gYieldM_0010_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 0 - 10 %
  double gYieldM_1020_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 10 - 20 %
  double gYieldM_2030_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 20 - 30 %
  double gYieldM_3060_npr_err[nphibins] = {0.05, 0.05, 0.05, 0.05}; // error of etHFm for 30 - 60 %

  const int nhf = 2;
  TGraphErrors *g[nsignals][ncentbins][nhf];
  double pr[nsignals][ncentbins][nhf], ep[nsignals][ncentbins][nhf], chi[nsignals][ncentbins][nhf], ndf[nsignals][ncentbins][nhf];
     
  double vpts1 = pts[pt_start];
  double vpts2 = pts[pt_end];
  
  double gYieldPlus[ncentbins];
  double gYieldMinus[ncentbins];
  double gYieldPlus_err[ncentbins];
  double gYieldMinus_err[ncentbins];

  for(int iSignal = signal_start; iSignal<signal_end; iSignal++)
      {
	const char* chosenSignal = signals[iSignal];
	for(int mcent = centrality_start; mcent < centrality_end; mcent++)
	  {
	   
	    TGraphErrors *pgTempP;
	    TGraphErrors *pgTempM;
	    switch (mcent) {
	    case 0:
	      if (iSignal == 0) 
		{
		  assignArray(gYieldP_0010_nsig,gYieldPlus);
		  assignArray(gYieldM_0010_nsig,gYieldMinus);

		  assignArray(gYieldP_0010_nsig_err,gYieldPlus_err);
		  assignArray(gYieldM_0010_nsig_err,gYieldMinus_err);
		  
		}
	      else 
		if (iSignal == 1)  
		  {
		    assignArray(gYieldP_0010_npr,gYieldPlus);
		    assignArray(gYieldM_0010_npr,gYieldMinus);
		 
		    assignArray(gYieldP_0010_npr_err,gYieldPlus_err);
		    assignArray(gYieldM_0010_npr_err,gYieldMinus_err);
		  }
		else 
		  {
		    assignArray(gYieldP_0010_nnp,gYieldPlus);
		    assignArray(gYieldM_0010_nnp,gYieldMinus);
		  
		    assignArray(gYieldP_0010_nnp_err,gYieldPlus_err);
		    assignArray(gYieldM_0010_nnp_err,gYieldMinus_err);
		  }
	      break;
	   
	    case 1:
	      if (iSignal == 0) 
		{
		  assignArray(gYieldP_1020_nsig,gYieldPlus);
		  assignArray(gYieldM_1020_nsig,gYieldMinus);
	
		  assignArray(gYieldP_1020_nsig_err,gYieldPlus_err);
		  assignArray(gYieldM_1020_nsig_err,gYieldMinus_err);
		}
	      else 
		if (iSignal == 1)  
		  {
		    assignArray(gYieldP_1020_npr,gYieldPlus);
		    assignArray(gYieldM_1020_npr,gYieldMinus);
		  
		    assignArray(gYieldP_1020_npr_err,gYieldPlus_err);
		    assignArray(gYieldM_1020_npr_err,gYieldMinus_err);
		  }
		else 
		  {
		    assignArray(gYieldP_1020_nnp,gYieldPlus);
		    assignArray(gYieldM_1020_nnp,gYieldMinus);
		  
		    assignArray(gYieldP_1020_nnp_err,gYieldPlus_err);
		    assignArray(gYieldM_1020_nnp_err,gYieldMinus_err);
		  }
	      break;
	    case 2:
	      if (iSignal == 0) 
		{
		  assignArray(gYieldP_2030_nsig,gYieldPlus);
		  assignArray(gYieldM_2030_nsig,gYieldMinus);
		
		  assignArray(gYieldP_2030_nsig_err,gYieldPlus_err);
		  assignArray(gYieldM_2030_nsig_err,gYieldMinus_err);
		}
	      else 
		if (iSignal == 1)  
		  {
		    assignArray(gYieldP_2030_npr,gYieldPlus);
		    assignArray(gYieldM_2030_npr,gYieldMinus);
		  
		    assignArray(gYieldP_2030_npr_err,gYieldPlus_err);
		    assignArray(gYieldM_2030_npr_err,gYieldMinus_err);
		  }
		else 
		  {
		    assignArray(gYieldP_2030_nnp,gYieldPlus);
		    assignArray(gYieldM_2030_nnp,gYieldMinus);
		 
		    assignArray(gYieldP_2030_nnp_err,gYieldPlus_err);
		    assignArray(gYieldM_2030_nnp_err,gYieldMinus_err);
		  }
	      break;
	      
	    case 3:
	      if (iSignal == 0) 
		{
		  assignArray(gYieldP_3060_nsig,gYieldPlus);
		  assignArray(gYieldM_3060_nsig,gYieldMinus);
		
		  assignArray(gYieldP_3060_nsig_err,gYieldPlus_err);
		  assignArray(gYieldM_3060_nsig_err,gYieldMinus_err);
		}
	      else 
		if (iSignal == 1)  
		  {
		    assignArray(gYieldP_3060_npr,gYieldPlus);
		    assignArray(gYieldM_3060_npr,gYieldMinus);
		  
		    assignArray(gYieldP_3060_npr_err,gYieldPlus_err);
		    assignArray(gYieldM_3060_npr_err,gYieldMinus_err);
		  }
		else 
		  {
		    assignArray(gYieldP_3060_nnp,gYieldPlus);
		    assignArray(gYieldM_3060_nnp,gYieldMinus);
		  
		    assignArray(gYieldP_3060_nnp_err,gYieldPlus_err);
		    assignArray(gYieldM_3060_nnp_err,gYieldMinus_err);
		  }
	      break;
	    default:
	      cout<< "Wrong centrality, please pick something else!"<<endl;
	    }//mcent switch

	    pgTempP = new TGraphErrors(nphibins, dPhiBinCenters, gYieldPlus, dPhiBinCentersErr, gYieldPlus_err);
	    pgTempM = new TGraphErrors(nphibins, dPhiBinCenters, gYieldMinus, dPhiBinCentersErr, gYieldMinus_err);
	    g[iSignal][mcent][0] = pgTempP;
	    g[iSignal][mcent][1] = pgTempM;
	    double cp[4] = {0.0, 0.0, 0.0, 0.0};
	    double cm[4] = {0.0, 0.0, 0.0, 0.0};
	    GetV2(g[iSignal][mcent][0], cp);
	    GetV2(g[iSignal][mcent][1], cm);
	    
	    pr[iSignal][mcent][0]  = cp[0];
	    ep[iSignal][mcent][0]  = cp[1];
	    chi[iSignal][mcent][0] = cp[2];
	    ndf[iSignal][mcent][0] = cp[3];

	    pr[iSignal][mcent][1]  = cm[0];
	    ep[iSignal][mcent][1]  = cm[1];
	    chi[iSignal][mcent][1] = cm[2];
	    ndf[iSignal][mcent][1] = cm[3];

	  }//centrality loop
	
	// #### Drawing: 
	TLatex *lt1 = new TLatex();
	lt1->SetNDC();
	lt1->SetTextSize(0.05);
	
	// frawing: 
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
	    if(!g[iSignal][mcent][0]|| !g[iSignal][mcent][1]) 
	      {
		cout<<"No graph! continued !!!!"<<endl;
		continue;
	      }
	    int vcts1 = cts[mcent]; 
	    int vcts2 = cts[mcent+1];
	    cout<<"Centrlality bin: "<<vcts1<<"-"<<vcts2<<endl;
	    pc1->cd(ind+1);
	    pp->SetMaximum(1.2);
	    pp->SetMinimum(0.2);		    
	    pp->Draw();

	    TGraphErrors *pgP = (TGraphErrors *) g[iSignal][mcent][0];
	    TGraphErrors *pgM = (TGraphErrors *) g[iSignal][mcent][1];

	    pgP->SetMarkerStyle(20);
	    pgM->SetMarkerStyle(24);
	    pgP->Draw("pz");
	    pgM->Draw("pz");
	    
	    if(ind==0 )
	      {
		 // drapwing v2 value and chi2
		 lt1->SetTextSize(0.045);
		 lt1->DrawLatex(0.24,0.23,Form("v_{2}^{etHFp} = %.2f #pm %.2f (#chi^{2} / ndf = %.2f / %.0f)", pr[iSignal][mcent][0], ep[iSignal][mcent][0], chi[iSignal][mcent][0], ndf[iSignal][mcent][0]));
		 lt1->DrawLatex(0.24,0.18,Form("v_{2}^{etHFm} = %.2f #pm %.2f (#chi^{2} / ndf = %.2f / %.0f)", pr[iSignal][mcent][1], ep[iSignal][mcent][1], chi[iSignal][mcent][1], ndf[iSignal][mcent][1]));
		 lt1->SetTextSize(0.055);
		 lt1->DrawLatex(0.24,0.3,Form("Cent. %d - %d %%",vcts1, vcts2)); // centrality 
	       }
	     else
	       {
		 lt1->SetTextSize(0.05);
		 lt1->DrawLatex(0.05,0.23,Form("v_{2}^{etHFp} = %.2f #pm %.2f (#chi^{2} / ndf = %.2f / %.0f)", pr[iSignal][mcent][0], ep[iSignal][mcent][0], chi[iSignal][mcent][0], ndf[iSignal][mcent][0]));
		 lt1->DrawLatex(0.05,0.18,Form("v_{2}^{etHFm} = %.2f #pm %.2f (#chi^{2} / ndf = %.2f / %.0f)", pr[iSignal][mcent][1], ep[iSignal][mcent][1], chi[iSignal][mcent][1], ndf[iSignal][mcent][1]));

		 lt1->SetTextSize(0.06);
		 if(ind == 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
		 if(ind > 1) lt1->DrawLatex(0.05,0.3,Form("Cent. %d - %d %%",vcts1, vcts2));
		 
		 if(ind == 3)
		   {
		     lt1->SetTextSize(0.06);
		     lt1->DrawLatex(0.45,0.8,Form("|y| < %.1f",raps[1]));       // rapidity
		     lt1->DrawLatex(0.45,0.73,Form("%.1f < p_{T} < %.0f GeV/c", vpts1, vpts2)); // pt
		   }
		 if(ind==1)
		   {
		     TLegend *legCent = new TLegend(0.7,0.8,0.95,0.95);
		     legCent->SetFillColor(0);
		     legCent->SetBorderSize(0);
		     legCent->SetTextSize(0.07);
		     legCent->AddEntry(pgP,"etHFp","P");
		     legCent->AddEntry(pgM,"etHFm","P");
		     legCent->Draw("same");
		   }
	       }
	     ind++;
	  }//mcent
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
	lt1->DrawLatex(0.05,0.89,Form("%s",legend[iSignal]));  // what signal is
	
	//_______
	if(bSavePlots)
	  {
	    pc1->SaveAs(Form("%s_Jpsi_%d_%d_a1_newErr_Nominal.png",chosenSignal,pt_start,pt_end));
 	    pc1->SaveAs(Form("%s_Jpsi_%d_%d_a1_newErr_Nominal.pdf",chosenSignal,pt_start,pt_end));
	  }
	  
	//	pc1->Clear();
    
	
	cout<<""<<endl;
	cout<<"%%%%% Category etHFp: "<<signals[iSignal]<<" %%%%%"<<endl;
	cout<<"|  Cent.  |  0 - 10 %  |  10 - 20 %  |  20 - 30 %  |  30 - 60 %  |"<<endl;
	cout<<"|  v2  |"<<"  "<<pr[iSignal][0][0]<<"  |  "<<pr[iSignal][1][0]<<"  |  "<<pr[iSignal][2][0]<<"  |  "<<pr[iSignal][3][0]<<"  |"<<endl;
	cout<<"|  err  |"<<"  "<<ep[iSignal][0][0]<<"  |  "<<ep[iSignal][1][0]<<"  |  "<<ep[iSignal][2][0]<<"  |  "<<ep[iSignal][3][0]<<"  |"<<endl;

	cout<<"%%%%% Category etHFm: "<<signals[iSignal]<<" %%%%%"<<endl;
	cout<<"|  Cent.  |  0 - 10 %  |  10 - 20 %  |  20 - 30 %  |  30 - 60 %  |"<<endl;
	cout<<"|  v2  |"<<"  "<<pr[iSignal][0][1]<<"  |  "<<pr[iSignal][1][1]<<"  |  "<<pr[iSignal][2][1]<<"  |  "<<pr[iSignal][3][1]<<"  |"<<endl;
	cout<<"|  err  |"<<"  "<<ep[iSignal][0][1]<<"  |  "<<ep[iSignal][1][1]<<"  |  "<<ep[iSignal][2][1]<<"  |  "<<ep[iSignal][3][1]<<"  |"<<endl;

	
      }// chose signal


}


//____________________________________________________________________________
void assignArray(double adIn[4], double adOut[4])
{
  // divide by dPhi
  for(int i=0; i<4;i++)
    {
      adOut[i] = adIn[i];
    }
  //cout<<"Finished copying"<<endl;
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
