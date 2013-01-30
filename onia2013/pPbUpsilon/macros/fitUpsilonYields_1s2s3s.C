#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooKeysPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

// Root stuff
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"

#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

// miscellaneous  
#include <fstream>
#include <iostream>

//#include "TMath"

/*INPUTs*/

double mass_l =  7.0;
double mass_h = 14.0;
double binw   = 0.1;    //bin width of the histogram

int choseFitParams             = 1; //0: 1s, 2s, 3s separatelly 1: 1s, 2s/1s; 3s/1s; 2: 1S, (2s+3s)/1s
int choseCase                  = 4; 
const int nChoices             = 11;
const char* fitCases[nChoices] = {"upsY_minM24_maxP146_cent0100",//0: defautl case
				  "upsY_minM24_maxP146_cent010",//1: integrated over rapidity and 0-10%
				  "upsY_minM24_maxP146_cent60100",//2: integrated over rapidity and 60-100%
				  "upsY_minM24_maxM047_cent0100", //3: 0-100%, bkwd integrated
				  "upsY_minM047_maxP146_cent0100",//4: 0-100%, fwd integrated
				  "upsY_minM24_maxM047_cent020",  //5: 0-20%, bkwd integrated
				  "upsY_minM047_maxP146_cent020", //6: 0-20%, fwd integrated
				  "upsY_minM24_maxM047_cent2050", //7: 20-50, bkwd integrated
				  "upsY_minM047_maxP146_cent2050", //8: 20-50, fwd integrated
				  "upsY_minM24_maxM047_cent50100", //9: 50-100, bkwd integrated
				  "upsY_minM047_maxP146_cent50100"//10: 50-100, fwd integrated
};
int choseData =  1;      //Input data sample.  0: pA data;   1: pp@7TeV data ; 2: pp@2.76TeV; 3: PbPb@276
const int nData   = 4;
const char* choseDataCases[nData] = {"pA_5teV",
				     "pp_7tev",
				     "pp_2p76tev",
				     "pbpb_2p76tev"
};

const char* choseDataLegend[nData] = {"pPb #sqrt{s_{NN}} = 5 TeV",
				      "pp #sqrt{s} = 7 TeV",
				      "pp #sqrt{s} = 2.76 TeV",
				      "PbPb #sqrt{s_{NN}} = 2.76 TeV"
};

const char* choseDataLumi[nData] = {"L_{int} = 5.5 nb^{-1}",
				    "L_{int} = xxx nb^{-1}",
				    "L_{int} = 231 nb^{-1}",
				    "L_{int} = 150 \mub^{-1}"
};

// do not touch these
int centrality_min    = 0;
int centrality_max    = 100;

double upsYCut_min    = -2.4;
double upsYCut_max    = 1.46;

double muonpTcut      = 4;      //single muon pT cut
double muonEtaCut_min =-2.4;
double muonEtaCut_max = 1.46; 

// for symmetric collisions, use the ***_coll corresponding to the ***_lab
const int nEtaBins = 2;

double muonEtaCut_coll[nEtaBins+1] = {-1.93, 0.0, 1.93};
double upsYCut_coll[nEtaBins+1]    = {-1.93, 0.0, 1.93};

double muonEtaCut_lab[nEtaBins+1] = {-2.4, -0.47, 1.46};
double upsYCut_lab[nEtaBins+1]    = {-2.4, -0.47, 1.46};

double muonEtaCut[nEtaBins+1];
double upsYCut[nEtaBins+1];


// -------- up to here for the moment

bool fitMB        = 1;      //1: fit the the Minbias sample(0-100%);   0: fit the each centrality bin
double width_     = 0.0782; //new resolution for 2011 HI data; taken from where?!!!!! ====== will have to change !!!!!!!!!!!!!!!!!!!!!
 
bool plotpars     = 1;      //1: plot parameters;   0: plot CMS label
bool doMinos      = 0;      //kFALSE;

int  bkgdModel    = 3;     //Background Model.  1: LS erf*exp + polynomial; 2: LS RookeyPdf + polynomial; 3: erf*exp; 4: polynomial; 5: erf*exp+polynomial
bool plotLikeSign = 1;     //0: hide likesign or trkRot; 1: plot likesign or trkRot data points and fit lines;
bool TRKROT       = 0;     //0: use likesign;   1: use track rotation
bool LS_constrain = 0;     //1: use constrain method

// replace these values
RooRealVar *f2Svs1S_pp = new RooRealVar("N_{2S}/N_{1S}pp","Y(2S)/Y(1S) yields pp ratio",0.59,-1,5);
RooRealVar *f3Svs1S_pp = new RooRealVar("N_{3S}/N_{1S}pp","Y(3S)/Y(1S) yields pp ratio",0.41,-1,5);
ofstream outfile("fitresults.out", ios_base::app);

TString finput;
TString figName_;

TString figs_("newplots/"); //output fig location
const TString dirname_("");//"upsilonYield"); tree location in input file 
TString paramOn_("");
TString suffix_, cut_;

using namespace RooFit;
using namespace RooStats;
void fitpeaks(int choseCase);
void sidebandFit(bool,        RooRealVar*,   RooRealVar*, RooRealVar*, RooAbsPdf*, RooAbsPdf*, RooDataSet*, RooPlot*);
void setUpLimits(float xMin, float xMax,     RooAbsPdf*,  RooDataSet*, RooRealVar*,float baseNll);
pair<double, double> ConfidencInterval(float,RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf);

void fitUpsilonYields_1s2s3s()
{
  // gROOT->LoadMacro("setTDRStyle_modified.C");
  //setTDRStyle();
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  
  switch (choseData)
    {
    case 0: // pPb @ 5TeV
      finput   = "outfiles/dimuonTree_upsiMiniTree_pA5tev_prompt5p5nb_RunR210498_210737_partial_trigBit1_allTriggers0_pt4.root";
	//"outfiles/dimuonTree_upsiMiniTree_pA5tev_prompt2p65nb_Runr210614_210638_trigBit1_allTriggers0_pt4.root";
      break;
    case 1: // pp @ 7TeV
      finput   = "outfiles/dimuonTree_upsiMiniTree_7tev_Runa_trigBit1_allTriggers0_pt4.root";
      mass_l = 8.5;
      mass_h = 11.5;
      break;
    case 2://pp @ 2.76TeV
      finput   = "outfiles/dimuonTree_upsiMiniTree_pp_276tev_trigBit1_pt4.root";
      break;
    case 3://PbPb @ 2.76TeV
      finput   = "outfiles/dimuonTree_upsiMiniTree_pbpb_276tev_Runa_trigBit1_allTriggers0_pt4.root";
      break;
    default:
      cout<<"You don't know what you are doing! Pick one of the available datasets in the choseDataCases[] array"<<endl;
      break;
  }

  figName_ = Form("masspeak_");
  //different binning:
  fitpeaks(choseCase);
}


//____________________________________________________________________________________________
void fitpeaks(int choseCase)
{
  if(choseData>0)// if symmetric collision: pp, PbPb
    { 
      upsYCut_min    = upsYCut_coll[0];
      upsYCut_max    = upsYCut_coll[nEtaBins];
      muonEtaCut_min = muonEtaCut_coll[0];
      muonEtaCut_max = muonEtaCut_coll[nEtaBins];
      cout<<"!!!! I'm in this case!!"<<endl;
    for (int ibin=0; ibin<nEtaBins; ibin++)
      {
	upsYCut[ibin]    = upsYCut_coll[ibin];
	muonEtaCut[ibin] = muonEtaCut_coll[ibin];
      }
    }
  else
    {
      for (int ibin=0; ibin<nEtaBins; ibin++)
	{
	  upsYCut[ibin]    = upsYCut_lab[ibin];
	  muonEtaCut[ibin] = muonEtaCut_lab[ibin];
	}
    }

  switch (choseCase)
    {
    case 1:// "upsY_minM24_maxP146_cent010": 0-10, all
      centrality_max = 10;
      break;
    case 2:// "upsY_minM24_maxP146_cent60100": 60-100, all
      centrality_min = 60;
      break;
    case 3:// "upsY_minM24_maxM047_cent0100": 0-100, bkwd integrated
      upsYCut_max    = upsYCut[1];
      break;
    case 4://"upsY_minM047_maxP146_cent0100": 0-100, fwd integrated
      upsYCut_min    =  upsYCut[1];
      break;
     case 5://"upsY_minM24_maxM047_cent020": 0-20%, bkwd integrated
      upsYCut_max    = upsYCut[1];
      centrality_max = 20;
      break;
    case 6://"upsY_minM047_maxP146_cent020": 0-20%, fwd integrated
      upsYCut_min    = upsYCut[1];
      centrality_max = 20;
      break;
    case 7:// "upsY_minM24_maxM047_cent2050": 20-50, bkwd integrated
      upsYCut_max    = upsYCut[1];
      centrality_min = 20;
      centrality_max = 50;
      break;
    case 8://"upsY_minM047_maxP146_cent2050": 20-50, fwd integrated
      upsYCut_min    = upsYCut[1];
      centrality_min = 20;
      centrality_max = 50;
      break;
    case 9: //"upsY_minM24_maxM047_cent50100": 50-100, bkwd integrated
      upsYCut_max    = upsYCut[1];
      centrality_min = 50;
      break;
    case 10:// "upsY_minM047_maxP146_cent50100": 50-100, fwd integrated
      upsYCut_min    = upsYCut[1];
      centrality_min = 50;
      break;
    default:
      cout<<"Fitting integrated in rapidity and centrlaity!"<<endl;
      break;
    }
   cut_   = Form("(%.2f<=upsRapidity && upsRapidity<%.2f) && (%d<=Centrality && Centrality<%d)", 
		 upsYCut_min, upsYCut_max, centrality_min, centrality_max);
   suffix_= Form("%s_%s",choseDataCases[choseData],fitCases[choseCase]);

   cout<<"Fitting: y["<< upsYCut_min <<","<<upsYCut_max<<"] and centrality ["<<centrality_min<<","<<centrality_max<<"]!!!!"<<endl;

  //f2Svs1S_pp->setVal(0.5569);
  //f3Svs1S_pp->setVal(0.4140);
  f2Svs1S_pp->setVal(0);
  f3Svs1S_pp->setVal(0);
 
  cout << "oniafitter processing"
       << "\n\tInput:  \t" << finput
       << "\n\tresults:\t" << figs_
       << endl;
  outfile<<endl<<"**********"<<suffix_<<"**********"<<endl<<endl;
  
  // ###### read the data
  TFile f(finput,"read");
  gDirectory->Cd(finput+":/"+dirname_);
  TTree* theTree       = (TTree*)gROOT->FindObject("UpsilonTree"); // OS --- all mass
  TTree* allsignTree   = (TTree*)gROOT->FindObject("UpsilonTree_allsign");//all sign and all mass
  TTree* trkRotTree=0;
  if (TRKROT) trkRotTree = (TTree*)gROOT->FindObject("UpsilonTree_trkRot");
  
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");
  RooRealVar* upsPt      = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,60,"GeV");
  RooRealVar* upsEta     = new RooRealVar("upsEta",  "upsEta"  ,-10,10);
  RooRealVar* upsRapidity= new RooRealVar("upsRapidity",  "upsRapidity",upsYCut_min, upsYCut_max);
  RooRealVar* vProb      = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);
  RooRealVar* QQsign     = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
  RooRealVar* Centrality = new RooRealVar("Centrality",  "Centrality"  ,0,1000);
  
  RooRealVar* muPlusPt   = new RooRealVar("muPlusPt","muPlusPt",muonpTcut,50);
  RooRealVar* muMinusPt  = new RooRealVar("muMinusPt","muMinusPt",muonpTcut,50);
  RooRealVar* muPlusEta  = new RooRealVar("muPlusEta","muPlusEta", muonEtaCut_min, muonEtaCut_max);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",muonEtaCut_min, muonEtaCut_max);
  
  // *************************************************** importing
  //##### import unlike-sign data set
  RooDataSet* data0, *data, *likesignData0, *likesignData, *TrkRotData0, *TrkRotData;
  data0 = new RooDataSet("data","data",theTree,RooArgSet(*mass,*upsRapidity,*vProb,*upsPt,*Centrality,*muPlusPt,*muMinusPt));
  data0->Print();

  cout<<"########3 CUT "<<cut_<<endl;
  data = ( RooDataSet*) data0->reduce(Cut(cut_));
  data->Print();
  
  // ###### import like-sign data set
  likesignData0 = new RooDataSet("likesignData","likesignData",allsignTree,
				 RooArgSet(*mass,*upsRapidity,*vProb,*upsPt,*Centrality,*muPlusPt,*muMinusPt,*QQsign));

  likesignData0->Print();
  likesignData = ( RooDataSet*) likesignData0->reduce(Cut(cut_+" && QQsign != 0"));
  likesignData->Print();
  
  //import track-rotation data set
  if (TRKROT) 
    {
      TrkRotData0 = new RooDataSet("TrkRotData","TrkRotData",trkRotTree,
				   RooArgSet(*mass,*upsRapidity,*vProb,*upsPt,*Centrality,*muPlusPt,*muMinusPt,*QQsign));
      TrkRotData0->Print();
      TrkRotData = ( RooDataSet*) TrkRotData0->reduce(Cut(cut_+" && QQsign != 0"));
      TrkRotData->Print();
    }
  
  // *************************************************** signal PDF
  const double M1S = 9.46;   //upsilon 1S pgd mass value
  const double M2S = 10.023;  //upsilon 2S pgd mass value
  const double M3S = 10.355;  //upsilon 3S pgd mass value
  
  RooRealVar  *mean = new RooRealVar("#mu_{#Upsilon(1S)}","#Upsilon mean",M1S,M1S-0.1,M1S+0.1);
  RooConstVar *rat2 = new RooConstVar("rat2", "rat2", M2S/M1S);
  RooConstVar *rat3 = new RooConstVar("rat3", "rat3", M3S/M1S);

  // scale mean and resolution by mass ratio
  RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
  RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
  RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));

 //detector resolution ?? where is this coming from?
  RooRealVar    *sigma1  = new RooRealVar("sigma1","sigma1",0.092,0.01,0.3);
  RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
  RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
  RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
  
  /// to describe final state radiation tail on the left of the peaks
  RooRealVar *alpha  = new RooRealVar("alpha","tail shift",5.6,0.1,10);
  RooRealVar *npow   = new RooRealVar("npow","power order",4.3,0.1,10);       // MB case

  if (choseCase!=0) 
    {
      //  alpha->setConstant(kTRUE);
      //  npow ->setConstant(kTRUE);
    }
 

  // relative fraction of the two Gaussians components for each CB
  RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.3,0.,1.);
  sigmaFraction->setVal(0);
  sigmaFraction->setConstant(kTRUE);
  
  /// Upsilon 1S
  RooCBShape  *gauss1S_1 = new RooCBShape ("gauss1S_1", "FSR cb 1s",
					  *mass,*mean1S,*sigma1,*alpha,*npow);
  RooCBShape  *gauss1S_2 = new RooCBShape ("gauss1S_2", "FSR cb 1s",
					  *mass,*mean1S,*sigma1S,*alpha,*npow);
  RooAddPdf      *sig1S  = new RooAddPdf  ("sig1S","1S mass pdf",
					   RooArgList(*gauss1S_1,*gauss1S_2),*sigmaFraction);
  // sig1S is then jsut gauss1S2: c*pdf_1+(1-c)*pdf_2
  
  /// Upsilon 2S
  RooCBShape  *gauss2S_1 = new RooCBShape ("gauss2S_1", "FSR cb 2s", 
					  *mass,*mean2S,*sigma1,*alpha,*npow); 
  RooCBShape  *gauss2S_2 = new RooCBShape ("gauss2S_2", "FSR cb 2s", 
					  *mass,*mean2S,*sigma2S,*alpha,*npow); 
  RooAddPdf      *sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
					  RooArgList(*gauss2S_1,*gauss2S_2),*sigmaFraction);
  
  /// Upsilon 3S
  RooCBShape  *gauss3S_1 = new RooCBShape ("gauss3S_1", "FSR cb 3s", 
					  *mass,*mean3S,*sigma1,*alpha,*npow); 
  RooCBShape  *gauss3S_2 = new RooCBShape ("gauss3S_2", "FSR cb 3s", 
	 				  *mass,*mean3S,*sigma3S,*alpha,*npow); 
  RooAddPdf      *sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
					  RooArgList(*gauss3S_1,*gauss3S_2),*sigmaFraction); // = gauss3S1*sigmaFrac + gauss3S2*(1-sigmaFrac)
  
  // *************************************************** free param in the fit
  int nt = 100000;
  RooRealVar *nsig1f   = new RooRealVar("N_{1S}","nsig1S",nt*0.25,0,10*nt);
  switch (choseFitParams)
    {
    case 0://use the YIELDs of 2S and 3S as free parameters
      RooRealVar *nsig2f  = new RooRealVar("N_{2S}","nsig2S",   nt*0.25,-1*nt,10*nt);
      RooRealVar *nsig3f  = new RooRealVar("N_{3S}","nsig3S",   nt*0.25,-1*nt,10*nt);
      break;
    case 1:  //use the RATIOs of 2S and 3S as free parameters
      RooRealVar *f2Svs1S   = new RooRealVar("N_{2S}/N_{1S}","f2Svs1S",0.21,-0.1,1.0);
      RooRealVar *f3Svs1S   = new RooRealVar("N_{3S}/N_{1S}","f3Svs1S",0.0,-0.1,1.0);
      RooFormulaVar *nsig2f = new RooFormulaVar("nsig2S","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
      RooFormulaVar *nsig3f = new RooFormulaVar("nsig3S","@0*@1", RooArgList(*nsig1f,*f3Svs1S));
      f2Svs1S->setConstant(kFALSE);
      f3Svs1S->setConstant(kFALSE);
      break;
    case 2:// do (2s+3s)/1s
      RooRealVar *f23vs1S   = new RooRealVar("N_{2S+3S}/N_{1S}","f23vs1S",0.45,-0.1,1);
      RooFormulaVar *nsig3f = new RooFormulaVar("nsig3S","@0*@2-@0*@1", 
						RooArgList(*nsig1f,*f2Svs1S,*f23vs1S));
      break;
    default:
      cout<<"Make a pick from choseFitPArams!!!"<<endl;
      break;

    }  


  // bkg Chebychev
  RooRealVar *nbkgd   = new RooRealVar("N_{bkg}","nbkgd",          nt*0.75,-100,10*nt);
  RooRealVar *bkg_a1  = new RooRealVar("bkg_{a1}", "bkg_{a1}", 0, -2, 2);
  RooRealVar *bkg_a2  = new RooRealVar("bkg_{a2}", "bkg_{a2}", 0, -2, 2);
  RooRealVar *bkg_a3  = new RooRealVar("bkg_{a3}", "bkg_{a3}", 0, -2, 2);

  RooAbsPdf  *bkgPdf  = new RooChebychev("bkgPdf","bkgPdf",
					 *mass, RooArgList(*bkg_a1,*bkg_a2));

  //  likesign
  RooRealVar *nLikesignbkgd = new RooRealVar("NLikesign_{bkg}","nlikesignbkgd",nt*0.75,0,10*nt);
  if (TRKROT) 
    {
      nLikesignbkgd->setVal(TrkRotData->sumEntries());
      nLikesignbkgd->setError(sqrt(TrkRotData->sumEntries()));
    }
  else 
    {
      nLikesignbkgd->setVal(likesignData->sumEntries());
      nLikesignbkgd->setError(sqrt(likesignData->sumEntries()));
    }
  
  if (LS_constrain) 
    {
      RooGaussian* nLikesignbkgd_constr = new RooGaussian("nLikesignbkgd_constr","nLikesignbkgd_constr",
							  *nLikesignbkgd,
							  RooConst(nLikesignbkgd->getVal()),    //mean
							  RooConst(nLikesignbkgd->getError())); //sigma
    }
  else nLikesignbkgd->setConstant(kTRUE);
  
  RooFormulaVar *nResidualbkgd = new RooFormulaVar("NResidual_{bkg}","@0-@1",
						   RooArgList(*nbkgd,*nLikesignbkgd));
 

 // *************************************************** bkgModel
  RooRealVar turnOn("turnOn","turnOn", 6., 0., 20.);
  turnOn.setConstant(false);
  RooRealVar width("width","width", 1., 0., 20.);
  width.setConstant(false);
  RooRealVar decay("decay","decay", 7., 0., 20.);
  decay.setConstant(false);
  RooGaussian* turnOn_constr;
  RooGaussian* width_constr;
  RooGaussian* decay_constr;

  switch (bkgdModel) 
    {
    case 1 :  //(err+pol2 ) to fit the SS, then fix the shape and fit OS, in case of LS_constrain option
      RooGenericPdf *thisPdf = new  RooGenericPdf("thisPdf","thisPdf",
						      "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
						      RooArgList(*mass,turnOn,width,decay));
      RooFitResult* fit_1st      = thisPdf->fitTo(*likesignData,Save()) ; // likesign data
      if (TRKROT) fit_1st        = thisPdf->fitTo(*TrkRotData,Save()) ;
        
      if (LS_constrain) 
	{ // allow parameters to vary within cenral value from above fit + their sigma
	  turnOn_constr = new RooGaussian("turnOn_constr","turnOn_constr",
					   turnOn,
					   RooConst(turnOn.getVal()),
					   RooConst(turnOn.getError()));
	  width_constr   = new RooGaussian("width_constr","width_constr",
					   width,
					   RooConst(width.getVal()),
					   RooConst(width.getError()));
	  decay_constr    = new RooGaussian("decay_constr","decay_constr",
					   decay,
					   RooConst(decay.getVal()),
					   RooConst(decay.getError()));
	}
      else 
	{
	  turnOn.setConstant(kTRUE);
 	  width.setConstant(kTRUE);
 	  decay.setConstant(kTRUE);
	}
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						      RooArgList(*bkgPdf,*thisPdf),
						      RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    case 2 : //use RooKeysPdf to smooth the SS, then fit OS with pol+keys
      RooKeysPdf *thisPdf        = new RooKeysPdf("thisPdf","thisPdf",*mass,*likesignData,
						  RooKeysPdf::MirrorBoth, 1.4);
      if (TRKROT) *thisPdf       = new RooKeysPdf("thisPdf","thisPdf",*mass,*TrkRotData,
						  RooKeysPdf::MirrorBoth, 1.5);
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						      RooArgList(*bkgPdf,*thisPdf),
						      RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    case 3 : //use error function to fit the OS directly
      RooGenericPdf *thisPdf     = new  RooGenericPdf("thisPdf","thisPdf",
						      "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
						      RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf("pdf_combinedbkgd","total combined background pdf",
						     RooArgList(*thisPdf),
						     RooArgList(*nbkgd));
      break;
      
    case 4 : //use pol 2 to fit the OS directly
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						      RooArgList(*bkgPdf),
						      RooArgList(*nbkgd));
      break;
    case 5 : //use ( error function + polynomial 2) to fit the OS directly
      RooGenericPdf *thisPdf     = new  RooGenericPdf("thisPdf","thisPdf",
						      "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
							  RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("pdf_combinedbkgd","total combined background pdf",
						      RooArgList(*bkgPdf,*thisPdf),
						      RooArgList(*nResidualbkgd,*nLikesignbkgd));
      break;
    default :
      cout<<"Donno what you are talking about! Pick anothe fit option!"<<endl;
      break;
    }
  

  //###### the nominal fit with default pdf 
  if (LS_constrain) 
    {
      RooAbsPdf  *pdf_unconstr   = new RooAddPdf ("pdf_unconstr","total signal+background pdf",
						  RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
						  RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
      RooProdPdf *pdf            = new RooProdPdf ("pdf","total constr pdf",
						   RooArgSet(*pdf_unconstr,*turnOn_constr,*width_constr,*decay_constr,*nLikesignbkgd_constr));
      RooFitResult* fit_2nd      = pdf->fitTo(*data,Constrained(),Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }
  else 
    {
      RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total signal+background pdf",
						   RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
						   RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
      RooFitResult* fit_2nd       = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),Minos(doMinos));
    }
  
  
  // *************************************************** plotting
  TCanvas c; c.cd();
  int nbins = ceil((mass_h-mass_l)/binw); 
  RooPlot* frame = mass->frame(Bins(nbins),Range(mass_l,mass_h));
 
  data->plotOn(frame,Name("theData"),MarkerSize(0.8));
  pdf->plotOn(frame,Name("thePdf"));
  if (plotLikeSign) {
    if (TRKROT) TrkRotData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kMagenta),MarkerStyle(22));
    else likesignData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kRed),MarkerStyle(24));
      // thisPdf->plotOn(frame,Name("theLikeSign"),VisualizeError(*fit_1st,1),FillColor(kOrange));
      //thisPdf->plotOn(frame,Name("theLikeSign"),LineColor(kRed));
  }
  RooArgSet * pars = pdf->getParameters(data);
  //RooArgSet * pars = thisPdf->getParameters(likesignData);
  
  //calculate chi2 in a mass range
  float bin_Min  = (8.2-mass_l)/binw;
  float bin_Max  = (10.8-mass_l)/binw;
  int binMin     = ceil(bin_Min);
  int binMax     = ceil(bin_Max);
  int nfloatpars = pars->selectByAttrib("Constant",kFALSE)->getSize();
  float myndof   = ceil((10.8-8.2)/binw) - nfloatpars;
  cout<<binMin<<" "<<binMax<<" "<<nfloatpars<<" "<<myndof<<endl;
  double mychsq  = frame->mychiSquare("thePdf","theData",nfloatpars,true,binMin,binMax)*myndof;
  //double mychsq = frame->mychiSquare("theLikeSign","theLikeSignData",nfloatpars,true,binMin,binMax)*myndof;
  
  /*
    int nfloatpars = pars->selectByAttrib("Constant",kFALSE)->getSize();
    float myndof = frame->GetNbinsX() - nfloatpars;
    double mychsq = frame->chiSquare("theLikeSign","theLikeSignData",nfloatpars)*myndof;
  */
  //plot parameters
  if(plotpars) {
    paramOn_ = "_paramOn_";
    pdf->paramOn(frame,Layout(0.15,0.6,0.4),Layout(0.5,0.935,0.97),Label(Form("#chi^{2}/ndf = %2.1f/%2.0f", mychsq,myndof)));
  }
  
  outfile<<"Y(1S) yield  : = "<<nsig1f->getVal()<<" +/- "<<nsig1f->getError()<<endl<<endl;
  outfile<<"free parameter = "<< nfloatpars << ", mychi2 = " << mychsq << ", ndof = " << myndof  << endl << endl;
  
  //draw the fit lines and save plots
  data->plotOn(frame,Name("theData"),MarkerSize(0.8));
  pdf->plotOn(frame,Components("bkgPdf"),Name("theBkg"),LineStyle(5),LineColor(kGreen));
  pdf->plotOn(frame,Components("pdf_combinedbkgd"),LineStyle(kDashed));
  if (plotLikeSign) 
    {
      if (TRKROT) pdf->plotOn(frame,Components("thisPdf"),Name("theLikeSign"),LineStyle(9),LineColor(kMagenta));
      else  pdf->plotOn(frame,Components("thisPdf"),Name("theLikeSign"),LineStyle(9),LineColor(kRed));
    }
  pdf->plotOn(frame,Name("thePdf"));
  data->plotOn(frame,MarkerSize(0.8));
  if (plotLikeSign) 
    {
      if (TRKROT) TrkRotData->plotOn(frame,Name("theTrkRotData"),MarkerSize(0.8),MarkerColor(kMagenta),MarkerStyle(22));
      else likesignData->plotOn(frame,Name("theLikeSignData"),MarkerSize(0.8),MarkerColor(kRed),MarkerStyle(24));
    }   
  
  
  frame->SetTitle( "" );
  frame->GetXaxis()->SetTitle("#mu^{+}#mu^{-} (GeV/c^{2})");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->Draw();
  
  
  //plot parameters
  TLatex latex1;
  latex1.SetNDC();
  if(!plotpars) 
    {
      paramOn_ = ""; 
      latex1.DrawLatex(0.55,1.-0.05*3,"CMS Preliminary");
      latex1.DrawLatex(0.55,1.-0.05*4,Form("%s",choseDataLegend[choseData]));
      latex1.DrawLatex(0.55,1.-0.05*5.5,Form("%s",choseDataLumi[choseData])); 
    }  

  latex1.SetTextSize(0.035);
  if(choseData!=1 && choseData!=2) latex1.DrawLatex(0.2,1.-0.05*3,Form("Cent. %d-%d%%",centrality_min,centrality_max));
  latex1.DrawLatex(0.2,1.-0.05*4,Form("%.2f< y <%.2f",upsYCut_min,upsYCut_max)); 
  latex1.DrawLatex(0.2,1.-0.05*5.5,Form("p_{T}^{#mu} > %.0f GeV/c",muonpTcut));
  
   c.SaveAs(figs_+figName_+paramOn_+suffix_+".gif");
   c.SaveAs(figs_+figName_+paramOn_+suffix_+".pdf");
  
  //setup the limits
  float baseNll = fit_2nd->minNll();
  
  // Print fit results 
  cout << endl << "figure name: "<< figName_ << endl << endl;
  cout << "the nominal fit with the default pdf " << endl ;
  fit_2nd->Print() ;
  cout << "  free parameter = "<< nfloatpars << ", mychi2 = " << mychsq << ", ndof = " << myndof << endl << endl;
}


//____________________________________________________________________________________
pair<double, double> ConfidencInterval(float CI, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf)  
{  
  //calculate the confidence interval with RooStats
  ProfileLikelihoodCalculator pl(*data,*pdf,*fnll);
  pl.SetConfidenceLevel(CI); // 95% interval
  LikelihoodInterval* interval = pl.GetInterval();
  LikelihoodIntervalPlot plot(interval);
  TCanvas c4; c4.cd(); 
  plot.Draw();
  TString intrvlName = fnll->GetTitle();
  c4.SaveAs(figs_+"nll_"+intrvlName+"_pt4.gif");
  c4.SaveAs(figs_+"nll_"+intrvlName+"_pt4.pdf");
  // print out the iterval on the Parameter of Interest
  cout <<endl<< CI <<"\% interval on " <<fnll->GetName()<<" is : ["<<
    interval->LowerLimit(*fnll) << ", "<<
    interval->UpperLimit(*fnll) << "] "<<endl;
  pair<double, double> CnfdncIntrvl;
  CnfdncIntrvl.first  = interval->LowerLimit(*fnll);
  CnfdncIntrvl.second = interval->UpperLimit(*fnll);
  return CnfdncIntrvl;
}

 
//____________________________________________________________________________________
void setUpLimits(float xMin, float xMax, RooAbsPdf *pdf, RooDataSet *data, RooRealVar *param, float baseNll)
{
  //setting up upper limits
  TCanvas c2; c2.cd();
  int totalBins = 100;
  float BinWidth = (xMax - xMin)/totalBins;
  TH1F *h1 = new TH1F("h1","h1",totalBins,xMin,xMax);
  h1->GetXaxis()->SetTitle(param->getTitle());
  //h1->GetXaxis()->SetRangeUser(-0.1,0.9);
  h1->GetYaxis()->SetTitle("Maximum likelihood");
  TH1F *h2 = new TH1F("h2","h2",totalBins,xMin,xMax);
  gStyle->SetOptStat(kFALSE);
  RooFitResult* fit_nll;
  double MinNll, L, cl;
  for (int i=1; i<=totalBins; i++) 
    {
      param->setVal(xMin + BinWidth*(i-1));
      param->setConstant(kTRUE);
      fit_nll = pdf->fitTo(*data,Save(kTRUE));
      MinNll = fit_nll->minNll()-baseNll;
      L = TMath::Exp(-MinNll);
      cout<<"x = "<<param->getVal()<<", MinNll = "<<MinNll<<", L = "<<L<<endl;
      h1->SetBinContent(i,L);
    }
  for (int i=1; i<=totalBins; i++) 
    {
      cout<<"bin "<<i<<" = "<< h1->GetBinContent(i)<<endl;
    }
  h1->Draw();
  cout<<endl<<"integral = "<<h1->Integral(1,totalBins)<<endl<<endl;
  h1->Scale(1.0/h1->Integral(1,totalBins));
  h1->Draw();
  c2.Update();
  
  //convoluted with systematic gaussian
  TH1F *h3 = new TH1F("h3","h3",totalBins,xMin,xMax);
  TH1F *h4 = new TH1F("h4","h4",totalBins,xMin,xMax);
  for (int i=1; i<=totalBins; i++) {
    float gausmean = xMin + BinWidth*(i-1);
    cout<<"gausmean = "<<gausmean;
    float gaussigma = 0.00001;
    TH1F *h_syst = new TH1F("h_syst","syst histogram",totalBins,xMin,xMax);
    for (Int_t k=0;k<10000;k++) {h_syst->Fill(gRandom->Gaus(gausmean,gaussigma));}
    double new_nll = 0;
    for (int j=1; j<=totalBins; j++){
      new_nll += (h_syst->GetBinContent(j) * h1->GetBinContent(j));
    }
    cout<<", new value after convolution in bin "<<i<<" = "<<new_nll<<endl;
    h3->SetBinContent(i,new_nll);
    delete h_syst;
  }
  cout<<endl<<"integral = "<<h3->Integral(1,totalBins)<<endl<<endl;
  h3->Scale(1.0/h3->Integral(1,totalBins));
  h3->SetLineColor(kMagenta);
  h3->SetLineWidth(2);
  h3->SetLineStyle(2);
  h3->Draw("same");
  
  //find out the upper limit
  c2.Update();
  float UpperLimit;
  double CI = 0.842; //0.842; //0.158;
  for (int i=1; i<=totalBins; i++) 
    {
      cl = h1->Integral(1,i);
      cout<<"x = "<< xMin + BinWidth*(i-1) << ", y = " << h1->GetBinContent(i) <<", cl = "<<cl<<endl;
      if (cl<CI) UpperLimit = xMin + BinWidth*(i-1+1);
      h4->SetBinContent(i,cl);
    }   
  float rightmax = 1.1*h4->GetMaximum();
  float y_scale = gPad->GetUymax()/rightmax;
  cout<<gPad->GetUymax()<<" rightmax = "<<rightmax<<", scale = "<<y_scale<<endl;
  h4->SetLineWidth(2);
  h4->SetLineColor(kRed);
  h4->Scale(y_scale);
  h4->Draw("same");
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed); axis->SetLabelColor(kRed); 
  axis->SetTitle("confidence level"); axis->SetTitleColor(kRed);
  axis->Draw();
  
  TLine L1(xMin,y_scale*CI,xMax,y_scale*CI);
  L1.SetLineColor(kBlue);
  c2.SaveAs(figs_+"UpperLimits_"+paramName+"_hi_syst_Extend1.pdf");
  delete h1; delete h2; delete h3; delete h4;
  param->setConstant(kFALSE);
}


//____________________________________________________________________________________
Double_t RooPlot::mychiSquare(const char* curvename, const char* histname, Int_t nFitParam, bool chisqRange, int startbin, int stopbin) const
{
  // Calculate and return reduced chi-squared of curve with given name with respect
  // to histogram with given name. If nFitParam is non-zero, it is used to reduce the
  // number of degrees of freedom for a chi^2 for a curve that was fitted to the
  // data with that number of floating parameters
  
  // Find curve object
  RooCurve* curve = (RooCurve*) findObject(curvename,RooCurve::Class()) ;
  if (!curve) {
    coutE(InputArguments) << "RooPlot::chiSquare(" << GetName() << ") cannot find curve" << endl ;
    return -1. ;
  }
  
  // Find histogram object
  RooHist* hist = (RooHist*) findObject(histname,RooHist::Class()) ;
  if (!hist) {
    coutE(InputArguments) << "RooPlot::chiSquare(" << GetName() << ") cannot find histogram" << endl ;
    return -1. ;
  }
  return curve->mychiSquare(*hist,nFitParam,chisqRange,startbin,stopbin) ;
}


//____________________________________________________________________________________
Double_t RooCurve::mychiSquare(const RooHist& hist, Int_t nFitParam, bool chisqRange, int startbin, int stopbin) const 
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit
  
  if (chisqRange) {
    int i = startbin;
    int np = stopbin;
  }
  else {
    Int_t i=0;
    int np = hist.GetN() ;
  }
  
  Double_t x,y,eyl,eyh,exl,exh ;
  
  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;
  
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  GetPoint(0,xstart,y) ;
  GetPoint(GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(this)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(this)->GetPoint(GetN()-1,xstop,y) ;
#endif
  
  Int_t nbin(0) ;
  
  Double_t chisq(0) ;
  for (i ; i<np ; i++) 
    {   
      
      // Retrieve histogram contents
      ((RooHist&)hist).GetPoint(i,x,y) ;
      
      // Check if point is in range of curve
      if (x<xstart || x>xstop) continue ;
      
      nbin++ ;
      eyl = hist.GetEYlow()[i] ;
      eyh = hist.GetEYhigh()[i] ;
      exl = hist.GetEXlow()[i] ;
      exh = hist.GetEXhigh()[i] ;
      
      // Integrate function over this bin
      Double_t avg = average(x-exl,x+exh) ;
      
      // Add pull^2 to chisq
      if (y!=0) {      
	Double_t pull = (y>avg) ? ((y-avg)/eyl) : ((y-avg)/eyh) ;
	chisq += pull*pull ;
      }
    }
  
  // Return chisq/nDOF 
  return chisq / (nbin-nFitParam) ;
}


//____________________________________________________________________________________
