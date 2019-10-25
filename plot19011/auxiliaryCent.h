#ifndef AUXILIARYCENT
#define AUXILIARYCENT

#include "TColor.h"

//binning
const unsigned int nBinsLow   = 2;
const unsigned int nBinsHigh  = 1;

double boxWidth=5;
// arrays of yields, errors, etc
// Labeling
  const char* bType[2] = {"B^{+}","B_{s}^{0}"};

  const char* xAxName[] = {"<N_{part}>"};
  const char* yAxName[] = {"1/T_{AA}dN/dp_{T} (pb c/GeV)", "Yields ratio"};


// numbers
// center of bins
double binLow[nBinsLow];
double binHigh[nBinsHigh];
  
// yield
double bs_low[nBinsLow], bs_high[nBinsHigh];
double bpl_low[nBinsLow],bpl_high[nBinsHigh];
  
// bins width along x
double bs_low_xErrL[nBinsLow], bs_low_xErrH[nBinsHigh];
double bpl_low_xErrL[nBinsLow],bpl_low_xErrH[nBinsHigh];
  
    double bs_high_xErrL[nBinsHigh], bs_high_xErrH[nBinsHigh];
    double bpl_high_xErrL[nBinsHigh],bpl_high_xErrH[nBinsHigh];
    
// stat uncert
    double bs_low_yStatL[nBinsLow],   bs_low_yStatH[nBinsHigh];
    double bs_high_yStatL[nBinsHigh], bs_high_yStatH[nBinsHigh];
    
    double bpl_low_yStatL[nBinsLow],  bpl_low_yStatH[nBinsHigh];
    double bpl_high_yStatL[nBinsHigh],bpl_high_yStatH[nBinsHigh];
    
// systematic uncertainties
    double bs_low_ySystL[nBinsLow],  bs_low_ySystH[nBinsHigh];
    double bs_high_ySystL[nBinsHigh],bs_high_ySystH[nBinsHigh];
    
    double bpl_low_ySystL[nBinsLow],  bpl_low_ySystH[nBinsHigh];
    double bpl_high_ySystL[nBinsHigh],bpl_high_ySystH[nBinsHigh];

// global uncertainty
    double glbSystUp   = 0.;
    double glbSystDown = 0.;


//-------------------------------------------------------------------
// ***** //Drawing
    
//++++++++++++++ coloring scheme
// first Bs, then Bpl
Int_t colorLow[2]   = {TColor::GetColor("#4F42B5"), TColor::GetColor("#299617")};
Int_t colorHigh[2]  = {TColor::GetColor("#4F42B5"), TColor::GetColor("#299617")};

Int_t markerHigh[2]  = {24, 25};
Int_t markerLow[2] = {20, 21};

double markerSizeLow[2]  = {1.2, 1.2};
double markerSizeHigh[2] = {1.2, 1.2};

//first is cent, 2nd is integrated
Int_t colorRatio[2] = {TColor::GetColor("#C32148"),TColor::GetColor("#C32148")};

Int_t markerRatio[2]     = {30, 29};
double markerSizeRatio[2] = {1.9, 1.9};


//---------------------------- text

double ltxSetTextSize1 = 0.06;//legend entries
double ltxSetTextSize2 = 0.04;//legend entry

// xsec specific
double legXsec_xLowStart= 0.8;
double legXsec_y= 0.8;
double legXsec_xLowEnd = 0.9;

double xsec_ltxText1_xStart=0.2; // kinematic
double xsec_ltxText1_yStart=0.83;

// ratio specifics
double ratio_ltxText1_xStart=0.22; // ratio mesons
double ratio_ltxText1_yStart=0.8;

double ratio_ltxText2_xStart=0.7; //
double ratio_ltxText2_yStart=0.85;

double legRatio_xLowStart= 0.65;
double legRatio_y= 0.65;
double legRatio_xLowEnd = 0.85;

double legRatioRef_xLowStart= 0.35;
double legRatioRef_y        = 0.65;
double legRatioRef_xLowEnd  = 0.55;
#endif
