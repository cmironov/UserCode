#ifndef DATA_2013_REF
#define DATA_2013_REF

// ################################ reference points
const int nRefBin = 1;
double refBins[]      = {0.5};
double refBins7[]     = {0.7};
double refBinsError[] = {0.10};

// trigger= DoubleMuOpen 2013 , pol2 with fixed n-power
double pp276_yPM1_2s1s[]           = {0.355};
double pp276_yPM1_2s1s_statError[] = {0.022};
double pp276_yPM1_2s1s_systError[] = {0.005};

double pp276_yPM1_3s1s[]           = {0.169};
double pp276_yPM1_3s1s_statError[] = {0.016};
double pp276_yPM1_3s1s_systError[] = {0.005};

// tracker tracker
// trigger= Dimuon0 
double pp7_yPM1_2s1s[]           = {0.327};
double pp7_yPM1_2s1s_statError[] = {0.009};
double pp7_yPM1_2s1s_systError[] = {0.002};

double pp7_yPM1_3s1s[]           = {0.193};
double pp7_yPM1_3s1s_statError[] = {0.007};
double pp7_yPM1_3s1s_systError[] = {0.002};


//trigger= DoubleMu3
// systm is difference between nomFit for Dimuon0 and 'best' fit for this case (free param all)
double dblmu3_pp7_yPM1_2s1s[]           = {0.314};
double dblmu3_pp7_yPM1_2s1s_statError[] = {0.006};
double dblmu3_pp7_yPM1_2s1s_systError[] = {0.004};

double dblmu3_pp7_yPM1_3s1s[]           = {0.192};
double dblmu3_pp7_yPM1_3s1s_statError[] = {0.005};
double dblmu3_pp7_yPM1_3s1s_systError[] = {0.004};


// --------------------------------------------------
// ###### PbPb
double refBins_pbpb_2s1s[] = {0.5};
double refBins_pbpb_3s1s[] = {0.6};

// 50-100% pp reco
double ppreco_yPM1_cent50100_2s1s[]           = {0.102};
double ppreco_yPM1_cent50100_2s1s_statError[] = {0.077};
double ppreco_yPM1_cent50100_2s1s_systError[] = {0.027};

// with 2011 nominal fit settings (erf*exp, alpha=1.4, npow=2.3,sigma=92MeV)
double gg_zhen2011_yPM1_cent50100_2s1s[]           = {0.143};
double gg_zhen2011_yPM1_cent50100_2s1s_statError[] = {0.077};
double gg_zhen2011_yPM1_cent50100_2s1s_systError[] = {0.01}; // difference with AN

// #### regit
double gg_regit_yPM1_cent50100_2s1s[]           = {0.195};
double gg_regit_yPM1_cent50100_2s1s_statError[] = {0.088};
double gg_regit_yPM1_cent50100_2s1s_systError[] = {0.011};

double gg_regit_yPM1_cent50100_3s1s[]           = {0.137};
double gg_regit_yPM1_cent50100_3s1s_statError[] = {0.077};
double gg_regit_yPM1_cent50100_3s1s_systError[] = {0.011};

// global reco
double ppreco_gg_ypm24_cent50100_2s1s[]           = {0.153};
double ppreco_gg_ypm24_cent50100_2s1s_statError[] = {0.074};
double ppreco_gg_ypm24_cent50100_2s1s_systError[] = {0};

double ppreco_gg_ypm1_cent50100_2s1s[]           = {0.17};
double ppreco_gg_ypm1_cent50100_2s1s_statError[] = {0.11};
double ppreco_gg_ypm1_cent50100_2s1s_systError[] = {0};

// trk-trk
double ppreco_ypm24_cent50100_2s1s[]           = {0.116};
double ppreco_ypm24_cent50100_2s1s_statError[] = {0.067};
double ppreco_ypm24_cent50100_2s1s_systError[] = {0};

double ppreco_ypm1_cent50100_2s1s[]           = {0.135};
double ppreco_ypm1_cent50100_2s1s_statError[] = {0.080};
double ppreco_ypm1_cent50100_2s1s_systError[] = {0};

// zhen’s tree: vtxProb 0.01, new fit
// 50-100% 
double zhen_gg_ypm24_cent50100_2s1s[]           = {0.147};
double zhen_gg_ypm24_cent50100_2s1s_statError[] = {0.081};
double zhen_gg_ypm24_cent50100_2s1s_systError[] = {0};

double zhen_gg_ypm1_cent50100_2s1s[]           = {0.178};
double zhen_gg_ypm1_cent50100_2s1s_statError[] = {0.092};
double zhen_gg_ypm1_cent50100_2s1s_systError[] = {0};

// regit tree
// 50-100% 
double regit_gg_ypm24_cent50100_2s1s[]           = {0.119};
double regit_gg_ypm24_cent50100_2s1s_statError[] = {0.065};
double regit_gg_ypm24_cent50100_2s1s_systError[] = {0};

double regit_gg_ypm1_cent50100_2s1s[]           = {0.138};
double regit_gg_ypm1_cent50100_2s1s_statError[] = {0.08};
double regit_gg_ypm1_cent50100_2s1s_systError[] = {0};

// published
// 50-100% 
double pub_gg_ypm24_cent50100_2s1s[]           = {0.152};
double pub_gg_ypm24_cent50100_2s1s_statError[] = {0.077};
double pub_gg_ypm24_cent50100_2s1s_systError[] = {0};


#endif
