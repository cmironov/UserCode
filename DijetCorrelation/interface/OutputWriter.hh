#ifndef  __OUTPUTWRITER_H__
#define  __OUTPUTWRITER_H__

#include "Correlation/DijetCorrelation/interface/Selector.hh"

namespace CLHEP {
  class RandFlat;
}


namespace LanlJetCorrelation{

  class OutputWriter {

  public:
    OutputWriter(Selector *selectorp);
    ~OutputWriter(){;}
    
    void CreateGeneric();
    void CreateAzimuth();
    void CreateParentRecon();
    
    void FillGlobal(float cent, float zvert);
    void FillSingleHistos(std::deque<GenericParticle*> PartList, int cBin, corrType type);
    void FillParentNtuple(fillType type, int cBin, std::deque<GenericParticle*> ParentList);
    void FillAzimuth(fillType typ, int cBin, GenericParticle*, GenericParticle*);
    
    float DeltaPhi(float phi1, float phi2); // REMMEMBER TO USE THE CMS RANDOM GENERATOR
    void  PrintStats();
    
  protected:
    Selector *selector;
    TString  hname, htit;
    //   RandomEngine* randEngine;
    CLHEP::RandFlat*        fRandomGenerator; 

    int nReal[MAXCENTRBIN], nMix[MAXCENTRBIN];
    int nTrigg[MAXCENTRBIN], nAssoc[MAXCENTRBIN];

    // Generic histograms: store chosen binning for analysis macro automatization
    TH1D *hBinsCent, *hBinsTrigg, *hBinsAssoc;
    
    // QA Histograms:
    TH1D *hCentr, *hZVert, *hDphiMixQA;
    
    // Reconstructed Parent Histograms: Invariant Mass/pT Spectra
    TNtuple *ParentNtu[2];
    
    // Correlations Histograms
    TH2D *hCorr[2][MAXCENTRBIN][MAXTRIGGBIN][MAXASSOCBIN];  
    TH1D *hTriggPt[MAXCENTRBIN][MAXTRIGGBIN], *hAssocPt[MAXCENTRBIN][MAXASSOCBIN];
    TH2D *hTriggAcceptance[MAXCENTRBIN][MAXTRIGGBIN], *hAssocAcceptance[MAXCENTRBIN][MAXASSOCBIN];
    
  };
} // namespace

#endif /*  __OUTPUTWRITER_H__  */
