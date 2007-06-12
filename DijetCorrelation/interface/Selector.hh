#ifndef __SELECTOR_H__
#define __SELECTOR_H__

// CMS include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
// input CMS containers
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
// correlation containers
#include "Correlation/DijetCorrelation/interface/GenericParticle.hh"
#include "Correlation/DijetCorrelation/interface/RecoParent.hh"
#include "Correlation/DijetCorrelation/interface/Hadron.hh"
#include "Correlation/DijetCorrelation/interface/Muon.hh"

namespace LanlJetCorrelation{
  
  class Selector{
    
  public:
    Selector(TString filename); 
    ~Selector(){;}
    void Print();
    
    //=== Build the deques of good particles ===
    std::deque<GenericParticle*> GetAssocList(edm::Handle<reco::TrackCollection> &inputlist);
    std::deque<GenericParticle*> GetTriggList(edm::Handle<reco::TrackCollection> &inputlist);
    std::deque<GenericParticle*> DecayRecon(std::deque<GenericParticle*> childlist);
    
    //=== Generic selections ============
    bool MixTwoHadrons()        {return ((int)CorrelType==0);}
    bool MixPi0Hadrons()        {return ((int)CorrelType==1);}
    bool MixPi0Muons()          {return ((int)CorrelType==2);}
    bool MixDielectronHadrons() {return ((int)CorrelType==3);}
    bool MixDimuonHadrons()     {return ((int)CorrelType==4);}
    bool UseHadrons()    {return (MixTwoHadrons()||MixPi0Hadrons()||MixDielectronHadrons()||MixDimuonHadrons());}
    bool UsePizeroes()   {return (MixPi0Hadrons()||MixPi0Muons());}   bool UsePhotons()   {return UsePizeroes();}
    bool UseDielectrons(){return MixDielectronHadrons();}    bool UseElectrons() {return MixDielectronHadrons();}
    bool UseDimuons()    {return MixDimuonHadrons();}               bool UseMuons()       {return MixPi0Muons();}
    bool SimilarCentrality(float cent1, float cent2) { return fabs(cent1-cent2)<maxDCent; }
    bool SimilarVertZ(float z1, float z2) { return fabs(z1-z2)<maxVertDZ; }
    float GetVtxRange() {return maxVertZRange;}

    // Hadron selections:
    float GetMinHadronPt() {return float(MinHadronPt);}
    bool TrackPairCuts(GenericParticle *Trigger, GenericParticle *Associated);

    // Dimuon selections:
    float GetMinMuonPt() {return float(MinMuonPt);}
    bool GoodParent(RecoParent* piZero, GenericParticle *p1, GenericParticle *p2);
    bool GoodMuon();
    bool InParentMassAndPtBin(RecoParent *parent);
   
    //=== Binning methods =========
    int GetNoOfBins (corrType ctype){return NoOfBins[ctype]-1;}  
    float GetBinBorder(corrType ctype, int ii);
    int GetBin(corrType ctype, float val);  
    
  protected:
    //=== Read the input card ======
    void DefineContent();
    void ReadIn();
    void AddLine(char *name, unsigned int addr, int *dim);
    ifstream list;
    TString fname;
    TString *var[MAXDIM];
    int readIn[MAXDIM];
    unsigned int storage[MAXDIM];
    int *storeDimP[MAXDIM], storeDim[MAXDIM], nullp;
    int counter;
    
    // Generic selections:
    float CorrelType, maxVertZRange, maxVertDZ, maxDCent;
    static const int mdim=20;
    int NoOfBins[4];    float Borders[4][mdim];
    
    // Hadron selection:
    float MinHadronPt;

    // Dimuon selections:
    float MinMuonPt;
  };
} // namespace
#endif 
