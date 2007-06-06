// -*- C++ -*-
//
// Package:    DijetCorrelation
// Class:      DijetCorrelation
// 
/**\class DijetCorrelation DijetCorrelation.cc Correlation/DijetCorrelation/src/DijetCorrelation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Camelia Mironov
//         Created:  Tue May 22 17:49:41 EDT 2007
// $Id$
//

#include "Correlation/DijetCorrelation/interface/Selector.hh"
#include "Correlation/DijetCorrelation/interface/OutputWriter.hh"
#include "Correlation/DijetCorrelation/interface/Mixer.hh"

namespace LanlJetCorrelation{

  class DijetCorrelation : public edm::EDAnalyzer {
  public:
    explicit DijetCorrelation(const edm::ParameterSet&);
    ~DijetCorrelation();
    
  private:
    virtual void beginJob(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void CrossCorrelate(std::deque<GenericParticle*>,std::deque<GenericParticle*>,int);
   
    int nTotalEvts;
    TH1D *DPhi;
    TFile *outFile;
    TString outFileName;
    
    std::string cardname_, ofilename_;
    
    std::deque<GenericParticle*> TriggList, AssocList;
    edm::InputTag trackTags_; //used to select what tracks to read from configuration file  
    Selector *selector;
    OutputWriter *outputer;
    Mixer *mixer;
  };
} // namespace
