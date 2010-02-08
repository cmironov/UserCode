// -*- C++ -*-
//
// Package:    DimuonAnalyzer
// Class:      DimuonAnalyzer
// 
/**\class DimuonAnalyzer DimuonAnalyzer.cc UserCode/DimuonAnalyzer/src/DimuonAnalyzer.cc

 Description: make dimuon from reco::Tracks and recoGenPArticles efficiency analysis

*/
//
// Original Author:  Camelia Mironov,40 1-A32,+41227679747,
//         Created:  Sun Feb  7 15:25:05 CET 2010
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <fstream>

// user include files
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/HeavyIon.h"


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

// data formats
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// root include files
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"


using namespace reco;
using namespace edm;
using namespace std;
using namespace HepMC;

const double m_mu = .105658;

//_______________________________________________________________________
class DimuonAnalyzer : public edm::EDAnalyzer 
{
// class declaration
public:
  explicit DimuonAnalyzer(const edm::ParameterSet&);
  ~DimuonAnalyzer();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  TH2D *phPtY_genDimuonAll; 
  TH1D *phPt_genDimuon;   
  TH1D *phY_genDimuon; 
  TH1D *phM_genDimuon;

  TH2D *phPtY_recoDimuonAll;
  TH1D *phPt_recoDimuon;
  TH1D *phY_recoDimuon;
  TH1D *phM_recoDimuon;

  TH2D *phPtEta_recoTrack;

  TNtuple *pnEventInfo;
 
  edm::InputTag     genparticletag;
  edm::InputTag     muontracktag; 
  edm::InputTag     tracktag;

  // cuts to be set/read in/from the configuration file 
  double                  massMax_dimuon;
  double                  massMin_dimuon;
  int                     pdg_dimuon;
  double                  ptMin_dimuon;
     
  double                  etaMax_muon;
  double                  etaMin_muon;
  double                  ptMin_muon;

  double                  etaMax_track;
  double                  etaMin_track;
  double                  ptMin_track;
  
  double                       b;

  Bool_t acceptDimuon(Float_t pt, Float_t m);
  Bool_t acceptMuon(Float_t pt, Float_t eta);
  Bool_t acceptTrack(Float_t pt, Float_t eta);

};


//_________________________________________________________________
DimuonAnalyzer::DimuonAnalyzer(const edm::ParameterSet& pset):
  genparticletag(pset.getParameter<edm::InputTag>("genParticle") ),
  muontracktag(pset.getUntrackedParameter<edm::InputTag>("muonTracks") ),
  tracktag(pset.getUntrackedParameter<edm::InputTag>("trackTracks") ),
  massMax_dimuon(pset.getParameter<double>("massMaxDimuon") ),
  massMin_dimuon(pset.getParameter<double>("massMinDimuon") ),
  pdg_dimuon(pset.getParameter<double>("pdgDimuon") ),
  ptMin_dimuon(pset.getParameter<double>("ptMinDimuon") ),
  etaMax_muon(pset.getParameter<double>("etaMaxMuon") ),
  etaMin_muon(pset.getParameter<double>("etaMinMuon") ),
  ptMin_muon(pset.getParameter<double>("ptMinMuon") ),
  etaMax_track(pset.getParameter<double>("etaMaxTrack") ),
  etaMin_track(pset.getParameter<double>("etaMinTrack") ),
  ptMin_track(pset.getParameter<double>("ptMinTrack") )
{
  // constructor
 
}


//_________________________________________________________________
DimuonAnalyzer::~DimuonAnalyzer()
{
  // destructor
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//_________________________________________________________________
void DimuonAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  // method called to for each event 
  edm::LogInfo("DimuonAnalyzer")<<"Start analyzing event ...";

  //------- reco::Track loop
  edm::Handle<edm::View<reco::Track> >  trackCollection;
  ev.getByLabel(tracktag, trackCollection); 
  if( trackCollection.isValid() ) 
    {
      for(edm::View<reco::Track> ::size_type i=0; i < trackCollection.product()->size(); i++)
	{
	  edm::RefToBase<reco::Track> recTrack(trackCollection, i);
	  if ( recTrack.isNull() ) continue;
	  if( !( acceptTrack(recTrack->pt(),recTrack->eta()) ) ) continue;
	  phPtEta_recoTrack->Fill(recTrack->pt(),recTrack->eta());
	}
    }
  else LogDebug("DimuonAnalyzer::analyze()") << "##### NO TrackCollection found for reconstructed tracks!";
 

  // ----- reco::Track->MUON loop
  edm::Handle<edm::View<reco::Track> >muonCollection;
  ev.getByLabel(muontracktag,muonCollection);

  if(muonCollection.isValid() && muonCollection.product()->size() >1 )
    {
      for(edm::View<reco::Track> ::size_type i=0; i < muonCollection.product()->size()-1; i++)
	{
	  edm::RefToBase<reco::Track> muTrack1(muonCollection, i);
	  if ( muTrack1.isNull() ) continue;
	  // single muon cuts
	  if( !( acceptMuon(muTrack1->pt(),muTrack1->eta()) ) ) continue;
	  
	  for(edm::View<reco::Track> ::size_type j=i+1; j < muonCollection.product()->size(); j++)
	    {
	      edm::RefToBase<reco::Track> muTrack2(muonCollection, j);
	      if ( muTrack2.isNull() ) continue;
	      // single muon cuts
	      if( !( acceptMuon(muTrack2->pt(),muTrack2->eta()) ) ) continue;
	      
	      // only opposite sign dimuons
	      if(muTrack1->charge() * muTrack2->charge()>0) continue;
	      
	      TLorentzVector kid1, kid2;
	      double en1 = sqrt(muTrack1->p()*muTrack1->p()+m_mu*m_mu);
	      double en2 = sqrt(muTrack2->p()*muTrack2->p()+m_mu*m_mu);
	      
	      kid1.SetPxPyPzE(muTrack1->px(),muTrack1->py(),muTrack1->pz(),en1);
	      kid2.SetPxPyPzE(muTrack2->px(),muTrack2->py(),muTrack2->pz(),en2);
	      
	      TLorentzVector dimuon;
	      dimuon = kid1 + kid2;
	      
	      phPtY_recoDimuonAll->Fill(dimuon.Pt(),dimuon.Rapidity());
	      if ( !(acceptDimuon(dimuon.Pt(),dimuon.M())) ) continue;
	      
	      phPt_recoDimuon->Fill(dimuon.Pt());
	      phY_recoDimuon->Fill(dimuon.Rapidity());
	      phM_recoDimuon->Fill(dimuon.M());
	      
	    }//2nd muon loop
      
	}// 1st muon loop
    }// muonCollectino is valid and has more than 2 muons reconstructed
  else  LogDebug("DimuonAnalyzer::analyze()") << "##### NO TrackCollection found for recoMuons or size<2!";

  //------ recoGenParticle loop
  edm::Handle<reco::GenParticleCollection> genCollection;
  ev.getByLabel(genparticletag,genCollection);
  if( genCollection.isValid() )
    {
      for(reco::GenParticleCollection::size_type i=0; i < genCollection.product()->size(); i++)
	{
	  const GenParticleRef genPart(genCollection,i);
	  if( genPart.isNull() ) continue;
	  
	  if (genPart->pdgId() == pdg_dimuon ) cout <<"\t status = "<< genPart->status()<<"\t pt="<<genPart->pt()<<endl;
	  if( genPart->pdgId() != pdg_dimuon || genPart->status() != 3) continue;
	  phPtY_genDimuonAll->Fill(genPart->pt(),genPart->y());
	  
	  //dimuon gen cuts
	  if( !acceptDimuon(genPart->pt(),genPart->mass()) ) continue;
	  const Candidate *kid0  = genPart->daughter(0);
	  const Candidate *kid1  = genPart->daughter(1);
	  
	  if( !(acceptMuon(kid0->pt(),kid0->eta())) || !(acceptMuon(kid1->pt(),kid1->eta())) ) continue;
	  phPt_genDimuon->Fill(genPart->pt());
	  phY_genDimuon->Fill(genPart->y());
	  phM_genDimuon->Fill(genPart->mass());
	}//genPArticle loop
    }
  else  LogDebug("DimuonAnalyzer::analyze()") << "##### NO genPArticleCollectino found!";

  // ---------- event information:
  edm::Handle<edm::HepMCProduct> hepEv;
  ev.getByLabel("generator",hepEv);

  const HepMC::HeavyIon* hi = hepEv->GetEvent()->heavy_ion();
  if(hi!=NULL)
    b = hi->impact_parameter();
  else
    b = -99.;

  // fill event info tuple
  pnEventInfo->Fill(trackCollection.product()->size(),muonCollection.product()->size(),b);
  

  edm::LogInfo("MCMatchAnalyzer")<<"Finished analyzing event ...";  
}


//_________________________________________________________________
void DimuonAnalyzer::beginJob()
{
  //method called once each job just before starting event loop
  edm::LogInfo("HiEtAnalyzer::beginJob()")<<"Begin initialization in beginJob()";
  edm::Service<TFileService> fs;

  // all dimuons, without any cuts
  phPtY_genDimuonAll  = fs->make<TH2D>("phPtY_genDimuonAll",";p_{T}[GeV/c];rapidity (y)",200,0,100,10,-5.,5.) ;

  // gen dimuons that passed single muon and dimuon cuts
  phPt_genDimuon      = fs->make<TH1D>("phPt_genDimuon",";p_{T}[GeV/c]",200,0,100);           phPt_genDimuon->Sumw2();
  phY_genDimuon       = fs->make<TH1D>("phY_genDimuon",";rapidity (y)",50,-2.5,2.5);          phY_genDimuon->Sumw2();
  phM_genDimuon       = fs->make<TH1D>("phM_genDimuon",";M[GeV/c^2]",400,0.,200);             phM_genDimuon->Sumw2();

  phPtY_recoDimuonAll = fs->make<TH2D>("phPtY_recoDimuonAll",";p_{T}[GeV/c];rapidity (y)",200,0,100,10,-5.,5.) ;
  phPt_recoDimuon     = fs->make<TH1D>("phPt_recoDimuon",";p_{T}[GeV/c]",200,0,100);          phPt_recoDimuon->Sumw2();
  phY_recoDimuon      = fs->make<TH1D>("phY_recoDimuon",";rapidity (y)",50,-2.5,2.5);         phY_recoDimuon->Sumw2();
  phM_recoDimuon      = fs->make<TH1D>("phM_recoDimuon",";M[GeV/c^2]",400,0.,200);            phM_recoDimuon->Sumw2();

  phPtEta_recoTrack   = fs->make<TH2D>("phPtEta_recoTrack",";p_{T}[GeV/c];#eta",200,0,100,50,-2.5,2.5) ;    

  pnEventInfo         = fs->make<TNtuple>("pnEventInfo","pnEventInfo","ntrk:nmu:b");
 
  edm::LogInfo("HiEtAnalyzer::beginJob()")<<"Done beginJob()";
  return ;

}

//_________________________________________________________________
void DimuonAnalyzer::endJob() 
{
  //method called once each job just after ending the event loop

}


//_________________________________________________________________________
Bool_t DimuonAnalyzer::acceptDimuon(Float_t pt, Float_t m)
{
  // cuts on the dimuon; pt and mass
  Bool_t ok = false;
  Bool_t okmass = ( massMin_dimuon<=m ) && ( m<=massMax_dimuon );
  Bool_t okpt = pt>ptMin_dimuon;

  if( okmass && okpt ) ok = true;

  return ok;
}


//_______________________________________________________________________
Bool_t DimuonAnalyzer::acceptMuon(Float_t pt, Float_t eta)
{
  // cuts on single muon: pt and eta
  Bool_t ok = false;
  Bool_t oketa = ( etaMin_muon<=eta ) && ( eta<=etaMax_muon );
  Bool_t okpt = pt>ptMin_muon;
  
  if ( oketa && okpt ) ok = true;

  return ok;
}


//_________________________________________________________________
Bool_t DimuonAnalyzer::acceptTrack(Float_t pt, Float_t eta)
{
  // cuts on reco::Track: pt and eta
  Bool_t ok = false;
  Bool_t oketa = ( etaMin_track<=eta ) && ( eta<=etaMax_track );
  Bool_t okpt = pt>ptMin_track;
  
  if ( oketa && okpt ) ok = true;

  return ok;
}


//_________________________________________________________________
//define this as a plug-in
DEFINE_FWK_MODULE(DimuonAnalyzer);
