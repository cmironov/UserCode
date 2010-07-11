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
// $Id: DimuonAnalyzer.cc,v 1.2 2010/02/08 17:20:39 mironov Exp $
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
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

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
  TH2D *phPtPhi_recoTrack;

  // By moon
  TH1D *nTMuTrk;

  TH1D *phPt_globalMu;
  TH1D *phEta_globalMu;
  TH1D *phPhi_globalMu;
  TH1D *nGLB;

  TH1D *phPt_staMu;
  TH1D *phEta_staMu;
  TH1D *phPhi_staMu;
  TH2D *phPtPhi_staMu;
  TH1D *nSTA;

  TH1D *phPt_otherMu;
  TH1D *phEta_otherMu;
  TH1D *phPhi_otherMu;
  TH1D *nOTH;
  
  TH2D *phPtEta_recoTrack;

  TNtuple *pnEventInfo;
 
  edm::InputTag     genparticletag;
  edm::InputTag     muontracktag; 
  edm::InputTag     muontag; 
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
  
  Bool_t doMC;
  Bool_t doSingMu;

  const  CentralityBins *cbins_;
  double bmean, bsigma;
  double npartmean, npartsigma;
  double ncollmean, ncollsigma;
  double hf;
  int    bin;

  //for vertex
  const reco::VertexCollection *vertices;
  edm::InputTag                vtxtag;
  int                          nvtx;
  double                       vx,vy,vz;
  
  TNtuple                      *pnSTAmuInfo;
  TNtuple                      *pnGLBmuInfo;
  TNtuple                      *pnDimuInfo;

};


//_________________________________________________________________
DimuonAnalyzer::DimuonAnalyzer(const edm::ParameterSet& pset):
  genparticletag(pset.getParameter<edm::InputTag>("genParticle") ),
  muontracktag(pset.getUntrackedParameter<edm::InputTag>("muonTracks") ),
  muontag(pset.getUntrackedParameter<edm::InputTag>("muons") ),
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
  ptMin_track(pset.getParameter<double>("ptMinTrack") ),
  doMC(pset.getParameter<bool>("doMC") ),
  doSingMu(pset.getParameter<bool>("doSingleMuon")),
  cbins_(0),
  vtxtag(pset.getUntrackedParameter<edm::InputTag>("vertices") )
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
  using namespace std;
  using namespace edm;

  // method called to for each event 
  edm::LogInfo("DimuonAnalyzer")<<"Start analyzing event ...";

  // Get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  ev.getByLabel(vtxtag,vertexCollection);
  vertices = vertexCollection.product();
  edm::LogInfo("DimuonAnalyzer::analyze()")<<"Number of vertices in the event = "<< vertexCollection.product()->size();
  bool vertexAvailable =  ev.getByLabel(vtxtag,vertexCollection);
  if (!vertexAvailable)
  {
    ev.getByLabel("pixelVertices",vertexCollection);
    edm::LogInfo("DimuonAnalyzer::analyze()")<<"Using the pp vertex collection: pixelVertices ";
  }
  nvtx = vertices->size();
  vx = vertices->begin()->position().x();
  vy = vertices->begin()->position().y();
  vz = vertices->begin()->position().z();


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
      phPtPhi_recoTrack->Fill(recTrack->pt(),recTrack->phi());
	  }
  }
  else LogDebug("DimuonAnalyzer::analyze()") << "##### NO TrackCollection found for reconstructed tracks!";
 
  // ---- Centrality
  cbins_ = 0;
  if (!cbins_) cbins_ = getCentralityBinsFromDB(iSetup);
  edm::Handle<reco::Centrality> cent;
  ev.getByLabel(InputTag("hiCentrality"),cent);
  hf = cent->EtHFhitSum();
  bin = cbins_->getBin(hf);
  bmean = cbins_->bMeanOfBin(bin);
  bsigma = cbins_->bSigmaOfBin(bin);
  npartmean = cbins_->NpartMeanOfBin(bin);
  npartsigma = cbins_->NpartSigmaOfBin(bin);
  ncollmean = cbins_->NcollMeanOfBin(bin);
  ncollsigma = cbins_->NcollSigmaOfBin(bin);

  // ----- reco::Track->MUON loop
  // for all muons info by Moon
  if (doSingMu)
  {
    //const reco::Vertex *vtx= vertices[0];
    const reco::Vertex *vtx = &(vertices->at(0));
    //if (!vtx) cout <<"NO VERTEX in the EVENTS" << endl;

   // edm::Handle<CandidateView> muons;
    edm::Handle<edm::View<reco::Muon> >muons;
    ev.getByLabel(muontag,muons);
    int nTMu = 0; int ngMu = 0; int nsMu = 0; int noMu = 0;
    for(unsigned int j = 0; j < muons->size(); ++j){
     //   CandidateBaseRef muCandRef = muons->refAt(j);
        edm::RefToBase<reco::Muon> muCandRef(muons, j);
        if ( muCandRef.isNull() ) continue;
        
        if(muCandRef->isGlobalMuon()){
            phPt_globalMu->Fill(muCandRef->globalTrack()->pt());
            phEta_globalMu->Fill(muCandRef->globalTrack()->eta());
            phPhi_globalMu->Fill(muCandRef->globalTrack()->phi());
            double dxy      = muCandRef->globalTrack()->dxy(vtx->position());
            double sigmaDxy = sqrt(muCandRef->globalTrack()->dxyError()*muCandRef->globalTrack()->dxyError() + vtx->yError()*vtx->yError()+vtx->xError()*vtx->xError());
            double dz      = muCandRef->globalTrack()->dz(vtx->position());
            double sigmaDz = sqrt(muCandRef->globalTrack()->dzError()*muCandRef->globalTrack()->dzError()+vtx->zError()*vtx->zError());
            pnGLBmuInfo->Fill(muCandRef->globalTrack()->pt(),muCandRef->globalTrack()->eta(),muCandRef->globalTrack()->phi(),dxy,dz,sigmaDxy,sigmaDz); 
            ngMu++;
        }
        if(muCandRef->isStandAloneMuon()){
            phPt_staMu->Fill(muCandRef->standAloneMuon()->pt());
            phEta_staMu->Fill(muCandRef->standAloneMuon()->eta());
            phPhi_staMu->Fill(muCandRef->standAloneMuon()->phi());
            phPtPhi_staMu->Fill(muCandRef->standAloneMuon()->pt(),muCandRef->standAloneMuon()->phi());
            double dxy      = muCandRef->standAloneMuon()->dxy(vtx->position());
            double sigmaDxy = sqrt(muCandRef->standAloneMuon()->dxyError()*muCandRef->standAloneMuon()->dxyError() + vtx->yError()*vtx->yError()+vtx->xError()*vtx->xError());
            double dz      = muCandRef->standAloneMuon()->dz(vtx->position());
            double sigmaDz = sqrt(muCandRef->standAloneMuon()->dzError()*muCandRef->standAloneMuon()->dzError()+vtx->zError()*vtx->zError());
            pnSTAmuInfo->Fill(muCandRef->standAloneMuon()->pt(),muCandRef->standAloneMuon()->eta(),muCandRef->standAloneMuon()->phi(),dxy,dz,sigmaDxy,sigmaDz);
            nsMu++;
        }
        if(muCandRef->isTrackerMuon()){
            phPt_otherMu->Fill(muCandRef->track()->pt());
            phEta_otherMu->Fill(muCandRef->track()->eta());
            phPhi_otherMu->Fill(muCandRef->track()->phi());
            noMu++;
        }
        nTMu++;
    }
    nTMuTrk->Fill(nTMu);nGLB->Fill(ngMu);nSTA->Fill(nsMu);nOTH->Fill(noMu);
  }

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
        pnDimuInfo->Fill(dimuon.Pt(),dimuon.Rapidity(),dimuon.M());
	      
        if ( !(acceptDimuon(dimuon.Pt(),dimuon.M())) ) continue;

	      LogInfo("DimuonAnalyzer:RecoMuonLoop:") << "RecoKid0: pt="<<kid1.Pt()<<"\t eta = "<< kid1.Eta();
	      LogInfo("DimuonAnalyzer:RecoMuonLoop")  << "RecoKid1: pt="<<kid2.Pt()<<"\t eta = "<< kid2.Eta();	   
	      
	      phPt_recoDimuon->Fill(dimuon.Pt());
	      phY_recoDimuon->Fill(dimuon.Rapidity());
	      phM_recoDimuon->Fill(dimuon.M());

	    }//2nd muon loop
  	}// 1st muon loop
  }// muonCollectino is valid and has more than 2 muons reconstructed
  else  LogDebug("DimuonAnalyzer::analyze()") << "##### NO TrackCollection found for recoMuons or size<2!";
 
  if(doMC){
      //------ recoGenParticle loop
      edm::Handle<reco::GenParticleCollection> genCollection;
      ev.getByLabel(genparticletag,genCollection);
      if( genCollection.isValid() )
      {
          for(reco::GenParticleCollection::size_type i=0; i < genCollection.product()->size(); i++)
          {
              const GenParticleRef genPart(genCollection,i);
              if( genPart.isNull() ) continue;

              if (genPart->pdgId() == pdg_dimuon ) 
                  LogInfo("DimuonAnalyzer")<<"Parent: status = "<< genPart->status()<<"\t pt ="<<genPart->pt()<<"\t y = "<< genPart->y()<<endl;

              if( genPart->pdgId() != pdg_dimuon || genPart->status() != 2) continue;
              phPtY_genDimuonAll->Fill(genPart->pt(),genPart->y());

              //dimuon gen cuts
              if( !acceptDimuon(genPart->pt(),genPart->mass()) ) continue;
              const Candidate *kid0  = genPart->daughter(0);
              const Candidate *kid1  = genPart->daughter(1);

              if( !(acceptMuon(kid0->pt(),kid0->eta())) || !(acceptMuon(kid1->pt(),kid1->eta())) ) continue;

              LogInfo("DimuonAnalyzer:GenParticleLoop:")<<"GenKid0: status = "<< kid0->status()<<"\t pt="<<kid0->pt()<<"\t eta = "<< kid0->eta();
              LogInfo("DimuonAnalyzer:GenParticleLoop")<<"GenKid1: status = "<< kid1->status()<<"\t pt="<<kid1->pt()<<"\t eta = "<< kid1->eta();

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
  }
  // fill event info tuple
  pnEventInfo->Fill(trackCollection.product()->size(),muonCollection.product()->size(),bmean,bsigma,npartmean,npartsigma,ncollmean,ncollsigma,bin,hf,nvtx,vx,vy,vz);

  edm::LogInfo("DimuonAnalyzer")<<"Finished analyzing event ...";  
}

//_________________________________________________________________
void DimuonAnalyzer::beginJob()
{
  //method called once each job just before starting event loop
  edm::LogInfo("DimuonAnalyzer::beginJob()")<<"Begin initialization in beginJob()";
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

  phPt_globalMu       = fs->make<TH1D>("phPt_globalMu",";p_{T}[GeV/c]",200,0,100);            phPt_globalMu->Sumw2();
  phPhi_globalMu      = fs->make<TH1D>("phPhi_globalMu",";#phi [rad]",100,-4,4);              phPhi_globalMu->Sumw2();
  phEta_globalMu      = fs->make<TH1D>("phEta_globalMu",";#eta [rad]",100,-4,4);              phEta_globalMu->Sumw2();
  
  phPt_staMu          = fs->make<TH1D>("phPt_staMu",";p_{T}[GeV/c]",200,0,100);               phPt_staMu->Sumw2();
  phPhi_staMu         = fs->make<TH1D>("phPhi_staMu",";#phi [rad]",100,-4,4);                 phPhi_staMu->Sumw2();
  phEta_staMu         = fs->make<TH1D>("phEta_staMu",";#eta [rad]",100,-4,4);                 phEta_staMu->Sumw2();
  phPtPhi_staMu       = fs->make<TH2D>("phPtPhi_staMu",";p_{T};#phi",100,0.,100,100,-4.,4.);   phPtPhi_staMu->Sumw2();

  phPt_otherMu        = fs->make<TH1D>("phPt_otherMu",";p_{T}[GeV/c]",200,0,100);             phPt_otherMu->Sumw2();
  phPhi_otherMu       = fs->make<TH1D>("phPhi_otherMu",";#phi [rad]",100,-4,4);               phPhi_otherMu->Sumw2();
  phEta_otherMu       = fs->make<TH1D>("phEta_otherMu",";#eta [rad]",100,-4,4);               phEta_otherMu->Sumw2();
  
  nTMuTrk             = fs->make<TH1D>("nTMuTrk",";number of all Muons",100,0,5);
  nGLB                = fs->make<TH1D>("nGLB",";number of global Muons",100,0,5);
  nSTA                = fs->make<TH1D>("nSTA",";number of sta Muons",100,0,5);
  nOTH                = fs->make<TH1D>("nOTH",";number of other Muons",100,0,5);

  phPtEta_recoTrack   = fs->make<TH2D>("phPtEta_recoTrack",";p_{T}[GeV/c];#eta",200,0,100,50,-2.5,2.5) ;    
  phPtPhi_recoTrack   = fs->make<TH2D>("phPtPhi_recoTrack",";p_{T}[GeV/c];#phi",200,0,100,100,-4.,4.) ;    

  pnEventInfo         = fs->make<TNtuple>("pnEventInfo","pnEventInfo","ntrk:nmu:bmean:bsigma:npartmean:npartsigma:ncollmean:ncollsigma:bin:hf:nvtx:vx:vy:vz");
  pnSTAmuInfo         = fs->make<TNtuple>("pnSTAmuInfo","pnSTAmuInfo","pt:eta:phi:dxy:dz:sigmaDxy:sigmaDz");
  pnGLBmuInfo         = fs->make<TNtuple>("pnGLBmuInfo","pnGLBmuInfo","pt:eta:phi:dxy:dz:sigmaDxy:sigmaDz");
  pnDimuInfo          = fs->make<TNtuple>("pnDimuInfo","pnDimuInfo","pt:rapidity:mass");
 
  edm::LogInfo("DimuonAnalyzer::beginJob()")<<"Done beginJob()";
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
