/*  Implementation: 
0. Analyzer that gets as input 2 recoTrack collection and their respectiv sim2reco and reco2sim maps
1. makes single and pair track analysis
2. output: TNtuple

*/


// framework includes
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

// data formats
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
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

// sim data formats
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

//tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"

// miscellaneous  
#include <fstream>
using namespace std;
using namespace reco;
const double m_mu    = 0.105658;
const double epsilon = 0.001;
//__________________________________________________________________________
class McMatchTrackAnalyzer : public edm::EDAnalyzer
{

public:
  explicit McMatchTrackAnalyzer(const edm::ParameterSet& pset);
  ~McMatchTrackAnalyzer();
  
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
  virtual void beginJob(const edm::EventSetup& es);
  virtual void endJob();
  
private:
  void fillPairTrackTuple(const edm::RefToBase<Track>& kidRef1,
			  const edm::RefToBase<Track>& kidRef2,
			  vector<float>& result);
  void fillRecoPairTrackTuple(const edm::RefToBase<Track>& trkRef,
			      reco::RecoToSimCollection& reco2Sim,
			      vector<float>& result);
  void fillMatchedRecoTrackTuple(const edm::RefToBase<Track>& trkRef,
				 Int_t& nmatches,
				 Float_t fFrac,
				 vector<float>& result);
  void fillMatchedSimTrackTuple(const TrackingParticleRef& trkRef,
				Int_t& nmatches,
				Float_t fFrac,
				vector<float>& result);
  void fillRecoTrackTuple(const edm::RefToBase<Track>& trkRef,
			  vector<float>& result);
  void fillSimTrackTuple(const TrackingParticleRef& trkRef,
			 vector<float>& result);
  
  
  edm::RefToBase<Track> findRecoTrackMatch(const TrackingParticleRef& trk,
					   reco::SimToRecoCollection& reco2Sim,
					   Int_t& nmatches, 
					   Float_t& fFracShared);
  TrackingParticleRef   findSimTrackMatch(const edm::RefToBase<Track>& trkRef,
					  reco::RecoToSimCollection& reco2Sim,
					  Int_t& nmatches,
					  Float_t& fFracShared);
  
  const reco::Vertex* getClosestVertex(edm::RefToBase<reco::Track>);
  Int_t  getSimParentId(const TrackingParticleRef& trk);
  
  void matchRecoPairTracks(edm::Handle<edm::View<Track> >& trackCollection,
			   reco::RecoToSimCollection& p,	
			   TNtuple* tuple);
  void matchSimPairTracks(edm::Handle<TrackingParticleCollection>& simCollection,	
			  edm::Handle<edm::HepMCProduct>& hepEvt,	 
			  reco::SimToRecoCollection& trks, 
			  TNtuple* tuple); 
  
  void matchRecoTracks(edm::Handle<edm::View<Track> >& trackCollection,
		       reco::RecoToSimCollection& p,
		       TNtuple* tuple);
  void matchSimTracks(edm::Handle<TrackingParticleCollection>& simCollection,
		      reco::SimToRecoCollection& q,
		      TNtuple* tuple); 
  
  void dummVectorEntry(vector<float>& result, Int_t entries);


  // ----- member data -----
  bool                         dohlt;
  bool                         dohiembedding;
  bool                         doparticlegun;
  bool                         doreco2sim;
  bool                         dosim2reco;
  bool                         matchpair;
  bool                         matchsingle;
  Int_t                        proc,procsgn,ntrk,nvtx;
  double                       b;
  Int_t                        npart,ncoll,nhard;
  int                          pdgpair;
  int                          pdgsingle;

  TNtuple                     *pnRecoSimPairType1;
  TNtuple                     *pnRecoSimPairType2;
  TNtuple                     *pnRecoSimTrkType1;
  TNtuple                     *pnRecoSimTrkType2;
  TNtuple                     *pnSimRecoPairType1;
  TNtuple                     *pnSimRecoPairType2;
  TNtuple                     *pnSimRecoTrkType1;
  TNtuple                     *pnSimRecoTrkType2;
  TNtuple                     *pnEventInfo;
  edm::Service<TFileService>   fs;
  const reco::VertexCollection *vertices;

  edm::InputTag                type1trktag;
  edm::InputTag                type1maptag;
  edm::InputTag                type2trktag;
  edm::InputTag                type2maptag;

  edm::InputTag                simtracktag;
  edm::InputTag                trktracktag; 
  edm::InputTag                trkmaptag;
  edm::InputTag                vtxtag; 
};


//___________________________________________________________________________
McMatchTrackAnalyzer::McMatchTrackAnalyzer(const edm::ParameterSet& pset):
  dohlt(pset.getParameter<bool>("doHLT")),
  dohiembedding(pset.getParameter<bool>("doHiEmbedding") ),
  doparticlegun(pset.getParameter<bool>("doParticleGun")),
  doreco2sim(pset.getParameter<bool>("doReco2Sim") ),
  dosim2reco(pset.getParameter<bool>("doSim2Reco") ),
  matchpair(pset.getParameter<bool>("matchPair") ),
  matchsingle(pset.getParameter<bool>("matchSingle") ),
  pdgpair(pset.getParameter<int>("pdgPair") ),
  pdgsingle(pset.getParameter<int>("pdgSingle") ),
  type1trktag(pset.getUntrackedParameter<edm::InputTag>("type1Tracks") ),
  type1maptag(pset.getUntrackedParameter<edm::InputTag>("type1MapTag") ),
  type2trktag(pset.getUntrackedParameter<edm::InputTag>("type2Tracks") ),
  type2maptag(pset.getUntrackedParameter<edm::InputTag>("type2MapTag") ),
  simtracktag(pset.getUntrackedParameter<edm::InputTag>("simTracks") ),
  vtxtag(pset.getUntrackedParameter<edm::InputTag>("vertices") )
{
  // constructor
  //from reco2sim maps
  pnRecoSimTrkType1   = fs->make<TNtuple>("pnRecoSimTrkType1","pnRecoSimTrkType1",
					  "charge:chi2:chi2ndof:dxy:dxyerr:dz:dzerr:nvalidhits:eta:phi:pt:" // reco track
					  "nmatch:frachitmatch:"
					  "chargesim:etasim:idsim:idparentsim:phisim:ptsim"//matched sim sta
					  );//matched sim track

  pnRecoSimTrkType2   = fs->make<TNtuple>("pnRecoSimTrkType2","pnRecoSimTrkType2",
					  "charge:chi2:chi2ndof:dxy:dxyerr:dz:dzerr:nvalidhits:eta:phi:pt:" // reco track
					  "nmatch:frachitmatch:"
					  "chargesim:etasim:idsim:idparentsim:phisim:ptsim"//matched sim sta
					  );//matched sim track

  pnRecoSimPairType1 = fs->make<TNtuple>("pnRecoSimPairType1","pnRecoSimPairType1",
					 "y:minv:phi:pt:" // pair global  --4
					 "charge1:chi21:chi2ndof1:dxy1:dxyerr1:dz1:dzerr1:nvalidhits1:eta1:phi1:pt1:"//kid1 --11
					 "nmatch1:frachitmatch1:"//matched sim kid 1  --8
					 "chargesim1:etasim1:idsim1:idparentsim1:phisim1:ptsim1:"
					 "chage2:chi22:chi2ndof2:dxy2:dxyerr2:dz2:dzerr2:nvalidhits2:eta2:phi2:pt2:"//kid2 --11
					 "nmatch2:frachitmatch2:"// matched sim kid2  --8
					 "chargesim2:etasim2:idsim2:idparentsim2:phisim2:ptsim2"
					 );	
  
  pnRecoSimPairType2 = fs->make<TNtuple>("pnRecoSimPairType2","pnRecoSimPairType2",
					 "y:minv:phi:pt:" // pair global  --4
					 "charge1:chi21:chi2ndof1:dxy1:dxyerr1:dz1:dzerr1:nvalidhits1:eta1:phi1:pt1:"//kid1 --11
					 "nmatch1:frachitmatch1:"//matched sim kid 1  --8
					 "chargesim1:etasim1:idsim1:idparentsim1:phisim1:ptsim1:"
					 "chage2:chi22:chi2ndof2:dxy2:dxyerr2:dz2:dzerr2:nvalidhits2:eta2:phi2:pt2:"//kid2 --11
					 "nmatch2:frachitmatch2:"// matched sim kid2  --8
					 "chargesim2:etasim2:idsim2:idparentsim2:phisim2:ptsim2"
					 );


  // from sim2reco maps
  pnSimRecoTrkType1  = fs->make<TNtuple>("pnSimRecoTrkType1","pnSimRecoTrkType1",
					 "charge:dca:vtxx:vtxy:eta:id:idparent:npixellayers:nsilhits:phi:pt:status:"  //sim track --12
					 "nmatch:frachitmatch:" //# of reco tracks matched to one sim track
					 "chargereco:chi2reco:chi2ndofreco:dxyreco:dxyerrreco:dzreco:dzerrreco:nvalidhitsreco:"
					 "etareco:phireco:ptreco"); // matched reco track
  pnSimRecoTrkType2  = fs->make<TNtuple>("pnSimRecoTrkType2","pnSimRecoTrkType2",
					 "charge:dca:vtxx:vtxy:eta:id:idparent:npixellayers:nsilhits:phi:pt:status:"  //sim track --12
					 "nmatch:frachitmatch:" //# of reco tracks matched to one sim track
					 "chargereco:chi2reco:chi2ndofreco:dxyreco:dxyerrreco:dzreco:dzerrreco:nvalidhitsreco:"
					 "etareco:phireco:ptreco"); // matched reco track
  
  pnSimRecoPairType1 = fs->make<TNtuple>("pnSimRecoPairType1","pnSimRecoPairType1",
					 "y:minv:phi:pt:" // sim Z0 --4 
					 "charge1:dca1:vtxx1:vtxy1:eta1:id1:idparent1:npixellayers1:nsilhits1:phi1:pt1:status1:"  //sim mu1 --12
					 "nmatch1:frachitmatch1:" // reco muon1 -- 13
					 "chargereco1:chi2reco1:chi2ndofreco1:dxyreco1:dxyerrreco1:dzreco1:dzerrreco1:nvalidhitsreco1:etareco1:phireco1:ptreco1:"
					 "charge2:dca2:vtx2:vtxy2:eta2:id2:idparent2:npixellayers2:nsilhits2:phi2:pt2:status2:"  //sim mu2 --12
					 "nmatch2:frachitmatch2:" // reco muon2--13
					 "chargereco2:chi2reco2:chi2ndofreco2:dxyreco2:dxyerrreco2:dzreco2:dzerrreco2:nvalidhitsreco2:etareco2:phireco2:ptreco2:"
					 "yreco:minvreco:phireco:ptreco"// reco  Z0 
					 ); // reco Z0 from trk
  
  pnSimRecoPairType2 = fs->make<TNtuple>("pnSimRecoPairType2","pnSimRecoPairType2",
					 "y:minv:phi:pt:" // sim Z0 --4 
					 "charge1:dca1:vtxx1:vtxy1:eta1:id1:idparent1:npixellayers1:nsilhits1:phi1:pt1:status1:"  //sim mu1 --12
					 "nmatch1:frachitmatch1:" // reco muon1 -- 13
					 "chargereco1:chi2reco1:chi2ndofreco1:dxyreco1:dxyerrreco1:dzreco1:dzerrreco1:nvalidhitsreco1:etareco1:phireco1:ptreco1:"
					 "charge2:dca2:vtx2:vtxy2:eta2:id2:idparent2:npixellayers2:nsilhits2:phi2:pt2:status2:"  //sim mu2 --12
					 "nmatch2:frachitmatch2:" // reco muon2--13
					 "chargereco2:chi2reco2:chi2ndofreco2:dxyreco2:dxyerrreco2:dzreco2:dzerrreco2:nvalidhitsreco2:etareco2:phireco2:ptreco2:"
					 "yreco:minvreco:phireco:ptreco"// reco  Z0 
					 ); // reco Z0 from trk
   
  pnEventInfo      = fs->make<TNtuple>("pnEventInfo","pnEventInfo","proc:procsgn:nsimtrk:ntrk1:ntrk2:nvtx:ncoll:nhard:npart:b");
  
}


//_____________________________________________________________________________
McMatchTrackAnalyzer::~McMatchTrackAnalyzer()
{
  // destructor

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//_____________________________________________________________________
void McMatchTrackAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  // method called each event  
 
  edm::LogInfo("MCMatchAnalyzer")<<"Start analyzing each event ...";

  // ---------Get generated event (signal or background)
  edm::Handle<edm::HepMCProduct> hepEvSgn, hepEv;
  ev.getByLabel("generator",hepEv);
  proc = hepEv->GetEvent()->signal_process_id();
  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<<"Process ID= " <<proc;

  // -------- embedded signal
  procsgn=0;
  if(dohiembedding)
    {
      ev.getByLabel("hiSignal",hepEvSgn);
      procsgn = hepEvSgn->GetEvent()->signal_process_id();
      edm::LogInfo("McMatchTrackAnalyzer::analyze()")<<"Process ID= " << procsgn;
    }

 const HepMC::HeavyIon* hi = hepEv->GetEvent()->heavy_ion();
 if(hi!=NULL)
   {
     ncoll = hi->Ncoll();
     nhard = hi->Ncoll_hard();
     if( hi->Npart_proj() + hi->Npart_targ() > 0)
       {
	 npart =  hi->Npart_proj() + hi->Npart_targ();
	 b = hi->impact_parameter();
       }
   }
 else
   {
     ncoll = 0;
     nhard = 0;
     npart = 0;
     b = -99.;  
   }

  // --------- sim tracks
  edm::Handle<TrackingParticleCollection> simCollection;
  ev.getByLabel(simtracktag,simCollection);
  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<< "Size of simCollection = " << simCollection.product()->size();  

  // ----------- vertex collection
  nvtx = 0;
  edm::Handle<reco::VertexCollection> vertexCollection;
  edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"Getting reco::VertexCollection - "<<vtxtag;
  bool vertexAvailable =  ev.getByLabel(vtxtag,vertexCollection);
  if (vertexAvailable) nvtx = vertexCollection.product()->size();
  else edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"\t VERTEX COLLECTION NOT FOUND.";

  edm::Handle<edm::View<Track> >  type1TrkCollection;
  ev.getByLabel(type1trktag, type1TrkCollection); 
  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<< "Size of type 1 Collection = " << type1TrkCollection.product()->size();
 
  edm::Handle<edm::View<Track> >  type2TrkCollection;
  ev.getByLabel(type2trktag, type2TrkCollection); 
  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<< "Size of type 2 Collection = " << type2TrkCollection.product()->size();
  

  // --------- event info TNtuple 
  Int_t ntrk1 = type1TrkCollection.product()->size();
  Int_t ntrk2 = type2TrkCollection.product()->size();
  Int_t nsimtrk = simCollection.product()->size();
  vector<float> result;
  result.push_back(proc); // proc
  result.push_back(procsgn); // proc signal
  result.push_back(nsimtrk); // ntrkr
  result.push_back(ntrk1); // ntrkr1
  result.push_back(ntrk2); // ntrkr2
  result.push_back(nvtx); // nvtxr
  result.push_back(ncoll);//ncoll
  result.push_back(nhard);//ncoll_hard
  result.push_back(npart);//npart
  result.push_back(b);//b
  pnEventInfo->Fill(&result[0]);
  result.clear();
 
  if(doreco2sim)
    {
      edm::Handle<reco::RecoToSimCollection> type1RecoSimHandle;
      edm::Handle<reco::RecoToSimCollection> type2RecoSimHandle;
      ev.getByLabel(type1maptag,type1RecoSimHandle);
      ev.getByLabel(type2maptag,type2RecoSimHandle);
      
      if(type1RecoSimHandle.isValid() && type2RecoSimHandle.isValid())
	{
	  reco::RecoToSimCollection type1Reco2Sim = *(type1RecoSimHandle.product());
	  reco::RecoToSimCollection type2Reco2Sim = *(type2RecoSimHandle.product());
	  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<<"Got reco2sim maps for type1/type2! Size ="<<type1Reco2Sim.size()<<"\t"<<type2Reco2Sim.size();
	  
	  if(matchsingle) 
	    {
	      edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"##### Matching single typeTracks reco2Sim !!!!! ";
	      matchRecoTracks(type1TrkCollection,type1Reco2Sim,pnRecoSimTrkType1);
	      matchRecoTracks(type2TrkCollection,type2Reco2Sim,pnRecoSimTrkType2);
	    }
	  // pair analysis
	  if( matchpair && type1TrkCollection.product()->size() > 1 )
	    matchRecoPairTracks(type1TrkCollection,type1Reco2Sim,pnRecoSimPairType1);
	  else 
	    edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"##### Size of type 1 Reco2Sim <2! No pairs!";
	  
	  if( matchpair && type2TrkCollection.product()->size() > 1 )
	    matchRecoPairTracks(type2TrkCollection,type2Reco2Sim,pnRecoSimPairType2);
	  else 
	    edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"##### Size of type 2 Reco2Sim <2! No pairs!";
	  
	} // valid tyeps handles
      else edm::LogInfo("McMatchTrackAnalyzer::analyze()") <<"##### No pre-existent reco2sim maps for type1/type2 tracks!!!!! ";
	 
    }//doreco2sim
   
  //-------------- SIM2RECO
  if(dosim2reco)
    {
      edm::Handle<reco::SimToRecoCollection> type1SimRecoHandle;
      edm::Handle<reco::SimToRecoCollection> type2SimRecoHandle;
      ev.getByLabel(type1maptag,type1SimRecoHandle);
      ev.getByLabel(type2maptag,type2SimRecoHandle);
      
      if(type1SimRecoHandle.isValid() && type2SimRecoHandle.isValid())
	{
	  reco::SimToRecoCollection type1Sim2Reco = *(type1SimRecoHandle.product());
	  reco::SimToRecoCollection type2Sim2Reco = *(type2SimRecoHandle.product());
	  edm::LogInfo("McMatchTrackAnalyzer::analyze()")<<"Got sim2reco maps! Size ="<<type1Sim2Reco.size()<<"\t"<<type1Sim2Reco.size();
	  
	  if(matchsingle)
	    {
	      LogDebug("McMatchTrackAnalyzer::analyze()") <<"##### Matching single typeTracks sim2reco !!!!! ";
	      matchSimTracks(simCollection,type1Sim2Reco,pnSimRecoTrkType1);
	      matchSimTracks(simCollection,type2Sim2Reco,pnSimRecoTrkType2);
	    }
	  
	  if(matchpair) 
	    {
	      if(dohiembedding)
		{ 
		  LogDebug("McMatchTrackAnalyzer::analyze()") <<"##### Matching pair typeTracks sim2reco !!!!! ";
		  matchSimPairTracks(simCollection,hepEvSgn,type1Sim2Reco,pnSimRecoPairType1);
		  matchSimPairTracks(simCollection,hepEvSgn,type2Sim2Reco,pnSimRecoPairType2);
		}
	      else 
		{
		  LogDebug("McMatchTrackAnalyzer::analyze()") <<"##### Matching pair typeTracks sim2reco !!!!! ";
		  matchSimPairTracks(simCollection,hepEv,type1Sim2Reco,pnSimRecoPairType1);
		  matchSimPairTracks(simCollection,hepEv,type2Sim2Reco,pnSimRecoPairType2);
		}
	    } 
	}
      else LogDebug("McMatchTrackAnalyzer::analyze()") <<"##### No pre-existent sim2reco maps for type1/type2 tracks!!!!! ";
    } // dosim2reco
      
  edm::LogInfo("McMatchTrackAnalyzer::analyze()") << "[McMatchTrackAnalyzer] done with event# " << ev.id();
}


//________________________________________________________________________
void McMatchTrackAnalyzer::matchRecoTracks(edm::Handle<edm::View<Track> >& trackCollection,
					   reco::RecoToSimCollection& p,
					   TNtuple* tuple)
{
  // fill the TNtuple pnRecoSimTracks, 
  // with reco tracks and their corresponding sim tracks information

  edm::LogInfo("MCMatchAnalyzer::matchRecoTracks()")<<"Start matching reco tracks ...";

  for(edm::View<Track> ::size_type i=0; i < trackCollection.product()->size(); ++i)
  {
    edm::RefToBase<Track> recTrack(trackCollection, i);
    if ( recTrack.isNull() ) continue;
    vector<float> result2;
   
    Int_t nSim  = 0;
    Float_t fFr = 0.;
    // reco tracker
    fillRecoTrackTuple(recTrack,result2);
    TrackingParticleRef matchedSimTrack = findSimTrackMatch(recTrack,p,nSim,fFr);
    fillMatchedSimTrackTuple(matchedSimTrack,nSim,fFr,result2);
      
    tuple->Fill(&result2[0]);
    result2.clear();
  }//i trackCollection
  edm::LogInfo("MCMatchAnalyzer::matchRecoTracks()")<<"Done matchRecoTracks().";
}


//_________________________________________________________________________
void McMatchTrackAnalyzer::matchSimTracks(edm::Handle<TrackingParticleCollection>& simCollection,
					  reco::SimToRecoCollection& q,
					  TNtuple* tuple)
{
  // fill TNtuple pnSimRecoTracks with sim track and corresponding reco tracks info
  // what are the trackingParticles with status = -99? are they from the decay in GEANT of the original particles? yes, they are the particles with no correspondence in the genPArticelCollection ... so if doing specific particle analysis, i should require only the 'signal' muons, not the decay one ...

  edm::LogInfo("MCMatchAnalyzer::matchSimTracks()")<<"Start matching simTracks ...";
 
  for(TrackingParticleCollection::size_type i=0; i < simCollection.product()->size(); i++)
    {
      const TrackingParticleRef simTrack(simCollection,i);
      if( simTrack.isNull() ) continue;

      // only muons
      if( (abs(simTrack->pdgId()) != pdgsingle) ) continue;
      if( simTrack->charge() == 0 ) continue;
  
      vector<float> result4;
      
      // sim track
      fillSimTrackTuple(simTrack,result4); // 12 fields
     
      // recoTracker matched 
      Int_t nRec = 0;
      Float_t fF = 0.;
      edm::RefToBase<Track> matchedRecTrack = findRecoTrackMatch(simTrack,q,nRec,fF);
      fillMatchedRecoTrackTuple(matchedRecTrack,nRec,fF,result4);
            
      // fill
      tuple->Fill(&result4[0]);

      result4.clear();
    }

  edm::LogInfo("MCMatchAnalyzer::matchSimTracks()")<<"Done matchSimTracks()!";
}


//______________________________________________
void McMatchTrackAnalyzer::matchRecoPairTracks(edm::Handle<edm::View<Track> >& trkTypeCollection,
					       reco::RecoToSimCollection& trkReco2Sim, 
					       TNtuple* tuple)
{
  // fill the TNtuple pnRecoSimPairTracks, 
  // with reco and corresponding sim tracks information

  edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"Start reconstructing PairTracks!";

  for(edm::View<Track> ::size_type i1=0; i1 < trkTypeCollection.product()->size()-1; i1++)
    {
      edm::RefToBase<Track> trk1(trkTypeCollection, i1);
      if( trk1.isNull() ) continue;
      for(edm::View<Track> ::size_type i2=i1; i2 < trkTypeCollection.product()->size(); i2++)
      	{
	  edm::RefToBase<Track> trk2(trkTypeCollection, i2);
	  if( trk2.isNull() ) continue;

	  if(trk1->charge()*trk2->charge()>0) continue; // just oppposite sign pairs
	  
	  vector<float> result;
	  // fill pair with global tracks as kids
	  fillPairTrackTuple(trk1,trk2,result);// 4 fields
		  
	  // first muon  
	  fillRecoTrackTuple(trk1,result); // 11 fields
	  Int_t nSim  = 0;
	  Float_t fFr = 0.;
	  // matched sim track
	  TrackingParticleRef matchedSimTrack1 = findSimTrackMatch(trk1,trkReco2Sim,nSim,fFr);
	  fillMatchedSimTrackTuple(matchedSimTrack1,nSim,fFr,result); // 8 fields
	  LogDebug("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"First kid field filled!";

	  // second muon
	  fillRecoTrackTuple(trk2,result); // 11 fields
	  LogDebug("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"Second kid field filled!";
	  TrackingParticleRef matchedSimTrack2 = findSimTrackMatch(trk2,trkReco2Sim,nSim,fFr);
	  fillMatchedSimTrackTuple(matchedSimTrack2,nSim,fFr,result);
	  LogDebug("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"First kid field filled!";

	  tuple->Fill(&result[0]);

	  edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"SIZE = "<<(int)result.size();

	  LogDebug("McMatchTrackAnalyzer::matchRecoPairTracks():result")<<"eta mu1 = "<<result[4]<<"\t eta mu2 = "<< result[13];
	  LogDebug("McMatchTrackAnalyzer::matchRecoPairTracks():direct")<<"eta mu1 = " << trk1->eta() <<"\t ta mu2 = " << trk2->eta();
	  result.clear();
	}// second muon loop
    }// first muon loop
  edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"Done matchRecoPairTracks()!";
}


//_________________________________________________________________________
void McMatchTrackAnalyzer::matchSimPairTracks(edm::Handle<TrackingParticleCollection>& simCollection,
					      edm::Handle<edm::HepMCProduct>& evt,
					      reco::SimToRecoCollection& trkSim2Reco,
					      TNtuple* tuple)
{
  // fill TNtuple pnSimRecoPairTracks with gen info for the gen Z0 and muon daughters and reco info for the reco muons

  // status == 1 : stable particle at pythia level, but which can be decayed by geant
  // status == 2 :
  // status == 3 : outgoing partons from hard scattering

  // The barcode is the particle's reference number, every vertex in the
  // event has a unique barcode. Particle barcodes are positive numbers,
  // vertex barcodes are negative numbers.

  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"Start matching pairs ...";
  Int_t statuscheck = 3;
  if(doparticlegun) statuscheck = 2;

  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  Bool_t foundmuons = false;
  vector<Int_t> ixmu;
  vector<float> result3;
  edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"SIZE_1_result3 = "<<result3.size();

  //====================================
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) 
    { 
      if ( ( abs((*p)->pdg_id()) == pdgpair ) && (*p)->status() == statuscheck ) 
	{ 
	  for( HepMC::GenVertex::particle_iterator aDaughter=(*p)->end_vertex()->particles_begin(HepMC::descendants); 
	       aDaughter !=(*p)->end_vertex()->particles_end(HepMC::descendants);
	       aDaughter++)
	    {
	      if ( abs((*aDaughter)->pdg_id())==pdgsingle && (*aDaughter)->status()==1 )
		{
		  ixmu.push_back((*aDaughter)->barcode());
		  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()") << "Stable muon from Z0" << "\tindex= "<< (*aDaughter)->barcode(); 
		  foundmuons = true;
		}//  z0 descendentes
	    }
	  if (foundmuons)
	    {
	      double rapidity = 0.5 * log(( (*p)->momentum().e()+(*p)->momentum().pz() )/( (*p)->momentum().e()-(*p)->momentum().pz() ));
	      result3.push_back(rapidity);
	      result3.push_back((*p)->momentum().m());
	      result3.push_back((*p)->momentum().phi());
	      result3.push_back((*p)->momentum().perp());
	    }
	}// Z0
    }// loop over genParticles
   //==================================
  edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"SIZE_2_result3 = "<<result3.size();
  if( foundmuons )
    {
      vector<edm::RefToBase<Track> > recokids;
      Int_t nrecokids    = 0;
      Int_t nsimtracks   = 0;
      for(TrackingParticleCollection::size_type i=0; i < simCollection.product()->size(); i++)
	{
	  const TrackingParticleRef simTrk(simCollection, i);
	  if( simTrk.isNull() || simTrk->status()!=1 ) continue;
	  const SimTrack *gTrack = &(*simTrk->g4Track_begin());
	  if( gTrack == NULL ) continue;

	  if(abs(gTrack->type())==pdgsingle && (gTrack->genpartIndex()==ixmu[0] || gTrack->genpartIndex()==ixmu[1] ))
	    {
	      edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"Got an xmuon!.";
	      nsimtracks++;
	      // sim track
	      fillSimTrackTuple(simTrk,result3); //12 fields
	   
	      // trackerMatched_muon		 
	      Int_t nTrk = 0;
	      Float_t fTrkHit = 0.; 
	      edm::RefToBase<Track> matchedRecTrk = findRecoTrackMatch(simTrk,trkSim2Reco,nTrk,fTrkHit);
	      fillMatchedRecoTrackTuple(matchedRecTrk,nTrk,fTrkHit,result3);  // 13 fields
	       edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"SIZE_3_result3 = "<<result3.size();
	      if(nTrk!=0)
		{
		  recokids.push_back(matchedRecTrk);
		  nrecokids++;
		}
	    }// muon id
	}// TrackCollection loop
      //---------------
      // have the sim kids in recokids vector
      edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"Found kids: "<<nrecokids<<"SIZE_4_result3 = "<<result3.size();
      if(nsimtracks==0) 
	{
	  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"No sim found for the gen kids";
	  dummVectorEntry(result3,54); //(2*25(sim+reco)+4 recoZ0)
	}
      if(nsimtracks==1)
	{
	  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"Only one sim out of the 2 kids were found";
	  dummVectorEntry(result3,29); //(1*25+4 recoZ0)
	  	  
	}
      if(nrecokids == 0 || nrecokids==1) // no or just one kid reconstructed
	{
	  LogDebug("MCMatchAnalyzer::matchSimPairTracks()")<<"No sim daughter reconstructed";
	  dummVectorEntry(result3,4); //(4 recoZ0)
	}
     
      // get the reco Z0 from global muons
      if(nrecokids == 2)
	{
	  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"########## Fill reco_pair and reco_kids info";
	  edm::RefToBase<Track> ch0 = recokids[0];
	  edm::RefToBase<Track> ch1 = recokids[1];
	  
	  // reconstructed Z0
	  fillPairTrackTuple(ch0,ch1,result3); // 4fields
	}

      edm::LogInfo("McMatchTrackAnalyzer::matchRecoPairTracks()")<<"SIZE_7_result3 = "<<result3.size();

      tuple->Fill(&result3[0]);
      recokids.clear();
      result3.clear();
    }//foundmuons
 
  result3.clear();
 

  edm::LogInfo("MCMatchAnalyzer::matchSimPairTracks()")<<"Done matchSimPairTracks()!.";
}


//_________________________________________________________________________
void  McMatchTrackAnalyzer::fillPairTrackTuple(const edm::RefToBase<Track>& kidRef1,
					       const edm::RefToBase<Track>& kidRef2,
					       vector<float>& result)
{
  // reconstruct the parent dilepton composed object: y,minv,phi,pT

  TLorentzVector child1, child2;
  double en1 = sqrt(kidRef1->p()*kidRef1->p()+m_mu*m_mu);
  double en2 = sqrt(kidRef2->p()*kidRef2->p()+m_mu*m_mu);
  
  child1.SetPxPyPzE(kidRef1->px(),kidRef1->py(),kidRef1->pz(),en1);
  child2.SetPxPyPzE(kidRef2->px(),kidRef2->py(),kidRef2->pz(),en2);
  
  TLorentzVector pair;
  pair = child1 + child2;
	 	  

  result.push_back(pair.Rapidity());
  result.push_back(pair.M());
  result.push_back(pair.Phi());
  result.push_back(pair.Pt());

  edm::LogInfo("MCMatchAnalyzer::fillPairTrackTuple()")<<"Done fillPairTrackTuple()!";
}
//_______________________________________________________________________
void McMatchTrackAnalyzer::fillRecoPairTrackTuple(const edm::RefToBase<Track>& trkRef,
						  reco::RecoToSimCollection& reco2Sim,
						  vector<float>& result)
{
  // get the match, 
  // fill the result with the id of match, id of parent, and and number of matches
  edm::LogInfo("MCMatchAnalyzer::fillRecoPairTrackTuple()")<<"Filling pair tuple ...";

  if(!trkRef.isNull())
    {
      Int_t nTrkSim = 0;
      Float_t fFraction = 0.;
      TrackingParticleRef matchedSim = findSimTrackMatch(trkRef,reco2Sim,nTrkSim,fFraction); 
      Int_t nId = -99;
      Int_t nIdParent = -99;
      
      if( !matchedSim.isNull() && nTrkSim!=0 )
	{
	  nId = matchedSim->pdgId();
	  nIdParent = getSimParentId(matchedSim);
	}
      
      result.push_back(nTrkSim);
      result.push_back(fFraction);
      result.push_back(nId);
      result.push_back(nIdParent);
    }
  else dummVectorEntry(result,4);
  edm::LogInfo("MCMatchAnalyzer::fillRecoPairTrackTuple()")<<"Done fillRecoPairTrackTuple()!";
 }


//____________________________________________________________________
void McMatchTrackAnalyzer::fillMatchedRecoTrackTuple(const edm::RefToBase<Track>& recoTrack,
						Int_t& nMatches, Float_t fFracMatch,vector<float>& result)
{
  // the info part of the matched sim track
  edm::LogInfo("MCMatchAnalyzer::fillMatchedTrackTuple()")<<"Filling matched reco track info ...";

  if(!recoTrack.isNull() && nMatches > 0)
    {
     
      double dt      = -99.;
      double sigmaDt = -99.;
      double dz      = -99.;
      double sigmaDz = -99. ;
   
      result.push_back(nMatches);            // number of reco tracks matched to one sim track (split tracks)
      result.push_back(fFracMatch);          // fraction of hits matched, for the last matche found, the one kept
      result.push_back(recoTrack->charge()); // chargereco
      result.push_back(recoTrack->chi2());                 // normalized chi2
      result.push_back(recoTrack->normalizedChi2());       // n && getClosestVertex(recoTrack)!=NULLormalized Chi2
      result.push_back(dt);                                // transverse DCA
      result.push_back(sigmaDt);                           // Sigma_Dca_t
      result.push_back(dz);                                // z DCA
      result.push_back(sigmaDz);                           // Sigma_Dca_z
      result.push_back(recoTrack->numberOfValidHits());    // number of hits found 
      result.push_back(recoTrack->eta());    // eta
      result.push_back(recoTrack->phi());    // phireco
      result.push_back(recoTrack->pt());     // ptreco
      
      edm::LogInfo("MCMatchAnalyzer::fillMatchedSimTrackTuple()")<<"pt= "<<recoTrack->pt();
    }
  else dummVectorEntry(result,13);
  edm::LogInfo("MCMatchAnalyzer::fillMatchedTrackTuple()")<<"Done fillMatchedTrackTuple()!"; 
}


//____________________________________________________________________
void McMatchTrackAnalyzer::fillMatchedSimTrackTuple(const TrackingParticleRef& simTrack,
					       Int_t& nMatches, Float_t fFrac,vector<float>& result)
{
  // the info part of the matched sim track
  edm::LogInfo("MCMatchAnalyzer::fillMatchedTrackTuple()")<<"Filling matched sim track info ...";

  if( !simTrack.isNull() && nMatches>0)
    {
      
      Int_t nId       = simTrack->pdgId();
      Int_t nIdParent = getSimParentId(simTrack);
      edm::LogInfo("MCMatchAnalyzer::fillMatchedSimTrackTuple()")<<"pt= "<<simTrack->pt()<<"\t particle id = "<<nId<<"\t parent id = "<<nIdParent<<"\tstatus= "<<simTrack->status();

      result.push_back(nMatches);                 // number of sim tracks matched to one reco track
      result.push_back(fFrac);                    // fraction of hits matched, for the last matche found, the one kept
      result.push_back(simTrack->charge());       // charge sim
      result.push_back(simTrack->eta());          // eta sim
      result.push_back(nId);                      // id
      result.push_back(nIdParent);                // id parent
      result.push_back(simTrack->phi());          // phisim
      result.push_back(simTrack->pt());           // pt sim
    }
  else dummVectorEntry(result,8);
  edm::LogInfo("MCMatchAnalyzer::fillMatchedSimTrackTuple()")<<"Done fillMatchedTrackTuple()!"; 
}


//_______________________________________________________________________
void McMatchTrackAnalyzer::fillRecoTrackTuple(const edm::RefToBase<Track>& trkRef, vector<float>& result)
{
  // fill reco track info to be fed later to the tuple
  edm::LogInfo("MCMatchAnalyzer::fillRecoTrackTuple()")<<"Filling reco track info ...";
 
   
  double dt       = -99.;
  double sigmaDt  = -99.;
  double dz       = -99.;
  double sigmaDz  = -99.;
  edm::LogInfo("MCMatchAnalyzer::fillRecoTrackTuple()")<<"Gate1 ...";
 
   result.push_back(trkRef->charge());                         // charge
   result.push_back(trkRef->chi2());                           // normalized chi2
   result.push_back(trkRef->normalizedChi2());                 // normalized Chi2
   result.push_back(dt);                                       // transverse DCA
   result.push_back(sigmaDt);                                  // Sigma_Dca_t
   result.push_back(dz);                                       // z DCA
   result.push_back(sigmaDz);                                  // Sigma_Dca_z
   result.push_back(trkRef->numberOfValidHits());              // number of hits found 
   result.push_back(trkRef->eta());                            // eta
   result.push_back(trkRef->phi());                            // phi
   result.push_back(trkRef->pt());                             // pt
   
   edm::LogInfo("MCMatchAnalyzer::fillRecoTrackTuple()")<<"Done fillRecoTrackTuple()!";
}


//____________________________________________________________________
void McMatchTrackAnalyzer::fillSimTrackTuple(const TrackingParticleRef& simTrack,
					     vector<float>& result)
{
  // the info part of the matched sim track
  edm::LogInfo("MCMatchAnalyzer::fillMatchedTrackTuple()")<<"Filling matched reco track info ...";

  Int_t nParentId    = getSimParentId(simTrack);
  result.push_back(simTrack->charge());  // chargereco
  result.push_back(simTrack->vertex().Rho()); // dca
  result.push_back(simTrack->vertex().x()); // x of the vtx
  result.push_back(simTrack->vertex().y()); // y of the vtx
  result.push_back(simTrack->eta());     // eta
  result.push_back(simTrack->pdgId());   // pdg id
  result.push_back(nParentId);           // parent id
  result.push_back(0);        // number of pixel layers with hits
  result.push_back(0);        // total number of sim hits
  result.push_back(simTrack->phi());     // phi
  result.push_back(simTrack->pt());      // pt
  result.push_back(simTrack->status());  // status
 
  edm::LogInfo("MCMatchAnalyzer::fillMatchedTrackTuple()")<<"Done fillMatchedTrackTuple()!"; 
}


//________________________________________________________________________
edm::RefToBase<Track> McMatchTrackAnalyzer::findRecoTrackMatch(const TrackingParticleRef& simTrack,
							       reco::SimToRecoCollection& sim2Reco, 
							       Int_t& nMatches, Float_t& fFracShared)
{
  // return the matched reco track
  edm::LogInfo("MCMatchAnalyzer::findRecoTrackMatch()")<<"Finding reco track match ...";
  edm::LogInfo("McMatchTrackAnalyzer::findRecoTrackMatch()")<<" for simTrack with ... Pdg_ID = "<<simTrack->pdgId()<<"\t status = "<<simTrack->status()<<"\t pt = "<<simTrack->pt();
 
  edm::RefToBase<Track> recoTrackMatch;
  vector<pair<edm::RefToBase<Track>, double> > recTracks;
  if(sim2Reco.find(simTrack) != sim2Reco.end()) 	 
    {
      recTracks = sim2Reco[simTrack];
      if(recTracks.size() !=0 )
	{
	  recoTrackMatch = recTracks.begin()->first;
	  fFracShared    = recTracks.begin()->second;
	  nMatches       = recTracks.size();
	}
    }
  edm::LogInfo("MCMatchAnalyzer::findRecoTrackMatch()")<<"Done findRecoTrackMatch(). Found "<<nMatches<<"\t matches";
  return recoTrackMatch;
}


//________________________________________________________________________
TrackingParticleRef McMatchTrackAnalyzer::findSimTrackMatch(const edm::RefToBase<Track>& recoTrack,
							    reco::RecoToSimCollection& reco2Sim, 
							    Int_t& nMatches,Float_t& fFracShared)
{
  // return the matched sim track
  edm::LogInfo("MCMatchAnalyzer::findSimTrackMatch()")<<"Finding sim track match ...";

  TrackingParticleRef simTrackMatch;
  vector<pair<TrackingParticleRef, double> > simTracks;

  if(reco2Sim.find(recoTrack) != reco2Sim.end())  
    {
      simTracks = reco2Sim[recoTrack];
      if(simTracks.size() !=0 )
	{
	  simTrackMatch = simTracks.begin()->first;
	  fFracShared   = simTracks.begin()->second;
	  nMatches      = simTracks.size();
	}
    }

  edm::LogInfo("MCMatchAnalyzer::findSimTrackMatch()")<<"Done finding sim track match! Found "<<nMatches<<"matches";

  return simTrackMatch;
}


//___________________________________________________________________________
const reco::Vertex* McMatchTrackAnalyzer::getClosestVertex(edm::RefToBase<reco::Track> recTrack)
{
  // get the xloseset event vertex to the track
  const reco::Vertex *closestVertex=0;
  if(vertices->size()>0 && !dohlt)
    {
      //look for th elcosest vertex in z

      Float_t dzmin = -1;
      for(reco::VertexCollection::const_iterator vertex = vertices->begin(); vertex != vertices->end(); vertex++)
	{
	  Float_t dz = fabs(recTrack->vertex().z() - vertex->position().z());
	  if(vertex == vertices->begin() || dz < dzmin)
	    {
	      dzmin = dz;
	      closestVertex = &(*vertex);
	    }
	}
    }
  else edm::LogInfo("MCMatchAnalyzer::getClossestVertex()")<<"Null vertex in the event";

  return closestVertex;

}


//____________________________________________________________________________
Int_t McMatchTrackAnalyzer::getSimParentId(const TrackingParticleRef& match)
{
  // return the particle ID of associated GEN track (trackingparticle = gen + sim(geant) )
  // it is a final, stable particle (status=1): after final state radiaiton

  // same particle before final state radiation (status=3); from this one get to the real parent;
  // this is the 'parent' of the final state particle, p; same pdg_id
  
  Int_t parentId = -99; 
  
  if( match.isNull() ) return parentId;
   
  edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<"Getting the parent Id for the sim particle with status..."<< match->status();
  TrackingParticle::genp_iterator b;
  TrackingParticle::genp_iterator in  = match->genParticle_begin();
  TrackingParticle::genp_iterator fin = match->genParticle_end();
  for (b = in; b != fin; ++b)
    {
      const HepMC::GenParticle *p = b->get();
      if( p==NULL ) 
	{
	  edm::LogInfo("MCMatchAnalyzer::getSimParentId()")<<"No gen particle associated with simTrack with status " << match->status();
	  continue;
	}
  
      edm::LogInfo("MCMatchAnalyzer::getSimParentId()")<<"Gen particle is : " << (*b)->pdg_id();
      const HepMC::GenParticle *mother = p->production_vertex() && 
	p->production_vertex()->particles_in_const_begin() != p->production_vertex()->particles_in_const_end() ?
	*(p->production_vertex()->particles_in_const_begin()) : 0;
    
      if( mother!=0 &&
	  !(isnan(mother->momentum().mag())) && !isnan(abs(mother->pdg_id())) )
	{
	  
	  if( abs(mother->pdg_id())<10 || (mother->pdg_id()<=100 && mother->pdg_id()>=80) || mother->pdg_id()==21)
	    { 
	      edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<"Id_parent = 0; comes from PV directly ";
	      return 0; // from PV (parton or cluster)
	    }
	  else
	    {
	      parentId = mother->pdg_id();
	      if( parentId != p->pdg_id()) 
		{
		  edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<"Id_parent = "<<parentId;
		  return parentId;
		}
	      edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<" 1 parentId = "<<parentId<<"\t parent_status= "<<mother->status();
	  
	      // real parent after final state radiation
	      const HepMC::GenParticle *motherMother = mother->production_vertex() &&
		mother->production_vertex()->particles_in_const_begin() != mother->production_vertex()->particles_in_const_end() ?
		*(mother->production_vertex()->particles_in_const_begin()) : 0 ;
	      
	      if( motherMother!=0  
		  && !(isnan(motherMother->momentum().perp()) ) && !(isnan(motherMother->pdg_id()) ) )
		{ 
		  if(abs(motherMother->pdg_id())<10 || 
		     (motherMother->pdg_id()<=100 && motherMother->pdg_id()>=80) || 
		     motherMother->pdg_id()==21)
		    {
		      edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<"Id_parent = 0; comes from PV indirectly ";
		      return 0;
		    }
		  else
		    {
		      parentId = motherMother->pdg_id();
		      edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<" 2 parentId = "<<parentId<<"\t parent_status= "<<motherMother->status();
		    }
		  motherMother = 0;
		}// valid motherMother
	    }//else: mother not from PV
	}//valid mother 
     mother = 0;
    }//loop over tracking particles

  edm::LogInfo("McMatchTrackAnalyzer::getSimParentId")<<"Done with getSimParentId!\t"<<parentId;
  return parentId;
}


//_______________________________________________________________________
void McMatchTrackAnalyzer::dummVectorEntry(vector<float>& result, Int_t nEntries)
{
  // add to vector nentries of '-99' 
  for(Int_t i = 1; i<= nEntries; i++)
    result.push_back(-99);
}


//____________________________________________________________________________
void McMatchTrackAnalyzer::beginJob(const edm::EventSetup& es)
{
 // method called once each job just before starting event loop
 
}


//_________________________________________________________________________
void McMatchTrackAnalyzer::endJob()
{
  //method called once each job just after ending the event loop 

}


//________________________________________________________________________

DEFINE_FWK_MODULE(McMatchTrackAnalyzer);

