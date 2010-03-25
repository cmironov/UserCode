/* Analyzer to plot reco2sim and sim2reco info:
-- efficiency
-- purities
---etc

-- it takes as input a file obtianed with McMatchAnalyzer, that containes 2 TNtuples: 
pnRecoSimMuons and pnSimRecoMuons
-- it also assumes that you got a bunch of files from the anlyzer, which you MERGED afterwords, so 
to this plotter comes just one file!
-- please ask, in order to provide the merger- macro
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include "TDirectory.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3F.h>
#include "TKey.h"
#include <TLeaf.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TLorentzVector.h>

// miscellaneous  
#include <fstream>
#include <map>
#include <string>

#include "headers/bins.h"
#include "headers/histos.h"

using namespace std;


//__________________________________________________________________________
class McMatchTrackPlotter : public edm::EDAnalyzer
{

public:
  explicit McMatchTrackPlotter(const edm::ParameterSet& pset);
  ~McMatchTrackPlotter();
  
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
  virtual void beginJob(const edm::EventSetup& es);
  virtual void endJob();
  

private:
  bool mergedrootfiles;
  bool dobarel;
  bool doendcap;
  bool doparentspecificcut;
  bool dopair;
  bool dosingle;
  std::string infilelabel;
  std::string typelabel;
  Double_t hitfraccut;
  Int_t  pdgpair;
  double masscutmax;
  double masscutmin;
  Double_t ptcutpair;
  Double_t ptcutsingle;
  bool reco2sim;
  bool sim2reco;
  
  TFile *pfInFile;
   
  //____________________________________
  void initHistograms();
 
  void analyzeReco2SimPair(TNtuple *pnt);
  void analyzeSim2RecoPair(TNtuple *pnt);

  void analyzeSim2RecoSingle(TNtuple *pnt);
  void analyzeReco2SimSingle(TNtuple *pnt);

  Bool_t acceptPair(Float_t pt, Float_t m);
  Bool_t acceptEta(Float_t e);
  Bool_t acceptSingle(Float_t pt, Float_t eta);
  Bool_t acceptMass(Float_t m);
  Bool_t acceptParent(Float_t id);
};


//___________________________________________________________________________
McMatchTrackPlotter::McMatchTrackPlotter(const edm::ParameterSet& pset):
  mergedrootfiles(pset.getUntrackedParameter<bool>("mergedRootFiles") ),
  dobarel(pset.getUntrackedParameter<bool>("doBarel") ),
  doendcap(pset.getUntrackedParameter<bool>("doEndcap") ),
  doparentspecificcut(pset.getUntrackedParameter<bool>("doParentSpecificCut") ),
  dopair(pset.getUntrackedParameter<bool>("doPair") ),
  dosingle(pset.getUntrackedParameter<bool>("doSingle") ),
  infilelabel(pset.getUntrackedParameter<string>("inFileName") ),
  typelabel(pset.getUntrackedParameter<string>("typeNumber") ),
  hitfraccut(pset.getUntrackedParameter<double>("hitFractionMatch") ),
  pdgpair(pset.getUntrackedParameter<int>("pdgPair") ),
  masscutmax(pset.getUntrackedParameter<double>("massCutMax") ),
  masscutmin(pset.getUntrackedParameter<double>("massCutMin") ),
  ptcutpair(pset.getUntrackedParameter<double>("ptCutPair") ),
  ptcutsingle(pset.getUntrackedParameter<double>("ptCutSingle") ),
  reco2sim(pset.getUntrackedParameter<bool>("reco2Sim") ),
  sim2reco(pset.getUntrackedParameter<bool>("sim2Reco") )
 
{
  // constructor
  cout<<"Input file name " << infilelabel.c_str()<<endl;
  pfInFile = new TFile(infilelabel.c_str(),"READ");
  if( pfInFile == NULL || !pfInFile->IsOpen() ) edm::LogWarning("Plotter")<<" Could not open file!!!";
  else initHistograms();
}


//_____________________________________________________________________________
McMatchTrackPlotter::~McMatchTrackPlotter()
{
  // destructor

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  pfInFile = 0;
}


//_____________________________________________________________________
void McMatchTrackPlotter::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
  // method called each event
  return;
}


//____________________________________________________________________
void McMatchTrackPlotter::analyzeReco2SimPair(TNtuple *pnTuple1)
{
// reco2sim analysis: purity,resolution, x_dimuReco,x_dimuRecoSim,x_dimuRecoFake, x_muRecoMulti
  edm::LogInfo("analyzeReco2SimPair")<<"Analyzing pair reco2sim tuple ...";

 //### reco dimuon
  TLeaf *plY                 = (TLeaf*)pnTuple1->GetLeaf("y");
  TLeaf *plMinv              = (TLeaf*)pnTuple1->GetLeaf("minv");
  TLeaf *plPt                = (TLeaf*)pnTuple1->GetLeaf("pt");
  
  //### kid1 1
  TLeaf *plEta_1           = (TLeaf*)pnTuple1->GetLeaf("eta1");
  TLeaf *plPt_1            = (TLeaf*)pnTuple1->GetLeaf("pt1");
  TLeaf *plNmatch_1       = (TLeaf*)pnTuple1->GetLeaf("nmatch1");
  TLeaf *plFrac_1          = (TLeaf*)pnTuple1->GetLeaf("frachitmatch1");


  TLeaf *plId_sim1        = (TLeaf*)pnTuple1->GetLeaf("idsim1");
  TLeaf *plIdParent_sim1  = (TLeaf*)pnTuple1->GetLeaf("idparentsim1");

  //### kid 2
  TLeaf *plEta_2           = (TLeaf*)pnTuple1->GetLeaf("eta2");
  TLeaf *plPt_2            = (TLeaf*)pnTuple1->GetLeaf("pt2");
  TLeaf *plNmatch_2       = (TLeaf*)pnTuple1->GetLeaf("nmatch2");
  TLeaf *plFrac_2          = (TLeaf*)pnTuple1->GetLeaf("frachitmatch2");

  TLeaf *plId_sim2        = (TLeaf*)pnTuple1->GetLeaf("idsim2");
  TLeaf *plIdParent_sim2  = (TLeaf*)pnTuple1->GetLeaf("idparentsim2");

  Int_t nEvents = pnTuple1->GetEntries();
  edm::LogInfo("analyzeReco2SimPair")<<"Processing "<<nEvents<<" events...";

  for(Long_t ie=0;ie<nEvents;ie++)
    {
      pnTuple1->GetEntry(ie);
      Double_t eta1 = plEta_1->GetValue();
      Double_t eta2 = plEta_2->GetValue();
      Double_t pt1 = plPt_1->GetValue();
      Double_t pt2 = plPt_2->GetValue();

      LogDebug("analyzeReco2SimPair")<<" Gate 1"<<"eta_1= "<<eta1<<"\t pt1= "<<pt1;
      LogDebug("analyzeReco2SimPair")<<" Gate 1"<<"eta_2= "<<eta2<<"\t pt2= "<<pt2;
      LogDebug("analyzeReco2SimPair")<<" Gate 1"<<"pt = "<<plPt->GetValue()<<"\t M= "<< plMinv->GetValue() ;

      if( !acceptPair(plPt->GetValue(),plMinv->GetValue()) || 
	  !acceptSingle(pt1,eta1) || 
	  !acceptSingle(pt2,eta2) ) continue;  // kinematic cuts

      // all reco pairs
      phPtY_pairReco->Fill(plPt->GetValue(),plY->GetValue());
      phPtMinv_pairReco->Fill(plPt->GetValue(),plMinv->GetValue());
   
      bool kid1Matched = ( plNmatch_1->GetValue()>0 && plFrac_1->GetValue()>hitfraccut &&
			   (abs(plId_sim1->GetValue())==13 && acceptParent(plIdParent_sim1->GetValue()) ) 
			   );
      bool kid2Matched = ( plNmatch_2->GetValue()>0  && plFrac_2->GetValue()>hitfraccut &&
			   (abs(plId_sim2->GetValue())==13 && acceptParent(plIdParent_sim2->GetValue()) ) 
			   );
     
      if(kid1Matched && kid2Matched)
	{
	  phPtY_pairRecoSim->Fill(plPt->GetValue(),plY->GetValue());
	  phPtMinv_pairRecoSim->Fill(plPt->GetValue(),plMinv->GetValue());
	}
      else 
	{
	  phPtY_pairRecoFake->Fill(plPt->GetValue(),plY->GetValue());
	  phPtMinv_pairRecoFake->Fill(plPt->GetValue(),plMinv->GetValue());
	}
    }
 }


//____________________________________________________________________
void McMatchTrackPlotter::analyzeSim2RecoPair(TNtuple *pnTuple4)
{
  // reco2sim analysis: purity,resolution, x_pairReco,x_pairRecoSim,x_pairRecoFake, x_RecoMulti
  edm::LogInfo("analyzeSim2RecoPair")<<"Analyzing dimuon sim2reco tuple ...";

 //### sim dimuon
  TLeaf *plY                 = (TLeaf*)pnTuple4->GetLeaf("y");
  TLeaf *plMinv              = (TLeaf*)pnTuple4->GetLeaf("minv");
  TLeaf *plPt                = (TLeaf*)pnTuple4->GetLeaf("pt");

  //### reco dimuon
  TLeaf *plY_reco            = (TLeaf*)pnTuple4->GetLeaf("yreco");
  TLeaf *plMinv_reco         = (TLeaf*)pnTuple4->GetLeaf("minvreco");
  TLeaf *plPt_reco           = (TLeaf*)pnTuple4->GetLeaf("ptreco");

  //### sim kid 1 1
  TLeaf *plEta_1           = (TLeaf*)pnTuple4->GetLeaf("eta1");
  TLeaf *plPt_1            = (TLeaf*)pnTuple4->GetLeaf("pt1");
  TLeaf *plNmatch_1       = (TLeaf*)pnTuple4->GetLeaf("nmatch1");
  TLeaf *plFrac_1          = (TLeaf*)pnTuple4->GetLeaf("frachitmatch1");
  //### sim kid 2
  TLeaf *plEta_2           = (TLeaf*)pnTuple4->GetLeaf("eta2");
  TLeaf *plPt_2            = (TLeaf*)pnTuple4->GetLeaf("pt2");
  TLeaf *plNmatch_2       = (TLeaf*)pnTuple4->GetLeaf("nmatch2");
  TLeaf *plFrac_2          = (TLeaf*)pnTuple4->GetLeaf("frachitmatch2");

  
  Int_t nEvents = pnTuple4->GetEntries();
  edm::LogInfo("analyzeSim2RecoPair")<<"Processing "<<nEvents<<" events...";

  for(Long_t ie=0;ie<nEvents;ie++)
    {
      pnTuple4->GetEntry(ie);

      Float_t y    = plY->GetValue();
      Float_t pt   = plPt->GetValue();

      Double_t eta1 = plEta_1->GetValue();
      Double_t eta2 = plEta_2->GetValue();

      //id+process cuts: done in the match track analyzer; on the sim2reco pair collection are only kept the muons coming from a specific process                                    
      if(!acceptEta(eta1) || !acceptEta(eta2) ) continue;   //-- look for sim muons only in a part of the det
      
      // all sim pairs
      phPtY_pairSim->Fill(pt,y); 
      phPtMinv_pairSim->Fill(pt,plMinv->GetValue());
 
      // the cut done on the recoDimuon; pt and mass
      // in this way, when you cut on these at the reco level, you do it to kick out bkg; but you do lose also signal;
      // this reflects in the efficiency plot
      if(!acceptPair(plPt_reco->GetValue(),plMinv_reco->GetValue()) ) continue;


      bool kid1Matched = (plNmatch_1->GetValue() > 0 && plFrac_1->GetValue()>hitfraccut) ;
      bool kid2Matched = (plNmatch_2->GetValue() > 0 && plFrac_2->GetValue()>hitfraccut);

      if(kid1Matched && kid2Matched) 
	{
	  Float_t y_r    = plY_reco->GetValue();
	  Float_t pt_r   = plPt_reco->GetValue();

	  Float_t res_y   = (y_r-y)/y;
	  Float_t res_pt  = (pt_r-pt)/pt;
	  Float_t res_m   = (plMinv_reco->GetValue()- plMinv->GetValue())/plMinv->GetValue();
	  // eta, phi, pt resolution vs pt
	  phPtY_pairResolution->Fill(pt,res_y); 
	  phPtPt_pairResolution->Fill(pt,res_pt);
	  phPtMinv_pairResolution->Fill(pt,res_m);

	  phPtY_pairSimReco->Fill(pt,y);
	  phPtMinv_pairSimReco->Fill(pt,plMinv->GetValue());
	} 
      else // fake/semifake/lost dimuon 
	{
	  phPtY_pairSimLost->Fill(pt,y);
	  phPtMinv_pairSimLost->Fill(pt,plMinv->GetValue());
	}
    }
}


//____________________________________________________________________
void McMatchTrackPlotter::analyzeReco2SimSingle(TNtuple *pnTuple3)
{
  // reco2sim analysis: purity,resolution
   edm::LogInfo("analyzeReco2SimMuons")<<"Analyzing reco2 sim track tuple ...";

  TLeaf *plCh           = (TLeaf*)pnTuple3->GetLeaf("charge");
  TLeaf *plEta          = (TLeaf*)pnTuple3->GetLeaf("eta");
  TLeaf *plPt           = (TLeaf*)pnTuple3->GetLeaf("pt");

  TLeaf *plNmatch       = (TLeaf*)pnTuple3->GetLeaf("nmatch");
  TLeaf *plFracMatch    = (TLeaf*)pnTuple3->GetLeaf("frachitmatch");
  
  TLeaf *plCh_sim       = (TLeaf*)pnTuple3->GetLeaf("chargesim");
  TLeaf *plId_sim       = (TLeaf*)pnTuple3->GetLeaf("idsim");
  TLeaf *plIdParent_sim = (TLeaf*)pnTuple3->GetLeaf("idparentsim");
 
  Int_t nEvents = pnTuple3->GetEntries();
  edm::LogInfo("analyzeReco2SimMuons")<<"Processing "<<nEvents<<" events...";

  for(Long_t ie=0;ie<nEvents;ie++)
    {
      LogDebug("McMatchTrackPlotter::analyzeReco2SimTracks()")<<"Event #= "<<ie;

      pnTuple3->GetEntry(ie);
   
      Float_t eta    = plEta->GetValue();
      Float_t pt     = plPt->GetValue();

      if( !acceptSingle(pt,eta) ) continue;
      	    
      phPtEta_trkReco->Fill(pt,eta);       
      // not enogh! select also the mother to see whcih one is 'fake' (even if real fake) and which one is from the process itself
      bool singleMatched = (plNmatch->GetValue()>0  && plFracMatch->GetValue()>hitfraccut &&
			   (abs(plId_sim->GetValue())==13 && acceptParent(plIdParent_sim->GetValue()) ) 
			   );    
      // if nmatch >0 but the frachitmatch<.9, it's crap;
      LogDebug("McMatchTrackPlotter::analyzeReco2SimTracks()")<<"Filled the all reco tracks";
      // 
      if(singleMatched) 
	{
	  phPtEta_trkRecoSim->Fill(pt,eta);
	  phPtCharge_trkRecoSim->Fill(pt,plCh->GetValue()*plCh_sim->GetValue());
	  phPtId_trkRecoSim->Fill(pt,plId_sim->GetValue());
	  phPtIdparent_trkRecoSim->Fill(pt,plIdParent_sim->GetValue());

	  if( plNmatch->GetValue()>1. ) // matched to more than a sim 
	    {
	      LogDebug("McMatchTrackPlotter::analyzeReco2SimTraks()")<<"Trackss with more matches";
	      phPtEta_trkRecoMany->Fill(pt,eta); 
	    }
	}
      else
	{
	  LogDebug("McMatchTrackPlotter::analyzeReco2SimTrackss()")<<"Fake track";
	  phPtEta_trkRecoFake->Fill(pt,eta); 
	}
    }//entries loop
 
   return;
}


//____________________________________________________________________
void McMatchTrackPlotter::analyzeSim2RecoSingle(TNtuple *pnTuple5)
{
  // reco2sim analysis: purity,resolution, x_pairReco,x_pairRecoSim,x_pairRecoFake, x_RecoMulti
  edm::LogInfo("analyzeSim2RecoTrackss")<<"Analyzing track sim2reco tuple ...";

 //### sim track
  TLeaf *plCharge            = (TLeaf*)pnTuple5->GetLeaf("charge");
  TLeaf *plEta               = (TLeaf*)pnTuple5->GetLeaf("eta");
  TLeaf *plPt                = (TLeaf*)pnTuple5->GetLeaf("pt");
  TLeaf *plId                = (TLeaf*)pnTuple5->GetLeaf("id");
  TLeaf *plIdParent          = (TLeaf*)pnTuple5->GetLeaf("idparent");

  TLeaf *plNmatch            = (TLeaf*)pnTuple5->GetLeaf("nmatch");
  TLeaf *plFracMatch         = (TLeaf*)pnTuple5->GetLeaf("frachitmatch");

  //### reco track
  TLeaf *plCharge_reco       = (TLeaf*)pnTuple5->GetLeaf("chargereco");
  TLeaf *plEta_reco          = (TLeaf*)pnTuple5->GetLeaf("etareco");
  TLeaf *plPt_reco           = (TLeaf*)pnTuple5->GetLeaf("ptreco");

  Int_t nEvents = pnTuple5->GetEntries();
  edm::LogInfo("analyzeSim2RecoTracks")<<"Processing "<<nEvents<<" events...";

  for(Long_t ie=0;ie<nEvents;ie++)
    {
      pnTuple5->GetEntry(ie);
      // for single sim track, the sim2reco info, contains all muons (should modify this, to have the same scheme for both single and double)
      if( (doparentspecificcut && !acceptParent(plIdParent->GetValue()))
	  || abs(plId->GetValue()) != 13 ) continue; // only for muons from a specific process/mother

      edm::LogInfo("analyzeSim2RecoTracks")<<"Id and Id parent "<<plId->GetValue()<<"\t"<<plIdParent->GetValue();      
      Float_t eta    = plEta->GetValue();
      Float_t pt     = plPt->GetValue();
      
      if( !acceptEta(eta) ) continue;

      // all sim tracks
      phPtEta_trkSim->Fill(pt,eta); 
      
      if( plNmatch->GetValue() >0 && plFracMatch->GetValue()>hitfraccut) // only matched tracks
	{
	
	  phPtEta_trkSimReco->Fill(pt,eta);
	  phPtCharge_trkSimReco->Fill(pt,plCharge->GetValue()*plCharge_reco->GetValue());   
	  phIdId_simTrkParent->Fill(plId->GetValue(),plIdParent->GetValue());

	  if(plNmatch->GetValue()>1) // multiple matches for a single sim track ... 'split' track
	    {
	      phPtEta_trkSimMany->Fill(pt,eta);
	    }
	  
	  // eta, phi, pt resolution vs pt
	  Float_t eta_r  = plEta_reco->GetValue();
	  Float_t pt_r   = plPt_reco->GetValue();
	  phPtEta_trkResolution->Fill(pt,(eta_r-eta)/eta); 
	  phPtPt_trkResolution->Fill(pt,(pt_r-pt)/pt);
	
	} //matched track
      else // screwed track
	{
	  phPtEta_trkSimLost->Fill(pt,eta); // lost sim track
	}
    }


  edm::LogInfo("analyzeSim2RecoTrackss")<<"Done.";
}


//____________________________________________________________________________
void McMatchTrackPlotter::beginJob(const edm::EventSetup& es)
{
  // method called once each job just before starting event loop
 
  edm::LogInfo("MCMatchPlotter::beginJob()")<<"beginJob() ...";

  char szBuf[256];
  if(reco2sim && dopair)
    {
      edm::LogInfo("MCMatchPlotter::beginJob()")<<"Reco2Sim pair plots.";
      TNtuple *pnRecoSimPairTuple = NULL;
      
      if(mergedrootfiles) sprintf(szBuf,"pnRecoSimPair%s",typelabel.c_str());
      else sprintf(szBuf,"mcmatchanalysis/pnRecoSimPair%s",typelabel.c_str());
      pnRecoSimPairTuple= (TNtuple*)pfInFile->Get(szBuf);      

      if( pnRecoSimPairTuple !=NULL ) 
        analyzeReco2SimPair(pnRecoSimPairTuple);
      else 
	edm::LogWarning("McMatchTrackPlotter::beginJob()")<<"NO pnRecoSimPair TUPLE!!!";
    }

  if(sim2reco && dopair)
    {
      edm::LogInfo("MCMatchPlotter::beginJob()")<<"Sim2Reco Pair plots.";
      TNtuple *pnSimRecoPairTuple = NULL;
      if(mergedrootfiles) sprintf(szBuf,"pnSimRecoPair%s",typelabel.c_str());
      else sprintf(szBuf,"mcmatchanalysis/pnSimRecoPair%s",typelabel.c_str());
      pnSimRecoPairTuple = (TNtuple*)pfInFile->Get(szBuf);

      if( pnSimRecoPairTuple !=NULL ) 
        analyzeSim2RecoPair(pnSimRecoPairTuple);
      else 
	edm::LogWarning("McMatchTrackPlotter::beginJob()")<<"NO pnSim2RecoPair TUPLE!!!";
    }

  if(sim2reco && dosingle)
    {
      edm::LogInfo("MCMatchPlotter::beginJob()")<<"Sim2Reco Track plots.";
      TNtuple *pnSimRecoTrackTuple = NULL;
    
      if(mergedrootfiles) sprintf(szBuf,"pnSimRecoTrk%s",typelabel.c_str());
      else sprintf(szBuf,"mcmatchanalysis/pnSimRecoTrk%s",typelabel.c_str());
      pnSimRecoTrackTuple = (TNtuple*)pfInFile->Get(szBuf);

      if( pnSimRecoTrackTuple !=NULL ) 
        analyzeSim2RecoSingle(pnSimRecoTrackTuple);
      else 
	edm::LogWarning("McMatchTrackPlotter::beginJob()")<<"NO pnSim2RecoTracks TUPLE!!!";
    }

  if(reco2sim && dosingle)
    {
      edm::LogInfo("MCMatchPlotter::beginJob()")<<"Reco2Sim Track plots.";
      TNtuple *pnRecoSimTrackTuple = NULL;
   
      if(mergedrootfiles) sprintf(szBuf,"pnRecoSimTrk%s",typelabel.c_str());
      else sprintf(szBuf,"mcmatchanalysis/pnRecoSimTrk%s",typelabel.c_str());

      pnRecoSimTrackTuple = (TNtuple*)pfInFile->Get(szBuf);

      if( pnRecoSimTrackTuple !=NULL ) 
        analyzeReco2SimSingle(pnRecoSimTrackTuple);
      else 
	edm::LogWarning("McMatchTrackPlotter::beginJob()")<<"NO pnRecoSimTracks TUPLE!!!";
    }
}


//_________________________________________________________________________
void McMatchTrackPlotter::endJob()
{
  //method called once each job just after ending the event loop 
  // make the efficiency plots
 
}


//_______________________________________________________________________
void McMatchTrackPlotter::initHistograms()
{
  // define the histograms to be filled
  // bin limits
 
  edm::LogInfo("McMatchTrackPlotter::initHistograms")<<"start creating the histos ... ";
  edm::Service<TFileService> fs; 

  initBins();  
 
  if(!fs) 
    {
      throw edm::Exception(edm::errors::NullPointerError, "McMatchTrackPlotter::initHistograms()\n")
	<< "Could not get pointer to THistFileService.\n";
    }

  //______________________ PAIR HISTOS
  // reco2sim 
  Int_t nChargeBins    = 5;
  Float_t *afChargeBins = getBins(nChargeBins,-2.5,2.5);
 
  if(reco2sim && dopair)
    {
           
// FRACTION
      phPtMinv_pairReco          = fs->make<TH2D>("phPtMinv_pairReco",";p_{T}^{reco}[GeV/c];M[GeV/c^{2}]",
						  nPtBins,afPtBins,nMassBins,afMassBins); phPtMinv_pairReco->Sumw2(); 
      phPtMinv_pairRecoSim       = fs->make<TH2D>("phPtMinv_pairRecoSim",";p_{T}^{reco}[GeV/c];M[GeV/c^{2}]",
						  nPtBins,afPtBins,nMassBins,afMassBins); phPtMinv_pairRecoSim->Sumw2();
      phPtMinv_pairRecoFake      = fs->make<TH2D>("phPtMinv_pairRecoFake",";p_{T}^{reco}[GeV/c];M[GeV/c^{2}]",
						  nPtBins,afPtBins,nMassBins,afMassBins); phPtMinv_pairRecoFake->Sumw2();
            

      phPtY_pairReco      = fs->make<TH2D>("phPtY_pairReco",";p_{T}^{reco}[GeV/c];y^{reco}",
					   nPtBins,afPtBins,nYBins,afYBins);    phPtY_pairReco->Sumw2(); 
      phPtY_pairRecoSim   = fs->make<TH2D>("phPtY_pairRecoSim",";p_{T}^{reco}[GeV/c];y^{reco}",
					   nPtBins,afPtBins, nYBins,afYBins);  phPtY_pairRecoSim->Sumw2();
      phPtY_pairRecoFake  = fs->make<TH2D>("phPtY_pairRecoFake",";p_{T}^{reco}[GeV/c];y^{reco};",
					   nPtBins,afPtBins,nYBins,afYBins); 
    }

  if(sim2reco && dopair)
    {
      phPtY_pairSim         = fs->make<TH2D>("phPtY_pairSim",";p_{T}^{sim}[GeV/c];y^{sim}",
					     nPtBins,afPtBins,nYBins,afYBins); phPtY_pairSim->Sumw2(); 
      phPtY_pairSimReco     = fs->make<TH2D>("phPtY_pairSimReco",";p_{T}^{sim}[GeV/c];y^{sim}",
					     nPtBins,afPtBins,nYBins,afYBins); phPtY_pairSimReco->Sumw2(); 
      phPtY_pairSimLost     = fs->make<TH2D>("phPtY_pairSimLost",";p_{T}^{sim}[GeV/c];y^{sim}",
					     nPtBins,afPtBins,nYBins,afYBins); phPtY_pairSimLost->Sumw2();

      phPtMinv_pairSim      = fs->make<TH2D>("phPtMinv_pairSim",";p_{T}^{sim}[GeV/c];M[GeV/c^{2}]",
					       nPtBins,afPtBins,nMassBins,afMassBins); 
      phPtMinv_pairSimReco  = fs->make<TH2D>("phPtMinv_pairSimReco",";p_{T}^{sim}[GeV/c];M[GeV/c^{2}]",
					       nPtBins,afPtBins,nMassBins,afMassBins); 
      phPtMinv_pairSimLost  = fs->make<TH2D>("phPtMinv_pairSimLost",";p_{T}^{sim}[GeV/c];M[GeV/c^{2}]",
					     nPtBins,afPtBins,nMassBins,afMassBins);

      // RESOLUTIONS
      phPtY_pairResolution  = fs->make<TH2D>("phPtY_pairResolution",
					      ";p_{T}^{sim}[GeV/c];(y^{reco}-y^{sim})/y^{sim}",
					     nPtBins,afPtBins,nNormBins,afNormBins);
      phPtMinv_pairResolution = fs->make<TH2D>("phPtMinv_pairResolution",
					       ";p_{T}^{sim}[GeV/c];(M^{reco}-M^{sim})/M^{sim}",
					       nPtBins,afPtBins,nNormBins,afNormBins); 
      phPtPt_pairResolution   = fs->make<TH2D>("phPtPt_pairResolution",
					       ";p_{T}^{sim}[GeV/c];(p_{T}^{reco}-p_{T}^{sim})/p_{T}^{sim}",
					       nPtBins,afPtBins,nNormBins,afNormBins); 
    }

  //_______________________________ SINGLE plots
  if(reco2sim && dosingle)
    {
      phPtCharge_trkRecoSim      = fs->make<TH2D>("phPtCharge_trkRecoSim",
						  ";p_{T}^{reco}[GeV/c];charge^{reco}*charge^{sim}", 
						  nPtBins,afPtBins,nChargeBins,afChargeBins);
     
      phPtEta_trkReco            = fs->make<TH2D>("phPtEta_trkReco",";p_{T}^{reco}[GeV/c];#eta^{reco}",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);  phPtEta_trkReco->Sumw2(); 
      phPtEta_trkRecoFake        = fs->make<TH2D>("phPtEta_trkRecoFake",";p_{T}^{reco}[GeV/c];#eta^{reco}",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);  phPtEta_trkRecoFake->Sumw2();
      phPtEta_trkRecoMany        = fs->make<TH2D>("phPtEta_trkRecoMany",";p_{T}^{reco}[GeV/c];#eta^{reco}",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);  phPtEta_trkRecoMany->Sumw2(); 
      phPtEta_trkRecoSim         = fs->make<TH2D>("phPtEta_trkRecoSim",";p_{T}^{reco}[GeV/c];#eta^{reco}",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);  phPtEta_trkRecoSim->Sumw2(); 

      Int_t nIdBins    = 5000;
      Float_t *afIdBins = getBins(nIdBins, -5000.,5000);
      phPtId_trkRecoSim          = fs->make<TH2D>("phPtId_trkRecoSim",";p_{T}^{reco}[GeV/c];Id",
						  nPtBins,afPtBins,nIdBins,afIdBins);
      phPtIdparent_trkRecoSim    = fs->make<TH2D>("phPtIdparent_trkRecoSim",";p_{T}^{reco}[GeV/c];Id_{parent}",
						  nPtBins,afPtBins,nIdBins,afIdBins);
    }

   
  if(sim2reco && dosingle)
    {
      phPtCharge_trkSimReco      = fs->make<TH2D>("phPtCharge_trkSimReco",
						  ";p_{T}^{sim}[GeV/c];charge^{sim}*charge^{reco}",
						  nPtBins,afPtBins,nChargeBins,afChargeBins);      
    
      phPtEta_trkSim             = fs->make<TH2D>("phPtEta_trkSim",";p_{T}^{sim}[GeV/c];#eta^{sim}",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);    phPtEta_trkSim->Sumw2(); 
      phPtEta_trkSimLost         = fs->make<TH2D>("phPtEta_trkSimLost",";p_{T}[GeV/c];#eta",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);    phPtEta_trkSimLost->Sumw2();
      phPtEta_trkSimMany         = fs->make<TH2D>("phPtEta_trkSimMany",";p_{T}[GeV/c];#eta",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);    phPtEta_trkSimMany->Sumw2();
      phPtEta_trkSimReco         = fs->make<TH2D>("phPtEta_trkSimReco",";p_{T}[GeV/c];#eta",
						  nPtBins,afPtBins,nEtaBins,afEtaBins);    phPtEta_trkSimReco->Sumw2();
       
      // RESOLUTIONS
      phPtPt_trkResolution       = fs->make<TH2D>("phPtPt_trkResolution",
						  ";p_{T}^{sim}[GeV/c];(p_{T}^{reco}-p_{T}^{sim})/p_{T}^{sim}",
						  nPtResBins,afPtResBins,nNormBins,afNormBins); 
      phPtEta_trkResolution      = fs->make<TH2D>("phPtEta_trkResolution",
						  ";p_{T}^{sim}[GeV/c];(#eta^{reco}-#eta^{sim})/#eta^{sim}",
						  nPtResBins,afPtResBins,nNormBins,afNormBins);
      
      Int_t nIntBins5000         = 5000;
      Float_t *afIntBins5000     = getBins(nIntBins5000, -5000., 5000.);
      phIdId_simTrkParent        = fs->make<TH2D>("phIdId_simTrkParent",";Id;Id_{parent}",
						  nIntBins5000,afIntBins5000,nIntBins5000,afIntBins5000);
         
    }
  edm::LogInfo("McMatchTrackPlotter::initHistogram()")<<"Finish initiating the histograms";
  return;
}


//_________________________________________________________________________
Bool_t McMatchTrackPlotter::acceptEta(Float_t eta)
{
  // cut on eta: barrel or endcap so far
  Bool_t etaOK = false;
  if( dobarel && doendcap && fabs(eta)<=2.4 ) etaOK = true;
  if( dobarel && fabs(eta)<=0.8 ) etaOK = true;  
  if(doendcap && (fabs(eta)>0.8 && fabs(eta)<2.4) ) etaOK = true;

  return etaOK;
    
}


//_________________________________________________________________________
Bool_t McMatchTrackPlotter::acceptMass(Float_t e)
{
  // cut on M, depending on the input pdgDilepton: 23 (Z0)
  Bool_t massOK = false;
  if (e<masscutmax && e>masscutmin) massOK = true;

  return massOK;

}


//_________________________________________________________________________
Bool_t McMatchTrackPlotter::acceptPair(Float_t pt, Float_t m)
{
  // cut on M, depending on the input pdgDilepton: 23 (Z0)
  Bool_t ok = false;
  if(acceptMass(m) && pt>ptcutpair) ok= true;

  return ok;
}


//_______________________________________________________________________
Bool_t McMatchTrackPlotter::acceptSingle(Float_t pt, Float_t eta)
{
  // cut on M, depending on the input pdgDilepton: 23 (Z0)
  Bool_t ok = false;
  if(acceptEta(eta) && pt>ptcutsingle) ok= true;

  return ok;
}


//_________________________________________________________________________
Bool_t McMatchTrackPlotter::acceptParent(Float_t idParent)
{
  // cut on parentId, depending on the input pdgpair
  Bool_t idOK = true;
 
  if(idParent != pdgpair) idOK = false;

  return idOK;

}


//________________________________________________________________________
DEFINE_FWK_MODULE(McMatchTrackPlotter);


