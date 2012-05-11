// -*- C++ -*-
//
// Package:    DiMuonOnia2DPlots
// Class:      DiMuonOnia2DPlots
// 
/**\class DiMuonOnia2DPlots DiMuonOnia2DPlots.cc DiMuonAna/DiMuonOnia2DPlots/src/DiMuonOnia2DPlots.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Dilep PING, Vineet Kumar, Prashant Shukla
// Editted by Dong Ho Moon
//         Created:  Wed May 12 13:45:14 CEST 2010
// $Id: DiMuonOnia2DPlots.cc,v 1.4 2011/11/01 08:22:57 kumarv Exp $
//
//
// system include files
#include <memory>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>



// user include files   
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TFile.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

// dmoon add
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"
#include "FWCore/Common/interface/TriggerNames.h"


using std::cout;
using std::endl;

using namespace pat;
using namespace reco;
using namespace edm;
using namespace std;
using namespace HepMC;
using namespace muon;

//
// class declaration
//

class DiMuonOnia2DPlots : public edm::EDAnalyzer {
    public:
        explicit DiMuonOnia2DPlots(const edm::ParameterSet&);
        ~DiMuonOnia2DPlots();

        const CentralityBins *cbins_;
        CentralityProvider *centrality_;    



        bool fisCuts;  
        std::string fOutputFileName;
        //std::string fGenLevel;
        //std::string fHLTPathName;

        std::string fIsGenInfo;
        std::string fIsPATInfo;
        std::string fHLTFilterName;  
        std::string fMotherID;  



        //edm::InputTag fHLTFilterName;


        TFile *In;

        TH1F *Centrality;
        TH1F *CentralityTest;
        TFile *fOutputFile ;

        TTree *SingleMuonTree;
        TTree *SingleGenMuonTree;

        int bin, gbin, rbin;   


    private:

        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        virtual bool matchPATMuon(const pat::Muon *pMuon);
        virtual void FillTree(const edm::Event&, const edm::EventSetup&);
        virtual void FillGenTree(const edm::Event&, const edm::EventSetup&);
        math::XYZPoint RefVtx;
        float nPV;
        double zVtx;
        //Tree variables defined here  
        int eventNb,runNb,lumiBlock;

        //1.) J/psi variables RECO                                                                                                                                     
        //init events                                                                                                                                             
        //int RecJPsiSize, RecMuPlusSize, RecMuMinusSize;

        double JpsiCharge,JpsiNo , JpsiMass , JpsiPt , JpsiRap ,JpsiVprob ;
        double JpsiPx , JpsiPy , JpsiPz ;

        // dmoon add
        double JpsiPhi, JpsiPsi[100], JpsiNfPsi[100];
        double JpsiGenPsi;
        double JpsiEta;
        int nEP, nNfEP;
        double rpAng[100], rpCos[100], rpSin[100], NfRpAng[100], NfRpCos[100], NfRpSin[100];
        double rgRpAng;
        int hbit1, hbit2, hbit3, hbit4, hbit5, hbit6, hCowboy;
        int ghbit1, ghbit2, ghbit3, ghbit4, ghbit5, ghbit6, ghCowboy;
        double gJpsiPt, gJpsiEta, gJpsiRap, gJpsiPsi, gJpsiPhi, gJpsiNo, gJpsiMass;
        double gJpsiPx, gJpsiPy, gJpsiPz;

        double RecoCtau,RecoCtauErr,RecoCtauTrue;  
        //TLorentzVector* JpsiP;

        //2.) muon variables RECO                                                                                                                                                                                  
        double muPosPx, muPosPy, muPosPz,  muPosEta,  muPosPhi;
        double muNegPx, muNegPy, muNegPz,  muNegEta,  muNegPhi;

        //3.) cut variables

        //(i). Positive Muon                                                                                                            
        double muPos_nchi2In, muPos_dxy, muPos_dz, muPos_nchi2Gl;
        int muPos_found, muPos_pixeLayers, muPos_nValidMuHits,muPos_arbitrated;
        bool muPos_matches, muPos_tracker, muPos_global;  
        //(ii).Negative Muon                                                                                                             
        double muNeg_nchi2In, muNeg_dxy, muNeg_dz, muNeg_nchi2Gl;
        int muNeg_found, muNeg_pixeLayers, muNeg_nValidMuHits,muNeg_arbitrated;
        bool muNeg_matches, muNeg_tracker,muNeg_global;  

        int GeventNb,GrunNb,GlumiBlock;
        //Gen JPsi Variables
        double GenJpsiMass, GenJpsiPt, GenJpsiRap;
        double GenJpsiPx, GenJpsiPy, GenJpsiPz;
        // dmoon add
        double GenJpsiPhi, GenJpsiPsi;
        double GenJpsiEta;
        double gRpAng;
        double hi_b, hi_npart, hi_ncoll, hi_nhard;
        double rhi_b, rhi_npart, rhi_ncoll, rhi_nhard;

        //2.) muon variables Gen                                         

        double GenmuPosPx, GenmuPosPy, GenmuPosPz,  GenmuPosEta, GenmuPosPhi;
        double GenmuNegPx, GenmuNegPy, GenmuNegPz,  GenmuNegEta, GenmuNegPhi;

        double rGenmuPosPx, rGenmuPosPy, rGenmuPosPz,  rGenmuPosEta, rGenmuPosPhi;
        double rGenmuNegPx, rGenmuNegPy, rGenmuNegPz,  rGenmuNegEta, rGenmuNegPhi;

        edm::InputTag triggerResults_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


DiMuonOnia2DPlots::DiMuonOnia2DPlots(const edm::ParameterSet& iConfig):
    centrality_(0),
    //cbins_(0),
    fisCuts(iConfig.getUntrackedParameter<bool>("IsCuts")),
    fOutputFileName(iConfig.getUntrackedParameter<string>("OutputFileName")),
    fIsGenInfo(iConfig.getUntrackedParameter<string>("IsGenInfo")), 
    fIsPATInfo(iConfig.getUntrackedParameter<string>("IsPATInfo")), 
    fHLTFilterName(iConfig.getUntrackedParameter<string>("HLTFilterName")),
    fMotherID(iConfig.getUntrackedParameter<string>("MotherID")),
    triggerResults_(iConfig.getParameter<edm::InputTag>("trgResults"))
    //fHLTFilterName(iConfig.getUntrackedParameter<edm::InputTag>("HLTFilterName"))
{ 

    //now do what ever initialization is needed

}


DiMuonOnia2DPlots::~DiMuonOnia2DPlots()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------

void DiMuonOnia2DPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //cout << " ------------------------------------------------------ " << endl;
    using namespace edm;
    using namespace std;
    nPV = 0 ;

    centrality_ = new CentralityProvider(iSetup);
    centrality_->newEvent(iEvent,iSetup);
    bin = centrality_->getBin();
    Centrality->Fill(bin);

    double etbins[] = { 0, 2.02978, 3.79592, 5.91489, 8.50671, 11.8798, 15.9046, 20.7224, 26.4392, 33.6236, 42.3893, 53.2046, 67.2965, 84.7366, 105.139, 128.815, 
        157.11, 190.711, 228.571, 270.976, 319.346, 373.91, 435.673, 504.613, 581.418, 665.433, 760.095, 863.427, 977.871, 1104.11, 1242.55, 1394.24, 1564.32, 
        1749.37, 1952.94, 2183.64, 2452.93, 2746.6, 3116.24, 3591.07 };

    double sumEt = centrality_->raw()->EtHFtowerSum();

    for (int i=0; i<40;++i) {
        if (sumEt >= etbins[i] && sumEt<etbins[i+1]) CentralityTest->Fill(i);
    }

    // Primary Vertex
    Handle<reco::VertexCollection> privtxs;
    iEvent.getByLabel("hiSelectedVertex", privtxs);
    VertexCollection::const_iterator privtx;
    nPV = privtxs->size();
    if(!nPV) return; 
    if ( privtxs->begin() != privtxs->end() ) {
        privtx=privtxs->begin();
        RefVtx = privtx->position();
        zVtx = RefVtx.Z();

    } else {
        RefVtx.SetXYZ(0.,0.,0.);
    }
    if(!strcmp(fIsPATInfo.c_str(),"TRUE")){FillTree(iEvent,iSetup);}
    if(!strcmp(fIsGenInfo.c_str(),"TRUE")){FillGenTree(iEvent,iSetup);}

}

// ------------ method called once each job just before starting event loop  ------------
    void 
DiMuonOnia2DPlots::beginJob()
{

    edm::Service<TFileService> fs;
    fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" );

    //cout<<"begin job"<<endl;

    SingleMuonTree = new TTree("SingleMuonTree","SingleMuonTree");

    // Event variables
    SingleMuonTree->Branch("eventNb",             &eventNb,             "eventNb/I");
    SingleMuonTree->Branch("runNb",               &runNb,               "runNb/I");
    SingleMuonTree->Branch("lumiBlock",           &lumiBlock,           "lumiBlock/I");
    SingleMuonTree->Branch("nPV",                 &nPV,                 "nPV/I");
    SingleMuonTree->Branch("zVtx",                &zVtx,                "zVtx/D");
    SingleMuonTree->Branch("rbin",                &rbin,                "rbin/I");

    // HLT bits
    SingleMuonTree->Branch("hbit1",               &hbit1,               "hbit1/I");
    SingleMuonTree->Branch("hbit2",               &hbit2,               "hbit2/I");
    SingleMuonTree->Branch("hbit3",               &hbit3,               "hbit3/I");
    SingleMuonTree->Branch("hbit4",               &hbit4,               "hbit4/I");
    SingleMuonTree->Branch("hbit5",               &hbit5,               "hbit5/I");
    SingleMuonTree->Branch("hbit6",               &hbit6,               "hbit6/I");
    SingleMuonTree->Branch("hCowboy",             &hCowboy,             "hCowboy/I");


    // dmoon add
    SingleMuonTree->Branch("nEP",   &nEP,   "nEP/I");
    SingleMuonTree->Branch("nNfEP", &nNfEP, "nNfEP/I");
    SingleMuonTree->Branch("rpAng", &rpAng, "rpAng[nEP]/D");
    SingleMuonTree->Branch("rpSin", &rpSin, "rpSin[nEP]/D");
    SingleMuonTree->Branch("rpCos", &rpCos, "rpCos[nEP]/D");
    SingleMuonTree->Branch("NfRpAng", &NfRpAng, "NfRpAng[nNfEP]/D");
    SingleMuonTree->Branch("NfRpSin", &NfRpSin, "NfRpSin[nNfEP]/D");
    SingleMuonTree->Branch("NfRpCos", &NfRpCos, "NfRpCos[nNfEP]/D");

    SingleMuonTree->Branch("rgRpAng", &rgRpAng, "rgRpAng/D");
    SingleMuonTree->Branch("rhi_b", &rhi_b, "rhi_b/D");
    SingleMuonTree->Branch("rhi_npart", &rhi_npart, "rhi_npart/D");
    SingleMuonTree->Branch("rhi_ncoll", &rhi_ncoll, "rhi_ncoll/D");
    SingleMuonTree->Branch("rhi_nhard", &rhi_nhard, "rhi_nhard/D");


    //SingleMuonTree->Branch("RecJPsiSize",   &RecJPsiSize,  "RecJPsiSize/I");  

    // Jpsi Variables dimuon variable.

    //RecJPsiSize  (if we want to put array
    //SingleMuonTree->Branch("JpsiCharge", JpsiCharge,  "JpsiCharge[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiMass",   JpsiMass,  "JpsiMass[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiPt",     JpsiPt,    "JpsiPt[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiRap",    JpsiRap,   "JpsiRap[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiPx",     JpsiPx,    "JpsiPx[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiPy",     JpsiPy,    "JpsiPy[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiPz",     JpsiPz,    "JpsiPz[RecJPsiSize]/D");
    //SingleMuonTree->Branch("JpsiVprob",  JpsiVprob, "JpsiVprob[RecJPsiSize]/D");



    //RecJPsiSize      
    SingleMuonTree->Branch("JpsiCharge", &JpsiCharge,  "JpsiCharge/D");
    SingleMuonTree->Branch("JpsiNo",     &JpsiNo,      "JpsiNo/D");
    SingleMuonTree->Branch("JpsiMass",   &JpsiMass,    "JpsiMass/D");
    SingleMuonTree->Branch("JpsiPt",     &JpsiPt,      "JpsiPt/D");
    SingleMuonTree->Branch("JpsiPhi",    &JpsiPhi,     "JpsiPhi/D"); // dmoon add
    SingleMuonTree->Branch("JpsiPsi",    &JpsiPsi,     "JpsiPsi[38]/D"); // dmoon add
    SingleMuonTree->Branch("JpsiGenPsi",    &JpsiGenPsi,     "JpsiGenPsi/D"); // dmoon add
    SingleMuonTree->Branch("JpsiNfPsi",  &JpsiNfPsi,   "JpsiNfPsi[38]/D"); // dmoon add
    SingleMuonTree->Branch("JpsiRap",    &JpsiRap,   "JpsiRap/D");
    SingleMuonTree->Branch("JpsiEta",    &JpsiEta,   "JpsiEta/D");
    SingleMuonTree->Branch("JpsiPx",     &JpsiPx,    "JpsiPx/D");
    SingleMuonTree->Branch("JpsiPy",     &JpsiPy,    "JpsiPy/D");
    SingleMuonTree->Branch("JpsiPz",     &JpsiPz,    "JpsiPz/D");
    SingleMuonTree->Branch("JpsiVprob",  &JpsiVprob, "JpsiVprob/D");
    SingleMuonTree->Branch("RecoCtau",   &RecoCtau,   "RecoCtau/D");
    SingleMuonTree->Branch("RecoCtauErr", &RecoCtauErr,   "RecoCtauErr/D");
    SingleMuonTree->Branch("RecoCtauTrue", &RecoCtauTrue,   "RecoCtauTrue/D");

    SingleMuonTree->Branch("gJpsiNo",     &gJpsiNo,      "gJpsiNo/D");
    SingleMuonTree->Branch("gJpsiMass",   &gJpsiMass,    "gJpsiMass/D");
    SingleMuonTree->Branch("gJpsiPt",     &gJpsiPt,      "gJpsiPt/D");
    SingleMuonTree->Branch("gJpsiPhi",    &gJpsiPhi,     "gJpsiPhi/D"); // dmoon add
    SingleMuonTree->Branch("gJpsiPsi",    &gJpsiPsi,     "gJpsiPsi/D"); // dmoon add
    SingleMuonTree->Branch("gJpsiRap",    &gJpsiRap,     "gJpsiRap/D");
    SingleMuonTree->Branch("gJpsiEta",    &gJpsiEta,     "gJpsiEta/D");
    SingleMuonTree->Branch("gJpsiPx",     &gJpsiPx,      "gJpsiPx/D");
    SingleMuonTree->Branch("gJpsiPy",     &gJpsiPy,      "gJpsiPy/D");
    SingleMuonTree->Branch("gJpsiPz",     &gJpsiPz,      "gJpsiPz/D");

    //muon variable
    SingleMuonTree->Branch("muPosPx",    &muPosPx,   "muPosPx/D");
    SingleMuonTree->Branch("muPosPy",    &muPosPy,   "muPosPy/D");
    SingleMuonTree->Branch("muPosPz",    &muPosPz,   "muPosPz/D");
    SingleMuonTree->Branch("muPosEta",    &muPosEta,  "muPosEta/D");
    SingleMuonTree->Branch("muPosPhi",    &muPosPhi,  "muPosPhi/D");

    SingleMuonTree->Branch("muNegPx",    &muNegPx,   "muNegPx/D");
    SingleMuonTree->Branch("muNegPy",    &muNegPy,   "muNegPy/D");
    SingleMuonTree->Branch("muNegPz",    &muNegPz,   "muNegPz/D");
    SingleMuonTree->Branch("muNegEta",   &muNegEta,   "muNegEta/D");
    SingleMuonTree->Branch("muNegPhi",   &muNegPhi,   "muNegPhi/D");


    //1). Positive Muon 
    SingleMuonTree->Branch("muPos_nchi2In", &muPos_nchi2In, "muPos_nchi2In/D");
    SingleMuonTree->Branch("muPos_dxy", &muPos_dxy, "muPos_dxy/D");
    SingleMuonTree->Branch("muPos_dz", &muPos_dz, "muPos_dz/D");
    SingleMuonTree->Branch("muPos_nchi2Gl", &muPos_nchi2Gl, "muPos_nchi2Gl/D");
    SingleMuonTree->Branch("muPos_found", &muPos_found, "muPos_found/I");
    SingleMuonTree->Branch("muPos_pixeLayers", &muPos_pixeLayers, "muPos_pixeLayers/I");
    SingleMuonTree->Branch("muPos_nValidMuHits", &muPos_nValidMuHits, "muPos_nValidMuHits/I");
    SingleMuonTree->Branch("muPos_arbitrated", &muPos_arbitrated, "muPos_arbitrated/I");
    SingleMuonTree->Branch("muPos_matches", &muPos_matches, "muPos_matches/O");
    SingleMuonTree->Branch("muPos_tracker", &muPos_tracker, "muPos_tracker/O");
    SingleMuonTree->Branch("muPos_global", &muPos_global, "muPos_global/O");




    //2). Negative Muon
    SingleMuonTree->Branch("muNeg_nchi2In", &muNeg_nchi2In, "muNeg_nchi2In/D");
    SingleMuonTree->Branch("muNeg_dxy", &muNeg_dxy, "muNeg_dxy/D");
    SingleMuonTree->Branch("muNeg_dz", &muNeg_dz, "muNeg_dz/D");
    SingleMuonTree->Branch("muNeg_nchi2Gl", &muNeg_nchi2Gl, "muNeg_nchi2Gl/D");
    SingleMuonTree->Branch("muNeg_found", &muNeg_found, "muNeg_found/I");
    SingleMuonTree->Branch("muNeg_pixeLayers", &muNeg_pixeLayers, "muNeg_pixeLayers/I");
    SingleMuonTree->Branch("muNeg_nValidMuHits", &muNeg_nValidMuHits, "muNeg_nValidMuHits/I");
    SingleMuonTree->Branch("muNeg_arbitrated", &muNeg_arbitrated, "muNeg_arbitrated/I");
    SingleMuonTree->Branch("muNeg_matches", &muNeg_matches, "muNeg_matches/O");
    SingleMuonTree->Branch("muNeg_tracker", &muNeg_tracker, "muNeg_tracker/O");
    SingleMuonTree->Branch("muNeg_global", &muNeg_global, "muNeg_global/O");


    SingleGenMuonTree = new TTree("SingleGenMuonTree","SingleGenMuonTree");


    // Event variables                                                              

    SingleGenMuonTree->Branch("GeventNb",   &GeventNb,       "GeventNb/I");
    SingleGenMuonTree->Branch("GrunNb",     &GrunNb,         "GrunNb/I");
    SingleGenMuonTree->Branch("GlumiBlock", &GlumiBlock,     "GlumiBlock/I");

    // dmoon add
    SingleGenMuonTree->Branch("gRpAng", &gRpAng, "gRpAng/D");
    SingleGenMuonTree->Branch("hi_b", &hi_b, "hi_b/D");
    SingleGenMuonTree->Branch("hi_npart", &hi_npart, "hi_npart/D");
    SingleGenMuonTree->Branch("hi_ncoll", &hi_ncoll, "hi_ncoll/D");
    SingleGenMuonTree->Branch("hi_nhard", &hi_nhard, "hi_nhard/D");


    SingleGenMuonTree->Branch("ghbit1",            &ghbit1,               "ghbit1/I");
    SingleGenMuonTree->Branch("ghbit2",            &ghbit2,               "ghbit2/I");
    SingleGenMuonTree->Branch("ghbit3",            &ghbit3,               "ghbit3/I");
    SingleGenMuonTree->Branch("ghbit4",            &ghbit4,               "ghbit4/I");
    SingleGenMuonTree->Branch("ghbit5",            &ghbit5,               "ghbit5/I");
    SingleGenMuonTree->Branch("ghbit6",            &ghbit6,               "ghbit6/I");
    SingleGenMuonTree->Branch("ghCowboy",          &ghCowboy,             "ghCowboy/I");

    //Gen Jpsi Variables                                                                                                                                                        
    SingleGenMuonTree->Branch("GenJpsiMass",   &GenJpsiMass,  "GenJpsiMass/D");
    SingleGenMuonTree->Branch("GenJpsiPt",     &GenJpsiPt,    "GenJpsiPt/D");
    SingleGenMuonTree->Branch("GenJpsiPhi",    &GenJpsiPhi,    "GenJpsiPhi/D"); // dmoon add
    SingleGenMuonTree->Branch("GenJpsiPsi",    &GenJpsiPsi,    "GenJpsiPsi/D"); // dmoon add
    SingleGenMuonTree->Branch("GenJpsiRap",    &GenJpsiRap,   "GenJpsiRap/D");
    SingleGenMuonTree->Branch("GenJpsiEta",    &GenJpsiEta,   "GenJpsiEta/D");
    SingleGenMuonTree->Branch("GenJpsiPx",     &GenJpsiPx,    "GenJpsiPx/D");
    SingleGenMuonTree->Branch("GenJpsiPy",     &GenJpsiPy,    "GenJpsiPy/D");
    SingleGenMuonTree->Branch("GenJpsiPz",     &GenJpsiPz,    "GenJpsiPz/D");
    SingleGenMuonTree->Branch("gbin",          &gbin,             "gbin/I");

    //muon variable
    SingleGenMuonTree->Branch("GenmuPosPx",    &GenmuPosPx,   "GenmuPosPx/D");
    SingleGenMuonTree->Branch("GenmuPosPy",    &GenmuPosPy,   "GenmuPosPy/D");
    SingleGenMuonTree->Branch("GenmuPosPz",    &GenmuPosPz,   "GenmuPosPz/D");
    SingleGenMuonTree->Branch("GenmuPosEta",   &GenmuPosEta,   "GenmuPosEta/D");
    SingleGenMuonTree->Branch("GenmuPosPhi",   &GenmuPosPhi,   "GenmuPosPhi/D");

    SingleGenMuonTree->Branch("GenmuNegPx",    &GenmuNegPx,   "GenmuNegPx/D");
    SingleGenMuonTree->Branch("GenmuNegPy",    &GenmuNegPy,   "GenmuNegPy/D");
    SingleGenMuonTree->Branch("GenmuNegPz",    &GenmuNegPz,   "GenmuNegPz/D");
    SingleGenMuonTree->Branch("GenmuNegEta",   &GenmuNegEta,   "GenmuNegEta/D");
    SingleGenMuonTree->Branch("GenmuNegPhi",   &GenmuNegPhi,   "GenmuNegPhi/D");



    //cout<<"Tree booked "<<endl;


    //Histograms       
    Centrality = new TH1F("Centrality","Centrality", 60,-10,50);
    CentralityTest = new TH1F("CentralityTest","CentralityTest;centrality bin;Number of Events",40,0,40);

    //h_ZetaGen_ = genParticleDir.make<TH1D>("generatedZeta","#eta of generated Z",100,-5.,5.); 

    // Write comments in a file
}



// ------------ method called once each job just after ending the event loop  ------------
void DiMuonOnia2DPlots::endJob() 
{
    //cout<<"End Job"<<endl;
    fOutputFile->cd();
    SingleMuonTree->Write();
    SingleGenMuonTree->Write();
    Centrality->Write(); 
    CentralityTest->Write(); 
    fOutputFile->Close();

}



bool DiMuonOnia2DPlots::matchPATMuon(const pat::Muon *pMuon) 
{
    return(
            //to match with filter name
            //(!pMuon->triggerObjectMatchesByFilter(fHLTFilterName).empty())  
            //to match with trigger name
            //!pMuon->triggerObjectMatchesByPath("HLT_HIL1DoubleMuOpen").empty()
            !pMuon->triggerObjectMatchesByPath(fHLTFilterName).empty()
          );
}

void DiMuonOnia2DPlots::FillTree(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //init events
    //RecJPsiSize=0;
    //reset J/psi RECO variables
    JpsiNo=-9999.;
    JpsiCharge=-9999.;
    JpsiMass=-9999.;
    JpsiPt=-9999.;
    JpsiPhi=-9999.; // dmoon add
    JpsiEta=-9999.; // dmoon add
    JpsiRap=-9999.;
    JpsiVprob=-9999.;
    JpsiPx=-9999.;
    JpsiPy=-9999.;
    JpsiPz=-9999.;
    RecoCtau=-9999;
    RecoCtauErr=-9999;
    RecoCtauTrue=-9999;

    muPosPx=-9999.;
    muPosPy=-9999.;
    muPosPz=-9999.;
    muPosEta=-9999.;
    muPosPhi=-9999.;

    muNegPx=-9999.;
    muNegPy=-9999.;
    muNegPz=-9999.;
    muNegEta=-9999.;
    muNegPhi=-9999.;


    //1).Positive Muon 
    muPos_nchi2In=-9999.;
    muPos_dxy=-9999.;
    muPos_dz=-9999.;
    muPos_nchi2Gl=-9999.;
    muPos_found=-9999;
    muPos_pixeLayers=-9999;
    muPos_nValidMuHits=-9999;
    muPos_arbitrated=-9999;
    muPos_matches=0;
    muPos_tracker=0;
    muPos_global=0;

    //2).Negtive Muon 
    muNeg_nchi2In=-9999.;
    muNeg_dxy=-9999.;
    muNeg_dz=-9999.;
    muNeg_nchi2Gl=-9999.;
    muNeg_found=-9999;
    muNeg_pixeLayers=-9999;
    muNeg_nValidMuHits=-9999;
    muNeg_arbitrated=-9999;
    muNeg_matches=0;
    muNeg_tracker=0;
    muNeg_global=0;

    //reset EVENT information 
    eventNb= 0 ;
    runNb= 0 ;
    lumiBlock= 0 ;
    rbin=0;

    // Event related infos
    eventNb= iEvent.id().event();
    runNb=iEvent.id().run();
    lumiBlock= iEvent.luminosityBlock();

    centrality_ = new CentralityProvider(iSetup);
    centrality_->newEvent(iEvent,iSetup);
    rbin = centrality_->getBin();

    Handle<TriggerResults> trigResults;
    iEvent.getByLabel(triggerResults_,trigResults);

    TriggerNames triggerNames_;
    std::string bit1 = "HLT_HIL1DoubleMu0_HighQ_v1";
    std::string bit2 = "HLT_HIL2DoubleMu3_v1";
    std::string bit3 = "HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1";
    std::string bit4 = "HLT_HIL3DoubleMuOpen_v1";
    std::string bit5 = "HLT_HIL2Mu7_v1";
    std::string bit6 = "HLT_HIL2Mu15_v1";
    /*
    Trigger Names : HLT_HIL1DoubleMuOpen_v1
    Trigger Names : HLT_HIL1DoubleMu0_HighQ_v1
    Trigger Names : HLT_HIL2Mu3_v1
    Trigger Names : HLT_HIL2Mu3_NHitQ_v1
    Trigger Names : HLT_HIL2Mu7_v1
    Trigger Names : HLT_HIL2Mu15_v1
    Trigger Names : HLT_HIL2DoubleMu0_v1
    Trigger Names : HLT_HIL2DoubleMu0_NHitQ_v1
    Trigger Names : HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1
    Trigger Names : HLT_HIL2DoubleMu3_v1
    Trigger Names : HLT_HIL3Mu3_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_SS_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_OS_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1
    */
    int nbit1=0; int nbit2=0; int nbit3=0; int nbit4=0; int nbit5=0; int nbit6=0;

    if (trigResults.isValid()) {
        int ntrigs = trigResults->size();
        if (ntrigs==0){
            //std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;
        }
        edm::TriggerNames const& triggerNames = iEvent.triggerNames(*trigResults);
        for (int itrig = 0; itrig != ntrigs; ++itrig){
            std::string trigName=triggerNames.triggerName(itrig);
            bool accept = trigResults->accept(itrig);
            //cout<<"Trigger Names : "<<trigName<<" bit : "<<accept<<endl;
            if (accept) {
                if(trigName == bit1){ nbit1 = 1;}
                if(trigName == bit2){ nbit2 = 1;}
                if(trigName == bit3){ nbit3 = 1;}
                if(trigName == bit4){ nbit4 = 1;}
                if(trigName == bit5){ nbit5 = 1;}
                if(trigName == bit6){ nbit6 = 1;}
                /*
                if(trigName == bit1){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit1 = 1;}
                if(trigName == bit2){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit2 = 1;}
                if(trigName == bit3){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit3 = 1;}
                if(trigName == bit4){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit4 = 1;}
                if(trigName == bit5){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit5 = 1;}
                if(trigName == bit6){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit6 = 1;}
                */
            }
        }
    }else{cout<<"HLT Info is not ok !!!!!"<<endl;}
    //cout<<"nbit1 : "<<nbit1<<" nbit2 : "<<nbit2<<" nbit3 : "<<nbit3<<endl;
    //cout<<"nbit4 : "<<nbit4<<" nbit5 : "<<nbit5<<" nbit6 : "<<nbit6<<endl;
    hbit1 = nbit1; hbit2 = nbit2; hbit3 = nbit3;
    hbit4 = nbit4; hbit5 = nbit5; hbit6 = nbit6;
    //cout<<"hbit1 : "<<hbit1<<" hbit2 : "<<hbit2<<" hbit3 : "<<hbit3<<endl;
    //cout<<"hbit4 : "<<hbit4<<" hbit5 : "<<hbit5<<" hbit6 : "<<hbit6<<endl;

    // dmoon add 
    edm::Handle<reco::EvtPlaneCollection> FlatEvtPlanes;
    edm::Handle<reco::EvtPlaneCollection> NoFlatEvtPlanes;
    iEvent.getByLabel("hiEvtPlaneFlat",FlatEvtPlanes);
    iEvent.getByLabel("hiEvtPlane","recoLevel",NoFlatEvtPlanes);

    int EP = 0; int NfEP = 0; 
    if(!FlatEvtPlanes.isValid()){
        //cout << "Error! Can't get hiEvtPlane product!" << endl;
        return ;
    }
    for (reco::EvtPlaneCollection::const_iterator rp = FlatEvtPlanes->begin();rp !=FlatEvtPlanes->end(); rp++) {
        rpAng[EP] = rp->angle();
        rpSin[EP] = rp->sumSin();
        rpCos[EP] = rp->sumCos();
        //rpLabel[nEP]  = rp->label();
        //cout<<"Flatten Evt plane method: "<< rp->label() <<", rpAng["<<EP<<"] : "<<rpAng[EP]<<endl;
        EP++;

    }
    nEP = EP;
    if(!NoFlatEvtPlanes.isValid()){
        //cout << "Error! Can't get hiEvtPlane product!" << endl;
        return ;
    }
    for (reco::EvtPlaneCollection::const_iterator rp = NoFlatEvtPlanes->begin();rp !=NoFlatEvtPlanes->end(); rp++) {
        NfRpAng[NfEP] = rp->angle();
        NfRpSin[NfEP] = rp->sumSin();
        NfRpCos[NfEP] = rp->sumCos();
        //NfRpLabel[nNfEP] = rp->label();
        //cout<<"Non Flatten Evt plane method: "<< rp->label() <<", NfRpAng["<<NfEP<<"] : "<<NfRpAng[NfEP]<<endl;
        NfEP++;
    }
    nNfEP = NfEP;



    //cout<<" rbin "<<rbin<<endl;

    //--------------------------------------Reco DimuonGlobal ---------------------------------------------------------------------

    using namespace edm;
    using namespace std;
    using namespace pat;

    edm::Handle<edm::View<pat::CompositeCandidate> > diMuonsPATCand;
    iEvent.getByLabel("onia2MuMuPatGlbGlb", diMuonsPATCand);
    if(!(diMuonsPATCand.isValid())) return;
    edm::View<pat::CompositeCandidate>dimuonsPATColl= *diMuonsPATCand;
    JpsiNo=dimuonsPATColl.size();
    //cout<<" reco Jpsi size : "<<dimuonsPATColl.size()<<endl;

    // dmoon add gen reaction plane
    const GenEvent* evt;
    edm::Handle<HepMCProduct> mc;
    iEvent.getByLabel("generator",mc);
    evt = mc->GetEvent();
    const HeavyIon* hi = evt->heavy_ion();
    rhi_b = hi->impact_parameter();
    rhi_npart = hi->Npart_proj()+hi->Npart_targ();
    rhi_ncoll = hi->Ncoll();
    rhi_nhard = hi->Ncoll_hard();
    rgRpAng = hi->event_plane_angle();
    if(rgRpAng < -TMath::Pi()) rgRpAng += 2.*TMath::Pi();
    if(rgRpAng > TMath::Pi()) rgRpAng -= 2.*TMath::Pi();
    if(rgRpAng < -TMath::Pi()/2) rgRpAng += TMath::Pi();
    if(rgRpAng > TMath::Pi()/2) rgRpAng -= TMath::Pi();
    //cout<<" gRpAng : "<<gRpAng<<endl;

    float mumass =0.105658389;

    edm::Handle<edm::View<reco::GenParticle> >gPar;
    iEvent.getByLabel("hiGenParticles",gPar) ;

    if(!(gPar.isValid())) return;

    edm::View<reco::GenParticle> genParticles = *gPar ;
    TLorentzVector  genvector1, genvector2;

    double px1[10000], py1[10000], pz1[10000], px2[10000], py2[10000], pz2[10000];
    unsigned int nplus = 0, nminus =0;

    for(size_t i = 0; i < genParticles.size(); ++ i) {
        const reco::GenParticle& part = (*gPar)[i];
        const  Candidate *mom = part.mother();


        if (part.numberOfMothers()!=1) continue;
        int momId = mom->pdgId();


        if ((abs(part.pdgId()) == 13) && ( momId == 443 ) ){

            if(part.pdgId() == 13 ){
                px1[nplus] = part.px();
                py1[nplus] = part.py();
                pz1[nplus] = part.pz();
                nplus++;

                rGenmuNegPx=part.px();
                rGenmuNegPy=part.py();
                rGenmuNegPz=part.pz();
                rGenmuNegEta=part.eta();
                rGenmuNegPhi=part.phi();

                //cout<<"motherID "<<MID<<endl;

            }

            if(part.pdgId()== -13) {
                px2[nminus] = part.px();
                py2[nminus] = part.py();
                pz2[nminus] = part.pz();
                nminus++;

                rGenmuPosPx=part.px();
                rGenmuPosPy=part.py();
                rGenmuPosPz=part.pz();
                rGenmuPosEta=part.eta();
                rGenmuPosPhi=part.phi();

            }
        }
    }
    //cout<<" nplus : "<<nplus<<"  "<<nminus<<endl;

    for(size_t i = 0; i < nplus; i++) {
        double en1 = sqrt(px1[i]*px1[i] + py1[i]*py1[i] + pz1[i]*pz1[i] + mumass*mumass);

        for(size_t j = 0; j< nminus; j++) {
            double en2 = sqrt(px2[j]*px2[j] + py2[j]*py2[j] + pz2[j]*pz2[j] + mumass*mumass);
            TLorentzVector  genvector1,genvector2,genvector3;

            genvector1.SetPxPyPzE(px1[i], py1[i], pz1[i], en1);
            genvector2.SetPxPyPzE(px2[j], py2[j], pz2[j], en2);

            genvector3=genvector1+genvector2;

            double GenDiMuonY=genvector3.Rapidity();
            double GenDiMuonMinv=genvector3.M();
            double GenDiMuonPt =genvector3.Pt();
            double GenDiMuonPhi =genvector3.Phi();
            double GenDiMuonEta =genvector3.Eta();
            double GenDiMuonPsi =genvector3.Phi() - rgRpAng;
            if(GenDiMuonPsi < -TMath::Pi()) GenDiMuonPsi += 2.*TMath::Pi();
            if(GenDiMuonPsi > TMath::Pi()) GenDiMuonPsi -= 2.*TMath::Pi();
            if(GenDiMuonPsi < -TMath::Pi()/2) GenDiMuonPsi += TMath::Pi();
            if(GenDiMuonPsi > TMath::Pi()/2) GenDiMuonPsi -= TMath::Pi();
            double GenDiMuonPx=genvector3.Px();
            double GenDiMuonPy=genvector3.Py();
            double GenDiMuonPz=genvector3.Pz();

            //cout<<" gen mass "<< GenDiMuonMinv   <<" pT "<< GenDiMuonPt<<endl; 
            gJpsiMass=GenDiMuonMinv;
            gJpsiPt=GenDiMuonPt;
            gJpsiPhi=GenDiMuonPhi;
            gJpsiEta=GenDiMuonEta;
            gJpsiPsi=GenDiMuonPsi;
            gJpsiRap=GenDiMuonY;
            gJpsiPx=GenDiMuonPx;
            gJpsiPy=GenDiMuonPy;
            gJpsiPz=GenDiMuonPz;

        }
    }


    for(size_t ii = 0; ii <dimuonsPATColl.size(); ++ ii) 
    {
        const pat::CompositeCandidate &p = (dimuonsPATColl)[ii];
        const reco::Candidate *dau0 = p.daughter(0);
        const pat::Muon *mu0 = dynamic_cast<const pat::Muon *>(dau0);
        const reco::Candidate *dau1 = p.daughter(1); 
        const pat::Muon *mu1 = dynamic_cast<const pat::Muon *>(dau1);
        const pat::Muon *muonPos = 0, *muonNeg = 0;

        if(mu0->charge() > 0){ muonPos = mu0; muonNeg = mu1;}
        else if(mu0->charge() < 0){ muonPos = mu1; muonNeg = mu0;}

        //---------------------------------------- Trigger Matches -----------------------------------------//

        //to match with filter name                                                                                                                                              
        //(!pMuon->triggerObjectMatchesByFilter(fHLTFilterName).empty())                                                                                                         
        //to match with trigger name                                                                                                                                             
        //!pMuon->triggerObjectMatchesByPath("HLT_HIL1DoubleMuOpen").empty()                                                                                                     
        //!pMuon->triggerObjectMatchesByPath(fHLTFilterName).empty()

        //TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT")

        //cout<< !muonPos->triggerObjectMatchesByPath("HLT_HIL1DoubleMu0_HighQ_v1").empty()  <<" matches "<<!muonNeg->triggerObjectMatchesByPath("HLT_HIL1DoubleMu0_HighQ_v1").empty()<<endl;                                                                                                                                                                                                                  
        //Trigger matches        
        //cout<<matchPATMuon(muonPos)<<" matches "<<matchPATMuon(muonNeg)<<endl;                                                                                                                                                                                                                  
        if(!muonPos->triggerObjectMatchesByPath("HLT_HIL1DoubleMu0_HighQ_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL2DoubleMu3_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL3DoubleMuOpen_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL2Mu3_NHitQ_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL2Mu7_v1").empty() ||
           !muonPos->triggerObjectMatchesByPath("HLT_HIL2Mu15_v1").empty()) {
            muPos_matches=1;
        }
        if(!muonNeg->triggerObjectMatchesByPath("HLT_HIL1DoubleMu0_HighQ_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL2DoubleMu3_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL3DoubleMuOpen_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL2Mu3_NHitQ_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL2Mu7_v1").empty() ||
           !muonNeg->triggerObjectMatchesByPath("HLT_HIL2Mu15_v1").empty()) {
            muNeg_matches=1;
        }

        if(muonPos->isTrackerMuon()){muPos_tracker=1;}
        if(muonNeg->isTrackerMuon()){muNeg_tracker=1;}

        if(muonPos->isGlobalMuon()){muPos_global=1;}
        if(muonNeg->isGlobalMuon()){muNeg_global=1;}


        //cout<<muPos_tracker<<"    "<<muNeg_tracker<<endl;
        //cout<<"muPos_matches "<<muPos_matches<<" muNeg_matches "<<muNeg_matches<<endl; 
        //cout<<"JpsiCharge " <<p.charge()<<" JpsiMass  "<<p.mass()<<endl;


        // write out JPsi RECO information
        //cout<<" inside loop RecJPsiSize "<<RecJPsiSize<<endl;

        JpsiCharge  = p.charge();
        JpsiMass =p.mass();
        JpsiPt =p.pt();
        JpsiPhi =p.phi();
        JpsiEta =p.eta();
        JpsiRap =p.rapidity();
        JpsiPx =p.px();
        JpsiPy =p.py();
        JpsiPz =p.pz();

        for(int l = 0; l < 38; l++){
            JpsiPsi[l] = -999.0;
            double del_ang = p.phi() - rpAng[l];
            if(del_ang < -TMath::Pi()) del_ang += 2.*TMath::Pi();
            if(del_ang > TMath::Pi()) del_ang -= 2.*TMath::Pi();
            if(del_ang < -TMath::Pi()/2) del_ang +=TMath::Pi();
            if(del_ang > TMath::Pi()/2) del_ang -=TMath::Pi();
            JpsiPsi[l] = del_ang;
        }

        for(int l = 0; l < 38; l++){
            JpsiNfPsi[l] = -999.0;
            double del_ang = p.phi() - NfRpAng[l];
            if(del_ang < -TMath::Pi()) del_ang += 2.*TMath::Pi();
            if(del_ang > TMath::Pi()) del_ang -= 2.*TMath::Pi();
            if(del_ang < -TMath::Pi()/2) del_ang +=TMath::Pi();
            if(del_ang > TMath::Pi()/2) del_ang -=TMath::Pi();
            JpsiNfPsi[l] = del_ang;
        }

        JpsiGenPsi = -999.0;
        double gdel_ang = p.phi() - rgRpAng;
        if(gdel_ang < -TMath::Pi()) gdel_ang += 2.*TMath::Pi();
        if(gdel_ang > TMath::Pi()) gdel_ang -= 2.*TMath::Pi();
        if(gdel_ang < -TMath::Pi()/2) gdel_ang +=TMath::Pi();
        if(gdel_ang > TMath::Pi()/2) gdel_ang -=TMath::Pi();

        JpsiGenPsi = gdel_ang;

        JpsiVprob =p.userFloat("vProb");
        RecoCtau=10.0*p.userFloat("ppdlPV");
        RecoCtauErr=10.0*p.userFloat("ppdlErrPV");
        RecoCtauTrue=10.*p.userFloat("ppdlTrue");

        // write out Muon RECO information                                                                                                                                                                                                     
        float f_muPosPx, f_muPosPy, f_muPosPz, f_muPosEta, f_muPosPhi;
        float f_muNegPx, f_muNegPy, f_muNegPz, f_muNegEta, f_muNegPhi;

        f_muPosPx = muonPos->px();
        f_muPosPy = muonPos->py();
        f_muPosPz = muonPos->pz();
        f_muPosEta = muonPos->eta();
        f_muPosPhi= muonPos->phi();

        f_muNegPx = muonNeg->px();
        f_muNegPy = muonNeg->py();
        f_muNegPz = muonNeg->pz();
        f_muNegEta = muonNeg->eta();
        f_muNegPhi = muonNeg->phi();

        muPosPx= f_muPosPx ;
        muPosPy= f_muPosPy ;
        muPosPz= f_muPosPz ;
        muPosEta= f_muPosEta;
        muPosPhi= f_muPosPhi;

        muNegPx= f_muNegPx ;
        muNegPy= f_muNegPy ;
        muNegPz= f_muNegPz ;
        muNegEta= f_muNegEta;
        muNegPhi= f_muNegPhi;

        //-----------------------------------------------------                                                                                                                                                                                
        //-----------additional Reco Muon Variables------------                                                                                                                                                                                
        //----------------------------------------------------- 



        //1.Positive Muon                                                                                                                                                                                                                      
        if(muonPos->isTrackerMuon())
        {
            TrackRef iTrack =muonPos->innerTrack();
            const reco::HitPattern& p1=iTrack->hitPattern();
            muPos_found=iTrack->found();
            muPos_nchi2In=iTrack->chi2()/iTrack->ndof();
            muPos_arbitrated=muonPos->muonID("TrackerMuonArbitrated");
            muPos_pixeLayers=p1.pixelLayersWithMeasurement();
            muPos_dxy=iTrack->dxy(RefVtx);
            muPos_dz=iTrack->dz(RefVtx);
            if(muonPos->isGlobalMuon())
            {
                TrackRef gTrack =muonPos->globalTrack();
                const reco::HitPattern& q1=gTrack->hitPattern();
                muPos_nValidMuHits=q1.numberOfValidMuonHits();
                muPos_nchi2Gl=gTrack->chi2()/gTrack->ndof();
            }
        }



        //2.Negative Muobn                                                                                                                                                                                                                     
        if(muonNeg->isTrackerMuon())
        {
            TrackRef iTrack =muonNeg->innerTrack();
            const reco::HitPattern& p2=iTrack->hitPattern();
            muNeg_found=iTrack->found();
            muNeg_nchi2In=iTrack->chi2()/iTrack->ndof();
            muNeg_arbitrated=muonNeg->muonID("TrackerMuonArbitrated");
            muNeg_pixeLayers=p2.pixelLayersWithMeasurement();
            muNeg_dxy=iTrack->dxy(RefVtx);
            muNeg_dz=iTrack->dz(RefVtx);

            if(muonNeg->isGlobalMuon())
            {
                TrackRef gTrack =muonNeg->globalTrack();
                const reco::HitPattern& q2=gTrack->hitPattern();
                muNeg_nValidMuHits=q2.numberOfValidMuonHits();
                muNeg_nchi2Gl=gTrack->chi2()/gTrack->ndof();
            }
        }     

        // Cowboy 
        double dPhi2mu = muPosPhi - muNegPhi;
        while (dPhi2mu > TMath::Pi()) dPhi2mu -= 2*TMath::Pi();
        while (dPhi2mu <= -TMath::Pi()) dPhi2mu += 2*TMath::Pi();
        double chkCowboy = 1*dPhi2mu;

        int nbit7 = 0;
        if((chkCowboy > 0.)) {
            nbit7=1;
        }
        hCowboy = nbit7;


        SingleMuonTree->Fill();
        //RecJPsiSize++;
        //cout<<"RecJPsiSiZe " <<RecJPsiSize<<endl;
    }

    //cout<<" fill tree called "<<endl;

}

void DiMuonOnia2DPlots::FillGenTree(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    GenJpsiMass=-9999.;
    GenJpsiPt=-9999.;
    GenJpsiPhi=-9999.; // dmoon add
    GenJpsiPsi=-9999.; // dmoon add
    GenJpsiRap=-9999.;
    GenJpsiEta=-9999.;
    GenJpsiPx=-9999.;
    GenJpsiPy=-9999.;
    GenJpsiPz=-9999.;

    GenmuPosPx=-9999.;
    GenmuPosPy=-9999.;
    GenmuPosPz=-9999.;
    GenmuPosEta=-9999.;
    GenmuPosPhi=-9999.;

    GenmuNegPx=-9999.;
    GenmuNegPy=-9999.;
    GenmuNegPz=-9999.;
    GenmuNegEta=-9999.; 
    GenmuNegPhi=-9999.; 


    //reset EVENT information  
    GeventNb= 0 ;
    GrunNb= 0 ;
    GlumiBlock= 0 ;
    gbin=0;

    // dmoon add
    gRpAng = -9999.;
    hi_b = -9999.; hi_npart = -9999.; hi_ncoll = -9999.; hi_nhard = -9999.;

    // Event related infos 

    GeventNb= iEvent.id().event();
    GrunNb=iEvent.id().run();
    GlumiBlock= iEvent.luminosityBlock();


    centrality_ = new CentralityProvider(iSetup);
    centrality_->newEvent(iEvent,iSetup);
    gbin = centrality_->getBin();
    //cout<<" gbin "<<gbin<<endl;


    //-----------------------------------------------------hiGenParticle----------------------------------------------------------------------

    float mumass =0.105658389;
    using namespace edm;
    using namespace std;


    // dmoon add gen reaction plane
    const GenEvent* evt;
    edm::Handle<HepMCProduct> mc;
    iEvent.getByLabel("generator",mc);
    evt = mc->GetEvent();
    const HeavyIon* hi = evt->heavy_ion();
    hi_b = hi->impact_parameter();
    hi_npart = hi->Npart_proj()+hi->Npart_targ();
    hi_ncoll = hi->Ncoll();
    hi_nhard = hi->Ncoll_hard();
    gRpAng = hi->event_plane_angle();
    if(gRpAng < -TMath::Pi()) gRpAng += 2.*TMath::Pi();
    if(gRpAng > TMath::Pi()) gRpAng -= 2.*TMath::Pi();
    if(gRpAng < -TMath::Pi()/2) gRpAng += TMath::Pi();
    if(gRpAng > TMath::Pi()/2) gRpAng -= TMath::Pi();
    //cout<<" gRpAng : "<<gRpAng<<endl;

    Handle<TriggerResults> trigResults;
    iEvent.getByLabel(triggerResults_,trigResults);

    TriggerNames triggerNames_;
    std::string bit1 = "HLT_HIL1DoubleMu0_HighQ_v1";
    std::string bit2 = "HLT_HIL2DoubleMu3_v1";
    std::string bit3 = "HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1";
    std::string bit4 = "HLT_HIL3DoubleMuOpen_v1";
    std::string bit5 = "HLT_HIL2Mu7_v1";
    std::string bit6 = "HLT_HIL2Mu15_v1";
    /*
    Trigger Names : HLT_HIL1DoubleMuOpen_v1
    Trigger Names : HLT_HIL1DoubleMu0_HighQ_v1
    Trigger Names : HLT_HIL2Mu3_v1
    Trigger Names : HLT_HIL2Mu3_NHitQ_v1
    Trigger Names : HLT_HIL2Mu7_v1
    Trigger Names : HLT_HIL2Mu15_v1
    Trigger Names : HLT_HIL2DoubleMu0_v1
    Trigger Names : HLT_HIL2DoubleMu0_NHitQ_v1
    Trigger Names : HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1
    Trigger Names : HLT_HIL2DoubleMu3_v1
    Trigger Names : HLT_HIL3Mu3_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_SS_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_OS_v1
    Trigger Names : HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1
    */
    int nbit1=0; int nbit2=0; int nbit3=0; int nbit4=0; int nbit5=0; int nbit6=0;

    if (trigResults.isValid()) {
        int ntrigs = trigResults->size();
        if (ntrigs==0){
            //std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;
        }
        edm::TriggerNames const& triggerNames = iEvent.triggerNames(*trigResults);
        for (int itrig = 0; itrig != ntrigs; ++itrig){
            std::string trigName=triggerNames.triggerName(itrig);
            bool accept = trigResults->accept(itrig);
            //cout<<"Trigger Names : "<<trigName<<" bit : "<<accept<<endl;
            if (accept) {
                if(trigName == bit1){ nbit1 = 1;}
                if(trigName == bit2){ nbit2 = 1;}
                if(trigName == bit3){ nbit3 = 1;}
                if(trigName == bit4){ nbit4 = 1;}
                if(trigName == bit5){ nbit5 = 1;}
                if(trigName == bit6){ nbit6 = 1;}
                /*
                if(trigName == bit1){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit1 = 1;}
                if(trigName == bit2){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit2 = 1;}
                if(trigName == bit3){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit3 = 1;}
                if(trigName == bit4){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit4 = 1;}
                if(trigName == bit5){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit5 = 1;}
                if(trigName == bit6){ cout<<"%HLTInfo : "<<trigName<<" is fired !!!!"<<endl; nbit6 = 1;}
                */
            }
        }
    }else{cout<<"HLT Info is not ok !!!!!"<<endl;}

    ghbit1 = nbit1; ghbit2 = nbit2; ghbit3 = nbit3;
    ghbit4 = nbit4; ghbit5 = nbit5; ghbit6 = nbit6;

    mumass =0.105658389;

    edm::Handle<edm::View<reco::GenParticle> >genPar;
    iEvent.getByLabel("hiGenParticles",genPar) ;

    if(!(genPar.isValid())) return;

    edm::View<reco::GenParticle> genParticles = *genPar ;
    TLorentzVector  genvector1, genvector2;

    double px1[10000], py1[10000], pz1[10000], px2[10000], py2[10000], pz2[10000];
    unsigned int nplus = 0, nminus =0;

    for(size_t i = 0; i < genParticles.size(); ++ i) {
        const reco::GenParticle& part = (*genPar)[i];
        const  Candidate *mom = part.mother();


        if (part.numberOfMothers()!=1) continue;
        int momId = mom->pdgId();

        //!strcmp(fIsPATInfo.c_str(),"TRUE")

        char* MID=(char*)"AB";

        if(momId == 443) MID=(char *)"JPsi";
        if(momId == 553) MID=(char *)"Upsilon1s";
        if(momId == 100553)MID=(char *)"Upsilon2s";



        if ((abs(part.pdgId()) == 13) && ( !strcmp(fMotherID.c_str(),MID)) ){




            if(part.pdgId() == 13 ){
                px1[nplus] = part.px();
                py1[nplus] = part.py();
                pz1[nplus] = part.pz();
                nplus++;

                GenmuNegPx=part.px();
                GenmuNegPy=part.py();
                GenmuNegPz=part.pz();
                GenmuNegEta=part.eta();
                GenmuNegPhi=part.phi();

                //cout<<"motherID "<<MID<<endl;


            }

            if(part.pdgId()== -13) {
                px2[nminus] = part.px();
                py2[nminus] = part.py();
                pz2[nminus] = part.pz();
                nminus++;

                GenmuPosPx=part.px();
                GenmuPosPy=part.py();
                GenmuPosPz=part.pz();
                GenmuPosEta=part.eta();
                GenmuPosPhi=part.phi();

            }
        }
    }
    //cout<<" nplus : "<<nplus<<"  "<<nminus<<endl;

    for(size_t i = 0; i < nplus; i++) {
        double en1 = sqrt(px1[i]*px1[i] + py1[i]*py1[i] + pz1[i]*pz1[i] + mumass*mumass);

        for(size_t j = 0; j< nminus; j++) {
            double en2 = sqrt(px2[j]*px2[j] + py2[j]*py2[j] + pz2[j]*pz2[j] + mumass*mumass);
            TLorentzVector  genvector1,genvector2,genvector3;

            genvector1.SetPxPyPzE(px1[i], py1[i], pz1[i], en1);
            genvector2.SetPxPyPzE(px2[j], py2[j], pz2[j], en2);

            genvector3=genvector1+genvector2;

            double GenDiMuonY=genvector3.Rapidity();
            double GenDiMuonMinv=genvector3.M();
            double GenDiMuonPt =genvector3.Pt();
            double GenDiMuonPhi =genvector3.Phi();
            double GenDiMuonEta =genvector3.Eta();
            double GenDiMuonPsi =genvector3.Phi() - gRpAng;
            if(GenDiMuonPsi < -TMath::Pi()) GenDiMuonPsi += 2.*TMath::Pi();
            if(GenDiMuonPsi > TMath::Pi()) GenDiMuonPsi -= 2.*TMath::Pi();
            if(GenDiMuonPsi < -TMath::Pi()/2) GenDiMuonPsi += TMath::Pi();
            if(GenDiMuonPsi > TMath::Pi()/2) GenDiMuonPsi -= TMath::Pi();
            double GenDiMuonPx=genvector3.Px();
            double GenDiMuonPy=genvector3.Py();
            double GenDiMuonPz=genvector3.Pz();

            //cout<<" gen mass "<< GenDiMuonMinv   <<" pT "<< GenDiMuonPt<<endl; 
            GenJpsiMass=GenDiMuonMinv;
            GenJpsiPt=GenDiMuonPt;
            GenJpsiPhi=GenDiMuonPhi;
            GenJpsiEta=GenDiMuonEta;
            GenJpsiPsi=GenDiMuonPsi;
            GenJpsiRap=GenDiMuonY;
            GenJpsiPx=GenDiMuonPx;
            GenJpsiPy=GenDiMuonPy;
            GenJpsiPz=GenDiMuonPz;

            //SingleGenMuonTree->Fill();
            //cout<<"gen Tree Filled "<<endl;

        }
    }

    // Cowboy 
    double gdPhi2mu = GenmuPosPhi - GenmuNegPhi;
    while (gdPhi2mu > TMath::Pi()) gdPhi2mu -= 2*TMath::Pi();
    while (gdPhi2mu <= -TMath::Pi()) gdPhi2mu += 2*TMath::Pi();
    double gchkCowboy = 1*gdPhi2mu;

    int gnbit7 = 0;
    if((gchkCowboy > 0.)) {
        gnbit7=1;
    }
    ghCowboy = gnbit7;

    SingleGenMuonTree->Fill();
    //cout<<"gen Tree Filled "<<endl;

}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonOnia2DPlots);



