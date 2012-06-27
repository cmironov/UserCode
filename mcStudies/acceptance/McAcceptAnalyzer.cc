/* This is en example for an Analyzer of a MCatNLO HeoMCProduct
   The analyzer fills a histogram with the event weights,
   and looks for electron/positron pairs and fills a histogram
   with the invaraint mass of the two. 
   */

//
// Original Author:  Fabian Stoeckli
//         Created:  Tue Nov 14 13:43:02 CET 2006
// $Id: ZmmAnalyzerv4_yz20.cc,v 1.1 2012/06/16 08:27:59 hckim Exp $
//
//


// system include files
#include <memory>
#include <iostream>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/WeightContainer.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "TH1D.h"
#include "TH2D.h"

#include "TFile.h"

#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include <TSystem.h>
//
// class decleration
//

class McAcceptAnalyzer : public edm::EDAnalyzer {
    public:
        explicit McAcceptAnalyzer(const edm::ParameterSet&);
        ~McAcceptAnalyzer();


    private:
        std::string fOutputFile_;
        TFile *hOutputFile;
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------

        std::string outputFilename;
        TH1D* weight_histo;
        TH1D* invmass_histo;
        TH1D* Zpt;
        TH1D* hardpt;
        TH1D* softpt;
        TH1D* hardeta;
        TH1D* softeta;
        TH1D* hardphi;
        TH1D* softphi;

        TH1D* AllDimuM;
        TH1D* AllDimuPt;
        TH1D* AllDimuEta;
        TH1D* AllDimuY;
        TH1D* AllDimuPhi;
        TH1D* AllDimuMA;
        TH1D* AllDimuPtA;
        TH1D* AllDimuEtaA;
        TH1D* AllDimuYA;
        TH1D* AllDimuPhiA;


        TH1D* DenDimuM;
        TH1D* DenDimuPt;
        TH1D* DenDimuEta;
        TH1D* DenDimuY;
        TH1D* DenDimuPhi;
        TH1D* NeuDimuM;
        TH1D* NeuDimuPt;
        TH1D* NeuDimuEta;
        TH1D* NeuDimuY;
        TH1D* NeuDimuPhi;

          TH1D* AllSimuPt;
          TH1D* DenSimuPt;
       TH1D* NeuSimuPt;
          TH1D* AllSimuEta;
          TH1D* DenSimuEta;
       TH1D* NeuSimuEta;

         TH1D* Den2DimuM;
        TH1D* Den2DimuPt;
        TH1D* Den2DimuEta;
        TH1D* Den2DimuY;
        TH1D* Den2DimuPhi;
        TH1D* Neu2DimuM;
        TH1D* Neu2DimuPt;
        TH1D* Neu2DimuEta;
        TH1D* Neu2DimuY;
        TH1D* Neu2DimuPhi;
        TH1D* Den2DimuMA;
        TH1D* Den2DimuPtA;
        TH1D* Den2DimuEtaA;
        TH1D* Den2DimuYA;
        TH1D* Den2DimuPhiA;
        TH1D* Neu2DimuMA;
        TH1D* Neu2DimuPtA;
        TH1D* Neu2DimuEtaA;
        TH1D* Neu2DimuYA;
        TH1D* Neu2DimuPhiA;
        TH1D* DenDimuMv;
        TH1D* DenDimuPtv;
        TH1D* DenDimuEtav;
        TH1D* DenDimuYv;
        TH1D* DenDimuPhiv;
        TH1D* NeuDimuMv;
        TH1D* NeuDimuPtv;
        TH1D* NeuDimuEtav;
        TH1D* NeuDimuYv;
        TH1D* NeuDimuPhiv;
        TH1D* DenDimuMA;
        TH1D* DenDimuPtA;
        TH1D* DenDimuEtaA;
        TH1D* DenDimuYA;
        TH1D* DenDimuPhiA;
        TH1D* NeuDimuMA;
        TH1D* NeuDimuPtA;
        TH1D* NeuDimuEtaA;
        TH1D* NeuDimuYA;
        TH1D* NeuDimuPhiA;
        TH1D* AccDimuM;
        TH1D* AccDimuPt;
        TH1D* AccDimuEta;

        TH2D* TwoAllDimuYPt;
        TH2D* TwoAllDimuEtaPt;
        TH2D* TwoDenDimuYPt;
        TH2D* TwoDenDimuEtaPt;
        TH2D* TwoNeuDimuYPt;
        TH2D* TwoNeuDimuEtaPt;
        TH2D* TwoDen2DimuYPt;
        TH2D* TwoDen2DimuEtaPt;

        TH2D* TwoAllDimuPtM;
        TH2D* TwoAllDimuYM;
        TH2D* TwoDenDimuPtM;
        TH2D* TwoDenDimuYM;
        TH2D* TwoNeuDimuPtM;
        TH2D* TwoNeuDimuYM;
        TH2D* TwoDen2DimuPtM;
        TH2D* TwoDen2DimuYM;

        double sumWeights;

};


McAcceptAnalyzer::McAcceptAnalyzer(const edm::ParameterSet& iConfig)
{

  fOutputFile_ = iConfig.getUntrackedParameter<std::string>("hOutputFile");
  sumWeights=0.0;


    weight_histo  = new TH1D("weight_histo","weight_histo",800,-400,400);
    invmass_histo = new TH1D("invmass_histo","invmass_histo",200,0,200);
    Zpt = new TH1D("Zpt","Zpt",200,0,200);
    hardpt = new TH1D("hardpt","hardpt",200,0,200);
    softpt = new TH1D("softpt","softpt",200,0,200);
    hardeta = new TH1D("hardeta","hardeta",200,-10,10);
    softeta = new TH1D("softeta","softeta",200,-10,10);

    AllDimuM = new TH1D("AllDimuM","AllDimuM",200,0,200);
    AllDimuPt = new TH1D("AllDimuPt","AllDimuPt",200,0,200);
    AllDimuEta = new TH1D("AllDimuEta","AllDimuEta",200,-10,10);
    AllDimuY = new TH1D("AllDimuY","AllDimuY",200,-10,10);
    AllDimuPhi = new TH1D("AllDimuPhi","AllDimuPhi",70,-3.5,3.5);

    AllDimuMA = new TH1D("AllDimuMA","AllDimuMA",200,0,200);
    AllDimuPtA = new TH1D("AllDimuPtA","AllDimuPtA",200,0,200);
    AllDimuEtaA = new TH1D("AllDimuEtaA","AllDimuEtaA",200,-10,10);
    AllDimuYA = new TH1D("AllDimuYA","AllDimuYA",200,-10,10);
    AllDimuPhiA = new TH1D("AllDimuPhiA","AllDimuPhiA",70,-3.5,3.5);

    DenDimuM = new TH1D("DenDimuM","DenDimuM",200,0,200);
    DenDimuPt = new TH1D("DenDimuPt","DenDimuPt",200,0,200);
    DenDimuEta = new TH1D("DenDimuEta","DenDimuEta",200,-10,10);
    DenDimuY = new TH1D("DenDimuY","DenDimuY",200,-10,10);
    DenDimuPhi = new TH1D("DenDimuPhi","DenDimuPhi",70,-3.5,3.5);

      AllSimuPt = new TH1D("AllSimuPt","AllSimuPt",200,0,200);
       DenSimuPt = new TH1D("DenSimuPt","DenSimuPt",200,0,200);
  NeuSimuPt = new TH1D("NeuSimuPt","NeuSimuPt",200,0,200);
      AllSimuEta = new TH1D("AllSimuEta","AllSimuEta",200,0,200);
       DenSimuEta = new TH1D("DenSimuEta","DenSimuEta",200,0,200);
  NeuSimuEta = new TH1D("NeuSimuEta","NeuSimuEta",200,0,200);

     NeuDimuM = new TH1D("NeuDimuM","NeuDimuM",200,0,200);
    NeuDimuPt = new TH1D("NeuDimuPt","NeuDimuPt",200,0,200);
    NeuDimuEta = new TH1D("NeuDimuEta","NeuDimuEta",200,-10,10);
    NeuDimuY = new TH1D("NeuDimuY","NeuDimuY",200,-10,10);
    NeuDimuPhi = new TH1D("NeuDimuPhi","NeuDimuPhi",70,-3.5,3.5);

    Den2DimuM = new TH1D("Den2DimuM","Den2DimuM",200,0,200);
    Den2DimuPt = new TH1D("Den2DimuPt","Den2DimuPt",200,0,200);
    Den2DimuEta = new TH1D("Den2DimuEta","Den2DimuEta",200,-10,10);
    Den2DimuY = new TH1D("Den2DimuY","Den2DimuY",200,-10,10);
    Den2DimuPhi = new TH1D("Den2DimuPhi","Den2DimuPhi",70,-3.5,3.5);

    Neu2DimuM = new TH1D("Neu2DimuM","Neu2DimuM",200,0,200);
    Neu2DimuPt = new TH1D("Neu2DimuPt","Neu2DimuPt",200,0,200);
    Neu2DimuEta = new TH1D("Neu2DimuEta","Neu2DimuEta",200,-10,10);
    Neu2DimuY = new TH1D("Neu2DimuY","Neu2DimuY",200,-10,10);
    Neu2DimuPhi = new TH1D("Neu2DimuPhi","Neu2DimuPhi",70,-3.5,3.5);

    Den2DimuMA = new TH1D("Den2DimuMA","Den2DimuMA",200,0,200);
    Den2DimuPtA = new TH1D("Den2DimuPtA","Den2DimuPtA",200,0,200);
    Den2DimuEtaA = new TH1D("Den2DimuEtaA","Den2DimuEtaA",200,-10,10);
    Den2DimuYA = new TH1D("Den2DimuYA","Den2DimuYA",200,-10,10);
    Den2DimuPhiA = new TH1D("Den2DimuPhiA","Den2DimuPhiA",70,-3.5,3.5);

    Neu2DimuMA = new TH1D("Neu2DimuMA","Neu2DimuMA",200,0,200);
    Neu2DimuPtA = new TH1D("Neu2DimuPtA","Neu2DimuPtA",200,0,200);
    Neu2DimuEtaA = new TH1D("Neu2DimuEtaA","Neu2DimuEtaA",200,-10,10);
    Neu2DimuYA = new TH1D("Neu2DimuYA","Neu2DimuYA",200,-10,10);
    Neu2DimuPhiA = new TH1D("Neu2DimuPhiA","Neu2DimuPhiA",70,-3.5,3.5);

    DenDimuMv = new TH1D("DenDimuMv","DenDimuMv",200,0,200);
    DenDimuPtv = new TH1D("DenDimuPtv","DenDimuPtv",200,0,200);
    DenDimuEtav = new TH1D("DenDimuEtav","DenDimuEtav",200,-10,10);
    DenDimuYv = new TH1D("DenDimuYv","DenDimuYv",200,-10,10);
    DenDimuPhiv = new TH1D("DenDimuPhiv","DenDimuPhiv",70,-3.5,3.5);
    NeuDimuMv = new TH1D("NeuDimuM","NeuDimuM",200,0,200);
    NeuDimuPtv = new TH1D("NeuDimuPt","NeuDimuPt",200,0,200);
    NeuDimuEtav = new TH1D("NeuDimuEta","NeuDimuEta",200,-10,10);
    NeuDimuYv = new TH1D("NeuDimuY","NeuDimuY",200,-10,10);
    NeuDimuPhiv = new TH1D("NeuDimuPhi","NeuDimuPhi",70,-3.5,3.5);
    DenDimuMA = new TH1D("DenDimuM","DenDimuM",200,0,200);
    DenDimuPtA = new TH1D("DenDimuPtA","DenDimuPtA",200,0,200);
    DenDimuEtaA = new TH1D("DenDimuEtaA","DenDimuEtaA",200,-10,10);
    DenDimuYA = new TH1D("DenDimuYA","DenDimuYA",200,-10,10);
    DenDimuPhiA = new TH1D("DenDimuPhiA","DenDimuPhiA",70,-3.5,3.5);
    NeuDimuMA = new TH1D("NeuDimuMA","NeuDimuMA",200,0,200);
    NeuDimuPtA = new TH1D("NeuDimuPtA","NeuDimuPtA",200,0,200);
    NeuDimuEtaA = new TH1D("NeuDimuEtaA","NeuDimuMA",200,-10,10);
    NeuDimuYA = new TH1D("NeuDimuYA","NeuDimuYA",200,-10,10);
    NeuDimuPhiA = new TH1D("NeuDimuPhiA","NeuDimuPhiA",70,-3.5,3.5);

    TwoAllDimuYPt = new TH2D("TwoAllDimuYPt","TwoAllDimuYPt",200,-10,10,200,0,200);
    TwoAllDimuEtaPt = new TH2D("TwoAllDimuEtaPt","TwoAllDimuEtaPt",200,-10,10,200,0,200);
    TwoDenDimuYPt = new TH2D("TwoDenDimuYPt","TwoDenDimuYPt",200,-10,10,200,0,200);
    TwoDenDimuEtaPt = new TH2D("TwoDenDimuEtaPt","TwoDenDimuEtaPt",200,-10,10,200,0,200);
    TwoNeuDimuYPt = new TH2D("TwoNeuDimuYPt","TwoNeuDimuYPt",200,-10,10,200,0,200);
    TwoNeuDimuEtaPt = new TH2D("TwoNeuDimuEtaPt","TwoNeuDimuEtaPt",200,-10,10,200,0,200);
    TwoDen2DimuYPt = new TH2D("TwoDen2DimuYPt","TwoDen2DimuYPt",200,-10,10,200,0,200);
    TwoDen2DimuEtaPt = new TH2D("TwoDen2DimuEtaPt","TwoDen2DimuEtaPt",200,-10,10,200,0,200);

    TwoAllDimuYM = new TH2D("TwoAllDimuYM","TwoAllDimuYM",200,-10,10,200,0,200);
    TwoAllDimuPtM = new TH2D("TwoAllDimuPtM","TwoAllDimuPtM",200,0,200,200,0,200);
    TwoDenDimuYM = new TH2D("TwoDenDimuYM","TwoDenDimuYM",200,-10,10,200,0,200);
    TwoDenDimuPtM = new TH2D("TwoDenDimuPtM","TwoDenDimuPtM",200,0,200,200,0,200);
    TwoNeuDimuYM = new TH2D("TwoNeuDimuYM","TwoNeuDimuYM",200,-10,10,200,0,200);
    TwoNeuDimuPtM = new TH2D("TwoNeuDimuPtM","TwoNeuDimuPtM",200,0,200,200,0,200);
    TwoDen2DimuYM = new TH2D("TwoDen2DimuYM","TwoDen2DimuYM",200,-10,10,200,0,200);
    TwoDen2DimuPtM = new TH2D("TwoDen2DimuPtM","TwoDen2DimuPtM",200,0,200,200,0,200);

}


McAcceptAnalyzer::~McAcceptAnalyzer()
{

}


// ------------ method called to for each event  ------------
void McAcceptAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;


    // get HepMC::GenEvent ...
    Handle<HepMCProduct> evt_h;
    iEvent.getByType(evt_h);
    HepMC::GenEvent* evt = new  HepMC::GenEvent(*(evt_h->GetEvent()));

    Handle<GenEventInfoProduct> evt_info;
    iEvent.getByType(evt_info);


    double weight = evt_info->weight();
    if(weight) weight_histo->Fill(weight);

    // look for stable muons/positrons
    std::vector<HepMC::GenParticle*> muons;   
    for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) 
      {
    
        if(abs((*it)->pdg_id())==13 && (*it)->status()==1)
	  muons.push_back(*it);
      }

    // if there are at least two muons/positrons, 
    // calculate invarant mass of first two and fill it into histogram

    double inv_mass = 0.0;
    double Zpt_ = 0.0;
    double Z_phi = 999.0;
    double Z_eta = -999.0;
    double Z_pz = 0.0;
    double Z_e = 0.0;
    double Z_y = 0.0;
    if(muons.size()>1) {
        math::XYZTLorentzVector tot_momentum(muons[0]->momentum());
        math::XYZTLorentzVector mom2(muons[1]->momentum());
        tot_momentum += mom2;
        inv_mass = sqrt(tot_momentum.mass2());
        Zpt_=tot_momentum.pt();
        Z_phi=tot_momentum.phi();
        Z_eta=tot_momentum.eta();
        Z_pz=tot_momentum.pz();
        Z_e=tot_momentum.e();
        Z_y=(0.5)*log((Z_e+Z_pz)/(Z_e-Z_pz));
        // IMPORTANT: use the weight of the event ...

        double weight_sign = (weight > 0) ? 1. : -1.;
        std::cout << "weight: " << weight << std::endl;
        invmass_histo->Fill(inv_mass,weight_sign);
        Zpt->Fill(Zpt_,weight_sign);
        AllDimuM->Fill(inv_mass,weight_sign);
        AllDimuPt->Fill(Zpt_,weight_sign);
        AllDimuEta->Fill(Z_eta,weight_sign);
        AllDimuY->Fill(Z_y,weight_sign);
        AllDimuPhi->Fill(Z_phi,weight_sign);

        TwoAllDimuYPt->Fill(Z_y,Zpt_,weight_sign);
        TwoAllDimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
        TwoAllDimuPtM->Fill(Zpt_,inv_mass,weight_sign);
        TwoAllDimuYM->Fill(Z_y,inv_mass,weight_sign);

        AllDimuMA->Fill(inv_mass);
        AllDimuPtA->Fill(Zpt_);
        AllDimuEtaA->Fill(Z_eta);
        AllDimuYA->Fill(Z_y);
        AllDimuPhiA->Fill(Z_phi);
        if(60 < inv_mass && inv_mass < 120)
	  {
            AllSimuPt->Fill(muons[0]->momentum().perp(),weight_sign);
	    AllSimuPt->Fill(muons[1]->momentum().perp(),weight_sign);
	    AllSimuEta->Fill(muons[0]->momentum().eta(),weight_sign);
            AllSimuEta->Fill(muons[1]->momentum().eta(),weight_sign);
	  }
        if(muons[0]->momentum().perp()>muons[1]->momentum().perp()) 
	  {
            hardpt->Fill(muons[0]->momentum().perp(),weight_sign);
            softpt->Fill(muons[1]->momentum().perp(),weight_sign);
            hardeta->Fill(muons[0]->momentum().eta(),weight_sign);
            softeta->Fill(muons[1]->momentum().eta(),weight_sign);
	  } else 
	  {
            hardpt->Fill(muons[1]->momentum().perp(),weight_sign);
            softpt->Fill(muons[0]->momentum().perp(),weight_sign);       
            hardeta->Fill(muons[1]->momentum().eta(),weight_sign);
            softeta->Fill(muons[0]->momentum().eta(),weight_sign);
	  }

        //if(60 < inv_mass && inv_mass < 120 && (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1) {
        //if(60 < inv_mass && inv_mass < 120 && (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1 && fabs(Z_y)<2.4) {   
         if(60 < inv_mass && inv_mass < 120 && (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1 && fabs(Z_y)<2.0) 
	   {   
	     DenDimuM->Fill(inv_mass,weight_sign);
	     DenDimuPt->Fill(Zpt_,weight_sign);
	     DenDimuEta->Fill(Z_eta,weight_sign);
	     DenDimuY->Fill(Z_y,weight_sign);
	     DenDimuPhi->Fill(Z_phi,weight_sign);
	     TwoDenDimuYPt->Fill(Z_y,Zpt_,weight_sign);
	     TwoDenDimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
	     TwoDenDimuPtM->Fill(Zpt_,inv_mass,weight_sign);
	     TwoDenDimuYM->Fill(Z_y,inv_mass,weight_sign);
	     
	     DenDimuMA->Fill(inv_mass);
	     DenDimuPtA->Fill(Zpt_);
	     DenDimuEtaA->Fill(Z_eta);
	     DenDimuYA->Fill(Z_y);
	     DenDimuPhiA->Fill(Z_phi);
	     
	     //DenSimuM->Fill(muons[0]inv_mass,weight_sign);
	     DenSimuPt->Fill(muons[0]->momentum().perp(),weight_sign);
	     DenSimuPt->Fill(muons[1]->momentum().perp(),weight_sign);
	     DenSimuEta->Fill(muons[0]->momentum().eta(),weight_sign);
	     DenSimuEta->Fill(muons[1]->momentum().eta(),weight_sign);
	     //mu1_e=muon[0]->momentum().e();
	     //mu1_pz=muon[0]->momentum().pz();
	     //mu1_y=(0.5)*log((mu1_e+mu1_pz)/(mu1_e-mu1_pz));
	     //DenSimuY->Fill(Z_y,weight_sign);
             //mu2_e=muon[1]->momentum().e();
	     //mu2_pz=muon[1]->momentum().pz();
	     //mu2_y=(0.5)*log((mu2_e+mu2_pz)/(mu2_e-mu2_pz));
	     //DenSimuPhi->Fill(Z_phi,weight_sign);
	     //TwoDenSimuYPt->Fill(Z_y,Zpt_,weight_sign);
	     //TwoDenSimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
	     //TwoDenSimuPtM->Fill(Zpt_,inv_mass,weight_sign);
	     //TwoDenSimuYM->Fill(Z_y,inv_mass,weight_sign);
	     //TwoDenSimuYPt->Fill(Z_y,Zpt_,weight_sign);
	     //TwoDenSimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
	     //TwoDenSimuPtM->Fill(Zpt_,inv_mass,weight_sign);
	     //TwoDenSimuYM->Fill(Z_y,inv_mass,weight_sign);
	     
	     /*
               DenDimuMv->Fill(inv_mass,weight_sign);
               DenDimuPtv->Fill(Zpt_,weight_sign);
               DenDimuEtav->Fill(tot_momentum.eta(),weight_sign);
               */
        }
        //if(60 < inv_mass && inv_mass < 120 && (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1 && fabs(muons[0]->momentum().eta()) < 2.4 && fabs(muons[1]->momentum().eta()) < 2.4 && fabs(muons[0]->momentum().perp()) > 10 && fabs(muons[1]->momentum().perp()) > 10) {
         if(60 < inv_mass && inv_mass < 120 && 
	    (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1 && 
	    fabs(muons[0]->momentum().eta()) < 2.4 && fabs(muons[1]->momentum().eta()) < 2.4 && 
	    fabs(muons[0]->momentum().perp()) > 10 && fabs(muons[1]->momentum().perp()) > 10 && 
	    fabs(Z_y)<2.0) 
	   {
	     NeuDimuM->Fill(inv_mass,weight_sign);
	     NeuDimuPt->Fill(Zpt_,weight_sign);
	     NeuDimuEta->Fill(tot_momentum.eta(),weight_sign);
	     NeuDimuY->Fill(Z_y,weight_sign);
	     NeuDimuPhi->Fill(Z_phi,weight_sign);
	     TwoNeuDimuYPt->Fill(Z_y,Zpt_,weight_sign);
	     TwoNeuDimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
	     TwoNeuDimuPtM->Fill(Zpt_,inv_mass,weight_sign);
	     TwoNeuDimuYM->Fill(Z_y,inv_mass,weight_sign);
	     
	     NeuDimuMA->Fill(inv_mass);
	     NeuDimuPtA->Fill(Zpt_);
	     NeuDimuEtaA->Fill(Z_eta);
	     NeuDimuYA->Fill(Z_y);
	     NeuDimuPhiA->Fill(Z_phi);
	     
            NeuSimuPt->Fill(muons[0]->momentum().perp(),weight_sign);
	    NeuSimuPt->Fill(muons[1]->momentum().perp(),weight_sign);
	    NeuSimuEta->Fill(muons[0]->momentum().eta(),weight_sign);
            NeuSimuEta->Fill(muons[1]->momentum().eta(),weight_sign);
	    
            /*
               NeuDimuMv->Fill(inv_mass,weight_sign);
               NeuDimuPtv->Fill(Zpt_,weight_sign);
               NeuDimuEtav->Fill(tot_momentum.eta(),weight_sign);
               */
        }
        if(60 < inv_mass && inv_mass < 120 && (muons[0]->pdg_id())*(muons[1]->pdg_id())<-1) 
	  {   
            Den2DimuM->Fill(inv_mass,weight_sign);
            Den2DimuPt->Fill(Zpt_,weight_sign);
            Den2DimuEta->Fill(Z_eta,weight_sign);
            Den2DimuY->Fill(Z_y,weight_sign);
            Den2DimuPhi->Fill(Z_phi,weight_sign);
            TwoDen2DimuYPt->Fill(Z_y,Zpt_,weight_sign);
            TwoDen2DimuEtaPt->Fill(Z_eta,Zpt_,weight_sign);
            TwoDen2DimuPtM->Fill(Zpt_,inv_mass,weight_sign);
            TwoDen2DimuYM->Fill(Z_y,inv_mass,weight_sign);

            Den2DimuMA->Fill(inv_mass);
            Den2DimuPtA->Fill(Zpt_);
            Den2DimuEtaA->Fill(Z_eta);
            Den2DimuYA->Fill(Z_y);
            Den2DimuPhiA->Fill(Z_phi);

            /*
               DenDimuMv->Fill(inv_mass,weight_sign);
               DenDimuPtv->Fill(Zpt_,weight_sign);
               DenDimuEtav->Fill(tot_momentum.eta(),weight_sign);
               */
        }
        /*
           AccDimuM->Fill(DenDimuM/NeuDimuM,weight_sign);
           AccDimuPt->Fill(DenDimuPt/NeuDimuPt,weight_sign);
           AccDimuEta->Fill(DenDimuEta/NeuDimuEta,weight_sign);
           */
        sumWeights+=weight_sign;
    }

    delete evt;
    }


    // ------------ method called once each job just before starting event loop  ------------
    void 
        McAcceptAnalyzer::beginJob()
        {
            hOutputFile = new TFile( fOutputFile_.c_str(), "RECREATE" ) ;

        }

    // ------------ method called once each job just after ending the event loop  ------------
    void 
        McAcceptAnalyzer::endJob() {

            std::cout<<" total sum wieghts = "<<sumWeights<<std::endl;

            hOutputFile->SetCompressionLevel(2);
            hOutputFile->cd();

            // save histograms into file
            //  TFile file(outputFilename.c_str(),"RECREATE");
            weight_histo->Write();
            invmass_histo->Write();
            Zpt->Write();
            hardpt->Write();
            softpt->Write();
            hardeta->Write();
            softeta->Write();
            DenDimuM->Write();
            DenDimuPt->Write();
            DenDimuEta->Write();
            DenDimuY->Write();
            DenDimuPhi->Write();
            NeuDimuM->Write();
            NeuDimuPt->Write();
            NeuDimuEta->Write();
            NeuDimuY->Write();
            NeuDimuPhi->Write();
            DenDimuMA->Write();
            DenDimuPtA->Write();
            DenDimuEtaA->Write();
            DenDimuYA->Write();
            DenDimuPhiA->Write();
            NeuDimuMA->Write();
            NeuDimuPtA->Write();
            NeuDimuEtaA->Write();
            NeuDimuYA->Write();
            NeuDimuPhiA->Write();

            Den2DimuM->Write();
            Den2DimuPt->Write();
            Den2DimuEta->Write();
            Den2DimuY->Write();
            Den2DimuPhi->Write();
            Den2DimuMA->Write();
            Den2DimuPtA->Write();
            Den2DimuEtaA->Write();
            Den2DimuYA->Write();
            Den2DimuPhiA->Write();

            Neu2DimuM->Write();
            Neu2DimuPt->Write();
            Neu2DimuEta->Write();
            Neu2DimuY->Write();
            Neu2DimuPhi->Write();
            Neu2DimuMA->Write();
            Neu2DimuPtA->Write();
            Neu2DimuEtaA->Write();
            Neu2DimuYA->Write();
            Neu2DimuPhiA->Write();

            AllDimuM->Write();
            AllDimuPt->Write();
            AllDimuEta->Write();
            AllDimuY->Write();
            AllDimuPhi->Write();
            AllDimuMA->Write();
            AllDimuPtA->Write();
            AllDimuEtaA->Write();
            AllDimuYA->Write();
            AllDimuPhiA->Write();

            TwoAllDimuYPt->Write();
            TwoAllDimuEtaPt->Write();
            TwoDenDimuYPt->Write();
            TwoDenDimuEtaPt->Write();
            TwoNeuDimuYPt->Write();
            TwoNeuDimuEtaPt->Write();
            TwoDen2DimuYPt->Write();
            TwoDen2DimuEtaPt->Write();

            TwoAllDimuPtM->Write();
            TwoAllDimuYM->Write();
            TwoDenDimuPtM->Write();
            TwoDenDimuYM->Write();
            TwoNeuDimuPtM->Write();
            TwoNeuDimuYM->Write();
            TwoDen2DimuPtM->Write();
            TwoDen2DimuYM->Write();

             AllSimuPt->Write();
           DenSimuPt->Write();
            NeuSimuPt->Write();
             AllSimuEta->Write();
           DenSimuEta->Write();
            NeuSimuEta->Write();

            double EntriesDen = DenDimuPt->GetEntries();
            double EntriesNeu = NeuDimuPt->GetEntries();
            double Acc = EntriesNeu /EntriesDen;
            std::cout << "Den : " << EntriesDen << ", Neu : " << EntriesNeu << ", *** Acc : " << Acc << std::endl;

            DenDimuPt->Rebin(10);
            NeuDimuPt->Rebin(10);
            /*
               for (int x=0;x<200;++x){
               DenDimuPtv->SetBinContent(x,fabs(DenDimuPt->GetBinContent(x)));
               NeuDimuPtv->SetBinContent(x,fabs(NeuDimuPt->GetBinContent(x)));
               std::cout << "Pt - x: " << x << ", Den: " << DenDimuPtv->GetBinContent(x) << ", Neu: " << NeuDimuPtv->GetBinContent(x) << std::endl;
            //  DenDimuPtv->Sumw2();
            //  NeuDimuPtv->Sumw2();
            }

            for (int x=0;x<200;++x){
            DenDimuEtav->SetBinContent(x,fabs(DenDimuEta->GetBinContent(x)));
            NeuDimuEtav->SetBinContent(x,fabs(NeuDimuEta->GetBinContent(x)));
            std::cout << "Eta - x: " << x << ", Den: " << DenDimuEtav->GetBinContent(x) << ", Neu: " << NeuDimuEtav->GetBinContent(x) << std::endl;
            //  DenDimuEtav->Sumw2();
            //  NeuDimuEtav->Sumw2();
            }


            DenDimuPtv->Sumw2();
            NeuDimuPtv->Sumw2();
            DenDimuEtav->Sumw2();
            NeuDimuEtav->Sumw2();

            TGraphAsymmErrors *gEffPt = new TGraphAsymmErrors();
            gEffPt->BayesDivide(NeuDimuPtv,DenDimuPtv);

            TGraphAsymmErrors *gEffEta = new TGraphAsymmErrors();
            gEffEta->BayesDivide(NeuDimuEtav,DenDimuEtav);

            EntriesDen = DenDimuPtv->GetEntries();
            EntriesNeu = NeuDimuPtv->GetEntries();
            Acc = EntriesNeu /EntriesDen;
            std::cout << "*** New Den : " << EntriesDen << ", Neu : " << EntriesNeu << ", *** Acc : " << Acc << std::endl;

            DenDimuMv->Write();
            DenDimuPtv->Write();
            DenDimuEtav->Write();
            NeuDimuMv->Write();
            NeuDimuPtv->Write();
            NeuDimuEtav->Write();
            */

            EntriesDen = DenDimuPtA->GetEntries();
            EntriesNeu = NeuDimuPtA->GetEntries();
            Acc = EntriesNeu /EntriesDen;
            std::cout << "*** ABS Den : " << EntriesDen << ", Neu : " << EntriesNeu << ", *** Acc : " << Acc << std::endl;

            EntriesDen = DenDimuPt->GetEntries();
            EntriesNeu = NeuDimuPt->GetEntries();
            Acc = EntriesNeu /EntriesDen;
            std::cout << "*** Rebin Den : " << EntriesDen << ", Neu : " << EntriesNeu << ", *** Acc : " << Acc << std::endl;

            DenDimuPt->Sumw2();
            NeuDimuPt->Sumw2();

            hOutputFile->Write();
            hOutputFile->Close();

        }

    //define this as a plug-in
    DEFINE_FWK_MODULE(McAcceptAnalyzer);
