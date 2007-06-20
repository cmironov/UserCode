#ifndef GENERICPARTICLE_HH
#define GENERICPARTICLE_HH

#include "Correlation/DijetCorrelation/interface/Constants.hh"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TString.h>
#include <TRandom.h>

// GenericParticle is at the top of preprocessor includes,
// hence we add the general libraries here: 
#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <memory>

namespace LanlJetCorrelation {

   class GenericParticle {
      
   public:
      virtual ~GenericParticle(){;}		//destructor
      
      virtual partType GetPID() const {return unknown;}
      virtual float  GetPt() const {return -999.;}
      virtual float  GetTheta() const {return -999.;}
      virtual float  GetPhi() const {return -999.;}
      virtual float  GetMass() const {return -999.;}
      virtual float  GetAssym() const {return -999.;}
      virtual float  GetOpenAng() const {return -999.;}
      virtual int    GetCharge() const {return -999;}
      virtual float  GetEner() const {return -999.;}
      virtual float  GetPx() const {return -999.;}
      virtual float  GetPy() const {return -999.;}
      virtual float  GetPz() const {return -999.;}
      
      virtual void SetPt(float pt) {;}
      virtual void SetTheta(float th) {;}
      virtual void SetPhi(float ph) {;}
      virtual void SetMass(float m) {;}
      virtual void SetAssym(float a) {;}
      virtual void SetOpenAng(float o) {;}
      virtual void SetCharge(int ch) {;}
      
   };
   
} // namespace declaration

#endif

