#ifndef DIMUONMC_H
#define DIMUONMC_H


#include "Correlation/DijetCorrelation/interface/RecoParent.hh"

namespace LanlJetCorrelation{

  class DiMuonMC : public RecoParent {

  public:
    DiMuonMC(){;}			//constructor
    ~DiMuonMC(){;}		//destructor

  private:

    int PunchThru;        // =0 if no punch-thrus, =1 if one is punch-thru, =2 if both are punch-thrus
    partType FirstType;   // type of first muon; see Constants.hh
    partType SecondType;  // type of second muon; see Constants.hh

  public:

    int  IsPunchThru() const {return PunchThru;}
    void  SetPunchThru(int pt) {PunchThru=pt;}
    partType  FirstMuonType() const {return FirstType;}
    void  SetFirstType(partType pt) {FirstType=pt;}
    partType  SecondMuonType() const {return SecondType;}
    void  SetSecondType(partType pt) {SecondType=pt;}
    
  };
}//namespace
#endif
