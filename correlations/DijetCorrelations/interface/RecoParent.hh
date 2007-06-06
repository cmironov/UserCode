#ifndef RECOPARENT_H
#define RECOPARENT_H

#include "Correlation/DijetCorrelation/interface/GenericParticle.hh"

namespace LanlJetCorrelation{
   
  class RecoParent : public GenericParticle {
      
  public:
    RecoParent(){PartSum.SetPxPyPzE(0.,0.,0.,0.);}   //constructor 
    ~RecoParent(){;}		      //destructor
    
    partType GetPID()  const {return PID;}
    float GetPt()      const {return PartSum.Pt();}
    float GetPhi()     const {return PartSum.Phi();}
    float GetTheta()   const {return PartSum.Theta();}
    float GetMass()    const {return PartSum.M();}
    float GetAssym()   const {return Assym;}
    float GetOpenAng() const {return OpenAng;}
    int   GetCharge()  const {return Charge;}
      
    TLorentzVector GetLorentzVector() {return PartSum;}
    void SetLorentzVector(GenericParticle*, GenericParticle*);
    
  private:
    TLorentzVector PartSum;
    float Assym;
    float OpenAng;
    int Charge;
    partType PID;
      
  };
  
} // namespace declaration

#endif

