#ifndef MUON_H
#define MUON_H

#include "Correlation/DijetCorrelation/interface/GenericParticle.hh"

namespace LanlJetCorrelation{

  class Muon : public GenericParticle {

  public:
    Muon(){;}			//constructor
    ~Muon(){;}		//destructor
    
  private:
    
    float  Pt;                       // pt        
    float  Theta;                    // theta at vertex
    float  Phi;                      // phi at vertex 
    int    Charge;           
    
  public:
    
    partType GetPID() const {return muon;}
    float  GetPt() const {return Pt;}
    float  GetTheta() const {return Theta;}
    float  GetPhi() const {return Phi;}
    float  GetMass() const {return MuonMass;}
    float  GetAssym() const {return -999.;}
    float  GetOpenAng() const {return -999.;}
    int    GetCharge() const {return Charge;}
    
    float GetMom()  const {return Pt/sin(Theta);}
    float GetPx()   const {return Pt*cos(Phi);}
    float GetPy()   const {return Pt*sin(Phi);}
    float GetPz()   const {return Pt/tan(Theta);}
    float GetEner() const {return sqrt(GetMom()*GetMom()+GetMass()*GetMass());}
    
    void     SetPt(float pt) {Pt=pt;}
    void     SetTheta(float th) {Theta=th;}
    void     SetPhi(float ph) {Phi=ph;}
    void     SetCharge(int q) {Charge=q;}

  };
}// namespace declaration
#endif
