#include "Correlation/DijetCorrelation/interface/RecoParent.hh"

using namespace LanlJetCorrelation;
using namespace std;

void RecoParent::SetLorentzVector(GenericParticle *p1, GenericParticle *p2){
  TLorentzVector Part1, Part2;

  Part1.SetPxPyPzE(p1->GetPx(), p1->GetPy(), p1->GetPz(), p1->GetEner());
  Part2.SetPxPyPzE(p2->GetPx(), p2->GetPy(), p2->GetPz(), p2->GetEner());
 
  PartSum = Part1 + Part2;
  Assym = fabs((p1->GetEner()-p2->GetEner())/(p1->GetEner()+p2->GetEner()));
  OpenAng = Part1.Angle(Part2.Vect());
  if(p1->GetPID()==photon && p2->GetPID()==photon) {
    PID = diphoton;
    Charge = 0;
  } else if(p1->GetPID()==electron && p2->GetPID()==electron) {
    PID = dielectron;
    if(p1->GetCharge()>0&&p2->GetCharge()>0) Charge = 1;
    else if(p1->GetCharge()<0&&p2->GetCharge()<0) Charge = 0;
    else Charge = p1->GetCharge()*p2->GetCharge();
  } else if(p1->GetPID()==muon && p2->GetPID()==muon) {
    PID = dimuon;
    if(p1->GetCharge()>0&&p2->GetCharge()>0) Charge = 1;
    else if(p1->GetCharge()<0&&p2->GetCharge()<0) Charge = 0;
    else Charge = p1->GetCharge()*p2->GetCharge();
  } else {cout<<"RecoParent::SetLorentzVector - WARNING: unknown particle!"<<endl;}
}

