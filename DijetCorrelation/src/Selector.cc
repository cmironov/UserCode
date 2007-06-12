#include "Correlation/DijetCorrelation/interface/Selector.hh"

using namespace LanlJetCorrelation;
using namespace std;
using namespace reco;

Selector::Selector(TString filename){
  cout<<"***ENTERING THE SELECTOR CONSTRUCTOR***"<<endl;
  for(int i=0; i<4; i++){
    NoOfBins[i]=-1;
    for(int j=0; j<mdim; j++) Borders[i][j]=-1.;
  }
  list.open(filename,ios::in);
  if(!list) {
    std::cout<<"Selector::Selector - ERROR: file <"<<filename<<" not found!"<<std::endl;
    exit(1);
  }
  fname = filename;
  DefineContent();
  ReadIn();
}

///////////////////////////////////////////////////////////
//
// Methods for reading in the running card:
//
///////////////////////////////////////////////////////////
void Selector::ReadIn(){
  char bufmain[200], buf[200];
  while(!list.eof()){
    list>>bufmain;
    // std::cout<<bufmain<<std::endl;
    for(int iv=0; iv<counter; iv++)
      if(var[iv]->CompareTo(bufmain)==0){
        int index=0;
        bool notDelim = true;
        while(notDelim){
          list>>buf;
          notDelim = (strcmp(buf,";") != 0);
          if(notDelim) {
            *((float*)storage[iv]+index) = atof(buf);
            index++;
          }
        }
        if ((int*)storeDimP[iv] != NULL)
          *(int*)storeDimP[iv] = index;
        else
          storeDim[iv] = index;
        readIn[iv]=1;
      }
  }
  for(int iv=0;iv<counter;iv++) if(readIn[iv]==0) std::cout<<*var[iv]<<" not found "<<std::endl;
}

void Selector::AddLine(char *name, unsigned int addr, int *dim = NULL){
  var[counter] = new TString(name);
  readIn[counter] = 0;
  storage[counter] = addr;
  storeDimP[counter++] = dim;
  if(counter>=MAXDIM) {
    std::cout<<"Error: too many input lines. Enlarge MAXDIM"<<std::endl;
    exit(0);
  } 
}

void Selector::Print(){
  std::cout<<std::endl<<"====== RUNNING CARD: "<<fname<<" ========="<<std::endl;
  for(int iv=0;iv<counter;iv++) {
    std::cout<< *var[iv] << "\t" ;
    if(strlen(*var[iv])<8) std::cout<<"\t";
    int dim;
    if((int*)storeDimP[iv]==NULL) dim = storeDim[iv]; 
    else dim = *(int*)storeDimP[iv];
    std::cout << dim << "\t";
    for(int id=0;id<dim;id++) 
      std::cout<<*((float*)storage[iv]+id)<<" "; 
    std::cout<<std::endl;
  }
}

void Selector::DefineContent(){
  counter=0;
  // Generic selections:
  AddLine("CorrelType",    (unsigned int)& CorrelType );
  AddLine("CentBinBorders",(unsigned int)& Borders[centr],  & NoOfBins[centr]);
  AddLine("TriggPtBorders",(unsigned int)& Borders[trigg],  & NoOfBins[trigg]);
  AddLine("AssocPtBorders",(unsigned int)& Borders[assoc],  & NoOfBins[assoc]);
  AddLine("maxZVert",(unsigned int)& maxVertZRange );
  AddLine("maxMixDZ",(unsigned int)& maxVertDZ );
  AddLine("maxDCent",(unsigned int)& maxDCent );

  // Hadron selections:
  AddLine("MinHadronPt",(unsigned int)& MinHadronPt );

  // Muon selections:
  AddLine("MinMuonPt",(unsigned int)& MinMuonPt );

}


///////////////////////////////////////////////////////////
//
// Methods for building the deques of good particles:
//
///////////////////////////////////////////////////////////

std::deque<GenericParticle*> Selector::GetAssocList(edm::Handle<reco::TrackCollection> &hadlist){
  cout<<"***ENTERINGING Selector::GetAssocList***"<<endl;
  std::deque<GenericParticle*> AssocHadList;

  reco::TrackCollection tCa = *(hadlist.product());
  int NTracks = tCa.size();
  if(NTracks<1) return AssocHadList;
  cout<<"Selector::GetAssocList - NTracks = "<<NTracks<<endl;

  for(int i=0; i<NTracks; i++){
    Hadron *had = new Hadron;
    
    float pT = tCa[i].pt();
    had->SetPt(pT);
    had->SetTheta(tCa[i].theta());
    had->SetPhi(tCa[i].phi());
    had->SetCharge(tCa[i].charge());
    had->SetMass(-999.);
    
    AssocHadList.push_back(had);    
  } 
  
  return AssocHadList;
}

std::deque<GenericParticle*> Selector::GetTriggList(edm::Handle<reco::TrackCollection> &inputlist){ 
  using LanlJetCorrelation::Muon;
  std::deque<GenericParticle*> TriggList;

  reco::TrackCollection tC = *(inputlist.product());
  unsigned int NParticles = tC.size();
  if(NParticles<1) return TriggList;
  cout<<"Selector::GetTriggList - NParticles = "<<NParticles<<endl;

  if(UseDimuons()){
    std::deque<GenericParticle*> GoodMuonList;
    Muon *muon = new Muon; 
    for(unsigned int i=0; i<NParticles; i++){
      if(!GoodMuon()) continue; // define good muons
      //muon->SetKinematics(tC[i].pt(),tC[i].phi(),tC[i].theta(),tC[i].charge());
      muon->SetPt(tC[i].pt());       muon->SetPhi(tC[i].phi());
      muon->SetTheta(tC[i].theta()); muon->SetCharge(tC[i].charge());
      cout<<"Selector::GetTriggList - Muon("<<i<<")="<<muon->GetPt()<<" "<<muon->GetPhi()<<endl;
      GoodMuonList.push_back(muon);
      cout<<"Selector::GetTriggList - GoodMuons = "<<GoodMuonList.size()<<endl;
    }
    TriggList = DecayRecon(GoodMuonList);
  }

  return TriggList;
}


std::deque<GenericParticle*> Selector::DecayRecon(std::deque<GenericParticle*> childlist){

  std::deque<GenericParticle*> parentlist;
  if(childlist.size()<1) return parentlist;

  for(unsigned int i=0; i<childlist.size()-1; i++){
    GenericParticle *child1 = (GenericParticle*)childlist[i];
    for(unsigned int j=i+1; j<childlist.size(); j++){
      GenericParticle *child2 = (GenericParticle*)childlist[j];
      RecoParent *parent = new RecoParent();
      parent->SetLorentzVector(child1, child2);
      //if(InParentMassAndPtBin(parent)) 
	parentlist.push_back(parent);
	cout<<"Selector::DecayRecon i="<<i<<" j="<<j
	    <<" parent="<<parent->GetMass()<<" "<<parent->GetPt()<<endl;
    }
  }
  
  return parentlist;
}

bool Selector::GoodMuon(){
  bool GoodMuonCuts = false;
    
  if(true) GoodMuonCuts = true; // REPLACE THE "true" INSIDE THE CONDITIONAL WITH REAL CUTS ON MUON

  return GoodMuonCuts;
}

bool Selector::GoodParent(RecoParent* parent, GenericParticle *p1, GenericParticle *p2){
  bool PassPairCuts = false;
    
  if(true) PassPairCuts = true; // REPLACE THE "true" INSIDE THE CONDITIONAL WITH REAL CUTS ON PARENT

  return PassPairCuts;
}

bool Selector::TrackPairCuts(GenericParticle *Trigger, GenericParticle *Associated){
  bool IsGoodPair = false;

  if(true) IsGoodPair = true; // REPLACE THE DUMMY true INSIDE THE CONDITIONAL WITH A REAL PAIR CUT

  return IsGoodPair;
}

//////////////////////////////////////////////////////////
//
// Binning Methods
//
/////////////////////////////////////////////////////////

float Selector::GetBinBorder( corrType ctype, int ii){
  if(ii>=0 && ii<=NoOfBins[ctype]) return Borders[ctype][ii]; 
  else { std::cout<<"Error: bin "<<ii<<" out!"<<std::endl;  exit(0); }
}

int Selector::GetBin(corrType ctype, float val){
  int iBin=-1;
  for(int i=0;i<NoOfBins[ctype];i++)
    if(Borders[ctype][i]<=val && val<Borders[ctype][i+1]) iBin=i;
  return iBin;
}
  
bool Selector::InParentMassAndPtBin(RecoParent *parent){ // DO NOT HARDCODE CUTS! PASS THEM THRU THE CARD! 
  partType type = parent->GetPID();
  float mass = parent->GetMass();
  bool InMassBin = false;

  // YOU CAN GET THE FIRST TRIGGER PT BIN BORDER WITH GetBinBorder(trigg, 0)
  
  if(type==diphoton && 0.127<=mass && mass<=0.157) InMassBin = true;
  if(type==dimuon && 81.<=mass && mass<=101.) InMassBin = true;
  if(type!=diphoton && type!=dimuon)
    std::cout<<"Selector::GetParentMassBin - WARNING: unknown particle type: "<<type<<" "<<mass<<std::endl;
  
  return InMassBin;
}
