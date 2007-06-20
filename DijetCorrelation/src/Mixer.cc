#include "Correlation/DijetCorrelation/interface/Mixer.hh"

using namespace LanlJetCorrelation;
using namespace std;

Mixer::Mixer(Selector* selectorP, OutputWriter *hin, int evtPoolDepth){
  
  selector = selectorP;
  outputer = hin;
  MAXNOEVENT = evtPoolDepth;
  std::cout<<"Mixing pool depth = "<< MAXNOEVENT <<std::endl;  
  ncbins = selector->GetNoOfBins(centr);
  
  ChildPool.List = *(new Vector3D(ncbins));
  ChildPool.eventNumber = *(new std::deque<std::deque<int> >(ncbins));
  ChildPool.eventCentrality = *(new std::deque<std::deque<float> >(ncbins));
  ChildPool.eventVertex = *(new std::deque<std::deque<float> >(ncbins));
  TriggPool.List = *(new Vector3D(ncbins));
  TriggPool.eventNumber = *(new std::deque<std::deque<int> >(ncbins));
  TriggPool.eventCentrality = *(new std::deque<std::deque<float> >(ncbins));
  TriggPool.eventVertex = *(new std::deque<std::deque<float> >(ncbins));
  AssocPool.List = *(new Vector3D(ncbins));
  AssocPool.eventNumber = *(new std::deque<std::deque<int> >(ncbins));
  AssocPool.eventCentrality = *(new std::deque<std::deque<float> >(ncbins));
  AssocPool.eventVertex = *(new std::deque<std::deque<float> >(ncbins));



  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << " OutputWriter requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }

  CLHEP::HepRandomEngine& engine = rng->getEngine();
  fRandomGenerator = new CLHEP::RandFlat(engine) ;


}

void Mixer::CopyAssocPool(){  
  AssocPoolCopy.List = AssocPool.List;
  AssocPoolCopy.eventNumber = AssocPool.eventNumber;
  AssocPoolCopy.eventCentrality = AssocPool.eventCentrality;
  AssocPoolCopy.eventVertex = AssocPool.eventVertex;
}

void Mixer::FillPool(std::deque<GenericParticle*> partList, TString poolType,
		     int evt, float cent, float zvtx){
  poolType.ToUpper();
  int cBin = selector->GetBin(centr,cent);
  if(partList.size()>0){   // fill pool only with non-empty events
    if(poolType.CompareTo("CHILD")==0){
      ChildPool.List[cBin].push_back(partList);
      ChildPool.eventNumber[cBin].push_back(evt);
      ChildPool.eventCentrality[cBin].push_back(cent);
      ChildPool.eventVertex[cBin].push_back(zvtx);
    } else if(poolType.CompareTo("TRIGG")==0){
      TriggPool.List[cBin].push_back(partList);
      TriggPool.eventNumber[cBin].push_back(evt);
      TriggPool.eventCentrality[cBin].push_back(cent);
      TriggPool.eventVertex[cBin].push_back(zvtx);
    } else if(poolType.CompareTo("ASSOC")==0){
      AssocPool.List[cBin].push_back(partList);
      AssocPool.eventNumber[cBin].push_back(evt);
      AssocPool.eventCentrality[cBin].push_back(cent);
      AssocPool.eventVertex[cBin].push_back(zvtx);
    } else {
      std::cout<<"Mixer::FillPool - ERROR: unknown pool type!"<<std::endl; exit(-1);
    }
  }
}

void Mixer::CleanPool(TString poolType){
  poolType.ToUpper();
  for(int cBin=0; cBin<ncbins; cBin++){
    if(poolType.CompareTo("CHILD")==0){
      ChildPool.List[cBin].clear();
      ChildPool.eventNumber[cBin].clear();
      ChildPool.eventCentrality[cBin].clear();
      ChildPool.eventVertex[cBin].clear();
    } else if(poolType.CompareTo("TRIGG")==0){
      TriggPool.List[cBin].clear();
      TriggPool.eventNumber[cBin].clear();
      TriggPool.eventCentrality[cBin].clear();
      TriggPool.eventVertex[cBin].clear();
    } else if(poolType.CompareTo("ASSOC")==0){
      AssocPool.List[cBin].clear();
      AssocPool.eventNumber[cBin].clear();
      AssocPool.eventCentrality[cBin].clear();
      AssocPool.eventVertex[cBin].clear();
    } else if(poolType.CompareTo("COPY")==0){
      AssocPoolCopy.List[cBin].clear();
      AssocPoolCopy.eventNumber[cBin].clear();
      AssocPoolCopy.eventCentrality[cBin].clear();
      AssocPoolCopy.eventVertex[cBin].clear();
    } else {
      std::cout<<"Mixer::CleanPool - ERROR: unknown pool type!"<<std::endl; exit(-1);
    }
  }
}

void Mixer::PopPool(TString poolType, int cBin){
  poolType.ToUpper();
  if(poolType.CompareTo("CHILD")==0){
    ChildPool.List[cBin].pop_front();
    ChildPool.eventNumber[cBin].pop_front();
    ChildPool.eventCentrality[cBin].pop_front();
    ChildPool.eventVertex[cBin].pop_front();
  } else if(poolType.CompareTo("TRIGG")==0){
    TriggPool.List[cBin].pop_front();
    TriggPool.eventNumber[cBin].pop_front();
    TriggPool.eventCentrality[cBin].pop_front();
    TriggPool.eventVertex[cBin].pop_front();
  } else if(poolType.CompareTo("ASSOC")==0){
    AssocPool.List[cBin].pop_front();
    AssocPool.eventNumber[cBin].pop_front();
    AssocPool.eventCentrality[cBin].pop_front();
    AssocPool.eventVertex[cBin].pop_front();
  } else {
    std::cout<<"Mixer::PopPool - ERROR: unknown pool type!"<<std::endl; exit(-1);
  }
}

// method for two-photon and two-electrons mixed spectra:
void Mixer::MixChildren(){
  std::cout<<"Start mixing children..."<<std::endl;

  std::deque<GenericParticle*> parentlist;
  for(int cBin=0; cBin<ncbins; cBin++){
    int nChildEvnt = (int)ChildPool.List[cBin].size();

    for(int iFirstEvnt=0; iFirstEvnt<nChildEvnt-1; iFirstEvnt++){
      int FirstEvnt = ChildPool.eventNumber[cBin][iFirstEvnt];
      float FirstCent = ChildPool.eventCentrality[cBin][iFirstEvnt];
      float FirstVert = ChildPool.eventVertex[cBin][iFirstEvnt];
      
      for(int iSecondEvnt=iFirstEvnt+1; iSecondEvnt<nChildEvnt; iSecondEvnt++){
	int SecondEvnt = ChildPool.eventNumber[cBin][iSecondEvnt];
	float SecondCent = ChildPool.eventCentrality[cBin][iSecondEvnt];
	float SecondVert = ChildPool.eventVertex[cBin][iSecondEvnt];

	if(selector->SimilarCentrality(FirstCent,SecondCent) &&
	   selector->SimilarVertZ(FirstVert,SecondVert) && FirstEvnt!=SecondEvnt){
	  
	  for(unsigned int ii=0; ii<ChildPool.List[cBin][iFirstEvnt].size(); ii++){
	    GenericParticle *child1 = (GenericParticle*)ChildPool.List[cBin][iFirstEvnt][ii];
       	    int jj = int(fRandomGenerator->fire()*ChildPool.List[cBin][iSecondEvnt].size());

	    GenericParticle *child2 = (GenericParticle*)ChildPool.List[cBin][iSecondEvnt][jj];
	    	    
	    RecoParent *parent = new RecoParent;
	    parent->SetLorentzVector(child1, child2);
	    if(selector->InParentMassAndPtBin(parent)) parentlist.push_back(parent);
	  } // loop over children

	  outputer->FillParentNtuple(mixed, cBin, parentlist);

	} // if the event pair is good for mixing
      } // second loop over events in the pool
    } // first loop over events in the pool
  } // loop over centrality bins
}

// mixing method to be called once:
// loop over the events in the trigger pool and, for each trigger event:
// (1) make a copy of the associated pool
// (2) pick a random associated event
// (3) if the trigger and associated events are within same centrality and vertex interval
// (4) for each trigger particle pick randomly one associated particle and mix them
// (5) throw out the associated event from the copy pool
// (6) repeat steps 2-5 until either the copy pool gets empty or the maximum number of events is reached
// (7) clean the copy of the associated pool and go to next trigger event
void Mixer::Mix(){   
  std::cout<<"Start the mixing..."<<std::endl;
  for(int cBin=0; cBin<ncbins; cBin++){
    for(unsigned int iTriggEvnt=0; iTriggEvnt<TriggPool.List[cBin].size(); iTriggEvnt++){
      int nTriggers = (int)TriggPool.List[cBin][iTriggEvnt].size();
      int TriggEvnt = TriggPool.eventNumber[cBin][iTriggEvnt];
      float TriggCent = TriggPool.eventCentrality[cBin][iTriggEvnt];
      float TriggVert = TriggPool.eventVertex[cBin][iTriggEvnt];

      CopyAssocPool(); // make a copy of the associated pool
      int NumUsedEvents = 0;
      while(!AssocPool.List[cBin].empty() && NumUsedEvents<MAXNOEVENT){

	int iAssocEvnt = int(fRandomGenerator->fire()*(AssocPool.List[cBin].size())); // pick a random associated event

	int nAssociates = (int)AssocPool.List[cBin][iAssocEvnt].size();
	int AssocEvnt = AssocPool.eventNumber[cBin][iAssocEvnt];
	float AssocCent = AssocPool.eventCentrality[cBin][iAssocEvnt];
	float AssocVert = AssocPool.eventVertex[cBin][iAssocEvnt];
      
	if(selector->SimilarCentrality(TriggCent,AssocCent) &&
	   selector->SimilarVertZ(TriggVert,AssocVert) && TriggEvnt!=AssocEvnt){
	  
	  for(int ii=0; ii<nTriggers; ii++){
	    GenericParticle *trigger = (GenericParticle*)TriggPool.List[cBin][iTriggEvnt][ii];
	    int jj = int(fRandomGenerator->fire()*nAssociates); // pick ONE random associated particle
	    GenericParticle *associated = (GenericParticle*)AssocPool.List[cBin][iAssocEvnt][jj];
	    outputer->FillAzimuth(mixed, cBin, trigger, associated); // MIX'EM, trigger first!
	  } // loop over trigger particles
	} // if the event pair Trigger-Associated is good for mixing
	NumUsedEvents++;
	ThrowAssocEvnt(cBin,iAssocEvnt); // throw out used associated event
      } // loop over events in the associated pool
      CleanPool("copy");
    } // loop over events in the trigger pool
  } // loop over centrality bins
}

void Mixer::ThrowAssocEvnt(int cBin, int iAssocEvnt){
 
  int ListBefore = AssocPool.List[cBin].size();
  std::deque<Vector1D>::iterator listiter = AssocPool.List[cBin].begin();
  listiter += iAssocEvnt;
  AssocPool.List[cBin].erase(listiter);
  int ListAfter = AssocPool.List[cBin].size();
  
  int EvntBefore = AssocPool.eventNumber[cBin].size();
  std::deque<int>::iterator eventiter = AssocPool.eventNumber[cBin].begin();
  eventiter += iAssocEvnt;
  AssocPool.eventNumber[cBin].erase(eventiter);
  int EvntAfter = AssocPool.eventNumber[cBin].size();
  
  int CentBefore = AssocPool.eventCentrality[cBin].size();
  std::deque<float>::iterator centiter = AssocPool.eventCentrality[cBin].begin();
  centiter += iAssocEvnt;
  AssocPool.eventCentrality[cBin].erase(centiter);
  int CentAfter = AssocPool.eventCentrality[cBin].size();
  
  int VertBefore = AssocPool.eventVertex[cBin].size();
  std::deque<float>::iterator vertiter = AssocPool.eventVertex[cBin].begin();
  vertiter += iAssocEvnt;
  AssocPool.eventVertex[cBin].erase(vertiter);
  int VertAfter = AssocPool.eventVertex[cBin].size();

  // TESTS:
  if(abs(ListBefore-ListAfter-1)>0 || abs(EvntBefore-EvntAfter-1)>0 ||
     abs(CentBefore-CentAfter-1)>0 || abs(VertBefore-VertAfter-1)>0){
    std::cout<<"Mixer::ThrowAssocEvnt - ERROR: event was not thrown from pool!"<<std::endl;
    exit(-1);
  }
  if(ListAfter!=EvntAfter || ListAfter!=CentAfter || ListAfter!=VertAfter){
    std::cout<<"Mixer::ThrowAssocEvnt - ERROR: event pool mismatch!"<<std::endl;
    exit(-1);
  }
}

// mixing method to be called once
// loop over the events in the trigger pool and, for each trigger event:
// (1) loop over events in the associated pool
// (2) if the pair of trigger and associated events are within same centrality and vertex interval 
// (3) for each trigger particle pick randomly one associated particle and mix them
void Mixer::MixAll(){   
  std::cout<<"Start mixing all..."<<std::endl;

  for(int cBin=0; cBin<ncbins; cBin++){
    int nTriggEvnt = (int)TriggPool.List[cBin].size();
    int nAssocEvnt = (int)AssocPool.List[cBin].size();

    for(int iTriggEvnt=0; iTriggEvnt<nTriggEvnt; iTriggEvnt++){
      int nTriggers = (int)TriggPool.List[cBin][iTriggEvnt].size();
      int TriggEvnt = TriggPool.eventNumber[cBin][iTriggEvnt];
      float TriggCent = TriggPool.eventCentrality[cBin][iTriggEvnt];
      float TriggVert = TriggPool.eventVertex[cBin][iTriggEvnt];
      
      for(int iAssocEvnt=0; iAssocEvnt<nAssocEvnt; iAssocEvnt++){
	int nAssociates = (int)AssocPool.List[cBin][iAssocEvnt].size();
	int AssocEvnt = AssocPool.eventNumber[cBin][iAssocEvnt];
	float AssocCent = AssocPool.eventCentrality[cBin][iAssocEvnt];
	float AssocVert = AssocPool.eventVertex[cBin][iAssocEvnt];

	if(selector->SimilarCentrality(TriggCent,AssocCent) &&
	   selector->SimilarVertZ(TriggVert,AssocVert) && TriggEvnt!=AssocEvnt){
	  
	  for(int ii=0; ii<nTriggers; ii++){
	    GenericParticle *trigger = (GenericParticle*)TriggPool.List[cBin][iTriggEvnt][ii];
	    int jj = int(fRandomGenerator->fire()*nAssociates); // pick ONE random associated particle
	    GenericParticle *associated = (GenericParticle*)AssocPool.List[cBin][iAssocEvnt][jj];
	    outputer->FillAzimuth(mixed, cBin, trigger, associated); // MIX'EM, trigger first!
	  } // loop over trigger particles
	} // if the event pair Trigger-Associated is good for mixing
      } // loop over events in the associated pool
    } // loop over events in the trigger pool
  } // loop over centrality bins
}

// mixing method to be called every event
// (Push One Event to the Bottom- Mix Top Event With Pool - Pop Top Event Out)
void Mixer::RollingMix(){   
  
  for(int cBin=0;cBin<ncbins; cBin++){
    if((int)AssocPool.List[cBin].size()>=MAXNOEVENT){
      while((int)TriggPool.List[cBin].size()>0){
	int nTriggers = (int)TriggPool.List[cBin][0].size();
	float TriggCent = TriggPool.eventCentrality[cBin][0];
	float TriggVert = TriggPool.eventVertex[cBin][0];
	int TriggEvnt = TriggPool.eventNumber[cBin][0];
	
	for(int iAssocEvnt=0; iAssocEvnt<(int)AssocPool.List[cBin].size(); iAssocEvnt++){
	  int nAssociates = (int)AssocPool.List[cBin][iAssocEvnt].size();
	  float AssocCent = AssocPool.eventCentrality[cBin][iAssocEvnt];
	  float AssocVert = AssocPool.eventVertex[cBin][iAssocEvnt];
	  int AssocEvnt = AssocPool.eventNumber[cBin][iAssocEvnt];
	    
	  if(selector->SimilarCentrality(TriggCent,AssocCent) &&
	     selector->SimilarVertZ(TriggVert,AssocVert) && TriggEvnt!=AssocEvnt){
	    
	    for(int ii=0;ii<nTriggers;ii++){
	      GenericParticle *trigger = (GenericParticle*)TriggPool.List[cBin][0][ii];
	      int jj = int(fRandomGenerator->fire()*nAssociates);
	      GenericParticle *associated = (GenericParticle*)AssocPool.List[cBin][iAssocEvnt][jj];
	      outputer->FillAzimuth(mixed, cBin, trigger, associated);
	    } // particle loop
	    
	  } // good event pair
	} // loop over associated pool
	PopPool("trigg",cBin);
      } // if trigger pool is not empty
      PopPool("assoc",cBin);
      if((int)AssocPool.List[cBin].size()>(MAXNOEVENT+1)){
	std::cout<<"Mixer::RollingMix - ERROR: Pool size out of bounds!"<<std::endl; exit(-1);
      }
    } // if the associated pool for one centBin is full
  } // loop over centBins
}

void Mixer::PrintStats(){
  std::cout<<std::endl<<"Pools statistics before starting to mix:"<<std::endl;
  for(int cBin=0; cBin<ncbins; cBin++){
    std::cout<<"Cent bin:"<<cBin
	     <<"   Trigger Pool Size:"<<TriggPool.List[cBin].size()
	     <<"   Associated Pool Size:"<<AssocPool.List[cBin].size()<<std::endl;
  }
}
