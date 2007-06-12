#include "Correlation/DijetCorrelation/interface/DijetCorrelation.hh"

using namespace LanlJetCorrelation;
using namespace std;
using namespace edm;

DijetCorrelation::DijetCorrelation(const edm::ParameterSet& pset):
  cardname_(pset.getUntrackedParameter<string>("cardName",string("runningCard.txt"))),
  ofilename_(pset.getUntrackedParameter<string>("oFileName",string("outTest.root"))),
  trackTags_(pset.getUntrackedParameter<edm::InputTag>("tracks"))
{
   //now do what ever initialization is needed
  cout<<"***ENTERING THE DRIVER CONSTRUCTOR***"<<endl;
  cout<<"Input card name is "<<cardname_<<endl;
  cout<<"Output file name is "<<ofilename_<<endl;
 
  outFileName = (TString)ofilename_;
  selector = new Selector((TString)cardname_);
  outputer = new OutputWriter(selector);
  mixer = new Mixer(selector,outputer,2000);

  AssocList = *(new std::deque<GenericParticle*> ());
  TriggList = *(new std::deque<GenericParticle*> ());
}


DijetCorrelation::~DijetCorrelation()
{
  delete selector;
  delete outputer;
  delete mixer;
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called once each job just before starting event loop  ------------
void 
DijetCorrelation::beginJob(const edm::EventSetup&)
{
  selector->Print();

  outFile = new TFile(outFileName,"RECREATE");
  outputer->CreateGeneric();
  outputer->CreateAzimuth();
  if(selector->UsePizeroes() || selector->UseDielectrons() || selector->UseDimuons())
    outputer->CreateParentRecon();

  std::cout<<std::endl<<"This production uses the following particles:"<<std::endl
	   <<" hadrons:"<<selector->UseHadrons()<<" photons:"<<selector->UsePhotons()
	   <<" pizeroes:"<<selector->UsePizeroes()<<" electrons:"<<selector->UseElectrons()
	   <<" dielectrons:"<<selector->UseDielectrons()<<" muons:"<<selector->UseMuons()
	   <<" dimuons:"<<selector->UseDimuons()<<std::endl;
  
  nTotalEvts = 0;

}


// ------------ method called to for each event  ------------
void
DijetCorrelation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using reco::TrackCollection;
   nTotalEvts++;

   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(trackTags_, tracks);

   float cent = 0.1;
   float zvtx = 0.1;
   outputer->FillGlobal(cent,zvtx);

   int centralityBin = selector->GetBin(centr,cent);
  
   if(selector->MixDimuonHadrons()){
     cout<<"DijetCorrelation::analyze - GETTING THE MUON LIST"<<endl;
     Handle<reco::TrackCollection> muons;
     iEvent.getByLabel("globalMuons",muons);
     reco::TrackCollection mu = *(muons.product());

     if(mu.size()>1){
       cout<<"DijetCorrelation::analyze - INSIDE DIMUON-HADRONS"<<endl;
       TriggList = (std::deque<GenericParticle*>)selector->GetTriggList(muons);
       outputer->FillSingleHistos(TriggList, centralityBin, trigg);
       AssocList = (std::deque<GenericParticle*>)selector->GetAssocList(tracks);
       outputer->FillSingleHistos(AssocList, centralityBin, assoc);
       cout<<"DijetCorrelation::analyze - NumberDimuons="<<TriggList.size()<<" NumberHadrons="<<AssocList.size()<<endl;

       outputer->FillParentNtuple(real, centralityBin, TriggList); // fills the dimuon ntuple
       CrossCorrelate(TriggList, AssocList, centralityBin);
       cout<<"!!!!! Passed correlation step safe!!!!!"<<endl;

       mixer->FillPool(TriggList,"trigg",nTotalEvts,cent,zvtx); 
       mixer->FillPool(AssocList,"assoc",nTotalEvts,cent,zvtx);
       mixer->RollingMix();
     }
   }

   if(selector->MixTwoHadrons()){
     TriggList = (std::deque<GenericParticle*>)selector->GetTriggList(tracks);
     AssocList = (std::deque<GenericParticle*>)selector->GetAssocList(tracks);
   }

   TriggList.clear(); AssocList.clear(); // clear the lists
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DijetCorrelation::endJob() 
{
  // Do the mixing:
  if(selector->UsePizeroes() || selector->UseDielectrons() || selector->UseDimuons()){
    mixer->MixChildren();
    mixer->CleanPool("child");
  }
  mixer->PrintStats();
  // mixer->MixAll(); // we do rolling buffer mixing in analyze(), so this should not be used
  // Print out stats:
  mixer->CleanPool("trigg");
  mixer->CleanPool("assoc");
  outputer->PrintStats();
  // Write output file:
  outFile->Write();
  outFile->Close();
}

void DijetCorrelation::CrossCorrelate(std::deque<GenericParticle*> trigglist, std::deque<GenericParticle*> assoclist, int cBin){
  if(trigglist.size()<1 || assoclist.size()<1) return;

  for(unsigned int i=0; i<trigglist.size(); i++){
    GenericParticle *trigger = (GenericParticle*)trigglist[i];
    for(unsigned int j=0; j<assoclist.size(); j++){
      GenericParticle *associated = (GenericParticle*)assoclist[j];
      outputer->FillAzimuth(real, cBin, trigger, associated); // trigger first!
    }
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(DijetCorrelation);
