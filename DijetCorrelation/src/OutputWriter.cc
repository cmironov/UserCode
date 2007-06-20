#include "Correlation/DijetCorrelation/interface/OutputWriter.hh"

using namespace LanlJetCorrelation;
using namespace std;

OutputWriter::OutputWriter(Selector* selectorp){
  
    selector = selectorp;
    std::cout<<"Binning used in this production:"<<std::endl
	     <<" Centrality:"<<selector->GetNoOfBins(centr)
	     <<" TriggPt:"<<selector->GetNoOfBins(trigg)
	     <<" AssocPt:"<<selector->GetNoOfBins(assoc)<<std::endl;

    if(selector->GetNoOfBins(centr) > MAXCENTRBIN ){
      std::cout<<"ERROR: No of Centrality bins exceed max dim in OutputWriter.cc "<<std::endl;
      exit(0);
    }
    if(selector->GetNoOfBins(trigg) > MAXTRIGGBIN ){
      std::cout<<"ERROR: No of Trigger pT bins exceed max dim in OutputWriter.cc "<<std::endl;
      exit(0);
    }
    if(selector->GetNoOfBins(assoc) > MAXASSOCBIN ){
      std::cout<<"ERROR: No of Associated pT bins exceed max dim in OutputWriter.cc "<<std::endl;
      exit(0);
    }    

    for(int i=0; i<selector->GetNoOfBins(centr); i++) {nReal[i]=0; nMix[i]=0; nTrigg[i]=0; nAssoc[i]=0;}

    // set random number generation
    edm::Service<edm::RandomNumberGenerator> rng;
    if ( ! rng.isAvailable()) {
      throw cms::Exception("Configuration")
	<< " OutputWriter requires the RandomNumberGeneratorService\n"
	"which is not present in the configuration file.  You must add the service\n"
	"in the configuration file or remove the modules that require it.";
    }
    //    randEngine = new RandomEngine(&(*rng));
    CLHEP::HepRandomEngine& engine = rng->getEngine();
    fRandomGenerator = new CLHEP::RandFlat(engine) ;


}


void OutputWriter::CreateGeneric(){
  std::cout<<std::endl<<"Creating Generic Histograms"<<std::endl;
  float pt_bw=1./10; // pT bin width
  
  for(int ic=0; ic<selector->GetNoOfBins(centr); ic++){

    for(int ia=0; ia<selector->GetNoOfBins(assoc); ia++){
      float pTb1 = selector->GetBinBorder(assoc, ia);
      float pTb2 = selector->GetBinBorder(assoc, ia+1);
      int NPtBins = int(TMath::Ceil((pTb2-pTb1)/pt_bw));
      hname = "hAssocPt"; hname += ic; hname += ia;
      htit = "cent:"; htit +=  ic; htit += "   assoc:"; htit += ia; 
      hAssocPt[ic][ia] = new TH1D(hname, htit,NPtBins,pTb1,pTb2);
      hname = "hAssocAcceptance"; hname += ic; hname += ia;
      htit = "cent:"; htit +=  ic; htit += "   assoc:"; htit += ia; 
      hAssocAcceptance[ic][ia] = new TH2D(hname, htit,200,-2*pi,2*pi,200,-5.,5.);
    } // associated pT binning loop

    for(int it=0; it<selector->GetNoOfBins(trigg); it++){
      float pTb1 = selector->GetBinBorder(trigg, it);
      float pTb2 = selector->GetBinBorder(trigg, it+1);
      int NPtBins = int(TMath::Ceil((pTb2-pTb1)/pt_bw));
      hname = "hTriggPt"; hname += ic; hname += it;
      htit = "cent:"; htit +=  ic; htit += "   trigg:"; htit += it; 
      hTriggPt[ic][it] = new TH1D(hname, htit,NPtBins,pTb1,pTb2);
      hname = "hTriggAcceptance"; hname += ic; hname += it;
      htit = "cent:"; htit +=  ic; htit += "   trigg:"; htit += it; 
      hTriggAcceptance[ic][it] = new TH2D(hname, htit,200,-2*pi,2*pi,200,-5.,5.);
    } // trigger pT binning loop

  } // centrality binning loop
 
  hBinsCent   = new TH1D("hBinsCent","centrality binning",selector->GetNoOfBins(centr)+1, 0, 1); 
  for(int i=1;i<=selector->GetNoOfBins(centr)+1; i++)
    hBinsCent->SetBinContent(i,selector->GetBinBorder(centr,i-1));
  hBinsTrigg    = new TH1D("hBinsTrigg","trigger pT binning",   selector->GetNoOfBins(trigg)+1, 0, 1); 
  for(int i=1;i<=selector->GetNoOfBins(trigg)+1; i++)
    hBinsTrigg->SetBinContent(i,selector->GetBinBorder(trigg,i-1));
  hBinsAssoc    = new TH1D("hBinsAssoc","associated pT binning",   selector->GetNoOfBins(assoc)+1, 0, 1);
  for(int i=1;i<=selector->GetNoOfBins(assoc)+1; i++)
    hBinsAssoc->SetBinContent(i,selector->GetBinBorder(assoc,i-1));

  hCentr = new TH1D("hCentr","centrality distribution",100, 0., 100.);
  hZVert = new TH1D("hZVert","vertex distribution",100, -50., 50.);
}

void OutputWriter::CreateAzimuth(){
  std::cout<<std::endl<<"Creating Azimuth Histograms"<<std::endl;
  int bins=80; 
  float lr=-0.5*pi, ur=1.5*pi;

  // Since mixing containes detector-induced correlations only (acceptance, efficiency, etc.)
  // one can hunt down "bad" runs by comparing their mixed DeltaPhi distributions: choose a couple of
  // good large runs and pair them with all the other runs and compute the Chi^2 of their hDphiMixQA;
  // "bad" runs (with detector problems) will have large Chi^2 (>3) when compared with these "good" runs
  hDphiMixQA = new TH1D("hDphiMixQA", "Mixed DeltaPhi for QA", bins, lr, ur);
  for (int hic=0; hic<selector->GetNoOfBins(centr); hic++){ // centrality loop
    float cb1 = selector->GetBinBorder(centr, hic);
    float cb2 = selector->GetBinBorder(centr, hic + 1);
    
    for (int hit=0; hit<selector->GetNoOfBins(trigg); hit++){ // trigg pt loop
      float tb1 = selector->GetBinBorder(trigg, hit);
      float tb2 = selector->GetBinBorder(trigg, hit + 1);
      
      for(int ityp=0; ityp<2; ityp++){ // correl type loop (real/mixed)
	
	for(int hia=0; hia<selector->GetNoOfBins(assoc); hia++){ // assoc pt loop      
	  float ab1 = selector->GetBinBorder(assoc, hia);
	  float ab2 = selector->GetBinBorder(assoc, hia + 1);
	  htit = "cent: "; htit += cb1; htit += "-"; htit += cb2;
	  htit += "   trigg:"; htit += tb1; htit += "-"; htit += tb2;
	  htit += "   assoc:"; htit += ab1; htit += "-"; htit += ab2;
	  
	  hname = "hCorr"; hname += ityp; hname += hic; hname += hit; hname += hia;
	  hCorr[ityp][hic][hit][hia] = new TH2D(hname, htit, bins, lr, ur, 100, -5., 5.);
	  
	} // assoc pt loop
      } // correl type loop
    } // trigg pt loop
  } // centrality loop
  
}

void OutputWriter::CreateParentRecon(){
  std::cout<<std::endl<<"Creating Parent reconstruction histograms"<<std::endl;

  for(int ityp=0; ityp<2; ityp++){
    hname = "ParentNtu"; hname += ityp;
    htit  = "Parent Ntuple: "; if(ityp==0) htit += "Real"; else htit += "Mixed"; 
    ParentNtu[ityp] = new TNtuple(hname,htit,"q:m:pT:phi:theta:alpha:openang");
  } // correl type loop
}

// NOW THE FILLER METHODS:

void OutputWriter::FillGlobal(float cent, float zvert){
  hCentr->Fill(cent);
  hZVert->Fill(zvert);
}

void OutputWriter::FillSingleHistos(std::deque<GenericParticle*> PartList, int cBin, corrType type){
  int NumPart = (int)PartList.size();
  if(NumPart<1) return;
  if(type!=trigg && type!=assoc) 
    {cout<<"OutputWriter::FillTriggHistos - ERROR: wrong corrType received!"<<endl<<endl; exit(-1);}
  if(cBin<0)
    {cout<<"OutputWriter::FillTriggHistos - ERROR: negative centrality index!"<<endl<<endl; exit(-1);}

  for(int i=0; i<NumPart; i++){
    GenericParticle *Particle = (GenericParticle*)PartList[i];
    float pt = Particle->GetPt();
    float phi = Particle->GetPhi();
    float eta = log(tan(Particle->GetTheta()/2.0));
    int   pTBin  = selector->GetBin(type, pt);
    cout<<"OutputWriter::FillSingleHistos "<<pt<<" "<<phi<<" "<<pTBin<<endl;
    if(pTBin<0)
      {cout<<"OutputWriter::FillTriggHistos - ERROR: negative pT index!"<<endl<<endl; exit(-1);}
    if(type==trigg){
      hTriggPt[cBin][pTBin]->Fill(pt);
      hTriggAcceptance[cBin][pTBin]->Fill(phi,eta);
    } else if(type==assoc) {
      hAssocPt[cBin][pTBin]->Fill(pt);
      hAssocAcceptance[cBin][pTBin]->Fill(phi,eta);
    }
  }
}

void OutputWriter::FillParentNtuple(fillType type, int cBin, std::deque<GenericParticle*> ParentList){
  int NumParent = (int)ParentList.size();
  if(NumParent<1) return;
  float ParentVar[7]={0.};

  for(int i=0; i<NumParent; i++){
    GenericParticle *parent = (GenericParticle*) ParentList[i];
    ParentVar[0] = parent->GetCharge();
    ParentVar[1] = parent->GetMass();
    ParentVar[2] = parent->GetPt();
    ParentVar[3] = parent->GetPhi();
    ParentVar[4] = parent->GetTheta();
    ParentVar[5] = parent->GetAssym();
    ParentVar[6] = parent->GetOpenAng();
    ParentNtu[type]->Fill(ParentVar);
  }

}

void OutputWriter::FillAzimuth(fillType fTyp, int cBin, GenericParticle *Trigger, GenericParticle *Associated){
  if(cBin<0){cout<<"OutputWriter::FillAzimuth - ERROR: centrality index negative!"<<endl; exit(-1);}

  // With two-hadron correlations (or other identical particle) one gets two-particle detector proximity effects
  // (like cluster merging, ghost tracks, etc.). Special two-particle cuts have to be applied to remove them 
  // from the near-side (DeltaPhi~0). Again, they appear only at small DeltaPhi and for same particle type.
  // Apply charged hadron pair cuts:
  if(Trigger->GetPID()==hadron && Associated->GetPID()==hadron)
    if(selector->TrackPairCuts(Trigger,Associated)) return; // true if any track pair cut fails

  // trigger information (this is why the first particle has to be the trigger):
  int it = selector->GetBin(trigg, Trigger->GetPt());
  float ptt  = Trigger->GetPt();
  float phit = Trigger->GetPhi();
  float thet = Trigger->GetTheta();
  // associated information:
  int ia  = selector->GetBin(assoc, Associated->GetPt());
  float pta  = Associated->GetPt();
  float phia = Associated->GetPhi();
  float thea = Associated->GetTheta();

  if(it<0 || ia<0) return; // only selected pairs
  if(fabs(ptt-pta)<EPS && fabs(phit-phia)<EPS && fabs(thet-thea)<EPS) return; // don't auto-correlate
  
  // Compute basic quantities:
  float dphi = DeltaPhi(phit,phia);  //   dphi is in (-pi/2, 3*pi/2)
  float deta = log(tan(thet/2.))-log(tan(thea/2.));
  // Fill histograms:
  if(fTyp==mixed) hDphiMixQA->Fill(dphi);
  hCorr[fTyp][cBin][it][ia]->Fill(dphi,deta);
  
  if(fTyp==real) nReal[cBin]++; else nMix[cBin]++;
}

void OutputWriter::PrintStats(){
  for(int i=0; i<selector->GetNoOfBins(centr); i++)
    std::cout<<"Cent bin:"<<i<<std::endl
	     <<"   Real Pairs:"<<nReal[i]<<"  Mixed Pairs:"<<nMix[i]
	     <<"  Frac:"<<float(nReal[i])/nMix[i]*100.<<"%"<<std::endl
	     <<"   Triggers:"<<nTrigg[i]<<"  Associated:"<<nAssoc[i]
	     <<"  Frac:"<<float(nTrigg[i])/nAssoc[i]*100.<<"%"<<std::endl;
}

float OutputWriter::DeltaPhi(float  phi1, float phi2){

//random number seed

//  double num=randGun.flatShoot();

  double num=fRandomGenerator->fire();
  float res = phi2-phi1;
  if(num<0.5) res = phi1 - phi2;
  if(1.5*pi<res && res<=2*pi) res-=2*pi;
  else if(-2*pi<=res && res<-0.5*pi) res+=2*pi;
  return res;  
}

