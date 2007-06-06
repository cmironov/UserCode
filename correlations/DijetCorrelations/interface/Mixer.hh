#ifndef MIXER_HH
#define MIXER_HH

#include  "Correlation/DijetCorrelation/interface/Selector.hh"
#include  "Correlation/DijetCorrelation/interface/OutputWriter.hh"

typedef std::deque<LanlJetCorrelation::GenericParticle*> Vector1D;
typedef std::deque<Vector1D> Vector2D;
typedef std::deque<Vector2D> Vector3D;

namespace LanlJetCorrelation{
  class Mixer {
    
  public:
    Mixer(Selector* selectorP, OutputWriter *hin, int evtPoolDepth );
    virtual ~Mixer(){;}
    
    struct EventPool{
      Vector3D List; //Three dimensional array: outer - centrality, middle - event, inner - particle
      std::deque< std::deque<int> > eventNumber;
      std::deque< std::deque<float> > eventCentrality;
      std::deque< std::deque<float> > eventVertex;
    } TriggPool, AssocPool, AssocPoolCopy, ChildPool;
    
    // pool manipulation:
    void CopyAssocPool();
    void FillPool(std::deque<GenericParticle*>, TString, int, float, float);
    void CleanPool(TString);
    void PopPool(TString,int);
    // mixing methods:  
    void MixChildren(); // invariant mass/pT spectra for reconstructed parents
    void Mix();         // 1st mixing method for relative angular (DeltaPhi x DeltaTheta) spectra
    void ThrowAssocEvnt(int, int);  
    void MixAll();      // 2nd mixing method for relative angular (DeltaPhi x DeltaTheta) spectra  
    void RollingMix();  // 3rd mixing method for relative angular (DeltaPhi x DeltaTheta) spectra
    // stats printout:
    void PrintStats();
    
  protected:
    Selector* selector;
    OutputWriter* outputer;
    int MAXNOEVENT, ncbins;
  };
} // namespace

#endif

