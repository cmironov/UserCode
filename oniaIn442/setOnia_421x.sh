cvs co -r V08-00-02 MuonAnalysis/TagAndProbe
cvs co -r V01-06-05 RecoVertex/PrimaryVertexProducer

cvs co -r pPbProd_v04 DataFormats/HeavyIonEvent
cvs co -r pPbProd_v05 HeavyIonsAnalysis/Configuration
cvs co -r pPbProd_v05 RecoHI/HiCentralityAlgos

cvs co -r branch_pA -d HiSkim/HiOnia2MuMu UserCode/tdahms/HiSkim/HiOnia2MuMu
cvs co -r branch_pA -d HiAnalysis/HiOnia UserCode/tdahms/HiAnalysis/HiOnia

cvs co -r HEAD DataFormats/Math/interface/VDTMath.h
cvs co -r V01-15-01 MuonAnalysis/MuonAssociators

rm HeavyIonsAnalysis/Configuration/python/*Skim*

cvs co UserCode/CMironov/oniaIn442
cp UserCode/CMironov/oniaIn442/VDTMath.h DataFormats/Math/interface/VDTMath.h

cp UserCode/CMironov/oniaIn442/HiOnia2MuMuPAT.cc HiSkim/HiOnia2MuMu/src/HiOnia2MuMuPAT.cc
cp UserCode/CMironov/oniaIn442/onia2MuMuPAT_cff.py HiSkim/HiOnia2MuMu/python/onia2MuMuPAT_cff.py
cp UserCode/CMironov/oniaIn442/onia2MuMuPATData_pp2760_CRAB_cfg.py HiSkim/HiOnia2MuMu/test/onia2MuMuPATData_pp2760_CRAB_cfg.py
