# script to set the muon regit

# muon specific part
cvs co -r V00-02-04 RecoHI/HiMuonAlgos
cvs co -r $CMSSW_VERSION RecoHi/HiMuonAlgos/plugin
#Either work in CMSSW_4_4_4 or higher, or check out the following tag for older 4_4_X releases
cvs co -r V01-05-03 RecoHI/HiTracking


scram b -j4
