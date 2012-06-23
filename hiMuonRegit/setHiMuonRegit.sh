# script to set the muon regit

# muon specific part
#cvs co -r V00-02-04 RecoHI/HiMuonAlgos

#------------ for first path through real data
cvs co -r $CMSSW_VERSION RecoHI/HiMuonAlgos

# added on June 23 ------
cvs co -r $CMSSW_VERSION RecoMuon/MuonIdentification
cvs co -r $CMSSW_VERSION RecoMuon/GlobalTrackingTools
# end of adding

# --- from camelia's user's area
cvs co UserCode/CMironov/hiMuonRegit
cp UserCode/CMironov/hiMuonRegit/HiRegitMuon*.py RecoHI/HiMuonAlgos/python
cp UserCode/CMironov/hiMuonRegit/hiMuonIterativeTk_cff.py RecoHI/HiMuonAlgos/python
cp UserCode/CMironov/hiMuonRegit/HiReRecoMuon_cff.py RecoHI/HiMuonAlgos/python

# added on June 23rd
cp UserCode/CMironov/hiMuonRegit/HIMuonTrackingRegionProducer.h RecoHI/HiMuonAlgos/plugins
cp UserCode/CMironov/hiMuonRegit/RecoMuon/MuonIdProducer.cc RecoMuon/MuonIdentification/plugins
cp UserCode/CMironov/hiMuonRegit/RecoMuon/MuonIdProducer.h RecoMuon/MuonIdentification/plugins
cp UserCode/CMironov/hiMuonRegit/RecoMuon/MuonTrackingRegionBuilder.cc RecoMuon/GlobalTrackingTools/src
# end of adding


#Either work in CMSSW_4_4_4 or higher, or check out the following tag for older 4_4_X releases
cvs co -r V01-05-03 RecoHI/HiTracking


scram b -j4
