# script to set the muon regit

# muon specific part
#cvs co -r V00-02-04 RecoHI/HiMuonAlgos

#------------ for first path through real data
cvs co -r RecoHI/HiMuonAlgos


# --- from camelia's user's area, unti lmatt comes back and puts a tag
cvs co UserCode/CMironov/hiMuonRegit
cp UserCode/CMironov/hiMuonRegit/HiRegitMuon*.py RecoHI/HiMuonAlgos/python
cp UserCode/CMironov/hiMuonRegit/hiMuonIterativeTk_cff.py RecoHI/HiMuonAlgos/python
cp UserCode/CMironov/hiMuonRegit/HiReRecoMuon_cff.py RecoHI/HiMuonAlgos/python



#Either work in CMSSW_4_4_4 or higher, or check out the following tag for older 4_4_X releases
cvs co -r V01-05-03 RecoHI/HiTracking


scram b -j4
