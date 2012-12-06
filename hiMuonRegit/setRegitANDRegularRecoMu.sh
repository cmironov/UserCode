#
# !!!!! have to be in CMSSW_4_4_5_patch1 !!!!!!

cvs co -r $CMSSW_VERSION RecoHI/HiMuonAlgos                
cvs co -r $CMSSW_VERSION RecoHI/Configuration
cvs co -r $CMSSW_VERSION Configuration/StandardSequences

cvs co -d runBothRegitRegular UserCode/CMironov/hiMuonRegit/runBothRegitRegular

cp runBothRegitRegular/HiReRecoMuon_cff.py RecoHI/HiMuonAlgos/python
cp runBothRegitRegular/RecoHiMuon_EventContent_cff.py RecoHI/HiMuonAlgos/python

cp runBothRegitRegular/Reconstruction_HI_cff.py RecoHI/Configuration/python
cp runBothRegitRegular/ReconstructionHeavyIons_cff.py Configuration/StandardSequences/python

rm -rf runBothRegitRegular

scram b -j4

#modify the reconstruction step then: "process.reconstruction_step = cms.Path(process.reconstructionHeavyIons,process.reMuonRecoPbPb)"