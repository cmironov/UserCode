#
# !!!!! have to be in CMSSW_4_4_5_patch1 !!!!!!

cvs co -r V03-00-00 RecoHI/HiMuonAlgos                
cvs co -r V00-00-20 RecoHI/Configuration

scram b -j4

#modify the reconstruction step then: "process.reconstruction_step = cms.Path(process.reconstructionHeavyIons,process.reMuonRecoPbPb)"