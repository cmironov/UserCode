import FWCore.ParameterSet.Config as cms

process = cms.Process("AccAna")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  #  fileNames = cms.untracked.vstring("file:/tmp/camelia/pythia_Z0jet.root"
                            fileNames = cms.untracked.vstring("file:/tmp/camelia/powhegAlice.root"
        )
)

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.ana = cms.EDAnalyzer("McAcceptAnalyzer",
	genSource = cms.InputTag("genParticles"),    
	hOutputFile = cms.untracked.string(
        "acc_powheg.root"
        ) 
)
process.p = cms.Path(process.ana)
