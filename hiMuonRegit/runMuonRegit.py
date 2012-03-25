import FWCore.ParameterSet.Config as cms

process = cms.Process("REGRECO")

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_44_V10::All'
process.GlobalTag.globaltag      = 'GR_P_V27A::All'
#process.GlobalTag.globaltag = 'STARTHI44_V7::All'


##################################################################################
# setup 'standard'  options
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)
                                       )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:/tmp/camelia/hiReco_RAW2DIGI_RECO_479_1_Qkt.root"),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            skipEvents=cms.untracked.uint32(0),
                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop *_muons*_*_RECO',
                                                                  'drop *_globalMuons_*_RECO'
                                                                  )
                            )
# Output file
process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep recoTracks_hiRegitMuGeneralTracks_*_*'
                                                                         'keep *_remuons*_*_*',
                                                                         'drop *_reglobalMuons_*_*'
                                                                         ),
                                  fileName = cms.untracked.string('/tmp/camelia/regitStd_bjpsi69Gun_v2.root')
                                  )

##################################################################################
# Some Services
process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet(default = cms.untracked.PSet(limit = cms.untracked.int32(0)
                                                                                           )
                                                              ),
                                    destinations = cms.untracked.vstring('cout')
                                    )

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.Timing = cms.Service("Timing")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


##################################################################################
# Schedule definition
# need to rereco the tracker recHits
process.trackerRecHits = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits)

# global iterative tracking
process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.hiTrackReco    = cms.Path(process.heavyIonTracking*process.hiIterTracking)

#muon regit
process.load("RecoHI.HiMuonAlgos.HiReRecoMuon_cff")
process.muRegit        = cms.Path(process.reMuonRecoPbPb)

process.out_step       = cms.EndPath(process.output)
process.endjob_step    = cms.Path(process.endOfProcess)

process.schedule       = cms.Schedule(process.trackerRecHits,process.hiTrackReco,process.muRegit)
process.schedule.extend([process.endjob_step,process.out_step])


# replace everywhere the names of the new muon collections 
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
massSearchReplaceAnyInputTag(process.patMuonSequence, 'muons','remuons')
massSearchReplaceAnyInputTag(process.patMuonSequence, 'globalMuons','reglobalMuons')
massSearchReplaceAnyInputTag(process.outOnia2MuMu, 'hiSelectedTracks','hiRegitMuGeneralTracks')
massSearchReplaceAnyInputTag(process.outTnP, 'hiSelectedTracks','hiRegitMuGeneralTracks')

