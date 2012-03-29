import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_44_V10::All'
#process.GlobalTag.globaltag      = 'GR_P_V27A::All'
process.GlobalTag.globaltag = 'STARTHI44_V7::All'


##################################################################################
# setup 'standard'  options
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10)
                                       )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:/tmp/camelia/mix_jpsi69.root"),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            skipEvents=cms.untracked.uint32(0),
                            )
# Output file
process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = process.RECOEventContent.outputCommands,
                                 # outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
                                  fileName = cms.untracked.string('/tmp/camelia/regitTest_mc_outputTest_light1.root')
                                  )
process.output.outputCommands.extend(cms.untracked.vstring('drop *_hiSelectedTracks_*_RERECO',# default, 1 iteration, global tracking
                                                           'keep *_hiGeneralTracks_*_*', # new collection with improved global tracking
                                                           'keep *_hiGeneralAndRegitMuTracks*_*_*', # tracks from the muon regit on top
                                                           'keep *_remuons*_*_*',
                                                           'keep *_reglobalMuons*_*_*',
                                                           'keep *_retev*_*_*',
                                                           'keep *_recalomuons*_*_*',
                                                           'keep *_hiGenParticles_*_*'
                                                           ))
##################################################################################
# Some Services
'''
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
'''

##################################################################################
# Schedule definition
# need to rereco the tracker recHits
process.trackerRecHits = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits)

# global iterative tracking
process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.hiTrackReco    = cms.Path(process.heavyIonTracking*process.hiIterTracking)
process.heavyIonTracking.remove(process.hiPixelVertices);

#muon regit
process.load("RecoHI.HiMuonAlgos.HiReRecoMuon_cff")
process.muRegit        = cms.Path(process.reMuonRecoPbPb)

process.out_step       = cms.EndPath(process.output)
process.endjob_step    = cms.Path(process.endOfProcess)

process.schedule       = cms.Schedule(process.trackerRecHits,process.hiTrackReco,process.muRegit)
process.schedule.extend([process.endjob_step,process.out_step])


