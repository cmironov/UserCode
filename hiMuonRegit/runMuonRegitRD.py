import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_44_V10::All'
process.GlobalTag.globaltag      = 'GR_P_V27::All'

##################################################################################
# setup 'standard'  options
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Input source
process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring("file:/tmp/camelia/486BE09E-6115-E111-A0F4-003048F1110E.root"),#8021
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/m/mironov/regitbigprod/96FB8587-E817-E111-8C7B-001D09F24FBA.root"),
                          #  noEventSort = cms.untracked.bool(True),
                          #  duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                         #   skipEvents=cms.untracked.uint32(500),
                          #  eventsToProcess = cms.untracked.VEventRange('182099:15706562-182099:15706562')#8021
                          #  eventsToProcess = cms.untracked.VEventRange('182324:21724065-182324:21724065')
                           
                            )
# Output file
process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                                                         'keep *_remuons_*_*',
                                                                         'keep *_recalomuons_*_*',
                                                                         'keep *_retevMuons_*_*',
                                                                         'keep *_reglobalMuons_*_*',
                                                                         'keep *_hiGeneralTracks_*_*',
                                                                         'keep *_hiGeneralAndRegitMuTracks_*_*'
                                                                         ),
                                  fileName = cms.untracked.string('regitTest_rd_default.root')
                                  )

##################################################################################
# Some Services
process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet(default = cms.untracked.PSet(limit = cms.untracked.int32(100000000)
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

#---------------------------------------- ADD THIS FOR NEW REGIT
# need to rereco the tracker recHits
process.trackerRecHits = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits)

# global iterative tracking
process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.hiTrackReco    = cms.Path(process.heavyIonTracking)#*process.hiIterTracking)
process.heavyIonTracking.remove(process.hiPixelVertices);

#muon regit
process.load("RecoHI.HiMuonAlgos.HiReRecoMuon_cff")
process.muRegit        = cms.Path(process.reMuonRecoPbPb)

#----------------------------------------DONE ADDING FOR NEW REGIT


process.out_step       = cms.EndPath(process.output)
process.endjob_step    = cms.Path(process.endOfProcess)

process.schedule       = cms.Schedule(process.trackerRecHits,process.hiTrackReco,process.muRegit)
process.schedule.extend([process.endjob_step,process.out_step])


