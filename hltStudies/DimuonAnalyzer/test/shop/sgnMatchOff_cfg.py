import FWCore.ParameterSet.Config as cms

process = cms.Process("Match")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')

process.load('Configuration/StandardSequences/MixingNoPileUp_cff')

process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')

process.load('Configuration/StandardSequences/Reconstruction_cff')
#process.load('Configuration/StandardSequences/ReconstructionHeavyIons_cff')

process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContentHeavyIons_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_3XY_V15::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(JOB_MAXEVENTS)
)
                                  
# Input source
from Correlation.DijetCorrelation.JOB_INPUTFILES import *
process.source = cms.Source("PoolSource",
                            fileNames =  readFiles,
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                           , skipEvents = cms.untracked.uint32(0)
                            )


process.GlobalTag.globaltag = 'MC_3XY_V15::All'

#--------------- matching:
import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.staMuonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
process.staMuonAssociatorByHits.tracksTag                     = 'standAloneMuons:UpdatedAtVtx'
process.staMuonAssociatorByHits.UseTracker                    = False
process.staMuonAssociatorByHits.UseMuon                       = True

import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.glbMuonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
process.glbMuonAssociatorByHits.tracksTag                     = 'globalMuons'
process.glbMuonAssociatorByHits.UseTracker                    = True
process.glbMuonAssociatorByHits.UseMuon                       = True

# --------------
process.mcmatchanalysis = cms.EDAnalyzer("McMatchTrackAnalyzer",
                                         doHLT              = cms.bool(False),
                                         doHiEmbedding      = cms.bool(False),
                                         doParticleGun      = cms.bool(True),
                                         doSecondarySingle  = cms.bool(False),
                                         doReco2Sim         = cms.bool(True),
                                         doSim2Reco         = cms.bool(True),
                                         matchPair          = cms.bool(True),
                                         matchSingle        = cms.bool(True),
                                         pdgPair            = cms.int32(PDGPAIR),
                                         pdgSingle          = cms.int32(13),
                                         type1Tracks        = cms.untracked.InputTag("globalMuons"),
                                         type1MapTag        = cms.untracked.InputTag("glbMuonAssociatorByHits"),
                                         type2Tracks        = cms.untracked.InputTag("standAloneMuons","UpdatedAtVtx"),
                                         type2MapTag        = cms.untracked.InputTag("staMuonAssociatorByHits"),
                                         simTracks          = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
                                         vertices           = cms.untracked.InputTag("hiPixelAdaptiveVertex")
                                         )

#--------------------
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('JOB_OUTPUTFILENAME')
                                   )


process.p7 = cms.Path(process.staMuonAssociatorByHits+process.glbMuonAssociatorByHits)
process.p8 = cms.Path(process.mcmatchanalysis)

process.schedule = cms.Schedule(process.p7,process.p8)
