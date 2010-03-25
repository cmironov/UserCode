import FWCore.ParameterSet.Config as cms

process = cms.Process("Match")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/HiEventMixing_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('SimGeneral/MixingModule/himixGEN_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('SimGeneral/MixingModule/himixSIMExtended_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('SimGeneral/MixingModule/himixDIGI_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('Configuration/StandardSequences/RawToDigi_cff')
process.load('Configuration/StandardSequences/ReconstructionHeavyIons_cff')
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
#==============================
process.load("SimGeneral.TrackingAnalysis.HiTrackingParticles_cff")
process.GlobalTag.globaltag = 'MC_3XY_V15::All'

process.mcmatchanalysis = cms.EDAnalyzer("McMatchTrackAnalyzer",
                                         doHLT              = cms.bool(False),
                                         doHiEmbedding      = cms.bool(True),
                                         doParticleGun      = cms.bool(True),
                                         doSecondarySingle  = cms.bool(False),
                                         doReco2Sim         = cms.bool(True),
                                         doSim2Reco         = cms.bool(True),
                                         matchPair          = cms.bool(True),
                                         matchSingle        = cms.bool(True),
                                         pdgPair            = cms.int32(PDGPAIR),
                                         pdgSingle          = cms.int32(13),
                                         type1Tracks        = cms.untracked.InputTag("hltL3Muons"),
                                         type1MapTag        = cms.untracked.InputTag("tpToL3MuonAssociation"),
                                         type2Tracks        = cms.untracked.InputTag("hltL2Muons","UpdatedAtVtx"),
                                         type2MapTag        = cms.untracked.InputTag("tpToL2UpdMuonAssociation"),
                                         simTracks          = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
                                         vertices           = cms.untracked.InputTag("hiPixelAdaptiveVertex")
                                         )

#--------------------
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('JOB_OUTPUTFILENAME')
                                   )

process.load("Validation.RecoMuon.associators_cff")
# to check this is working
process.p7 = cms.Path(process.muonAssociationHLT_seq)
process.p8 = cms.Path(process.mcmatchanalysis)

process.schedule = cms.Schedule(process.p7,process.p8)
