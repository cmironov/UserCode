# Auto generated configuration file
# using: 
# Revision: 1.341.2.2 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s RECO --eventcontent RECO --conditions FrontierConditions_GlobalTag,GR_P_V27A::All --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('ppRECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring('file:/tmp/camelia/B46EE82F-9B14-E111-A4CC-E0CB4E4408D1.root')
                                                   )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool( False ))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.341.2.2 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('/tmp/camelia/reco_ppRECO_mbUpc.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('reconstruction_step')
    )
)
# Additional output definition
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*muons_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*Muons_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_muIsoDepositTk_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_muIsoDepositCalByAssociatorTowers_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*Jets_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_offlineBeamSpot_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*Cluster*_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_photon*_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_Castor*_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_ak7CastorJetID_*_RECO'))

process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_siPixelRecHits_*_ppRECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_siPixelClusters_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiCentrality_*_RECO'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiEvtPlane_*_RECO'))

# Other statements
process.GlobalTag.globaltag = 'GR_P_V27A::All'

# ___________________________________________________________________________________________
# #### set the centrality  selection filter
process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")

process.HeavyIonGlobalParameters = cms.PSet(centralityVariable = cms.string("HFtowers"),
                                            centralitySrc = cms.InputTag("hiCentrality")
                                            )
process.centralityFilter.selectedBins = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];

# #### HLT filters:
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasUPC          = process.hltHighLevel.clone()
process.hltMinBiasUPC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v*",
                                  "HLT_HIUPCNeuEG2Pixel_SingleTrack_v*",
                                  "HLT_HIUPCNeuEG5Pixel_SingleTrack_v*",
                                  "HLT_HIUPCNeuMuPixel_SingleTrack_v*",
                                  "HLT_HIMinBiasZDC_Calo_PlusOrMinus_v*",
                                  "HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v*",
                                  "HLT_HIZeroBias_v*"
                                  ]

process.GoodEventFilterSequence       = cms.Sequence(process.hltMinBiasUPC*process.centralityFilter)

# ___________________________________________________________________________________________
process.reconstruction_step       = cms.Path(process.reconstruction_fromRECO)
process.endjob_step               = cms.EndPath(process.endOfProcess)
process.RECOoutput_step           = cms.EndPath(process.RECOoutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step, process.endjob_step,process.RECOoutput_step)

# filter all path with the good event filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.GoodEventFilterSequence * getattr(process,path)._seq 
