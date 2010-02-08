# Auto generated configuration file
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECOMIX')

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('rfio:///castor/cern.ch/user/j/jrobles/HIon/cms340/Z0/root/HIzmumu_job_99.root'),
#'/store/relval/CMSSW_3_4_0/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V14-v1/0008/FCC672B1-99E9-DE11-9313-001731AF6BCD.root'),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            skipEvents = cms.untracked.uint32(0)
                            )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  fileName = cms.untracked.string('mixHI_reco.root'),
                                  dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW-RECO'),
        filterName = cms.untracked.string('')
        )
                                  )

# Other statements
process.GlobalTag.globaltag = 'MC_3XY_V15::All' 

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)


# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.out_step)


