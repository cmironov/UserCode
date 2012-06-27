# Auto generated configuration file
# using: 
# Revision: 1.177 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: ZMM_cfi -s GEN:ProductionFilterSequence -n 1000 --conditions FrontierConditions_GlobalTag,MC_37Y_V5::All --datatier GEN --eventcontent FEVTDEBUGHLT --no_exec --python_filename ZMM_Gen_cfg.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.177 $'),
    annotation = cms.untracked.string('ZMM_cfi nevts:1000'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)
process.options = cms.untracked.PSet(

)
# Input source
process.source = cms.Source("EmptySource")

#process.RandomNumberGeneratorService.generator.initialSeed = 65485432 
#process.RandomNumberGeneratorService.VtxSmeared.initialSeed = 76546234
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string(
        '/tmp/camelia/pythia_Z0photonjet.root'
      # '/tmp/camelia/pythia_Z0.root'
        ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
from Configuration.Generator.PyquenDefaultSettings_cff import *
process.GlobalTag.globaltag = 'STARTHI44_V7::All'
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(2760.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(pyquenPythiaDefaultBlock,
                                parameterSets = cms.vstring('pythiaUESettings',
                                                            'pythiaZmumu',
                                                            'kinematics'
                                                            ),
                                pythiaZmumu = cms.vstring('MSEL = 0 ! users defined processes only',
                                                          'MSUB(15)=1          !qq->Z0/gamma*+g',
                                                          'MSUB(30)=1          !qg->Z0/gamma*+q',  
                                                        #   'MSTP(43)=2  !Only Z0',
                                                          'MDME( 174,1) = 0    !Z decay into d dbar', 
                                                          'MDME( 175,1) = 0    !Z decay into u ubar', 
                                                          'MDME( 176,1) = 0    !Z decay into s sbar', 
                                                          'MDME( 177,1) = 0    !Z decay into c cbar', 
                                                          'MDME( 178,1) = 0    !Z decay into b bbar', 
                                                          'MDME( 179,1) = 0    !Z decay into t tbar', 
                                                          'MDME( 182,1) = 0    !Z decay into e- e+', 
                                                          'MDME( 183,1) = 0    !Z decay into nu_e nu_ebar', 
                                                          'MDME( 184,1) = 1    !Z decay into mu- mu+', 
                                                          'MDME( 185,1) = 0    !Z decay into nu_mu nu_mubar', 
                                                          'MDME( 186,1) = 0    !Z decay into tau- tau+', 
                                                          'MDME( 187,1) = 0    !Z decay into nu_tau nu_taubar', 
                                                          ),
                                kinematics = cms.vstring ('CKIN(3) = -1       !(D=0 GeV) lower lim pT_hat',
                                                          'CKIN(4) = 60.      !(D=-1 GeV) upper lim pT_hat, if < 0 innactive',
                                                          "CKIN(7)=-3.5",  #min rapidity
                                                          "CKIN(8)=3.5"    #max rapidity
                                                          )
                                )
                                 )
process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.endjob_step,process.out_step)

# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq
