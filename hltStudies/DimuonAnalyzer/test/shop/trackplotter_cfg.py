import FWCore.ParameterSet.Config as cms

process = cms.Process("Plot")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories = cms.untracked.vstring('plot',
                                                         'FwkJob', 'FwkReport', 'FwkSummary', 'Root_NoDictionary')

process.MessageLogger.cout = cms.untracked.PSet(
    noTimeStamps = cms.untracked.bool(True),
    threshold = cms.untracked.string('DEBUG'),
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(1000000)
        ),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
        ),
    plot = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
        ),
    
    FwkReport = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(0)
        ),
    FwkSummary = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        limit = cms.untracked.int32(10000000)
    ),
    FwkJob = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
        ),
    Root_NoDictionary = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
        )
    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# fake input file
process.source = cms.Source("EmptySource")

process.plot = cms.EDAnalyzer('McMatchTrackPlotter',
                              inFileName      = cms.untracked.string("JOB_INPUTFILENAME"),
                              typeNumber      = cms.untracked.string("TYPE"),
                              mergedRootFiles =cms.untracked.bool(False),
                              doBarel         = cms.untracked.bool(DOBARREL),
                              doEndcap        = cms.untracked.bool(DOENDCAP),
                              doParentSpecificCut = cms.untracked.bool(True),
                              doSingle        = cms.untracked.bool(True),
                              doPair          = cms.untracked.bool(True),
                              hitFractionMatch= cms.untracked.double(HITFRACMATCH),
                              pdgPair         = cms.untracked.int32(PDGPAIR),
                              ptCutSingle     = cms.untracked.double(PTCUTSINGLE),
                              ptCutPair       = cms.untracked.double(PTCUTPAIR),
                              massCutMax      = cms.untracked.double(MASSMAX),
                              massCutMin      = cms.untracked.double(MASSMIN),
                              reco2Sim        = cms.untracked.bool(True),
                              sim2Reco        = cms.untracked.bool(True)
                              )

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("JOB_OUTPUTFILENAME")
                                   )

process.p = cms.Path(process.plot)
