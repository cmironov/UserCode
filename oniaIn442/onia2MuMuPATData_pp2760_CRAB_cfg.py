# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("Onia2MuMuPAT")

process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT42_V24_AN1::All' 


# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1Tech_BSC_minBias_OR.v*', 'L1Tech_HCAL_HF_coincidence_PM.v*', 'L1Tech_BPTX_plus_AND_minus.v*')
    )
    
# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = ["HLT_DoubleMu3_v3",
                              "HLT_L1DoubleMu0_v1",
                              "HLT_L2DoubleMu0_v2",
                              "HLT_Mu0_v3",
                              "HLT_Mu3_v3",
                              "HLT_Mu5_v3"]
process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True

from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *

onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=True)

process.onia2MuMuPatTrkTrk.addMuonlessPrimaryVertex = False
process.onia2MuMuPatTrkTrk.resolvePileUpAmbiguity   = False

process.primaryVertexFilter.src = "offlinePrimaryVertices"

process.collisionEventSelection = cms.Sequence(process.hfCoincFilter *
                                               process.primaryVertexFilter
                                               #process.siPixelRecHits *
                                               #process.hltPixelClusterShapeFilter
                                               )

process.source.fileNames = cms.untracked.vstring(
    'file:/tmp/camelia/9C1EFD64-63B1-E011-8A2C-00E0817917AB.root'
    )

process.e = cms.EndPath(process.outOnia2MuMu)

process.outOnia2MuMu.fileName = cms.untracked.string( 'onia2MuMuPAT_2760.root' )

#this is a fix for messages like "two EventSetup Producers want to deliver type="CaloSubdetectorGeometry" label="TOWER"
process.es_prefer_CaloTowerGeometryFromDBEP = cms.ESPrefer( "CaloTowerGeometryFromDBEP", "" )
process.es_prefer_CastorGeometryFromDBEP    = cms.ESPrefer( "CastorGeometryFromDBEP", "" )
process.es_prefer_EcalGeometryFromDBEP      = cms.ESPrefer( "EcalBarrelGeometryFromDBEP", "" )
process.es_prefer_EcalEndcapGeometryFromDBEP= cms.ESPrefer( "EcalEndcapGeometryFromDBEP", "" )
process.es_prefer_EcalPreshowerGeometryFromDBEP = cms.ESPrefer( "EcalPreshowerGeometryFromDBEP", "" )
process.es_prefer_HcalGeometryFromDBEP          = cms.ESPrefer( "HcalGeometryFromDBEP", "" )
process.es_prefer_ZdcGeometryFromDBEP           = cms.ESPrefer( "ZdcGeometryFromDBEP", "" )
#------------------------------------------

process.schedule                   = cms.Schedule(process.L1Reco_step,process.Onia2MuMuPAT, process.e)
