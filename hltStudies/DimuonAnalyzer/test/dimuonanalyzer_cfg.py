import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.MessageLogger = cms.Service("MessageLogger",
                                    cout = cms.untracked.PSet(
       default = cms.untracked.PSet(
            limit = cms.untracked.int32(1) ## kill all messages in the log
            )
        ),
                                    destinations = cms.untracked.vstring('cout')
                                    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                              'rfio:/castor/cern.ch/user/d/dmoon/cms370/Hydjet_MinBias_2.76TeV_Z0_Emb_Reco/Hydjet_MinBias_2.76TeV_Z0Emb_Reco_e10_01_1.root'
                            )
                           )

process.demo = cms.EDAnalyzer('DimuonAnalyzer',
                              genParticle  = cms.InputTag("hiGenParticles"),#no mixing:genParticles
                              muonTracks   = cms.untracked.InputTag("globalMuons"),
                              trackTracks  = cms.untracked.InputTag("hiGlobalPrimTracks"), # pp reco: "generalTracks"),   
                              massMaxDimuon= cms.double(200),
                              massMinDimuon= cms.double(0),
                              pdgDimuon    = cms.double(23),
                              ptMinDimuon  = cms.double(0.),
                              etaMaxMuon   = cms.double(2.5),
                              etaMinMuon   = cms.double(-2.5),
                              ptMinMuon    = cms.double(3.5),
                              etaMaxTrack  = cms.double(2.5),
                              etaMinTrack  = cms.double(-2.5),
                              ptMinTrack   = cms.double(1.)
                              )

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("testJorgeMix.root")
                                   )

process.p = cms.Path(process.demo)
