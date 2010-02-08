import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:mixHI_sgn1.root'
        )
                            )

process.demo = cms.EDAnalyzer('DimuonAnalyzer',
                              genParticle  = cms.InputTag("hiGenParticles"),#no mixing:genParticles
                              muonTracks   = cms.untracked.InputTag("globalMuons"),
                              trackTracks  = cms.untracked.InputTag("hiGlobalPrimTracks"), # pp reco: "generalTracks"),   
                              massMaxDimuon= cms.double(120),
                              massMinDimuon= cms.double(70),
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
                                   fileName = cms.string("dimuonGenEff.root")
                                   )

process.p = cms.Path(process.demo)
