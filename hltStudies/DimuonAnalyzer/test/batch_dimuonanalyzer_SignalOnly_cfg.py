import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                           fileNames = cms.untracked.vstring(    'rfio:///castor/cern.ch/user/s/silvest/rootfiles/SignalOnly/PARTICLE/root/RECO/PARTICLE_SignalOnly_HIofflineReco_job1_RECO.root')
                            )

process.load("my_code.DimuonAnalyzer.dimuonanalyzer_cfi")
process.demo.massMaxDimuon= cms.double(MASS_MAX_DIMUON)
process.demo.massMinDimuon= cms.double(MASS_MIN_DIMUON)
process.demo.pdgDimuon    = cms.double(PDGID)
process.demo.ptMinDimuon  = cms.double(PT_MIN_DIMUON)
process.demo.etaMaxMuon   = cms.double(ETA_MAX_MU)
process.demo.etaMinMuon   = cms.double(ETA_MIN_MU)
process.demo.ptMinMuon    = cms.double(PT_MIN_MU)
process.demo.etaMaxTrack  = cms.double(ETA_MAX_TRACK)
process.demo.etaMinTrack  = cms.double(ETA_MIN_TRACK)
process.demo.ptMinTrack   = cms.double(PT_MIN_TRACK)
process.demo.genParticle = cms.InputTag("genParticles")
process.demo.muonTracks = cms.InputTag("MUON_TAG")


process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("OUTPUTFILE")
                                   )

process.p = cms.Path(process.demo)
