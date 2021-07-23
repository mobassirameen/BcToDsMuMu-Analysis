import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('BcToDsMuMu',
                          dimuons = cms.InputTag("slimmedMuonsWithTrigger"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secondaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          prescales     = cms.InputTag("patTrigger"),
                          objects       = cms.InputTag("slimmedPatTrigger"),
                          TriggerNames  = cms.vstring(
                              "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"),
                          PuInfoTag      = cms.InputTag("slimmedAddPileupInfo"),
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          )
