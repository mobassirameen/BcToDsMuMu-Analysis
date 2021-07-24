import FWCore.ParameterSet.Config as cms
import os
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('MiniAOD_MC_tifr.txt')
#mylist = FileUtils.loadListFromFile ('BcToDsMuMu_mc_miniAOD_v1.txt')
infiles = cms.untracked.vstring(*mylist)

process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag.globaltag = cms.string('102X_dataRun2_v12')
process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v4') #For  Private MC
#process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v21') #For official MC provided Bhai Chandi Bhai

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    #duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
        #'/store/data/Run2016C/Charmonium/MINIAOD/17Jul2018-v1/20000/9C03CBE2-4B8B-E811-9299-0CC47AC17678.root',
        #'/store/data/Run2016H/DoubleMuonLowMass/MINIAOD/17Jul2018-v1/50000/9A79D8BA-D38B-E811-BCD3-0090FAA57AE0.root',               
	mylist,
        '/store/mc/RunIIAutumn18MiniAOD/DsToPhiPi_ToMuMu_MuFilter_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/7DF46828-0E2C-684D-8305-FA95C167B3AC.root',
    )
)

process.load("slimmedMuonsTriggerMatcher_cfi")

process.load("bctodsmumu-analysis.BcToDsMuMuPAT.BcToDsMuMuRootupler_cfi")
process.rootuple.isMC = cms.bool(True)

process.TFileService = cms.Service("TFileService",

       #fileName = cms.string('BctoDsMuMu_oldMiniAOD_mc_ntuple_10event_test.root'),
       #fileName = cms.string('DsMuMu_new_MiniAOD_mc_ntuple_all_event.root'),
       #fileName = cms.string('DsMuMu_old_MiniAOD_mc_ntuple_geninf.root'),
       #fileName = cms.string('BctoDsMuMu__mc_miniaod-2018_ntuple_test_1_fortime.root'),
       #fileName = cms.string('DstoPhipi_official_mc_miniaod-2018.root'),
       #fileName = cms.string('Phimass_official_mc_miniaod-2018.root'),
       #fileName = cms.string('Dsmas_Phimass_Dspy_pionpt.root'),
       #fileName = cms.string('Dsmas_vertex_Phimass_Dspy_pionpt.root'),
       #fileName = cms.string('Dsmas_vertex_Phimass_Dspy_pionpt_test_for_500event.root'),
       fileName = cms.string('Dsmas_vertex_private_mc_test_for_100event.root'),
)

process.p = cms.Path(process.slimmedMuonsWithTriggerSequence *process.rootuple)
