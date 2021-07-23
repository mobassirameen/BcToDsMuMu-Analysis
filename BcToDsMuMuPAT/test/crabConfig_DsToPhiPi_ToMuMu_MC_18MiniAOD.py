from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'job_crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_18A_v1_Jhovanny'
config.General.requestName = 'job_crab_mc_DsToPhiPi_ToMuMu_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1'
#config.General.workArea = 'crab_data_Charmonium_finaljob_MINIAOD_CMSSW10218_18A_v1_Jhovanny'
config.General.workArea = 'crab_mc_DsToPhiPi_ToMuMu_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'Psiks0Rootupler.py'
config.JobType.psetName = 'BcToDsMuMuRootupler.py'

#config.Data.inputDataset = '/Charmonium/Run2018A-17Sep2018-v1/MINIAOD'
config.Data.inputDataset = '/DsToPhiPi_ToMuMu_MuFilter_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'

config.Data.outLFNDirBase = '/store/user/moameen/'
config.Data.publication = False

config.Site.storageSite = 'T2_IN_TIFR'
