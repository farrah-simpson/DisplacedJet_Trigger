from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DisplacedJets_Run2023C_PromptReco_v4_trigger'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Anl_data_HT_base_PF_cfg.py'

config.Data.inputDataset = '/DisplacedJet/Run2023C-PromptReco-v4/MINIAOD' 
config.Data.inputDBS='global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
#NJOBS = 1000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'#'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/DCSOnly_JSONS/dailyDCSOnlyJSON/Collisions23_13p6TeV_366403_367554_DCSOnly_TkPx.json' 
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/fsimpson/CRAB/Run2023Trigger'
config.Data.publication = False
config.Data.outputDatasetTag = 'DisplacedJets_Run2023C_miniAOD_v4_trigger'

config.Site.storageSite = 'T3_US_FNALLPC'
