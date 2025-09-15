from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DisplacedJets_Run2022C_PromptReco_v1_SM_trigger'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Anl_data_SM_base_cfg.py'

config.Data.inputDataset = '/SingleMuon/Run2022C-PromptReco-v1/AOD' 
config.Data.inputDBS='global'
#config.Data.outputPrimaryDataset = 'HiggsPortal'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
#NJOBS = 1000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_357900_Golden.json'#'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-324209_13TeV_PromptReco_Collisions18_JSON.txt'
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/jingyu/CRAB/Run2022Trigger'
config.Data.publication = False
config.Data.outputDatasetTag = 'DisplacedJets_Run2022C_PromptReco_v1_SM_trigger'

config.Site.storageSite = 'T3_US_FNALLPC'
