from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'LLGDV_Analysis_ttbar_PU20bx25' 
config.General.workArea = 'crabJobs'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_08042015.py'

config.Data.inputDataset = '/TT_Tune4C_13TeV-pythia8-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFN = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'LLGDV_Analysis_AnalysisNTuple_ttbar_PU20bx25'

config.section_("Site")
config.Site.storageSite = 'T2_EE_Estonia'

