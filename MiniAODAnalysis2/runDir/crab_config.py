from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'testRun_Analysis' 
config.General.workArea = 'crabtest'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFileForCrab_cfg.py'

config.Data.inputDataset = '/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFN = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_EE_Estonia'

