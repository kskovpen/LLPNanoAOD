output_base_path = "/nfs/dust/cms/user/lrygaard/ttalps_cms"

nEvents = 2000

filesPerJob = 1
maxJobs = 100

run_mini_and_nano = True
save_mini = False

dbs_instance = "prod/global"
# dbs_instance = "prod/phys03"

LLPcollections = (
  "DSAMuon",
  "BS",
  "GenPart",
  # "DGLMuon",
  "RefittedTracks",
)

# Run 3 years options for global tags:
# Data: 2022ReReco, 2022Prompt, 2023
# MC = 2022PreEE, 2022PostEE, 2023PreBPix, 2023PostBPix
# year = "2022PostEE"
year = "2022ReReco"

datasets = (
  ("backgrounds2022PostEE/TTtoLNu2Q", "/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),
  # ("backgrounds2022PostEE/TWminustoLNu2Q", "/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),
  # ("backgrounds2022PostEE/TbarWplustoLNu2Q", "/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),
  ("backgrounds2022PostEE/TTLL_MLL-4to50", "/TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),
  ("backgrounds2022PostEE/TTLL_MLL-50", "/TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),
  ("backgrounds2022PostEE/TTLNu-1Jets", "/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4/AODSIM"),
  # ("backgrounds2022PostEE/ttHto2B_M-125", "/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM"),

    # ("data2022/SingleMuonC", "/SingleMuon/Run2022C-27Jun2023-v1/AOD"),
)

dataset = ""
output_dir = ""
run_on_das = ""
run_on_data = ""
