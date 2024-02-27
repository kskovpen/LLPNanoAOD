home_path = "/afs/desy.de/user/l/lrygaard"

nEvents = 10

filesPerJob = 1
maxJobs = 1

run_mini_and_nano = True
save_mini = False

homePath = "/afs/desy.de/user/l/lrygaard"
output_base_path = "/nfs/dust/cms/user/lrygaard/ttalps_cms"

LLPcollections = (
  "DSAMuon",
  "BS",
  "GenPart",
  # "DGLMuon",
)

datasets = (
  # ("backgrounds2018/QCD_Pt-15To20", "/QCD_Pt-15To20_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/QCD_Pt-20To30", "/QCD_Pt-20To30_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/QCD_Pt-30To50", "/QCD_Pt-30To50_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/QCD_Pt-50To80", "/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/QCD_Pt-80To120", "/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/QCD_Pt-120To170", "/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-170To300", "/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-300To470", "/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-470To600", "/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-600To800", "/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-800To1000", "/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/QCD_Pt-1000", "/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),

  ("backgrounds2018/TTToSemiLeptonic", "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # # ("backgrounds2018/TTToHadronic", "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # # ("backgrounds2018/TTTo2L2Nu", "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/ST_tW_antitop", "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/ST_t-channel_antitop", "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/ST_tW_top", "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # # ("backgrounds2018/ST_t-channel_top", "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # # ("backgrounds2018/ttZJets", "/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/TTZToLLNuNu_M-10", "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/TTZToLL_M-1to10", "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/TTZZ", "/TTZZ_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/TTZH", "/TTZH_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v3/AODSIM"),
  # ("backgrounds2018/TTTT", "/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # # ("backgrounds2018/ttWJets", "/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # # ("backgrounds2018/TTWJetsToLNu", "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/ttHToMuMu", "/ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/ttHTobb", "/ttHTobb_ttToSemiLep_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/ttHToNonbb", "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # # ("backgrounds2018/WJetsToLNu", "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM"),
  # ("backgrounds2018/DYJetsToMuMu_M-50", "/DYJetsToMuMu_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),
  # ("backgrounds2018/DYJetsToMuMu_M-10to50", "/DYJetsToMuMu_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM"),

  # ("data2018/SingleMuonA", "/SingleMuon/Run2018A-12Nov2019_UL2018-v5/AOD"),
  # ("data2018/SingleMuonB", "/SingleMuon/Run2018B-12Nov2019_UL2018-v3/AOD"),
  # ("data2018/SingleMuonC", "/SingleMuon/Run2018C-12Nov2019_UL2018-v3/AOD"),
  # ("data2018/SingleMuonD", "/SingleMuon/Run2018D-12Nov2019_UL2018-v8/AOD"),

  # ("signals_test/tta_mAlp-0p35GeV_ctau-1e0mm", "/nfs/dust/cms/user/lrygaard/ttalps_cms/signals_RECO/tta_mAlp-0p35GeV_ctau-1e0mm_nEvents-100"),
  # ("signals/tta_mAlp-0p35GeV_ctau-1e1mm", "/nfs/dust/cms/user/lrygaard/ttalps_cms/signals_RECO/tta_mAlp-0p35GeV_ctau-1e1mm_nEvents-100"),
  # ("signals/tta_mAlp-0p35GeV_ctau-1e2mm", "/nfs/dust/cms/user/lrygaard/ttalps_cms/signals_RECO/tta_mAlp-0p35GeV_ctau-1e2mm_nEvents-100"),
  # ("signals/tta_mAlp-0p35GeV_ctau-1e3mm", "/nfs/dust/cms/user/lrygaard/ttalps_cms/signals_RECO/tta_mAlp-0p35GeV_ctau-1e3mm_nEvents-100"),
  # ("signals/tta_mAlp-0p35GeV_ctau-1e5mm", "/nfs/dust/cms/user/lrygaard/ttalps_cms/signals_RECO/tta_mAlp-0p35GeV_ctau-1e5mm_nEvents-1000"),
)

dataset = ""
output_dir = ""
run_on_das = ""
run_on_data = ""