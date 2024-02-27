import argparse
# import importlib.util
import imp
import uuid
import os
import math
# from SubmissionManager import SubmissionManager, SubmissionSystem

# from Logger import *

def get_args():
  parser = argparse.ArgumentParser(description="Submitter")

  parser.add_argument("--config", type=str, required=True, help="config with settings and datasets to be executed by the app")
    
  parser.add_argument("--local", action="store_true", default=False, help="run locally")
  parser.add_argument("--condor", action="store_true", default=False, help="run on condor")

  parser.add_argument("--app", type=str, default="run_LLPnanoAOD.py", help="Python or executable app to run. Default: run_LLPnanoAOD.py")
  
  parser.add_argument(
    "--job_flavour", 
    type=str, 
    default="longlunch",
    help="condor job flavour: espresso (20 min), microcentury (1h), longlunch (2h), workday (8h), tomorrow (1d), testmatch (3d), nextweek (1w)"
  )
  parser.add_argument(
    "--run_time", 
    type=int, 
    default=0,
    help="condor run time (seconds). Default: 0 = will use job_flavour"
  )
  parser.add_argument(
    "--resubmit_job", 
    type=int,
    default=None,
    help="use this option to resubmitt a specific job"
  )
  
  parser.add_argument("--dry", action="store_true", default=False, help="dry run, without submitting to condor")
  
  args = parser.parse_args()
  return args


def get_config(config_path):
  print("Reading config from path: ", config_path)
  config = imp.load_source('files_module', config_path)
  return config


def update_config(path, key, value):
  with open(path, "r") as f:
    lines = f.readlines()
  with open(path, "w") as f:
    for line in lines:
      if line.strip().startswith(key.strip()):
        line = "{0} {1}".format(key, value)
      f.write(line)


def prepare_tmp_files(args):
  hash_string = str(uuid.uuid4().hex[:6])
  tmp_config_path = "tmp/config_" + hash_string + ".py"
  
  print("Creating a temporary config: ", tmp_config_path)
  os.system("mkdir -p tmp")
  os.system("cp " + args.config + " " + tmp_config_path)
  
  return tmp_config_path
  

def get_dataset_files_list(dataset_name):
  das = True
  if dataset_name.startswith("/nfs/") or dataset_name.startswith("/afs/"):
    das = False
  if das:
    das_command = "dasgoclient -query='file dataset=" + dataset_name + "'"
    print("\n\nExecuting ", das_command)
    files = os.popen(das_command).read().splitlines()
  else:
    filenames = os.listdir(dataset_name)
    files = [filename for filename in filenames if filename.endswith(".root")]

  return files


def create_condor_directories():
  for path in ("error", "log", "output", "tmp"):
    if not os.path.exists(path):
      os.makedirs(path)


def make_output_paths(config):
  path = config.output_dir + "/LLPnanoAOD"
  if not os.path.exists(path):
      os.makedirs(path)
  if config.run_mini_and_nano:
    path = config.output_dir + "/LLPminiAOD"
    if not os.path.exists(path):
        os.makedirs(path)


def is_dataset_from_das(datasetpath):
  if datasetpath.startswith("/nfs/") or datasetpath.startswith("/afs/"):
    return False
  return True


def is_data(datasetpath):
  if is_dataset_from_das(datasetpath):
    if not datasetpath.endswith("SIM"):
      return True
  return False


def set_run_time(job_flavour, run_time):
  flavor_to_seconds = {
    'espresso': 20 * 60,            # 20 minutes
    'microcentury': 60 * 60,        # 1 hour
    'longlunch': 2 * 60 * 60,       # 2 hours
    'workday': 8 * 60 * 60,         # 8 hours
    'tomorrow': 24 * 60 * 60,       # 1 day
    'testmatch': 3 * 24 * 60 * 60,  # 3 days
    'nextweek': 7 * 24 * 60 * 60,   # 1 week
  }
  if(run_time == 0):
    if job_flavour not in flavor_to_seconds:
      return None
    return flavor_to_seconds[job_flavour]
  return run_time


def setup_condor_submit_files(config, condor_run_script_name, hash_string, run_time, resubmit_job, n_jobs):
  condor_config_name = "tmp/condor_config_%s.sub" % hash_string

  os.system("cp LLPNanoAOD/LLPnanoAOD/test/templates/condor_config.template.sub " + condor_config_name)
  print("Stored condor config at: " + condor_config_name)

  condor_run_script_name_escaped = condor_run_script_name.replace("/", "\/")
  os.system("sed -i 's/<executable>/{}/g' {}".format(condor_run_script_name_escaped, condor_config_name))
  os.system("sed -i 's/<run_time>/{}/g' {}".format(run_time, condor_config_name))

  if resubmit_job is not None:
    os.system("sed -i 's/$(ProcId)/{}/g' {}".format(resubmit_job_escaped, condor_config_name))
    os.system("sed -i 's/<n_jobs>/1/g' '{}'".format(escaped_condor_config_name))
  else:
    jobs = int(n_jobs)
    if config.maxJobs > 0 and config.maxJobs < n_jobs:
      jobs = config.maxJobs
    os.system("sed -i 's/<n_jobs>/{}/g' '{}'".format(jobs, condor_config_name))
  
  return condor_config_name


def setup_run_files(config):
  hash_string = str(uuid.uuid4().hex[:6])

  if config.run_mini_and_nano:
    condor_run_script_name = "tmp/run_LLPnanoAOD_from_AOD_%s.sh" % hash_string
  else:
    condor_run_script_name = "tmp/run_LLPnanoAOD_only_%s.sh" % hash_string
  input_files_list_file_name = "tmp/input_files_%s.txt" % hash_string

  
  if config.run_mini_and_nano:
    os.system("cp LLPNanoAOD/LLPnanoAOD/test/templates/run_LLPnanoAOD_from_AOD.template.sh " + condor_run_script_name)
  else:
    os.system("cp LLPNanoAOD/LLPnanoAOD/test/templates/run_LLPnanoAOD_only.template.sh " + condor_run_script_name)
  os.system("chmod 700 " + condor_run_script_name)
  print("Stored run shell script at: " + condor_run_script_name)
  
  dataset_files = get_dataset_files_list(config.dataset)
  n_jobs = math.ceil(len(dataset_files) / config.filesPerJob)

  with open(input_files_list_file_name, "w") as file:
    for input_file_path in dataset_files:
        file.write(input_file_path + "\n")

  dataset_name = os.path.basename(config.output_dir)
  output_file = config.output_dir + "/LLPnanoAOD/" + dataset_name
  output_file = output_file.replace("/", "\/")
    
  voms_proxy_path = os.popen("voms-proxy-info -path").read().strip().replace("/", "\/")    
  os.system("cp " + voms_proxy_path + " voms_proxy")
  
  home_path = config.home_path.replace("/", "\/")
  os.system("sed -i 's/<home_path>/{}/g' {}".format(home_path,condor_run_script_name))
  dataset = config.dataset.replace("/", "\/")
  os.system("sed -i 's/<dataset>/{}/g' {}".format(dataset,condor_run_script_name))

  input_files_list_path = input_files_list_file_name.replace("/", "\/")
  os.system("sed -i 's/<input_files_path>/{}/g' {}".format(input_files_list_path,condor_run_script_name))
  os.system("sed -i 's/<output_path>/{}/g' {}".format(output_file,condor_run_script_name))
  os.system("sed -i 's/<n_events>/{}/g' {}".format(config.nEvents,condor_run_script_name))
  os.system("sed -i 's/<n_files>/{}/g' {}".format(config.filesPerJob,condor_run_script_name))

  os.system("sed -i 's/<is_data>/{}/g' {}".format(config.run_on_data,condor_run_script_name))
  os.system("sed -i 's/<run_on_das>/{}/g' {}".format(config.run_on_das,condor_run_script_name))

  includeDSAMuon = True if "DSAMuon" in config.LLPcollections else False
  includeBS = True if "BS" in config.LLPcollections else False
  includeGenPart = True if "GenPart" in config.LLPcollections else False
  includeDGLMuon = True if "DGLMuon" in config.LLPcollections else False
  os.system("sed -i 's/<include_DSAMuon>/{}/g' {}".format(includeDSAMuon,condor_run_script_name))
  os.system("sed -i 's/<include_BS>/{}/g' {}".format(includeBS,condor_run_script_name))
  os.system("sed -i 's/<include_GenPart>/{}/g' {}".format(includeGenPart,condor_run_script_name))
  os.system("sed -i 's/<include_DGLMuon>/{}/g' {}".format(includeDGLMuon,condor_run_script_name))

  workDir = os.getcwd().replace("/", "\/")
  os.system("sed -i 's/<work_dir>/{}/g' {}".format(workDir,condor_run_script_name))
  
  if config.run_mini_and_nano:
    LLPminiAOD_file = config.output_dir + "/LLPminiAOD/" + dataset_name
    LLPminiAOD_file = LLPminiAOD_file.replace("/", "\/")
    os.system("sed -i 's/<LLPminiAOD_path>/{}/g' {}".format(LLPminiAOD_file,condor_run_script_name))
    os.system("sed -i 's/<save_LLPminiAOD>/{}/g' {}".format(config.save_mini,condor_run_script_name))
    print("save_mini: ", config.save_mini)

  return condor_run_script_name, hash_string, n_jobs


def run_locally(condor_run_script_name, process_idx):
  command = "./" + condor_run_script_name + " " + str(process_idx)
  os.system(command)


def run_condor(condor_config_name, dry):
  command = "condor_submit " + condor_config_name
  print("Submitting to condor: " + command)
  if not dry:  
    os.system(command)


def main():
  args = get_args()
  
  if not args.condor and not args.local:
    print("Please select either --local or --condor")
    exit()
  
  config = get_config(args.config)
  if not hasattr(config, "datasets"):
    print("datasets missing in config")
    exit()
  
  tmp_configs_paths = []
  datasets = config.datasets
    
  for datasetpath, dataset in datasets:
    tmp_config_path = prepare_tmp_files(args)
    update_config(tmp_config_path, "dataset = ", "\"%s\"\n" % dataset)
    update_config(tmp_config_path, "output_dir = ", "\"%s/%s\"\n" % (config.output_base_path, datasetpath))
    run_on_das = is_dataset_from_das(dataset)
    run_on_data = is_data(dataset)
    update_config(tmp_config_path, "run_on_data = ", "\"%s\"\n" % run_on_data)
    update_config(tmp_config_path, "run_on_das = ", "\"%s\"\n" % run_on_das)
    tmp_configs_paths.append(tmp_config_path)

  if args.condor:
    run_time = set_run_time(args.job_flavour, args.run_time)
    if run_time == None:
      print("Error: invalid job_flaovur {}.".format(args.job_flavour))
  
  for i, tmp_config_path in enumerate(tmp_configs_paths):
    if args.condor:
      create_condor_directories()

    tmp_config = get_config(tmp_config_path)

    make_output_paths(tmp_config)

    condor_run_script_name, hash_string, n_jobs = setup_run_files(tmp_config)

    if args.local:
      run_locally(condor_run_script_name, i)
    if args.condor:
      condor_config_name = setup_condor_submit_files(tmp_config, condor_run_script_name, hash_string, run_time, args.resubmit_job, n_jobs)
      run_condor(condor_config_name, args.dry)


if __name__=="__main__":
  main()
