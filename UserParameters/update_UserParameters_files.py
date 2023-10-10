import subprocess
import importlib.util
import os, sys
import glob
import fileinput

def get_userparams(userparams_input):
  spec = importlib.util.spec_from_file_location(userparams_input[:-3], userparams_input)
  userparams = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(userparams)
  return userparams

def get_file_count(input_format,dataset):
  file_count = 0
  if(input_format == "das"):
    try:
      subprocess.check_output("voms-proxy-info -exists", shell=True)
    except subprocess.CalledProcessError:
      print("No valid proxy certificate found. Please run voms-proxy-init.")
      sys.exit(1)

    file_count_cmd = f"dasgoclient -query='file dataset={dataset}' -limit=0 | wc -l"
    file_count = int(subprocess.check_output(file_count_cmd, shell=True, universal_newlines=True).strip())
  elif (input_format == "local"):
    file_count = len(glob.glob(os.path.join(dataset, "*.root")))
  return file_count

def update_userparams_files(userparams_input, dataset_type, dataset_name, file_count):
  updated_lines = []

  with fileinput.FileInput(userparams_input, inplace=True) as file:
    inside_dataset = False
    inside_dataset_type = False
    for line in file:
      print(line)
      if dataset_type in line:
        inside_dataset_type = True
      if inside_dataset_type and dataset_name in line:
        inside_dataset = True
      if inside_dataset and '"files": 0' in line:
          # Replace "files": 0 with "files": file_count
          line = line.replace('"files": 0', f'"files": {file_count}')
          inside_dataset = False
      updated_lines.append(line)

  # Rewrite the file with the updated lines
  with open(userparams_input, 'w') as file:
      file.writelines(updated_lines)

if __name__=="__main__":

  if len(sys.argv) != 2:
    print("Error missing path to UserParameter config")
    exit(1)

  userparams_input = sys.argv[1]
  userparams = get_userparams(userparams_input)

  update_files=False
  for dataset_type in userparams.datasets:
    for dataset_name in userparams.datasets[dataset_type]:
      if userparams.datasets[dataset_type][dataset_name]["files"] == 0:

        dataset = userparams.datasets[dataset_type][dataset_name]["dataset"]
        input_format = userparams.datasets[dataset_type][dataset_name]["input_format"]

        file_count = get_file_count(input_format,dataset)
        userparams.datasets[dataset_type][dataset_name]["files"] = file_count
        print(dataset_type, ": ", dataset_name,": ", file_count)
    
        update_userparams_files(userparams_input, dataset_type, dataset_name, file_count)