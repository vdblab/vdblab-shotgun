# Mar/17/2025
# This file prepares the `config.yaml` for shotgun pipeline from config_raw.yaml
#
#Usage
#step 1: Prepare config file from sample manifest file:
#python prepare_config_from_sample.py --raw_config config/config_raw.yaml --out_config config/config_toy2.yaml --resource_path /labs/microbiome/apollo/resources/ --directory tmpout/ --sample_manifest manifest_sample_toy2.csv 
#
#step 2: Run pipeline
#PROJECT_NAME="tmp_preprocess"
#nohup snakemake --config nshards=1 stage=preprocess dedup_platform="HiSeq" --directory tmpout_preprocess --configfile config/config_toy2.yaml --jobname "Project:{$PROJECT_NAME}_Sample:{wildcards.sample}_{rule}_{jobid}" --profile run_and_log_light  &>nohup_preprocess &
#
#
## nohup snakemake --config nshards=4 stage=biobakery dedup_platform="HiSeq" --directory tmpout --profile run_and_log_heavy &> nohup_test &


#sample_manifest_file
#3 columns csv file: sampleid, R1, R2


import argparse
import yaml
import csv

parser = argparse.ArgumentParser(description="It prepares config.yaml from config_raw.yaml.",
                                 epilog="Usage: python prepare_config_from_sample.py --raw_config config/config_raw.yaml --out_config MY_config.yaml --resource_path /net/nfs-irwrsrchnas01/labs/rjenq/apollo/resources/ --directory tmpout/ --sample_manifest sample_manifest_file")

parser.add_argument('-rc', '--raw_config',
                    default="config/config_raw.yaml")      # option that takes a value
parser.add_argument('-oc', '--out_config',
                    default="config/config.yaml")
parser.add_argument('-rp', '--resource_path',
                    default="/net/nfs-irwrsrchnas01/labs/rjenq/apollo/resources")
parser.add_argument('-d', '--directory',
                    default="tmpout/")
parser.add_argument('-sm', '--sample_manifest')

args = parser.parse_args()

config_raw_file = args.raw_config
out_config_file = args.out_config
old_prefix = "/data/brinkvd/resources/"
new_prefix = args.resource_path
if(new_prefix[-1]!="/"):
    new_prefix = new_prefix + "/"
# Read YAML file
with open(config_raw_file, 'r') as file:
    data = yaml.safe_load(file)

# Modify values that start with old_prefix
def update_values(d):
    if isinstance(d, dict):
        return {k: update_values(v) for k, v in d.items()}
    elif isinstance(d, list):
        return [update_values(v) for v in d]
    elif isinstance(d, str) and d.startswith(old_prefix):
        return d.replace(old_prefix, new_prefix, 1)
    return d

updated_data = update_values(data)

## Updating config file with sample and fastq fields.
# Adding fastq and sample manifest file
updated_data["sample"] = []  # Initialize as a list
updated_data["R1"] = {}  # Ensure 'R1' key exists
updated_data["R2"] = {}  # Ensure 'R2' key exists
sample_manifest_file = args.sample_manifest

#Read manifest CSV file
with open(sample_manifest_file, newline='') as file:
    reader = csv.DictReader(file)
    for row in reader:  # ✅ Loop is inside the with-block
        sample = row["sample"]
        updated_data["sample"].append(sample)

        updated_data["R1"][sample] = row["R1"].strip().split(",")
        updated_data["R2"][sample] = row["R2"].strip().split(",")

directory_file=args.directory
updated_data["directory"] = directory_file

# Write back to YAML file
with open(out_config_file, 'w') as file:
    yaml.safe_dump(updated_data, file, default_flow_style=False)





