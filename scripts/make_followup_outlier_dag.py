import yaml
from astropy.io import fits
from tqdm import tqdm

from paws.filepaths import PathManager
from paws.workflow.manager import WorkflowManager

# 1. Load Configs
CONFIG_FILE = "/home/hoitim.cheung/galacticCenter/config/config.yaml"
TARGET_FILE = "/home/hoitim.cheung/galacticCenter/config/gal.yaml"

with open(CONFIG_FILE, "r") as f:
    config = yaml.safe_load(f)

with open(TARGET_FILE, "r") as f:
    target = yaml.safe_load(f)

paths = PathManager(config, target)
manager = WorkflowManager(config, target)  # PathManager is initialized inside here too

# freq 298 - 310 not yet ready, so start from 310
sat_band_list = [299, 302, 303, 306, 307]
fmin, fmax = 20, 400
cluster = True

#################################################################
is_injection = False
prev_stage = "followup-1"
prev_tcoh = 10
prev_freq_deriv_order = 2

stage = "followup-2"
tcoh = 20
freq_deriv_order = 2

n_sky = 57

if is_injection:
    num_toplist = 1
else:
    num_toplist = 10

zero_threshold = True  # follow-up stage keeps all candidates instead of re-thresholding
max_workers = 32

request_memory = "8GB"
request_disk = "8GB"
request_cpu = 1

exe = "osdf:///igwn/cit/staging/hoitim.cheung/scripts/make_outlier_cli.py"
image = "osdf:///igwn/cit/staging/hoitim.cheung/images/paws.sif"
#################################################################

home_dir = config["home_dir"]
dag_list_path = f"{home_dir}dagFiles/{stage}_{target['name']}_dag{fmin}-{fmax}Hz.txt"

skipped_freqs = []

with open(dag_list_path, "w") as f_daglist:
    for freq in tqdm(range(fmin, fmax), total=(fmax - fmin)):
        if freq in sat_band_list:
            print(f"Skipping saturated band {freq} Hz")
            continue

        prev_taskname = f"{target['name']}_{prev_stage}_TCoh{prev_tcoh}_O{prev_freq_deriv_order}_{freq}Hz"
        taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

        prev_outlier_file = paths.outlier_file(
            freq, prev_taskname, prev_stage, cluster=cluster
        )

        try:
            n_jobs = fits.getdata(prev_outlier_file)["mean2F threshold"].size
        except FileNotFoundError as e:
            print(f"Error loading previous stage outlier file for {freq} Hz: {e}")
            skipped_freqs.append(freq)
            continue

        dag_file = manager.make_outlier_dag(
            CONFIG_FILE,
            TARGET_FILE,
            taskname,
            freq,
            stage,
            freq_deriv_order,
            n_jobs,
            prev_outlier_file,
            exe,
            num_toplist=num_toplist,
            n_sky=n_sky,
            zero_threshold=zero_threshold,
            cluster=cluster,
            separate_saturated=False,
            is_injection=is_injection,
            max_workers=max_workers,
            request_memory=request_memory,
            request_disk=request_disk,
            request_cpu=request_cpu,
            image=image,
        )

        f_daglist.write(f"{dag_file}\n")

if skipped_freqs:
    print(f"Skipped frequencies: {skipped_freqs}")
