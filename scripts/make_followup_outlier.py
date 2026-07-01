import numpy as np
import yaml
from astropy.io import fits
from tqdm import tqdm

from paws.analysis.outlier import ResultAnalysisManager
from paws.filepaths import PathManager
from paws.workflow.manager import WorkflowManager
from scripts.make_followup_dag import e

# 1. Load Configs
with open("/home/hoitim.cheung/galacticCenter/config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

with open("/home/hoitim.cheung/galacticCenter/config/gal.yaml", "r") as f:
    target = yaml.safe_load(f)

paths = PathManager(config, target)
manager = WorkflowManager(config, target)  # PathManager is initialized inside here too
result_manager = ResultAnalysisManager(config, target)

# freq 298 - 310 not yet ready, so start from 310
sat_band_list = [299, 302, 303, 306, 307]
THREADS = 32
fmin, fmax = 20, 400
f0_band = config["f0_band"]
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

#################################################################

for i, freq in tqdm(enumerate(range(fmin, fmax)), total=(fmax - fmin)):
    if freq in sat_band_list:
        print(f"Skipping saturated band {freq} Hz")
        continue

    prev_taskname = f"{target['name']}_{prev_stage}_TCoh{prev_tcoh}_O{prev_freq_deriv_order}_{freq}Hz"

    taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

    data_file = paths.outlier_file(freq, prev_taskname, prev_stage, cluster=cluster)

    mean2f_th = fits.getdata(data_file)["mean2F threshold"]
    mean2f_th = np.zeros_like(mean2f_th)
    n_jobs = mean2f_th.size

    result_file = result_manager.make_outlier(
        taskname,
        freq,
        mean2f_th,
        n_jobs,
        num_toplist=num_toplist,
        stage=stage,
        freq_deriv_order=freq_deriv_order,
        n_sky=n_sky,
        cluster=cluster,
        work_in_local_dir=False,
        separate_saturated=False,
        is_injection=is_injection,
        max_workers=THREADS,
    )
