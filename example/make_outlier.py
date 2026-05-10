import numpy as np
import yaml
from paws.workflow.manager import WorkflowManager
from paws.analysis.outlier import ResultAnalysisManager
from paws.filepaths import PathManager
from astropy.io import fits
from tqdm import tqdm

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
fmin, fmax = 186, 400
f0_band = config["f0_band"]
cluster = True

stage = "injections-3"
tcoh = 40
freq_deriv_order = 3
# n_jobs = 200

for i, freq in tqdm(enumerate(range(fmin, fmax)), total=(fmax - fmin)):
    if freq in sat_band_list:
        print(f"Skipping saturated band {freq} Hz")
        continue

    prev_taskname = f"GalacticCenter_injections-2_TCoh20_O3_{freq}Hz"

    taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

    data_file = paths.outlier_file(freq, prev_taskname, "injections-2", cluster=cluster)

    mean2f_th = fits.getdata(data_file)["mean2F threshold"]
    mean2f_th = np.zeros_like(mean2f_th)
    n_jobs = mean2f_th.size

    result_file = result_manager.make_outlier(
        taskname,
        freq,
        mean2f_th,
        n_jobs,
        num_toplist=1,
        stage=stage,
        freq_deriv_order=freq_deriv_order,
        cluster=cluster,
        work_in_local_dir=False,
        separate_saturated=False,
        is_injection=True,
        max_workers=THREADS,
    )

# for i, freq in tqdm(enumerate(range(fmin, fmax)), total=(fmax-fmin)):

#     prev_taskname = f'GalacticCenter_search-0_TCoh5_O2_{freq}Hz'

#     taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

#     data_file = paths.outlier_file(freq, prev_taskname, 'search-0', cluster=True)

#     mean2f_th = fits.getval(data_file, 'mean2F_th', 0)

#     result_file = result_manager.make_outlier(
#         taskname,
#         freq,
#         mean2f_th,
#         n_jobs,
#         num_toplist=1,
#         stage=stage,
#         freq_deriv_order=freq_deriv_order,
#         cluster=cluster,
#         work_in_local_dir=False,
#         separate_saturated=True,
#         max_workers=THREADS
#     )
