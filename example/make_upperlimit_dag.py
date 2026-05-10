# /home/hoitim.cheung/.conda/envs/paws/bin/python
import numpy as np
import yaml
from paws.workflow.manager import WorkflowManager
from paws.filepaths import PathManager
from astropy.io import fits
from paws.definitions import phase_param_name

from tqdm import tqdm

# 1. Load Configs
config_file = "/home/hoitim.cheung/galacticCenter/config/config.yaml"
with open(config_file, "r") as f:
    config = yaml.safe_load(f)

target_file = "/home/hoitim.cheung/galacticCenter/config/gal.yaml"
with open(target_file, "r") as f:
    target = yaml.safe_load(f)

manager = WorkflowManager(
    config=config, target=target
)  # PathManager is initialized inside here too
paths = PathManager(config=config, target=target)

obs_day = 630
coh_day = 5
stage = "upperlimit-v2"
freq_deriv_order = 2
inj_freq_deriv_order = 4

fmin, fmax = 20, 400
metric_file = "osdf:///igwn/cit/staging/hoitim.cheung/metricSetup/Start1368970000_TCoh432000_N107_Spin2.fts"

exe = "osdf:///igwn/cit/staging/hoitim.cheung/scripts/calUL.py"

image = "osdf:///igwn/cit/staging/hoitim.cheung/images/paws_gc.sif"

h0_file = "/home/hoitim.cheung/gc/results/upperLimit-rt-200inj/NGC6544/NGC6544_upperLimit-rt-200inj_20-475Hz_clustered.txt"
freq, h0_arr, _ = np.loadtxt(h0_file).T

_, freq_deriv_param_names = phase_param_name(freq_deriv_order)

dag_list_path = f"/home/hoitim.cheung/galacticCenter/dagFiles/{stage}_{target['name']}_dag{fmin}-{fmax}Hz.txt"

skipped_freqs = []

with open(dag_list_path, "w") as f_daglist:
    for freq, h0_est in tqdm(
        zip(range(fmin, fmax), h0_arr), desc="Generating DAGs", total=fmax - fmin
    ):
        # --- Collect SFT files ---
        sft_files = []
        files = paths.sft_ensemble(freq, use_osdf=True)
        sft_files.extend(files)  # Use extend, not append, to flatten the list

        data_taskname = f"GalacticCenter_search-0_TCoh5_O2_{freq}Hz"

        data_file = paths.outlier_file(freq, data_taskname, "search-0", cluster=True)

        df_grid = [
            fits.getval(data_file, param_name, 0)
            for param_name in freq_deriv_param_names
        ]

        mean2f_th = fits.getval(data_file, "mean2F_th", 0)

        non_sat_bands = fits.getdata(data_file, extname="non_sat_band")["non_sat_band"]

        if len(non_sat_bands) == 0:
            print(
                f"Warning: No non-saturated bands found for {freq} Hz. Skipping DAG generation."
            )
            skipped_freqs.append(freq)
            continue
        # --- Make DAG ---
        # Generate task name
        taskname = f"GalacticCenter_{stage}_TCoh5_O2_{freq}Hz"

        dagfile = manager.make_upperlimit_dag(
            config_file,
            target_file,
            taskname,
            freq,
            stage,
            freq_deriv_order,
            sft_files,
            metric_file,
            mean2f_th,
            non_sat_bands,
            exe,
            df_grid=df_grid,
            inj_freq_deriv_order=4,
            num_toplist=1,
            sky_uncertainty=1e-4,
            h0_est=h0_est,
            n_inj=200,
            request_memory="16GB",
            request_disk="6GB",
            request_cpu=8,
            cluster=True,
            work_in_local_dir=True,
            save_intermediate=False,
            image=image,
        )

        f_daglist.write(f"{dagfile}\n")

if skipped_freqs:
    print(f"Skipped frequencies due to no non-saturated bands: {skipped_freqs}")
