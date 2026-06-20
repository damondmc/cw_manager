import numpy as np
import yaml
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

from paws.definitions import phase_param_name
from paws.filepaths import PathManager
from paws.params.followup import FollowUpParamGenerator
from paws.params.models import PowerLawModel
from paws.workflow.manager import WorkflowManager

# 1. Load Configs
with open("/home/hoitim.cheung/galacticCenter/config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

home_dir = config["home_dir"]

with open(f"{home_dir}config/gal.yaml", "r") as f:
    target = yaml.safe_load(f)

paths = PathManager(config, target)
manager = WorkflowManager(config, target)

fmin = 20
fmax = 400
use_osg = True
use_osdf = True
cluster = True
is_injection = True  # True to carry injections from prev stage into DAG

################################################
prev_stage = "search-0"
prev_tcoh = 5
prev_freq_deriv_order = 2
################################################

################################################
stage = "followup-1"
tcoh = 10
freq_deriv_order = 2
################################################

################################################
inj_freq_deriv_order = 4
n_seg = 54
tasks_per_job = 1000
sky_radius = 0
spacing_alpha = None
spacing_delta = None
sky_grid_file = "tests/results_sky/gc_sky_grid.txt"  # actual Weave sky grid offsets (d_alpha, d_delta)
################################################

extra_stats = "coh2F_det,mean2F,coh2F_det,mean2F_det"
weave_exe = config["executables"]["weave"]
num_top_list = config["num_toplist"]
metric_file = f"osdf:///igwn/cit/staging/hoitim.cheung/metricSetup/o4ab_t{tcoh}_s{freq_deriv_order}.fts"

_, freq_deriv_param_names = phase_param_name(prev_freq_deriv_order)

skipped_freqs = []

dag_list_path = f"{home_dir}dagFiles/{stage}_{target['name']}_dag{fmin}-{fmax}Hz.txt"

with open(dag_list_path, "w") as f_daglist:
    for freq in tqdm(range(fmin, fmax), desc="Generating DAGs", total=fmax - fmin):
        sft_files = []
        files = paths.sft_ensemble(freq, use_osdf=use_osdf)
        sft_files.extend(files)

        data_taskname = f"{target['name']}_{prev_stage}_TCoh{prev_tcoh}_O{prev_freq_deriv_order}_{freq}Hz"

        data = []
        injection_data = None

        try:
            data_file = paths.outlier_file(
                freq, data_taskname, prev_stage, cluster=cluster
            )
            data = fits.getdata(data_file, ext=1)
            if is_injection:
                injection_data = fits.getdata(data_file, extname="injection")
        except FileNotFoundError as e:
            print(f"Error loading data for frequency {freq} Hz: {e}")

        if len(data) == 0:
            print(
                f"No outliers found for frequency {freq} Hz. Skipping follow-up generation."
            )
            skipped_freqs.append(freq)
            continue

        df_grid = [
            fits.getval(data_file, param_name, 0)
            for param_name in freq_deriv_param_names
        ]

        _, freq_deriv_names = phase_param_name(prev_freq_deriv_order)
        if len(df_grid) != len(freq_deriv_names):
            raise ValueError(
                f"Length of df_grid ({len(df_grid)}) does not match number of frequency derivative names ({len(freq_deriv_names)})"
            )

        df_grid = {name: df_grid[i] for i, name in enumerate(freq_deriv_names)}

        if freq < 200:
            tau = 86400 * 365.25 * 300
        else:
            tau = 86400 * 365.25 * (300 + (freq - 199) * 0.5)

        model = PowerLawModel(nc_min=config["nc_min"], nc_max=config["nc_max"], tau=tau)
        followup_generator = FollowUpParamGenerator(model)

        search_data = followup_generator.generate_parameter(
            alpha=target["alpha"],
            dalpha=target["dalpha"],
            delta=target["delta"],
            ddelta=target["ddelta"],
            data=data,
            old_freq_deriv_order=prev_freq_deriv_order,
            new_freq_deriv_order=freq_deriv_order,
            spacing=df_grid,
            n_spacing=config["followup_n_spacing"],
            sky_radius=sky_radius,
            spacing_alpha=spacing_alpha,
            spacing_delta=spacing_delta,
        )

        # Tile each outlier across the actual Weave sky grid: (d_alpha, d_delta)
        # offsets centered on the grid's mean, read from sky_grid_file. One of
        # these points is already nearly coincident with the outlier's own
        # (0,0) position, so it is not added separately.
        if len(search_data.data) > 0:
            d_alpha, d_delta = np.loadtxt(sky_grid_file, unpack=True)
            n_sky = len(d_alpha)
            n_out = len(search_data.data)
            idx = np.repeat(np.arange(n_out), n_sky)
            sky_idx = np.tile(np.arange(n_sky), n_out)

            table = Table(search_data.data)[idx]
            table["alpha"] = table["alpha"] + d_alpha[sky_idx]
            table["delta"] = table["delta"] + d_delta[sky_idx]
            search_data = fits.BinTableHDU(data=table)

            if injection_data is not None:
                injection_data = np.repeat(injection_data, n_sky)

        taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

        dag_file = manager.make_search_dag(
            taskname,
            freq,
            search_data.data,
            num_top_list=num_top_list,
            stage=stage,
            freq_deriv_order=freq_deriv_order,
            n_seg=n_seg,
            sft_files=sft_files,
            metric_file=metric_file,
            request_memory="8GB",
            request_disk="8GB",
            request_cpu=1,
            use_osg=use_osg,
            use_osdf=use_osdf,
            inj_params=injection_data,
            inj_freq_deriv_order=inj_freq_deriv_order if is_injection else None,
            tasks_per_job=tasks_per_job,
        )

        f_daglist.write(f"{dag_file}\n")

if skipped_freqs:
    print(f"Skipped frequencies: {skipped_freqs}")
