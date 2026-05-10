import numpy as np
import yaml
from astropy.io import fits
from tqdm import tqdm

from paws.filepaths import PathManager
from paws.workflow.manager import WorkflowManager
from paws.definitions import phase_param_name
from paws.params.models import PowerLawModel
from paws.params.injections import InjectionParamGenerator

# 1. Load Configs
with open("/home/hoitim.cheung/galacticCenter/config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

home_dir = config["home_dir"]

with open(f"{home_dir}config/gal.yaml", "r") as f:
    target = yaml.safe_load(f)

paths = PathManager(config, target)
manager = WorkflowManager(config, target)

# injection-2 308-318 dag was accidentally generated with new ver of paws, but results were done with old vers so keep it.
fmin = 20
fmax = 400
stage = "injections-0"
tcoh = 5
freq_deriv_order = 2
inj_freq_deriv_order = 4
sky_uncertainty = 1e-4
ref_time = 1395617439
n_seg = 107
n_inj = 200

extra_stats = "coh2F_det,mean2F,coh2F_det,mean2F_det"

weave_exe = config["executables"]["weave"]


num_top_list = config["num_toplist"]
metric_file = "osdf:///igwn/cit/staging/hoitim.cheung/metricSetup/Start1368970000_TCoh432000_N107_Spin2.fts"
use_osg = True
use_osdf = True

h0_file = f"{home_dir}postprocess/data/{target['name']}_upperlimit-v2_TCoh{tcoh}_O2_{fmin}-{fmax}Hz_clustered.txt"
freq, h0_arr, _ = np.loadtxt(h0_file).T

_, freq_deriv_param_names = phase_param_name(freq_deriv_order)

skipped_freqs = []

dag_list_path = f"{home_dir}dagFiles/{stage}_{target['name']}_dag{fmin}-{fmax}Hz.txt"

with open(dag_list_path, "w") as f_daglist:
    for freq, h0 in tqdm(
        zip(range(fmin, fmax), h0_arr), desc="Generating DAGs", total=fmax - fmin
    ):
        sft_files = []
        files = paths.sft_ensemble(freq, use_osdf=True)
        sft_files.extend(files)  # Use extend, not append, to flatten the list

        data_taskname = (
            f"{target['name']}_search-0_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"
        )

        data_file = paths.outlier_file(freq, data_taskname, "search-0", cluster=True)

        df_grid = [
            fits.getval(data_file, param_name, 0)
            for param_name in freq_deriv_param_names
        ]

        _, freq_deriv_names = phase_param_name(freq_deriv_order)
        if len(df_grid) != len(freq_deriv_names):
            raise ValueError(
                f"Length of df_grid ({len(df_grid)}) does not match number of frequency derivative names ({len(freq_deriv_names)})"
            )

        df_grid = {name: df_grid[i] for i, name in enumerate(freq_deriv_names)}

        mean2f_th = fits.getval(data_file, "mean2F_th", 0)

        non_sat_bands = fits.getdata(data_file, extname="non_sat_band")["non_sat_band"]

        if len(non_sat_bands) == 0:
            print(
                f"Warning: No non-saturated bands found for {freq} Hz. Skipping DAG generation."
            )
            skipped_freqs.append(freq)
            continue

        if freq < 200:
            tau = 86400 * 365.25 * 300
        else:
            tau = 86400 * 365.25 * (300 + (freq - 199) * 0.5)

        model = PowerLawModel(nc_min=config["nc_min"], nc_max=config["nc_max"], tau=tau)
        injection_generator = InjectionParamGenerator(
            model=model, ref_time=ref_time, f0_band=config["f0_band"]
        )

        search_data, injection_data = injection_generator.generate_parameters(
            alpha=target["alpha"],
            dalpha=target["dalpha"],
            delta=target["delta"],
            ddelta=target["ddelta"],
            non_sat_bands=non_sat_bands,
            spacing=df_grid,
            h0=h0,
            freq=freq,
            n_inj=n_inj,
            n_spacing=config["followup_n_spacing"],
            inj_freq_deriv_order=inj_freq_deriv_order,
            freq_deriv_order=freq_deriv_order,
            sky_uncertainty=sky_uncertainty,
        )

        # --- Make DAG ---
        taskname = f"{target['name']}_{stage}_TCoh{tcoh}_O{freq_deriv_order}_{freq}Hz"

        dag_file = manager.make_search_dag(
            taskname,
            freq,
            search_data[str(freq)].data,
            num_top_list=num_top_list,
            stage=stage,
            freq_deriv_order=freq_deriv_order,
            n_seg=n_seg,
            sft_files=sft_files,
            metric_file=metric_file,
            request_memory="8GB",
            request_disk="4GB",
            request_cpu=1,
            use_osg=use_osg,
            use_osdf=use_osdf,
            inj_params=injection_data[str(freq)].data,
            inj_freq_deriv_order=inj_freq_deriv_order,
            tasks_per_job=200,
        )

        f_daglist.write(f"{dag_file}\n")

if skipped_freqs:
    print(f"Skipped frequencies due to no non-saturated bands: {skipped_freqs}")
