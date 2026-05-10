import yaml
import os
from paws.params.search import SearchParamGenerator
from paws.params.models import UniformModel
from paws.workflow.manager import WorkflowManager
from paws.filepaths import PathManager
from paws.definitions import task_name
from pyfstat.utils import get_sft_as_arrays

from tools import generate_hexagonal_grid

# ==============================================================================
# 1. Load Configs & SFT Timestamps (Done once to save time)
# ==============================================================================
with open("/home/hoitim.cheung/o3_validate/config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

with open("/home/hoitim.cheung/o3_validate/config/g347.yaml", "r") as f:
    target = yaml.safe_load(f)

paths = PathManager(config, target)
manager = WorkflowManager(config, target)

# Extract timestamps for the exact data you are using
timestamps_dict = {}
try:
    sft_h1 = get_sft_as_arrays(
        "/osdf/igwn/cit/staging/hoitim.cheung/SFTs/o3_data/SFTs/H1/31/H-11025_H1_1800SFT_NBF0030Hz900W0002Hz0_O3_Gated_Sub60Hz-1238166483-31197072.sft"
    )
    timestamps_dict["H1"] = sft_h1[1]["H1"]
    sft_l1 = get_sft_as_arrays(
        "/osdf/igwn/cit/staging/hoitim.cheung/SFTs/o3_data/SFTs/L1/31/L-10133_L1_1800SFT_O3GatedSub60Hz_NBF0030Hz900W0002Hz0-1238184901-31178322.sft"
    )
    timestamps_dict["L1"] = sft_l1[1]["L1"]
except Exception as e:
    print(f"Warning: Could not load SFTs. Error: {e}")

# ==============================================================================
# 2. Search Setup & Invariant Math
# ==============================================================================

observation_time = "o3ab"
coh_day = 360
freq_deriv_order = 2
n_seg = 1
num_top_list = 200000
stage = f"narrow_{observation_time}_{coh_day}_sky_shift"
metric_file = "osdf:///igwn/cit/staging/hoitim.cheung/metricSetup/o3ab_t360_sky_s2.fts"

ref_time = 1246197626.5

span_df0 = 5e-8
span_df1 = 1e-14
span_df2 = 1.5e-21

# The base step sizes you defined
step_alpha = 6e-6
step_delta = 3e-5

os.makedirs("./dagFiles", exist_ok=True)
dag_list_file = "./dagFiles/o3_dag.txt"

# Clear the file if it already exists from a previous run
with open(dag_list_file, "w") as f:
    f.write("")

num_layers = 12

# Fetch the grouped grid from your tool
layer_points = generate_hexagonal_grid(num_layers=num_layers, filter_threshold=6)

# Flatten the dictionary into a single list to iterate over
all_grid_points = [pt for points in layer_points.values() for pt in points]

# Loop through the generated grid
for q, r, i, j in all_grid_points:
    print("\n=======================================================")
    print(f"Generating DAG for hex position: dq={q}, dr={r}")
    print("=======================================================")

    # Calculate the actual shift for this specific grid point
    current_shift_alpha = i * step_alpha
    current_shift_delta = j * step_delta

    new_alpha = target["alpha"] + current_shift_alpha
    new_delta = target["delta"] + current_shift_delta

    # Apply the widths around the newly computed center
    fmin = target["f0"] - span_df0
    fmax = target["f0"] + span_df0

    f1_lim = (target["f1dot"] - span_df1, target["f1dot"] + span_df1)
    f2_lim = (target["f2dot"] - span_df2, target["f2dot"] + span_df2)

    # Generate Parameters
    model = UniformModel(f1_lim, f2_lim)
    generator = SearchParamGenerator(model, freq_deriv_order)

    params = generator.generate_parameters(
        new_alpha,
        target["dalpha"],
        new_delta,
        target["ddelta"],
        fmin,
        fmax,
        df=1e-7,
        df1=1e-12,
        df2=5e-21,
    )

    freq = int(fmin)
    sftFiles = paths.sft_ensemble(freq, use_osdf=True)

    # Generate a UNIQUE task name so the 8 DAGs don't overwrite each other
    # Added grid coordinates _dx{i}_dy{j} to the filename
    current_taskname = (
        task_name(target["name"], stage, coh_day, freq_deriv_order)
        + f"_{freq}Hz_dq{q}_dr{r}"
    )

    dag_file = manager.make_search_dag(
        current_taskname,
        freq,
        params[freq].data,
        num_top_list=num_top_list,
        stage=stage,
        freq_deriv_order=freq_deriv_order,
        n_seg=n_seg,
        sft_files=sftFiles,
        metric_file=metric_file,
        request_memory="2GB",
        request_disk="3GB",
        request_cpu=1,
        use_osg=True,
        use_osdf=True,
        tasks_per_job=1,
    )

    print(f"Successfully created DAG: {dag_file}")

    with open(dag_list_file, "a") as f:
        f.write(f"{dag_file}\n")
