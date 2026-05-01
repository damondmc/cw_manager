import time
import numpy as np
from pathlib import Path
from tqdm import tqdm

from .writer import write_search_subfile, write_search_dagfile
from paws.io import make_dir
from paws.definitions import phase_param_name, task_name, ext_param_name
from paws.filepaths import PathManager

class WorkflowManager:
    """
    Central manager for creating HTCondor DAG and SUB files for all 
    analysis stages (Search, Follow-up, Upper Limit).
    """

    def __init__(self, config, target):
        """
        Initialize the WorkflowManager.

        Parameters:
            config (dict): Configuration dictionary (user, accounting, etc.).
            target (dict): Target object containing astronomical target info.
        """
        self.config = config
        self.target = target
        self.paths = PathManager(config, target)

        # Internal state for search parameter names
        self.freq_param_names = []
        self.freq_deriv_param_names = []
        self.num_top_list = 0

    # =========================================================================
    #  SECTION 1: SEARCH & FOLLOW-UP STAGE
    # =========================================================================

    def _get_execution_kwargs(self, n_seg):
        """Helper to generate common keyword arguments for the search executable."""
        extra_stats = "coh2F_det,mean2F,coh2F_det,mean2F_det"
        kwargs = {
            "semi-max-mismatch": self.config['semi_mm'],
            "toplist-limit": self.num_top_list,
            "extra-statistics": extra_stats
        }
        
        if n_seg != 1:
            kwargs["coh-max-mismatch"] = self.config['coh_mm']
            
        return kwargs

    def _format_injection_str(self, colnames, inj_param):
        """Formats injection parameters into a semicolon-separated string."""
        return ";".join([f"{col}={inj_param[col]}" for col in colnames])

    def _create_wrapper_script(self, stage):
        """Creates a bash script that runs multiple Weave commands sequentially."""
        wrapper_path = self.paths.dag_file(0, "wrapper", stage).parent / f"run_weave_batch_{stage}.sh"
        wrapper_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(wrapper_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write('TASK_FILE=$1\n')
            f.write('EXE=$2\n\n')
            f.write('while IFS= read -r weave_args; do\n')
            f.write('    $EXE $weave_args\n')
            f.write('    if [ $? -ne 0 ]; then\n')
            f.write('        echo "Error: Weave failed on args: $weave_args"\n')
            f.write('        exit 1\n')
            f.write('    fi\n')
            f.write('done < "$TASK_FILE"\n')
            
        return wrapper_path

    def _search_batch_args(self, freq, stage, chunk, taskname, n_seg, sft_files, 
                           node_index, use_osg, metric_file, exe, wrapper_path, tasks_per_job, inj_colnames=None):
        """Generates task files and DAG VARS for a chunk of grouped jobs."""
        
        # 1. Create a task file for this specific Condor node (saved in the DAG directory)
        task_dir = self.paths.dag_file(freq, taskname, stage).parent / "tasks"
        task_dir.mkdir(parents=True, exist_ok=True)
        task_file = task_dir / f"{taskname}_task_{node_index}.txt"
        
        output_files = []
        remap_strings = []
        task_lines = []
        
        kwargs = self._get_execution_kwargs(n_seg)
        
        # 2. Build literal command lines for the wrapper script to run
        for i, job_data in enumerate(chunk):
            job_index = (node_index - 1) * tasks_per_job + i + 1
            
            # Unpack based on injection mode
            if inj_colnames:
                search_param, inj_param = job_data
                inj_str = self._format_injection_str(inj_colnames, inj_param)
            else:
                search_param = job_data
                inj_str = None
                
            result_file = self.paths.weave_output_file(freq, taskname, job_index, stage)
            make_dir([result_file])
            
            # OSG expects local filenames (dynamically pulling the correct extension); Local execution expects absolute paths
            local_out_name = f"{freq}Hz_out.fts.{job_index}" if use_osg else str(result_file)
            sft_names = ";".join([Path(s).name if use_osg else str(s) for s in sft_files])
            metric_name = Path(metric_file).name if use_osg else str(metric_file)
            
            cmd_parts = [
                f"--output-file={local_out_name}",
                f'--sft-files="{sft_names}"',
                f"--setup-file={metric_name}"
            ]
            
            for key, value in kwargs.items():  
                cmd_parts.append(f"--{key}={value}")
            
            cmd_parts.append(f"--alpha={search_param['alpha']}/{search_param['dalpha']}")
            cmd_parts.append(f"--delta={search_param['delta']}/{search_param['ddelta']}")
            
            for key1, key2 in zip(self.freq_param_names, self.freq_deriv_param_names): 
                cmd_parts.append(f"--{key1}={search_param[key1]}/{search_param[key2]}")
                
            if inj_str:
                cmd_parts.append(f'--injections={{{inj_str}}}')
                
            task_lines.append(" ".join(cmd_parts))
            
            if use_osg:
                output_files.append(local_out_name)
                remap_strings.append(f"{local_out_name}={result_file}")
            
        # Write task file for the worker node
        with open(task_file, "w") as f:
            f.write("\n".join(task_lines) + "\n")
            
        # 3. Construct DAG VARS mapping for this node
        cmd_args = f"{task_file.name} {exe}" if use_osg else f"{task_file} {exe}"
        args_list = [f'CMD_ARGS="{cmd_args}"']

        if use_osg:
            args_list.append(f'OUTPUT_FILES="{", ".join(output_files)}"')
            args_list.append(f'REMAP_OUTPUT_FILES="{";".join(remap_strings)}"')
            
            transfer_files = [str(s) for s in sft_files] + [str(metric_file), str(task_file), str(wrapper_path)]
            args_list.append(f'TRANSFER_FILES="{", ".join(transfer_files)}"')
            
        return " ".join(args_list) + " "

    def make_search_dag(self, taskname, freq, params, num_top_list, stage, freq_deriv_order, n_seg,
                        sft_files, metric_file, request_memory='18GB', request_disk='5GB', request_cpu=1, 
                        use_osg=False, use_osdf=False, exe=None, image=None,
                        inj_params=None, inj_freq_deriv_order=None, tasks_per_job=10):
        """
        Creates the DAG and SUB files for the Search/Follow-up stage using grouped batching.
        """
        t0 = time.time()
        is_injection = inj_params is not None and inj_freq_deriv_order is not None
        print(f"Generating DAG for {taskname} (Mode: {stage}, Injections: {is_injection}, Batch Size: {tasks_per_job})...")

        if use_osdf and not use_osg:
            print('Warning: SFTs from OSDF requested but not using OSG resources.')

        self.freq_param_names, self.freq_deriv_param_names = phase_param_name(freq_deriv_order)
        self.num_top_list = num_top_list

        inj_colnames = None
        if is_injection:
            inj_freq_names, _ = phase_param_name(inj_freq_deriv_order)
            # Combine external params + "Freq" + remaining freq derivatives (f1dot, etc.)
            inj_colnames = ext_param_name() + ["Freq"] + inj_freq_names[1:]
            
        dag_file_path = self.paths.dag_file(freq, taskname, stage)
        dag_file_path.parent.mkdir(parents=True, exist_ok=True)
        dag_file_path.unlink(missing_ok=True)
        
        cr_files = self.paths.condor_record_files(freq, taskname, stage)
        make_dir(cr_files)

        wrapper_path = self._create_wrapper_script(stage)
        exe = exe if exe else self.paths.weave_executable
        
        sub_file_path = self.paths.condor_sub_file(freq, taskname, stage)
        sub_file_path.parent.mkdir(parents=True, exist_ok=True)
        sub_file_path.unlink(missing_ok=True)
        
        # Executable is /bin/bash; arguments point to the wrapper script and dynamic job tasks
        write_search_subfile(
            filename=str(sub_file_path), executable_path="/bin/bash", transfer_executable=False, 
            output_path=str(cr_files[0]), error_path=str(cr_files[1]), log_path=str(cr_files[2]), 
            arg_list_string=f"{wrapper_path.name if use_osg else wrapper_path} $(CMD_ARGS)", 
            accounting_group=self.config['acc_group'], user=self.config['user'],
            request_memory=request_memory, request_disk=request_disk, request_cpu=request_cpu,
            use_osg=use_osg, use_osdf=use_osdf, image=image
        )   
        
        # Combine jobs into chunks
        job_list = list(zip(params, inj_params) if is_injection else params)
        chunks = [job_list[i:i + tasks_per_job] for i in range(0, len(job_list), tasks_per_job)]
        
        for node_index, chunk in tqdm(enumerate(chunks, 1), total=len(chunks)):
            arg_list = self._search_batch_args(
                freq, stage, chunk, taskname, n_seg, sft_files, node_index, 
                use_osg, metric_file, exe, wrapper_path, tasks_per_job, inj_colnames
            )
            write_search_dagfile(str(dag_file_path), taskname, str(sub_file_path), node_index, arg_list)

        elapsed = time.time() - t0
        print(f'Finished writing {stage} dag files. Time: {elapsed:.2f}s')
        return dag_file_path

    # =========================================================================
    #  SECTION 2: UPPER LIMIT STAGE
    # =========================================================================

    def _ul_args(self, config_file, target_file, taskname, 
                 freq, stage, freq_deriv_order, df_grid, inj_freq_deriv_order,
                 sft_files_local, metric_file_local, n_inj, num_toplist, h0_est, 
                 mean2f_th, non_sat_bands, sky_uncertainty, 
                 cluster, work_in_local_dir, save_intermediate, request_cpu):
        """Generates command line arguments for the python upper limit script."""
 
        bands_str = ' '.join(map(str, non_sat_bands))
        df_grid_str = ' '.join(map(str, df_grid))
        
        arg_list_string = (
            f'--config_file {config_file} --target_file {target_file} --taskname {taskname} '
            f'--freq {freq} --stage {stage} --freq_deriv_order {freq_deriv_order} '
            f'--df_grid {df_grid_str} --inj_freq_deriv_order {inj_freq_deriv_order} '
            f'--sft_files {sft_files_local} --metric_file {metric_file_local} '
            f'--n_inj {n_inj} --num_toplist {num_toplist} --h0_est {h0_est} '
            f'--mean2f_th {mean2f_th} --non_sat_bands {bands_str} '
            f'--sky_uncertainty {sky_uncertainty} --n_cpus {request_cpu}'
        )

        if cluster:
            arg_list_string += " --cluster"
        if work_in_local_dir:
            arg_list_string += " --work_in_local_dir"
        if save_intermediate:
            arg_list_string += " --save_intermediate"
            
        return arg_list_string

    def _ul_transfer_args(self, config_file, target_file, taskname, freq, metric_file, stage, sft_files, cluster, exe, image):
        """Generates VARS for OSG file transfers for Upper Limits."""
        
        input_files_list = [str(exe), str(image), str(config_file), str(target_file)]
        input_files_list.extend([str(s) for s in sft_files])
        input_files_list.append(str(metric_file))
        
        input_files_str = ", ".join(input_files_list)

        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=cluster)
        make_dir([outlier_file_path])

        arg_list = (
            f'OUTPUT_FILES="{Path(outlier_file_path).name}" '
            f'REMAP_OUTPUT_FILES="{Path(outlier_file_path).name}={outlier_file_path}" '
            f'TRANSFERFILES="{input_files_str}" '
        )
        return arg_list

    def make_upperlimit_dag(self, config_file, target_file, taskname, freq, stage, freq_deriv_order,
                            sft_files, metric_file, mean2f_th, non_sat_bands, exe,
                            df_grid=[1e-6, 1e-13, 1e-20], inj_freq_deriv_order=4, num_toplist=1,
                            sky_uncertainty=1e-5, h0_est=6e-26, n_inj=64, 
                            request_memory='4GB', request_disk='4GB', request_cpu=32,
                            cluster=False, work_in_local_dir=False, save_intermediate=False,
                            image=None):
        """
        Creates the DAG and SUB files for the Upper Limit stage.
        """
        print(f"Generating UPPERLIMIT DAG for {taskname}...")

        dag_file_path = self.paths.dag_file(freq, taskname, stage)
        dag_file_path.parent.mkdir(parents=True, exist_ok=True)
        dag_file_path.unlink(missing_ok=True)
        
        cr_files = self.paths.condor_record_files(freq, taskname, stage)
        make_dir(cr_files)

        sub_file_path = self.paths.condor_sub_file(freq, taskname, stage)
        sub_file_path.unlink(missing_ok=True)
 
        sft_files_local = ';'.join([Path(s).name for s in sft_files])
        metric_file_local = Path(metric_file).name
        
        python_args = self._ul_args(
            config_file=Path(config_file).name, target_file=Path(target_file).name, 
            taskname=taskname, freq=freq, stage=stage, freq_deriv_order=freq_deriv_order, 
            df_grid=df_grid, inj_freq_deriv_order=inj_freq_deriv_order,
            sft_files_local=sft_files_local, metric_file_local=metric_file_local, 
            n_inj=n_inj, num_toplist=num_toplist, h0_est=h0_est, 
            mean2f_th=mean2f_th, non_sat_bands=non_sat_bands, sky_uncertainty=sky_uncertainty, 
            cluster=cluster, work_in_local_dir=work_in_local_dir, 
            save_intermediate=save_intermediate, request_cpu=request_cpu
        )

        full_arg_string = f"{Path(exe).name} {python_args}"

        write_search_subfile(
            filename=str(sub_file_path), executable_path="/opt/conda/bin/python3", transfer_executable=False,
            output_path=str(cr_files[0]), error_path=str(cr_files[1]), log_path=str(cr_files[2]),
            arg_list_string=full_arg_string, accounting_group=self.config['acc_group'], user=self.config['user'], 
            request_memory=request_memory, request_disk=request_disk, request_cpu=request_cpu,
            use_osg=True, use_osdf=True, image=Path(image).name,
        )

        arg_list = self._ul_transfer_args(
            config_file, target_file, taskname, freq, metric_file, stage, sft_files, cluster, exe, image
        )
        
        write_search_dagfile(str(dag_file_path), taskname, str(sub_file_path), 1, arg_list)
        return dag_file_path
