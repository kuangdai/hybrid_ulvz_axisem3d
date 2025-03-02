import json
import os
import sys
from pathlib import Path

import numpy as np

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_solve_3d')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'cp templates/solve/*.e {input_dir}/')
    os.system(f'cp templates/solve/*.yaml {input_dir}/')
    os.system(f'cp templates/solve/inparam.stream_wj {input_dir}/')
    if args["wave_extrapolation"]["save_outgoing"]:
        os.system(f'cp templates/solve/inparam.stream_3d_extrapolation {input_dir}/inparam.stream')
    elif args['ulvz_animation']['enabled']:
        os.system(f'cp templates/solve/inparam.stream_3d_animation {input_dir}/inparam.stream')
    else:
        os.system(f'cp templates/solve/inparam.stream_3d {input_dir}/inparam.stream')

    # inparam.stream
    replace_in_file(input_dir / 'inparam.stream',
                    {'__SAMPLE__': args['time_series']['sample_interval']})
    if args['ulvz_animation']['enabled']:
        replace_in_file(input_dir / 'inparam.stream',
                        {'__SAMPLE_ANIMATION__': args['ulvz_animation']['sample_interval']})

    # inparam.stream_wj
    quad_bound = np.loadtxt(f'outputs/{run_name}/@@_exodus/quad_bound.txt', dtype=int)
    args_mesh = json.load(open(f'outputs/{run_name}/@@_exodus/args.json'))
    replace_in_file(input_dir / 'inparam.stream_wj',
                    {'__NUM_IN_ELEM__': len(quad_bound),
                     '__IN_ELEM_TAGS__': ' '.join(quad_bound.astype(str)),
                     '__IN_ELEM_R0__': args_mesh['rb_min'],
                     '__IN_ELEM_R1__': args_mesh['rb_max'],
                     '__IN_ELEM_T1__': args_mesh['tb_max'],
                     '__INFO_ONLY__': 0})

    # job
    nodes = args['slurm']['3d']['nodes']
    nproc = args['slurm']['3d']['nproc']
    wct = args['slurm']['3d']['wct_solve']
    assert nproc % nodes == 0
    with open(f"../axisem3d/hpc/{args['slurm']['hpc']}/header_mpi.sh", 'r') as fs:
        header = fs.read()
    replace_in_file('templates/submit_solve_3d.sh',
                    {'__SBATCH__': header,
                     '__NODES__': nodes,
                     '__N_PER_NODE__': nproc // nodes,
                     '__N__': nproc,
                     '__T__': wct,
                     '__J__': 'solve_3d',
                     '__HPC__': args['slurm']['hpc']},
                    dest=out_dir / 'submit_solve_3d.sh')

    # truncate incident
    truncate_incident = args.get('time_series', {}).get('truncate_incident', False)
    if truncate_incident:
        # copy python
        stf_dir = Path(f'outputs/{run_name}/@@_stf_3d')
        replace_in_file('templates/truncate_taper_incident.py', {},
                        dest=stf_dir / 'truncate_taper_incident.py')
        # add to stf script
        t0 = args['time_series']['truncate_incident_t0']
        t1 = args['time_series']['truncate_incident_t1']
        with open(stf_dir / "submit_stf_3d.sh", 'r') as fs:
            content = fs.read()
        with open(stf_dir / "submit_stf_3d.sh", 'w') as fs:
            cmd = (f"\n# Truncate and taper incident"
                   f"\npython truncate_taper_incident.py --t0 {t0} --t1 {t1}")
            fs.write(content + cmd)
        # change inparams
        replace_in_file(out_dir / 'input/inparam.stream_wj',
                        {
                            "injection_SOLID.nc": "injection_SOLID_truncated.nc",
                            "injection_FLUID.nc": "injection_FLUID_truncated.nc"})
        replace_in_file(out_dir / 'input/inparam.stream',
                        {"-9999999": str(t0)})
