import json
import os
import sys
from pathlib import Path

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_solve_2d')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'cp templates/solve/*.e {input_dir}/')
    os.system(f'cp templates/solve/*.yaml {input_dir}/')
    os.system(f'cp templates/solve/inparam.stream_2d {input_dir}/inparam.stream')

    # inparam.stream
    replace_in_file(input_dir / 'inparam.stream',
                    {'__SAMPLE__': args['time_series']['sample_interval'],
                     '__MESH_DIR__': '@@_mesh_2d'})

    # job
    nodes = args['slurm']['2d']['nodes']
    nproc = args['slurm']['2d']['nproc']
    wct = args['slurm']['2d']['wct']
    assert nproc % nodes == 0
    with open(f"../axisem3d/hpc/{args['slurm']['hpc']}/header_mpi.sh", 'r') as fs:
        header = fs.read()
    replace_in_file('templates/submit_mesh_solve.sh',
                    {'__SBATCH__': header,
                     '__NODES__': nodes,
                     '__N_PER_NODE__': nproc // nodes,
                     '__N__': nproc,
                     '__T__': wct,
                     '__J__': 'mesh_solve_2d',
                     '__HPC__': args['slurm']['hpc'],
                     '__MESH_DIR__': '@@_mesh_2d',
                     '__SOLVE_DIR__': '@@_solve_2d'},
                    dest=out_dir / 'submit_mesh_solve.sh')
