import json
import os
import sys
from pathlib import Path

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_mesh_2d')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'cp templates/mesh/inparam.* {input_dir}/')

    # inparam.advanced
    replace_in_file(input_dir / 'inparam.advanced',
                    {'__DD_REGIONS_ENHANCE__': 'none'})

    # inparam.nu
    replace_in_file(input_dir / 'inparam.nu',
                    {'__NU__': 0 if args['event']['monopole'] else 2})

    # inparam.time_src_recv
    replace_in_file(input_dir / 'inparam.time_src_recv',
                    {'__DT__': args['time_series']['dt'],
                     '__LENGTH__': args['time_series']['length'],
                     '__SAMPLE__': args['time_series']['sample_interval'],
                     '__STATIONS__': 'STATIONS_ARRAY',
                     '__HDUR__': args['mesh']['period']
                     if args['time_series']['use_period_for_stf'] else "0.0"})

    # ULVZ.txt
    args_stations = json.load(open(f'outputs/{run_name}/@@_stations/args.json'))
    layers = args['ulvz']['layers']
    replace_in_file('templates/mesh/ULVZ_true.txt',
                    {'__N_LAYER__': len(layers),
                     '__LAYERS__': '\n'.join(layers),
                     '__TMIN__': args_stations['ulvz_dist_near'],
                     '__TMAX__': args_stations['ulvz_dist_far']},
                    dest=input_dir / 'ULVZ.txt')

    # CMTSOLUTION
    os.system(f'cp inputs/{run_name}/CMTSOLUTION {input_dir}/')
