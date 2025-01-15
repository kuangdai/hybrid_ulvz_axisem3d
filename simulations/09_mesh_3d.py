import json
import os
import sys
from pathlib import Path

import numpy as np

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_mesh_3d')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'cp templates/mesh/inparam.* {input_dir}/')

    # inparam.advanced
    replace_in_file(input_dir / 'inparam.advanced',
                    {'__DD_REGIONS_ENHANCE__': 'none'})

    # inparam.nu
    args_mesh = json.load(open(f'outputs/{run_name}/@@_exodus/args.json'))
    replace_in_file(input_dir / 'inparam.nu',
                    {'__NU__': args_mesh['nu_to_use']})

    # inparam.time_src_recv
    if args['array']['wave_extrapolation']:
        replace_in_file(input_dir / 'inparam.time_src_recv',
                        {'__DT__': args['time_series']['dt'],
                         '__LENGTH__': args['time_series']['length'],
                         '__SAMPLE__': args['time_series']['sample_interval'],
                         '__STATIONS__': 'STATIONS_ARRAY_EXTRAPOLATION',
                     '__HDUR__': args['mesh']['period']
                     if args['time_series']['use_period_for_stf'] else "0.0"})
    else:
        replace_in_file(input_dir / 'inparam.time_src_recv',
                        {'__DT__': args['time_series']['dt'],
                         '__LENGTH__': args['time_series']['length'],
                         '__SAMPLE__': args['time_series']['sample_interval'],
                         '__STATIONS__': 'STATIONS_ARRAY',
                     '__HDUR__': args['mesh']['period']
                     if args['time_series']['use_period_for_stf'] else "0.0"})

    # ULVZ.txt
    if args['ulvz']['nc_3d']:
        os.system(f'cp templates/mesh/ULVZ_false.txt {input_dir}/ULVZ.txt')
        os.system(f'cp {input_dir}/inparam.nc_3d_model {input_dir}/inparam.model')
        os.system(f'cp inputs/{run_name}/ulvz.nc {input_dir}/')
    else:
        layers = args['ulvz']['layers']
        replace_in_file('templates/mesh/ULVZ_true.txt',
                        {'__N_LAYER__': len(layers),
                         '__LAYERS__': '\n'.join(layers),
                         '__TMIN__': 0,
                         '__TMAX__': np.pi / args_mesh['NEX'] * args_mesh['NEX_U']},
                        dest=input_dir / 'ULVZ.txt')

    # CMTSOLUTION
    args_stations = json.load(open(f'outputs/{run_name}/@@_stations/args.json'))
    replace_in_file('templates/mesh/CMTSOLUTION_ULVZ',
                    {'__ULVZ_LAT__': args_stations['ulvz_lat'],
                     '__ULVZ_LON__': args_stations['ulvz_lon']},
                    dest=input_dir / 'CMTSOLUTION')
