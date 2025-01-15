import json
import os
import sys
from pathlib import Path

import numpy as np

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_mesh_prem')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)
    os.system(f'cp templates/mesh/inparam.* {input_dir}/')

    # inparam.advanced
    box_bot = (3480. - args['box']['height']) * 1e3
    box_top = (3480. + args['ulvz']['height'] + args['box']['height']) * 1e3
    grid_dist = np.loadtxt(f'outputs/{run_name}/@@_stations/grid_dist.txt')
    box_near = grid_dist[0]
    box_far = grid_dist[-1]
    replace_in_file(input_dir / 'inparam.advanced',
                    {'__DD_REGIONS_ENHANCE__': f"rt${box_bot}${box_top}${box_near}${box_far}$1000"})

    # inparam.nu
    replace_in_file(input_dir / 'inparam.nu',
                    {'__NU__': 0 if args['event']['monopole'] else 2})

    # inparam.time_src_recv
    replace_in_file(input_dir / 'inparam.time_src_recv',
                    {'__DT__': args['time_series']['dt'],
                     '__LENGTH__': args['time_series']['length'] + 1.,
                     '__SAMPLE__': args['time_series']['sample_interval'],
                     '__STATIONS__': 'STATIONS_ARRAY_INCIDENT',
                     '__HDUR__': args['mesh']['period']
                     if args['time_series']['use_period_for_stf'] else "0.0"})

    # ULVZ.txt
    os.system(f'cp templates/mesh/ULVZ_false.txt {input_dir}/ULVZ.txt')

    # CMTSOLUTION
    os.system(f'cp inputs/{run_name}/CMTSOLUTION {input_dir}/')
