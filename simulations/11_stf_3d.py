import json
import sys
from pathlib import Path

from utils import replace_in_file


def is_lucky_number(n, forceOdd):
    num = n
    if forceOdd and num % 2 == 0:
        return False
    if not forceOdd and num % 2 != 0 and num > 10:
        return False
    for i in range(2, num + 1):
        while num % i == 0:
            num //= i
            if i > 13:
                return False
    num = n
    e = 0
    while num % 11 == 0:
        num //= 11
        e += 1
    num = n
    f = 0
    while num % 13 == 0:
        num //= 13
        f += 1
    if e + f > 1:
        return False
    return True


def next_lucky_number(n, forceOdd=False):
    while True:
        if is_lucky_number(n, forceOdd):
            return n
        n += 1


if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_stf_3d')
    out_dir.mkdir(parents=True, exist_ok=True)

    # input folder
    input_dir = out_dir / 'input'
    input_dir.mkdir(parents=True, exist_ok=True)

    # job
    with open(f"../axisem3d/hpc/{args['slurm']['hpc']}/header_omp.sh", 'r') as fs:
        header = fs.read()
    args_mesh = json.load(open(f'outputs/{run_name}/@@_exodus/args.json'))
    args_stations = json.load(open(f'outputs/{run_name}/@@_stations/args.json'))
    nu_1d = 0 if args['event']['monopole'] else 2
    replace_in_file('templates/submit_stf_3d.sh',
                    {'__SBATCH__': header,
                     '__N__': args['slurm']['3d']['nproc_stf'],
                     '__T__': args['slurm']['3d']['wct_stf'],
                     '__J__': 'stf_3d',
                     '__HPC__': args['slurm']['hpc'],
                     '__E_LAT__': args_stations['event_lat'],
                     '__E_LON__': args_stations['event_lon'],
                     '__U_LAT__': args_stations['ulvz_lat'],
                     '__U_LON__': args_stations['ulvz_lon'],
                     '__NR_EVT__': nu_1d * 2 + 1,
                     '__NR_ULVZ__': next_lucky_number(args_mesh['nu_to_use'] * 2 + 1)},
                    dest=out_dir / 'submit_stf_3d.sh')
