import json
import os
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from obspy.core.trace import Stats, Trace, UTCDateTime
from scipy.signal import fftconvolve
from tqdm import tqdm


def read_process_sac(nc, T, sac_dir=None, channel=2, aggressive_conv_filter_factor=0.8,
                     sampling_rate=None):
    seis = {}
    with Dataset(nc, 'r') as nc:
        # header
        times = nc.variables['time_points'][:2]
        stats = Stats()
        stats.npts = len(times)
        stats.delta = times[1] - times[0]
        stats.starttime = UTCDateTime(times[0])

        # source-time function
        half_duration = T * aggressive_conv_filter_factor
        decay = 1.628
        n_stf = int(np.ceil(2.5 * half_duration / stats.delta))
        t_stf = np.arange(-n_stf, n_stf + 1) * stats.delta
        stf = np.exp(-np.power((decay / half_duration * t_stf), 2.)) * decay / (
                half_duration * np.sqrt(np.pi))
        stf /= stf.sum()

        # filter window
        tmin = T * aggressive_conv_filter_factor
        tmax = 200.

        # read, conv, filter, resample, trim
        for st_key in tqdm(nc.variables.keys(), total=len(nc.variables.keys())):
            if st_key == 'time_points':
                continue
            # read
            data = nc.variables[st_key][:, channel]
            # convolve with stf
            data = fftconvolve(data, stf, mode='same')
            # trace
            trace = Trace(data, header=stats)
            # filter
            trace = trace.filter('bandpass', freqmin=1 / tmax, freqmax=1 / tmin)
            # resample and trim
            if sampling_rate is None:
                sampling_rate = stats.sampling_rate
            trace.interpolate(sampling_rate, starttime=0)
            seis[st_key] = trace

        # save to sac
        if sac_dir is not None:
            Path(sac_dir).mkdir(exist_ok=True, parents=True)
            for key, trace in seis.items():
                trace.write(f'{sac_dir}/{key}.sac', format='SAC')
    return seis


if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'../simulations/inputs/{run_name}/args.json'))
    period = args['mesh']['period']
    out_dir = Path(f'outputs/sac/{run_name}')
    for case in ['prem', '1d', '2d']:
        read_process_sac(f'outputs/nc/{run_name}/{case}.nc', period,
                         sac_dir=out_dir / f'wave_{case}')
    # event info
    info_dir = Path(out_dir / 'info')
    info_dir.mkdir(exist_ok=True, parents=True)
    e_data = np.loadtxt(f'../simulations/inputs/{run_name}/CMTSOLUTION',
                        skiprows=1, dtype=str, delimiter=':')
    e_lat_lon_dep = e_data[3:6, 1]
    np.savetxt(info_dir / 'event.txt', e_lat_lon_dep, fmt='%s',
               header='latitude longitude and depth of event')

    # array grid
    if args['array']['use_grid']:
        st_dist_crds = np.arange(args['array']['dist0_from_event'],
                                 args['array']['dist1_from_event'] + 1e-4,
                                 args['array']['delta_dist'])
        st_azim_crds = np.arange(args['array']['azim0_from_ulvz'],
                                 args['array']['azim1_from_ulvz'] + 1e-4,
                                 args['array']['delta_azim'])
        grid_dir = Path(info_dir / 'array_grid')
        grid_dir.mkdir(exist_ok=True, parents=True)
        np.savetxt(grid_dir / 'dist.txt', st_dist_crds)
        np.savetxt(grid_dir / 'azim.txt', st_azim_crds)

        # copy python
        os.system(f'cp plot_dist_azim.py {out_dir}')
