import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from obspy.core.trace import UTCDateTime
from obspy.taup import TauPyModel

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot options.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simulation', type=str, required=True,
                        help='simulation name')
    parser.add_argument('-r', '--reference', type=str, default=None,
                        help="reference simulation name; "
                             "if given, plot scattered wave (simulation-reference); "
                             "otherwise, plot total wave")
    parser.add_argument('-d', '--distance-index', type=int, default=None,
                        help="index of distance; "
                             "if given, plot azimuthal record section at this distance")
    parser.add_argument('-a', '--azimuth-index', type=int, default=0,
                        help="index of azimuth; "
                             "if given, plot distance record section at this azimuth")
    parser.add_argument('-R', '--station-range', type=float, nargs=2, default=None,
                        help="range of distance or azimuth")
    parser.add_argument('-f', '--sac-filename-formatter', type=str, default='ARRAY.D%d_A%d.sac',
                        help='sac filename formatter using distance index and azimuth index')
    parser.add_argument('-M', '--amplitude-normalize-mode',
                        choices=['all', 'individual'], default='individual',
                        help="mode for amplitude normalization")
    parser.add_argument('-A', '--amplitude-zoom-factor', type=float, default=-4.0,
                        help="factor for amplitude zooming")
    parser.add_argument('-p', '--aligned-phase', type=str, default='ScP',
                        help="phase used as time origin")
    parser.add_argument('-t', '--travel-time-phases', type=str, nargs='+', default=['S'],
                        help="phases for plotting travel time curves")
    parser.add_argument('-w', '--window', type=float, nargs=2, default=[-10.0, 30.0],
                        help="window length before/after time origin")
    parser.add_argument('-F', '--figure-size', type=float, default=None,
                        help='figure size')
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        help='output filename')
    args = parser.parse_args()

    # taup
    needed_phases = list(np.unique([args.aligned_phase] + args.travel_time_phases))
    taup_model = TauPyModel(model='prem')

    # event
    event_depth = np.loadtxt('info/event.txt', skiprows=1)[-1]

    # read grid
    dist_grid = np.loadtxt('info/array_grid/dist.txt')
    azim_grid = np.loadtxt('info/array_grid/azim.txt')

    # set y-axis
    if args.station_range is None:
        args.station_range = [-10000, 10000]
    if args.distance_index is not None:
        filenames = []
        y_ticks = []
        for i_azim, azim in enumerate(azim_grid):
            if args.station_range[0] < azim < args.station_range[1]:
                filenames.append(
                    args.sac_filename_formatter % (args.distance_index, i_azim))
                y_ticks.append(azim)
        y_label = 'Azimuth (deg)'
        y_title = f'Distance={np.round(dist_grid[args.distance_index], decimals=3)} deg'
        # travel time
        tt = {ph: np.full_like(y_ticks, np.nan) for ph in needed_phases}
        arrivals = taup_model.get_travel_times(source_depth_in_km=event_depth,
                                               distance_in_degree=dist_grid[args.distance_index],
                                               phase_list=needed_phases)
        for arr in arrivals:
            tt[arr.name][:] = arr.time
    elif args.azimuth_index is not None:
        filenames = []
        y_ticks = []
        for i_dist, dist in enumerate(dist_grid):
            if args.station_range[0] < dist < args.station_range[1]:
                filenames.append(
                    args.sac_filename_formatter % (i_dist, args.azimuth_index))
                y_ticks.append(dist)
        y_label = 'Distance (deg)'
        y_title = f'Azimuth={np.round(azim_grid[args.azimuth_index], decimals=3)} deg'
        # travel time
        tt = {ph: np.full_like(y_ticks, np.nan) for ph in needed_phases}
        for i_dist, dist in enumerate(y_ticks):
            arrivals = taup_model.get_travel_times(source_depth_in_km=event_depth,
                                                   distance_in_degree=dist,
                                                   phase_list=needed_phases)
            for arr in arrivals:
                tt[arr.name][i_dist] = arr.time
    else:
        raise ValueError('Either --distance-index or --azimuth-index must be provided.')

    # read data
    traces = []
    global_max_amp = -1e10
    r_title = args.simulation \
        if args.reference is None else f'{args.simulation} - {args.reference}'
    for i, (y_tick, filename) in enumerate(zip(y_ticks, filenames)):
        # window
        t0 = UTCDateTime(0) + tt[args.aligned_phase][i] + args.window[0]
        t1 = UTCDateTime(0) + tt[args.aligned_phase][i] + args.window[1]
        # trace
        tr = read(f'{args.simulation}/{filename}')[0]
        tr.trim(starttime=t0, endtime=t1)
        if args.reference is not None:
            tr_ref = read(f'{args.reference}/{filename}')[0]
            tr_ref.trim(starttime=t0, endtime=t1)
            min_len = min(len(tr), len(tr_ref))
            tr.data[:min_len] -= tr_ref.data[:min_len]
        traces.append(tr)
        global_max_amp = max(global_max_amp, np.max(np.abs(tr.data)))

    # plot waveforms
    if args.figure_size is None:
        args.figure_size = [(args.window[1] - args.window[0]) / 5, len(y_ticks) / 4]
    plt.figure(dpi=200, figsize=args.figure_size)
    ax = plt.gca()
    tt_plot = {key: val.copy() for key, val in tt.items()}
    for i, (y_tick, filename) in enumerate(zip(y_ticks, filenames)):
        amp_factor = global_max_amp \
            if args.amplitude_normalize_mode == 'all' else np.max(np.abs(traces[i]))
        ax.plot(traces[i].times() + args.window[0],
                traces[i].data / amp_factor * args.amplitude_zoom_factor + i, c='k', lw=1)
        for key in tt_plot.keys():
            tt_plot[key][i] = tt_plot[key][i] - tt[args.aligned_phase][i]

    # plot travel times
    for ph in needed_phases:
        ax.plot(tt_plot[ph], range(len(y_ticks)), lw=2, label=ph)

    ax.set_xlabel(f'Time after {args.aligned_phase} phase (sec)')
    ax.set_ylabel(y_label)
    ax.set_yticks(range(len(y_ticks)), np.round(y_ticks, decimals=3))
    ax.legend(loc='upper right')
    ax.set_title('; '.join([r_title, y_title]))
    ax.set_xlim(args.window[0], + args.window[1])
    ax.set_ylim(-2, len(y_ticks) + 2)
    if args.output_file is not None:
        fig_dir = Path('figures')
        fig_dir.mkdir(exist_ok=True, parents=True)
        plt.savefig(fig_dir / args.output_file, bbox_inches='tight',
                    pad_inches=0.02, facecolor='w')
    else:
        plt.show()
