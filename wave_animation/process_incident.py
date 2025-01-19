import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from tqdm import trange

from geodetic import GeoPoints
from synthetics import AxiSEM3DSynthetics


def dist_azim_interp(data, dist_grid, target_dist_azim, batch_size=1024):
    n_depth, n_dist, n_azim, _ = target_dist_azim.shape
    target_dist_azim = target_dist_azim.reshape(target_dist_azim.shape[0], -1, 2)
    output = np.empty((data.shape[0], data.shape[1],
                       target_dist_azim.shape[0],
                       target_dist_azim.shape[1]), dtype=np.float32)
    nfft = data.shape[-1]
    k = (np.fft.fftfreq(nfft, d=(2 * np.pi / nfft)) * nfft)[:(nfft // 2 + 1)]
    # depth loop
    for i_depth in trange(n_depth, leave=False):
        # distance interpolator
        interp_dist = interp1d(dist_grid, data[:, :, i_depth, :, :], axis=2, kind="linear")

        # batch loop
        for start_ in range(0, target_dist_azim.shape[1], batch_size):
            dist = target_dist_azim[i_depth, start_:start_ + batch_size, 0]
            azim = target_dist_azim[i_depth, start_:start_ + batch_size, 1]

            # dist interpolation
            dist_data = interp_dist(dist)

            # azimuth interpolation
            if dist_data.shape[-1] > 1:
                fft_data = np.fft.rfft(dist_data, axis=3, n=nfft)
                exp_terms = np.exp(1j * np.outer(k, azim))
                azim_data = np.einsum("TCda,ad->TCd", fft_data, exp_terms).real
            else:
                azim_data = dist_data[..., 0]  # monopole

            # dist interp
            output[:, :, i_depth, start_:start_ + batch_size] = azim_data
    return output.reshape(data.shape[0], data.shape[1], n_depth, n_dist, n_azim)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Incident Wave.")
    parser.add_argument("--input-path", type=str, required=True,
                        help="Path of AxiSEM3D output.")
    parser.add_argument("--t0", type=float, required=True,
                        help="Starting time.")
    parser.add_argument("--t1", type=float, required=True,
                        help="Ending time.")
    parser.add_argument("--time-interval", type=int, default=1,
                        help="Time interval (step).")
    parser.add_argument("--medium", choices=["solid", "fluid"], required=True,
                        help="Medium.")
    parser.add_argument("--view", choices=["side", "top"], required=True,
                        help="View.")
    parser.add_argument("--batch-size", type=int, default=8192,
                        help="Batch size.")
    parser.add_argument("--output-path", type=str, default="./outputs",
                        help="Output path.")
    args = parser.parse_args()

    # Read meta
    in_path = Path(args.input_path)
    with open(in_path / "02_stations/args.json", "r") as fs:
        meta = json.load(fs)

    # Read grids
    grid_dist = np.loadtxt(in_path / "02_stations/grid_dist_anim_incident.txt")
    grid_azim = np.loadtxt(in_path / "02_stations/grid_azim_anim_incident.txt")
    if not hasattr(grid_azim, "__len__"):
        grid_azim = np.array([grid_azim])
    grid_depth = np.loadtxt(in_path / f"02_stations/grid_depth_anim_incident_{args.medium}.txt")
    grid_dist_anim = np.loadtxt(in_path / f"02_stations/grid_dist_anim_{args.view}.txt")
    grid_azim_anim = np.loadtxt(in_path / f"02_stations/grid_azim_anim_{args.view}.txt")
    grid_depth_top = np.array([2891.])

    # Read data
    print("Reading raw data...")
    ds = AxiSEM3DSynthetics(in_path / f"04_solve_prem/output/stations/ANIM_INCIDENT_{args.medium.upper()}",
                            in_path / f"02_stations/STATIONS_ANIM_INCIDENT_{args.medium.upper()}",
                            meta["event_lat"], meta["event_lon"], src_spz=False)
    t0_id = np.searchsorted(ds.times, args.t0)
    t1_id = np.searchsorted(ds.times, args.t1)
    raw_data = ds.get(start_time=t0_id, end_time=t1_id, time_interval=args.time_interval)
    raw_data = raw_data.reshape(raw_data.shape[0], raw_data.shape[1],
                                len(grid_dist),
                                len(grid_depth),
                                len(grid_azim)).swapaxes(2, 3)

    # Handle depth of top view
    if args.view == "top":
        if args.medium == "solid":
            raw_data = raw_data[:, :, -1:, :, :]
        else:
            raw_data = raw_data[:, :, :1, :, :]
        grid_depth = np.array([2891.])

    # Read stations
    print("Reading stations...")
    ex_st = np.loadtxt(in_path / f"02_stations/STATIONS_ANIM_{args.view.upper()}_{args.medium.upper()}", dtype=str)
    ex_lat = ex_st[:, 2].astype(float)
    ex_lon = ex_st[:, 3].astype(float)
    ex_dep = ex_st[:, 5].astype(float)
    points = GeoPoints(np.array([ex_lat, ex_lon, ex_dep / 1e3]).T)
    dist_azim = points.get_rtp_src_centered(meta["event_lat"], meta["event_lon"])[:, 1:]
    dist_azim = dist_azim.reshape(len(grid_dist_anim),
                                  len(grid_depth),
                                  len(grid_azim_anim), 2).swapaxes(0, 1)

    # Interpolation
    print("Spatial interpolation...")
    anim_data = dist_azim_interp(raw_data, grid_dist, dist_azim, args.batch_size)
    anim_data = anim_data.reshape(anim_data.shape[0],
                                  anim_data.shape[1],
                                  len(grid_depth),
                                  len(grid_dist_anim),
                                  len(grid_azim_anim))

    # Save
    print("Saving...")
    out_path = Path(args.output_path)
    out_path.mkdir(parents=True, exist_ok=True)
    np.savez(out_path / f"incident_{args.view}_{args.medium}.npz", anim_data)
    np.savetxt(out_path / f"incident_{args.view}_{args.medium}_times.txt", ds.times[t0_id:t1_id:args.time_interval])
