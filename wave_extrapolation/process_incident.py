import argparse
import json
from pathlib import Path

import numpy as np
import tqdm
from scipy.interpolate import interp1d

from geodetic import GeoPoints
from synthetics import AxiSEM3DSyntheticsLoader, rotate


def dist_azim_interp(data, dist_grid, target_dist_azim, batch_size=1024):
    n_depth, n_dist, n_azim, _ = target_dist_azim.shape
    target_dist_azim = target_dist_azim.reshape(target_dist_azim.shape[0], -1, 2)
    output = np.empty((data.shape[0], data.shape[1],
                       target_dist_azim.shape[0],
                       target_dist_azim.shape[1]), dtype=np.float32)
    nfft = data.shape[-1]
    freq = np.fft.fftfreq(nfft, d=1.0) * nfft
    # depth loop
    for i_depth in tqdm.trange(n_depth, leave=False):
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
                fft_data = np.fft.fft(dist_data, axis=3, n=nfft) / nfft
                exp_terms = np.exp(1j * np.outer(freq, azim))
                azim_data = np.einsum("TCda,ad->TCd", fft_data, exp_terms).real
            else:
                azim_data = dist_data[..., 0]  # monopole

            # dist interp
            output[:, :, i_depth, start_:start_ + batch_size] = azim_data
    return output.reshape(data.shape[0], data.shape[1], n_depth, n_dist, n_azim)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Incident Wave.")
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
    parser.add_argument("--batch-size", type=int, default=8192,
                        help="Batch size.")
    parser.add_argument("--device", type=str, default="cpu",
                        help="Device.")
    parser.add_argument("--output-path", type=str, default="./outputs",
                        help="Output path.")
    args = parser.parse_args()

    # Read meta
    in_path = Path(args.input_path)
    with open(in_path / "02_stations/args.json", "r") as fs:
        meta = json.load(fs)
    out_path = Path(args.output_path)
    out_path.mkdir(parents=True, exist_ok=True)
    with open(out_path / "args.json", "w") as fs:
        json.dump(meta, fs)

    # Read grids
    grid_dist_incident = np.loadtxt(in_path / "02_stations/grid_dist.txt")
    grid_azim_incident = np.loadtxt(in_path / "02_stations/grid_azim.txt")
    if np.ndim(grid_azim_incident) == 0:
        grid_azim_incident = np.array([grid_azim_incident])
    grid_depth_extrap = np.loadtxt(in_path / f"02_stations/grid_depth_{args.medium}.txt")
    grid_dist_extrap = np.loadtxt(in_path / f"02_stations/grid_dist_extrapolation.txt")
    grid_azim_extrap = np.loadtxt(in_path / f"02_stations/grid_azim_extrapolation.txt")

    # Read data
    print("Reading raw data...")
    medium_u = args.medium
    if args.medium == "fluid":
        medium_u = "fluid_u"
    ds = AxiSEM3DSyntheticsLoader(in_path / f"04_solve_prem/output/stations/INCIDENT_{medium_u.upper()}",
                                  in_path / f"02_stations/STATIONS_INCIDENT_{medium_u.upper()}")
    t0_id = np.searchsorted(ds.times, args.t0)
    t1_id = np.searchsorted(ds.times, args.t1)
    extrap_data = ds.get(start_time=t0_id, end_time=t1_id, time_interval=args.time_interval)
    extrap_data = extrap_data.reshape(extrap_data.shape[0], extrap_data.shape[1],
                                      len(grid_dist_incident),
                                      len(grid_depth_extrap),
                                      len(grid_azim_incident)).swapaxes(2, 3)

    # Read stations
    print("Reading stations...")
    st = np.loadtxt(in_path / f"02_stations/STATIONS_EXTRAPOLATION_{args.medium.upper()}",
                    dtype=str)
    lat = st[:, 2].astype(float)
    lon = st[:, 3].astype(float)
    dep = st[:, 5].astype(float)
    points = GeoPoints(np.array([lat, lon, dep / 1e3]).T)
    dist_azim = points.get_rtp_src_centered(meta["event_lat"], meta["event_lon"])[:, 1:]
    dist_azim = dist_azim.reshape(len(grid_dist_extrap),
                                  len(grid_depth_extrap),
                                  len(grid_azim_extrap), 2).swapaxes(0, 1)

    # Interpolation
    print("Spatial interpolation...")
    extrap_data = dist_azim_interp(extrap_data, grid_dist_incident, dist_azim, args.batch_size)

    # Rotate
    if args.medium == "solid" or True:  # Now fluid also uses U
        print("Rotating...")
        fr_frame = points.form_spz_frame(meta["event_lat"], meta["event_lon"])
        fr_frame = fr_frame.reshape(len(grid_dist_extrap),
                                    len(grid_depth_extrap),
                                    len(grid_azim_extrap), 3, 3).swapaxes(0, 1).reshape(-1, 3, 3)
        to_frame = points.form_RTZ_frame(meta["ulvz_lat"], meta["ulvz_lon"])
        to_frame = to_frame.reshape(len(grid_dist_extrap),
                                    len(grid_depth_extrap),
                                    len(grid_azim_extrap), 3, 3).swapaxes(0, 1).reshape(-1, 3, 3)
        extrap_data = extrap_data.reshape(extrap_data.shape[0], extrap_data.shape[1], -1)
        extrap_data = rotate(extrap_data, fr_frame, to_frame, args.batch_size, args.device)
        extrap_data = extrap_data.reshape(extrap_data.shape[0], extrap_data.shape[1],
                                          len(grid_depth_extrap),
                                          len(grid_dist_extrap),
                                          len(grid_azim_extrap))

    # Save
    print("Saving...")
    np.savez(out_path / f"incident_with_box_{args.medium}.npz", extrap_data)
    np.savetxt(out_path / f"incident_with_box_{args.medium}_times.txt", ds.times[t0_id:t1_id:args.time_interval])
