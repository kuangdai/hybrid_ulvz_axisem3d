import argparse
import json
from pathlib import Path

import numpy as np

from geodetic import GeoPoints
from synthetics import AxiSEM3DSyntheticsLoader, rotate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Outgoing Wave.")
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

    # Read grids
    grid_depth_extrap = np.loadtxt(in_path / f"02_stations/grid_depth_{args.medium}.txt")
    grid_dist_extrap = np.loadtxt(in_path / f"02_stations/grid_dist_extrapolation.txt")
    grid_azim_extrap = np.loadtxt(in_path / f"02_stations/grid_azim_extrapolation.txt")

    # Read data
    print("Reading raw data...")
    ds = AxiSEM3DSyntheticsLoader(
        in_path / f"12_solve_3d/output/stations/EXTRAPOLATION_{args.medium.upper()}",
        in_path / f"02_stations/STATIONS_EXTRAPOLATION_{args.medium.upper()}")
    t0_id = np.searchsorted(ds.times, args.t0)
    t1_id = np.searchsorted(ds.times, args.t1)
    extrap_data = ds.get(start_time=t0_id, end_time=t1_id, time_interval=args.time_interval)
    extrap_data = extrap_data.reshape(extrap_data.shape[0], extrap_data.shape[1],
                                      len(grid_dist_extrap),
                                      len(grid_depth_extrap),
                                      len(grid_azim_extrap)).swapaxes(2, 3)

    # Read stations
    print("Reading stations...")
    st = np.loadtxt(in_path / f"02_stations/STATIONS_EXTRAPOLATION_{args.medium.upper()}", dtype=str)
    lat = st[:, 2].astype(float)
    lon = st[:, 3].astype(float)
    dep = st[:, 5].astype(float)
    phi_under_ulvz = np.zeros((len(grid_dist_extrap), len(grid_depth_extrap), len(grid_azim_extrap)))
    phi_under_ulvz[:, :, :] = grid_azim_extrap[None, None, :]
    points = GeoPoints(np.array([lat, lon, dep / 1e3]).T, phi_under_ulvz.reshape(-1))

    # Rotate
    if args.medium == "solid" or True:
        print("Rotating...")
        fr_frame = points.form_spz_frame(meta["ulvz_lat"], meta["ulvz_lon"])
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

    # Add incident
    print("Handling boundary...")
    out_path = Path(args.output_path)
    incident_data = np.load(out_path / f"incident_with_box_{args.medium}.npz")["arr_0"]
    extrap_data += incident_data

    if args.medium == "fluid":
        extrap_data[:, :, :-1, :-1, :] -= incident_data[:, :, :-1, :-1, :]
        extrap_data[:, :, -2, :, :] = (extrap_data[:, :, -2, :, :] + extrap_data[:, :, -1, :, :]) / 2
        extrap_data[:, :, :, -2, :] = (extrap_data[:, :, :, -2, :] + extrap_data[:, :, :, -1, :]) / 2
        extrap_data = extrap_data[:, :, :-1, :-1, :]
        incident_data[:, :, -2, :, :] = (incident_data[:, :, -2, :, :] + incident_data[:, :, -1, :, :]) / 2
        incident_data[:, :, :, -2, :] = (incident_data[:, :, :, -2, :] + incident_data[:, :, :, -1, :]) / 2
        incident_data = incident_data[:, :, :-1, :-1, :]
        grid_depth_extrap[-2] = (grid_depth_extrap[-2] + grid_depth_extrap[-1]) / 2
        grid_depth_extrap = grid_depth_extrap[:-1]
    else:
        extrap_data[:, :, 1:, :-1, :] -= incident_data[:, :, 1:, :-1, :]
        extrap_data[:, :, 1, :, :] = (extrap_data[:, :, 1, :, :] + extrap_data[:, :, 0, :, :]) / 2
        extrap_data[:, :, :, -2, :] = (extrap_data[:, :, :, -2, :] + extrap_data[:, :, :, -1, :]) / 2
        extrap_data = extrap_data[:, :, 1:, :-1, :]
        incident_data[:, :, 1, :, :] = (incident_data[:, :, 1, :, :] + incident_data[:, :, 0, :, :]) / 2
        incident_data[:, :, :, -2, :] = (incident_data[:, :, :, -2, :] + incident_data[:, :, :, -1, :]) / 2
        incident_data = incident_data[:, :, 1:, :-1, :]
        grid_depth_extrap[1] = (grid_depth_extrap[0] + grid_depth_extrap[1]) / 2
        grid_depth_extrap = grid_depth_extrap[1:]
    grid_dist_extrap[-2] = (grid_dist_extrap[-2] + grid_dist_extrap[-1]) / 2
    grid_dist_extrap = grid_dist_extrap[:-1]

    # Save
    print("Saving...")
    np.savez(out_path / f"outgoing_ml_ready_{args.medium}.npz", extrap_data)
    np.savetxt(out_path / f"outgoing_ml_ready_{args.medium}_times.txt", ds.times[t0_id:t1_id:args.time_interval])
    np.savez(out_path / f"incident_ml_ready_{args.medium}.npz", incident_data)
    np.savetxt(out_path / f"incident_ml_ready_{args.medium}_times.txt", ds.times[t0_id:t1_id:args.time_interval])
    np.savetxt(out_path / f"grid_depth_ml_ready_{args.medium}.txt", grid_depth_extrap)
    np.savetxt(out_path / f"grid_dist_ml_ready.txt", grid_dist_extrap)
    np.savetxt(out_path / f"grid_azim_ml_ready.txt", grid_azim_extrap)
