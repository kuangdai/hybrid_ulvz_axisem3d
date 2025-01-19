import argparse
import json
from pathlib import Path

import numpy as np
import torch
import tqdm

from geodetic import GeoPoints
from loader import AxiSEM3DSyntheticsLoader

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
    parser.add_argument("--view", choices=["side", "top"], required=True,
                        help="View.")
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
    if args.view == "top":
        grid_depth_anim = np.array([2891.])
    else:
        grid_depth_anim = np.loadtxt(in_path / f"02_stations/grid_depth_anim_side_{args.medium}.txt")
    grid_dist_anim = np.loadtxt(in_path / f"02_stations/grid_dist_anim_{args.view}.txt")
    grid_azim_anim = np.loadtxt(in_path / f"02_stations/grid_azim_anim_{args.view}.txt")

    # Read data
    print("Reading raw data...")
    ds = AxiSEM3DSyntheticsLoader(
        in_path / f"12_solve_3d/output/stations/ANIM_{args.view.upper()}_{args.medium.upper()}",
        in_path / f"02_stations/STATIONS_ANIM_{args.view.upper()}_{args.medium.upper()}",
        meta["ulvz_lat"], meta["ulvz_lon"], src_spz=True)
    t0_id = np.searchsorted(ds.times, args.t0)
    t1_id = np.searchsorted(ds.times, args.t1)
    anim_data = ds.get(start_time=t0_id, end_time=t1_id, time_interval=args.time_interval)
    anim_data = anim_data.reshape(anim_data.shape[0], anim_data.shape[1],
                                  len(grid_dist_anim),
                                  len(grid_depth_anim),
                                  len(grid_azim_anim)).swapaxes(2, 3)

    # Read stations
    print("Reading stations...")
    st = np.loadtxt(in_path / f"02_stations/STATIONS_ANIM_{args.view.upper()}_{args.medium.upper()}", dtype=str)
    lat = st[:, 2].astype(float)
    lon = st[:, 3].astype(float)
    dep = st[:, 5].astype(float)
    points = GeoPoints(np.array([lat, lon, dep / 1e3]).T)

    # Rotate
    print("Rotating...")
    to_frame = points.form_RTZ_frame(meta["event_lat"], meta["event_lon"])
    to_frame = to_frame.reshape(len(grid_dist_anim),
                                len(grid_depth_anim),
                                len(grid_azim_anim), 3, 3).swapaxes(0, 1).reshape(-1, 3, 3)
    rot = np.einsum('nij,njk->nik', ds.frame, to_frame.swapaxes(1, 2))
    anim_data = anim_data.reshape(anim_data.shape[0], anim_data.shape[1], -1)
    for start in tqdm.trange(0, anim_data.shape[-1], args.batch_size, leave=False):
        u = torch.from_numpy(anim_data[:, :, start:start + args.batch_size]).to(args.device)
        m = torch.from_numpy(rot[start:start + args.batch_size]).to(args.device, torch.float32)
        u1 = torch.einsum('nij,tjn->tin', m, u)
        anim_data[:, :, start:start + args.batch_size] = u.cpu().numpy()
    anim_data = anim_data.reshape(anim_data.shape[0], anim_data.shape[1],
                                  len(grid_depth_anim),
                                  len(grid_dist_anim),
                                  len(grid_azim_anim))

    # Save
    print("Saving...")
    out_path = Path(args.output_path)
    out_path.mkdir(parents=True, exist_ok=True)
    np.savez(out_path / f"outgoing_{args.view}_{args.medium}.npz", anim_data)
    np.savetxt(out_path / f"outgoing_{args.view}_{args.medium}_times.txt", ds.times[t0_id:t1_id:args.time_interval])
