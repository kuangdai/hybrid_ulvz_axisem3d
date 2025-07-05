import argparse
import json
import numpy as np
import os
from pathlib import Path

from geodetic import GeoPoints


def replace_in_file(source, replace_dict, dest=None):
    if dest is None:
        dest = source
    with open(source, 'r') as fs:
        text = fs.read()
    for key, val in replace_dict.items():
        text = text.replace(key, str(val))
    with open(dest, 'w') as fs:
        fs.write(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate reciprocal wave.")
    parser.add_argument("--run-name", type=str, required=True,
                        help="Name of simulation.")
    parser.add_argument("--station-name", type=str, required=True,
                        help="Station Name.")
    parser.add_argument("--station-lat", type=float, required=True,
                        help="Station latitude.")
    parser.add_argument("--station-lon", type=float, required=True,
                        help="Station longitude.")
    parser.add_argument("--component", choices=["R", "T", "Z"], required=True,
                        help="Station component to use.")
    parser.add_argument("--dt", type=float, default=None,
                        help="Time interval.")
    parser.add_argument("--length", type=float, default=1000.,
                        help="Record length.")
    args = parser.parse_args()

    # Prepare output
    out_dir = Path(f"./reciprocal_simulations/{args.station_name}/input")
    out_dir.mkdir(parents=True, exist_ok=True)
    os.system(f"cp template/*.yaml {out_dir}")

    args_input = json.load(open(f'../simulations/inputs/{args.run_name}/args.json'))
    in_path = Path(f'../simulations/outputs/{args.run_name}')

    # Copy mesh
    os.system(f"cp {in_path / '01_exodus/mesh.e'} {out_dir}")

    # Compute direction of vector force
    with open(in_path / "02_stations/args.json", "r") as fs:
        meta = json.load(fs)
    e_lat, e_lon = meta["event_lat"], meta["event_lon"]
    s_lat, s_lon = args.station_lat, args.station_lon
    if s_lon > 180.:
        s_lon = s_lon - 360.
    s_pnt = GeoPoints(lld=np.array([[s_lat, s_lon, 0.]]))
    e_rtz = s_pnt.form_RTZ_frame(e_lat, e_lon)[0]
    force_vec = e_rtz["RTZ".index(args.component)]
    s_rtz = s_pnt.form_RTZ_frame(min(s_lat + 1., 89.9999), s_lon)[0]
    if args.component == "Z":
        amp_r = 0.
        amp_t = 0.
        amp_z = 1.
    else:
        amp_r = np.dot(force_vec, s_rtz[0])
        amp_t = np.dot(force_vec, s_rtz[1])
        amp_z = 0.

    # Prepare source
    dt = args.dt
    if dt is None:
        dt = args_input["time_series"]["dt"]
        if dt == 0.0:
            raise ValueError("Cannot determine DT from config. Specify --dt.")

    replace_in_file(f"{out_dir}/inparam.source.yaml",
                    {
                        "__DT__": dt,
                        "__LAT__": s_lat,
                        "__LON__": s_lon,
                        "__R__": amp_r,
                        "__T__": amp_t,
                        "__Z__": amp_z,
                        "__LENGTH__": args.length,
                    })

    # Prepare stations
    replace_in_file("../wave_convolve/surface_nodes_cache/reciprocal_stations_solid.txt",
                    replace_dict={},
                    dest=f"{out_dir}/reciprocal_stations_solid.txt")
    replace_in_file("../wave_convolve/surface_nodes_cache/reciprocal_stations_fluid.txt",
                    replace_dict={},
                    dest=f"{out_dir}/reciprocal_stations_fluid.txt")

    # Also save args
    with open(out_dir / "args.json", "w") as fs:
        json.dump({
            "station_name": args.station_name,
            "station_lat": args.station_lat,
            "station_lon": args.station_lon
        }, fs)

# Example
"""
python xgen.py --run-name test --station-name US.ERP --station-lat 42.12 --station-lon 280.01 --component T
"""
