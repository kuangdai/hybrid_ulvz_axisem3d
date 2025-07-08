import argparse
import json
import numpy as np

from geodetic import GeoPoints
from synthetics import AxiSEM3DSyntheticsLoader, rotate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Outgoing Wave.")
    parser.add_argument("--station-name", type=str, required=True,
                        help="Station name.")
    parser.add_argument("--medium", choices=["solid", "fluid"], required=True,
                        help="Medium.")
    parser.add_argument("--batch-size", type=int, default=8192,
                        help="Batch size.")
    parser.add_argument("--device", type=str, default="cpu",
                        help="Device.")
    args = parser.parse_args()

    # Read meta
    with open("../wave_extrapolation/outputs/args.json", "r") as fs:
        meta = json.load(fs)
        ulvz_lat, ulvz_lon = meta["ulvz_lat"], meta["ulvz_lon"]
    with open(f"reciprocal_simulations/{args.station_name}/input/args.json", "r") as fs:
        meta = json.load(fs)
        st_lat, st_lon = meta["station_lat"], meta["station_lon"]

    # Read data
    print("Reading raw data...")
    ds = AxiSEM3DSyntheticsLoader(
        f"reciprocal_simulations/{args.station_name}/output/stations/{args.medium}",
        f"reciprocal_simulations/{args.station_name}/input/reciprocal_stations_{args.medium}.txt", )
    data = ds.get()

    # Read stations
    print("Reading reciprocal stations...")
    st = np.loadtxt(f"reciprocal_simulations/{args.station_name}/input/reciprocal_stations_{args.medium}.txt",
                    dtype=str)
    lat = st[:, 2].astype(float)
    lon = st[:, 3].astype(float)
    dep = st[:, 5].astype(float)
    points = GeoPoints(np.array([lat, lon, dep / 1e3]).T)

    # Rotate
    if args.medium == "solid" or True:  # Fluid also uses U
        print("Rotating...")
        fr_frame = points.form_RTZ_frame(st_lat, st_lon)
        to_frame = points.form_RTZ_frame(ulvz_lat, ulvz_lon)
        data = rotate(data, fr_frame, to_frame, args.batch_size, args.device)

    # Save
    print("Saving...")
    np.savez(f"reciprocal_simulations/{args.station_name}/output/rotated_{args.medium}.npz", data)
