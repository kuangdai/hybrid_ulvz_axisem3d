import xarray as xr
import numpy as np
import argparse


def find_nearest_left(time_array, target_time):
    """Find the nearest left value in a sorted time array."""
    return time_array[time_array <= target_time][-1]  # Last occurrence before or equal to target_time


def truncate_and_taper(file_path, t0, t1):
    """
    Truncate early time steps and apply tapering between t0 and t1.
    Output is saved as input_name_truncated.nc.

    Parameters:
        file_path (str): Path to the NetCDF file.
        t0 (float): The user-specified start time for tapering.
        t1 (float): The user-specified end time for tapering.
    """
    print(f"Truncating and tapering file: {file_path}")

    # Load dataset
    ds = xr.open_dataset(file_path)

    # Identify the time variable and its dimension
    time_var = "time_points"
    if time_var not in ds:
        raise ValueError(f"{time_var} not found in dataset.")
    time_dim = ds[time_var].dims[0]  # Detect the actual time dimension name

    # Extract time values
    time_values = ds[time_var].values

    # Find actual t0 and t1 using nearest left search
    actual_t0 = find_nearest_left(time_values, t0)
    actual_t1 = find_nearest_left(time_values, t1)
    print(f"Using t0: {actual_t0}, t1: {actual_t1}")

    # Compute tapering weights
    taper_weights = np.ones_like(time_values)
    taper_weights[time_values < actual_t0] = 0  # Fully remove before actual_t0
    taper_weights[time_values > actual_t1] = 1  # Fully preserve after actual_t1

    # Apply smooth tapering between actual_t0 and actual_t1 using a cosine function
    mask = (time_values >= actual_t0) & (time_values <= actual_t1)
    taper_weights[mask] = 0.5 * (1 - np.cos(np.pi * (time_values[mask] - actual_t0) / (actual_t1 - actual_t0)))

    # Apply tapering to all time-dependent variables except time_points
    for var in ds.data_vars:
        if time_dim in ds[var].dims and var != time_var:
            ds[var] = ds[var] * taper_weights[:, np.newaxis]  # Broadcasting over second dimension

    # Remove truncated time steps (keep only t0 and beyond)
    ds = ds.sel({time_dim: ds[time_var] >= actual_t0})

    # Define output filename
    output_path = file_path.replace(".nc", "_truncated.nc")
    ds.to_netcdf(output_path, mode="w")
    print(f"Saved truncated file: {output_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Truncate and taper NetCDF files.")
    parser.add_argument("--t0", type=float, required=True, help="Start time for truncation (before this is removed).")
    parser.add_argument("--t1", type=float, required=True, help="End time for tapering (smooth transition).")
    parser.add_argument("--files", nargs="+",
                        default=["output/injection_FLUID.nc",
                                 "output/injection_SOLID.nc"], help="NetCDF file(s) to process.")
    args = parser.parse_args()

    # Process all provided files
    for file in args.files:
        truncate_and_taper(file, args.t0, args.t1)
