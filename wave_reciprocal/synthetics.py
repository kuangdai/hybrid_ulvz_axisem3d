from pathlib import Path

import numpy as np
import torch
import tqdm
import xarray as xr


class AxiSEM3DSyntheticsLoader:
    def __init__(self, nc_path, station_file, in_memory=False):
        self.nc_path = Path(nc_path)
        self.station_file = station_file

        # Read station data
        station_file_data = np.loadtxt(station_file, dtype=str)
        self.st_keys = np.char.add(np.char.add(station_file_data[:, 1], "."),
                                   station_file_data[:, 0])

        # Build dictionaries to map between keys to ranks
        rank_station_data = np.loadtxt(self.nc_path / "rank_station.info", dtype=str, skiprows=1)
        self.station_to_rank = {}
        self.rank_to_station = {}
        for row in rank_station_data:
            rank, key, index = int(row[0]), row[1], int(row[2])
            self.station_to_rank[key] = (rank, index)
            if rank not in self.rank_to_station:
                self.rank_to_station[rank] = []
            self.rank_to_station[rank].append(key)

        # Determine dimensions (components and time)
        first_rank = self.station_to_rank[self.st_keys[0]]
        first_file = self.nc_path / f"axisem3d_synthetics.nc.rank{first_rank[0]}"
        with xr.open_dataset(first_file) as ds:
            self.n_dim = ds["data_wave"].shape[1]
            self.n_time = ds["data_wave"].shape[0]
            self.times = np.array(ds["data_time"].values)

        # Open files for read
        self.data_on_rank = {}
        for rank in self.rank_to_station.keys():
            nc_file = self.nc_path / f"axisem3d_synthetics.nc.rank{rank}"
            if in_memory:
                with xr.open_dataset(nc_file) as ds:
                    self.data_on_rank[rank] = np.array(ds["data_wave"].values)
            else:
                self.data_on_rank[rank] = xr.open_dataset(nc_file)["data_wave"]

    def get(self, st_keys=None, start_time=None, end_time=None, time_interval=1):
        if st_keys is None:
            st_keys = self.st_keys  # get all stations

        # Read
        rank_to_local_global = {}
        for i_st, key in enumerate(st_keys):
            rank, index = self.station_to_rank[key]
            if rank not in rank_to_local_global:
                rank_to_local_global[rank] = ([], [])
            rank_to_local_global[rank][0].append(index)
            rank_to_local_global[rank][1].append(i_st)

        # Time
        if start_time is None:
            start_time = 0
        if end_time is None:
            end_time = self.n_time
        n_time_out = len(np.arange(start_time, end_time, time_interval))

        # Create
        u = np.empty((n_time_out, self.n_dim, len(st_keys)), dtype=np.float32)

        # Read
        for rank, (local_ids, global_ids) in rank_to_local_global.items():
            u[:, :, global_ids] = self.data_on_rank[rank][start_time:end_time:time_interval, :, local_ids]
        return u


def rotate(data, fr_frame, to_frame, batch_size, device):
    data_new = data.copy()
    rot = np.einsum('nij,njk->nik', fr_frame, to_frame.swapaxes(1, 2))
    for start in tqdm.trange(0, data.shape[-1], batch_size, leave=False):
        u = torch.from_numpy(data[:, :, start:start + batch_size]).to(device)
        m = torch.from_numpy(rot[start:start + batch_size]).to(device, torch.float32)
        u1 = torch.einsum('tin,nij->tjn', u, m)
        data_new[:, :, start:start + batch_size] = u1.cpu().numpy()
    return data_new
