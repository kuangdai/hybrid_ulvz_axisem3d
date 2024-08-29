import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm

from geodetic import GeoPoints, rotate_vector_frame


def merge_nc(nc_dir):
    # rank-station map
    rank_station = np.genfromtxt(f'{nc_dir}/rank_station.info',
                                 dtype=str, skip_header=1)
    ranks_unique = np.unique(rank_station[:, 0])
    rank_station_dict = {rank: [] for rank in ranks_unique}
    for rank, st_key, _ in rank_station:
        rank_station_dict[rank].append(st_key)

    # get dimensions
    with Dataset(f'{nc_dir}/axisem3d_synthetics.nc.rank{ranks_unique[0]}',
                 'r') as nc_in:
        times = nc_in.variables['time_points'][:]

    # create
    with Dataset(f'{nc_dir}/axisem3d_synthetics.nc', 'w') as nc_out:
        nc_out.createDimension('dim_time', len(times))
        nc_out.createDimension('dim_channel', 3)
        nc_out.createVariable('time_points', float, ('dim_time',))
        nc_out.variables['time_points'][:] = times
        for rank in tqdm(ranks_unique, total=len(ranks_unique)):
            with Dataset(f'{nc_dir}/axisem3d_synthetics.nc.rank{rank}',
                         'r') as nc_in:
                data = nc_in.variables['data'][:, :, :]
                for i_st, st_key in enumerate(rank_station_dict[rank]):
                    nc_out.createVariable(st_key, 'float32',
                                          ('dim_time', 'dim_channel'))
                    nc_out[st_key][:, :] = data[:, :, i_st]


def rotate_nc(nc_dir, st_file, e_lat, e_lon, u_lat, u_lon):
    # station file
    txt = np.genfromtxt(st_file, dtype=str)
    # keys
    st_keys = np.core.defchararray.add(
        txt[:, 1], np.core.defchararray.add(".", txt[:, 0]))
    key_index = {st_key: i_st for i_st, st_key in enumerate(st_keys)}
    # frames
    lat = txt[:, 2].astype(float)
    lon = txt[:, 3].astype(float)
    dep = txt[:, 5].astype(float)
    st_points = GeoPoints(np.array([lat, lon, dep / 1e3]).T)
    x_spz_ULVZ = st_points.form_spz_frame(u_lat, u_lon)
    x_RTZ_EVENT = st_points.form_RTZ_frame(e_lat, e_lon)

    # rank-station map
    rank_station = np.genfromtxt(f'{nc_dir}/rank_station.info',
                                 dtype=str, skip_header=1)
    ranks_unique = np.unique(rank_station[:, 0])
    rank_station_dict = {rank: [] for rank in ranks_unique}
    for rank, st_key, _ in rank_station:
        rank_station_dict[rank].append((st_key, key_index[st_key]))

    # get dimensions
    with Dataset(f'{nc_dir}/axisem3d_synthetics.nc.rank{ranks_unique[0]}',
                 'r') as nc_in:
        times = nc_in.variables['time_points'][:]

    # create
    with Dataset(f'{nc_dir}/axisem3d_synthetics.nc', 'w') as nc_out:
        nc_out.createDimension('dim_time', len(times))
        nc_out.createDimension('dim_channel', 3)
        nc_out.createVariable('time_points', float, ('dim_time',))
        nc_out.variables['time_points'][:] = times
        for rank in tqdm(ranks_unique, total=len(ranks_unique)):
            with Dataset(f'{nc_dir}/axisem3d_synthetics.nc.rank{rank}',
                         'r') as nc_in:
                data = nc_in.variables['data'][:, :, :]
                for i_st, (st_key, st_index) in enumerate(rank_station_dict[rank]):
                    u_spz_U = data[:, :, i_st]
                    u_RTZ_E = rotate_vector_frame(u_spz_U,
                                                  x_spz_ULVZ[st_index],
                                                  x_RTZ_EVENT[st_index])
                    nc_out.createVariable(st_key, 'float32',
                                          ('dim_time', 'dim_channel'))
                    nc_out[st_key][:, :] = u_RTZ_E
