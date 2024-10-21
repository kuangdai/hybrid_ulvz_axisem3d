import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

from geodetic import GeoPoints


def coords_to_grids(*coords):
    grids = np.meshgrid(*coords, indexing='ij')
    return [g.reshape(-1) for g in grids], [len(c) for c in coords]


def reshape(*grids, grid_shape):
    return [g.reshape(grid_shape) for g in grids]


def to_station_file(grid_depth, media, grid_dist_, grid_azim_, cen_lat, cen_lon, wave_name):
    # compute lld
    (g_dist, g_depth, g_azim), grid_shape = coords_to_grids(
        grid_dist_, grid_depth, grid_azim_)
    st_points = GeoPoints.create_rtp_src_centered(
        np.array([6371. - g_depth, g_dist, g_azim]).T, cen_lat, cen_lon)
    st_lat_ = st_points.lld[:, 0]
    st_lon_ = st_points.lld[:, 1]
    st_dep_ = st_points.lld[:, 2]
    st_lat_, st_lon_, st_dep_ = reshape(st_lat_, st_lon_, st_dep_,
                                        grid_shape=grid_shape)
    # write
    with open(out_dir / f'STATIONS_{wave_name}_{media.upper()}', 'w') as fstream:
        for i_dist_, _ in enumerate(grid_dist_):
            for i_depth_, _ in enumerate(grid_depth):
                for i_azim_, _ in enumerate(grid_azim_):
                    fstream.write('%s %s %.18f %.18f 0 %.18f force_%s\n' % (
                        'ULVZ_%d_%d_%d' % (i_dist_, i_depth_, i_azim_),
                        f'INCIDENT_{media.upper()}',
                        st_lat_[i_dist_, i_depth_, i_azim_],
                        st_lon_[i_dist_, i_depth_, i_azim_],
                        st_dep_[i_dist_, i_depth_, i_azim_] * 1e3,
                        media.lower()))


if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_stations')
    out_dir.mkdir(parents=True, exist_ok=True)

    #########
    # event #
    #########
    e_data = np.loadtxt(f'inputs/{run_name}/CMTSOLUTION',
                        skiprows=1, dtype=str, delimiter=':')
    e_lat = float(e_data[3, 1])
    e_lon = float(e_data[4, 1])
    e_dep = float(e_data[5, 1])
    if e_lon < 0:
        e_lon += 360.

    ########
    # ULVZ #
    ########
    if args['ulvz']['center_loc_by_lat_lon']:
        u_lat = args['ulvz']['lat']
        u_lon = args['ulvz']['lon']
        u_pnt = GeoPoints(np.array([[u_lat, u_lon, 2891.]]))
        u_rtp = u_pnt.get_rtp_src_centered(e_lat, e_lon)
        u_dist, u_azim = u_rtp[0, 1:3]
    else:
        u_dist = np.radians(args['ulvz']['dist'])
        u_azim = np.radians(args['ulvz']['azim'])
        u_pnt = GeoPoints.create_rtp_src_centered(
            np.array([[3480., u_dist, u_azim]]), e_lat, e_lon)
        u_lat, u_lon = u_pnt.lld[0, 0:2]

    # furthest edge point
    r_ulvz = np.radians(args['ulvz']['radius_deg'])
    u_pnt = GeoPoints.create_rtp_src_centered(
        np.array([[3480., u_dist + r_ulvz, u_azim]]), e_lat, e_lon)
    u_lat_far, u_lon_far = u_pnt.lld[0, 0:2]

    # nearest edge point
    u_pnt = GeoPoints.create_rtp_src_centered(
        np.array([[3480., u_dist - r_ulvz, u_azim]]), e_lat, e_lon)
    u_lat_near, u_lon_near = u_pnt.lld[0, 0:2]

    #########
    # ARRAY #
    #########
    if args['array']['use_grid']:
        st_dist_crds = np.radians(np.arange(args['array']['dist0_from_event'],
                                            args['array']['dist1_from_event'] + 1e-4,
                                            args['array']['delta_dist']))
        st_azim_crds = np.radians(np.arange(args['array']['azim0_from_ulvz'],
                                            args['array']['azim1_from_ulvz'] + 1e-4,
                                            args['array']['delta_azim']))
        st_azim_crds += u_azim
        (st_dist, st_azim), shape = coords_to_grids(st_dist_crds, st_azim_crds)
        nst = len(st_dist)
        st_rtp = np.array([np.ones(nst) * 6371., st_dist, st_azim]).T
        st_pnt = GeoPoints.create_rtp_src_centered(st_rtp, e_lat, e_lon)
        st_lat = st_pnt.lld[:, 0]
        st_lon = st_pnt.lld[:, 1]
        st_lat, st_lon = reshape(st_lat, st_lon, grid_shape=shape)
        with open(out_dir / 'STATIONS_ARRAY', 'w') as fs:
            for i_dist in range(shape[0]):
                for i_azim in range(shape[1]):
                    name = f'D{i_dist}_A{i_azim}'
                    fs.write(f'{name} ARRAY {st_lat[i_dist, i_azim]:.18f} '
                             f'{st_lon[i_dist, i_azim]:.18f} 0 0\n')
    else:
        st = np.loadtxt(f'inputs/{run_name}/STATIONS', dtype=str)
        st_out = st.copy()
        st_out[:, 0] = np.char.add(st[:, 1], np.full_like(st[:, 1], '.'))
        st_out[:, 0] = np.char.add(st_out[:, 0], st[:, 0])
        st_out[:, 1] = 'ARRAY'
        np.savetxt(out_dir / 'STATIONS_ARRAY', st_out, fmt='%s')

    ############
    # plot map #
    ############
    plt.figure(dpi=200)
    # map
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='ivory', lake_color='lightblue')
    m.drawmapboundary(fill_color='lightblue')

    # earthquake
    m.scatter(e_lon, e_lat, c='red', s=50, marker='*')

    # ulvz
    m.scatter(u_lon, u_lat, c='orange', s=200)
    m.scatter(u_lon_far, u_lat_far, c='blue', s=10, marker='+')
    m.scatter(u_lon_near, u_lat_near, c='blue', s=10, marker='+')

    # stations
    s_data = np.loadtxt(out_dir / 'STATIONS_ARRAY', dtype=str)
    s_lat = s_data[:, 2].astype(float)
    s_lon = s_data[:, 3].astype(float)
    m.scatter(s_lon, s_lat, c='green', s=1, marker='v', linewidth=1)
    plt.savefig(out_dir / 'map_array.png', bbox_inches='tight', pad_inches=0.0)

    #####################
    # incident stations #
    #####################
    # dist sampling
    args_mesh = json.load(open(f'outputs/{run_name}/@@_exodus/args.json'))
    dist_delta = np.pi / (args_mesh['NEX'] * 4)
    dist_half_range = (args_mesh['NEX_U'] + 2) * 4 * dist_delta
    n_dist = int(np.ceil(dist_half_range * 2 / dist_delta))
    grid_dist = np.linspace(u_dist - dist_half_range,
                            u_dist + dist_half_range, n_dist)

    # depth sampling
    depth_delta = args['ulvz']['height'] / (args_mesh['n_elem_layer_ulvz'] * 4)
    height_solid = args['ulvz']['height'] + args['box']['height']
    n_solid = int(np.ceil(height_solid / depth_delta))
    grid_depth_solid = np.linspace(2891. - height_solid, 2891., n_solid)
    height_fluid = args['box']['height']
    n_fluid = int(np.ceil(height_fluid / depth_delta))
    grid_depth_fluid = np.linspace(2891., 2891. + height_fluid, n_fluid)

    # azimuth sampling
    if args['event']['monopole']:
        grid_azim = np.array([u_azim])
    else:
        grid_azim = np.radians(np.array([0., 72., 144., 216., 288.]))

    # store grid
    np.savetxt(out_dir / 'grid_dist.txt', grid_dist)
    np.savetxt(out_dir / 'grid_depth_SOLID.txt', grid_depth_solid)
    np.savetxt(out_dir / 'grid_depth_FLUID.txt', grid_depth_fluid)
    np.savetxt(out_dir / 'grid_azim.txt', grid_azim)

    # write stations
    to_station_file(grid_depth_solid, 'solid', grid_dist, grid_azim, e_lat, e_lon, "INCIDENT")
    to_station_file(grid_depth_fluid, 'fluid', grid_dist, grid_azim, e_lat, e_lon, "INCIDENT")

    # final combination
    fout = open(out_dir / 'STATIONS_ARRAY_INCIDENT', 'w')
    fin = open(out_dir / 'STATIONS_ARRAY')
    fout.write(fin.read())
    fin = open(out_dir / 'STATIONS_INCIDENT_SOLID')
    fout.write(fin.read())
    fin = open(out_dir / 'STATIONS_INCIDENT_FLUID')
    fout.write(fin.read())
    fout.close()

    # write derived info for later use
    args_stations = {'ulvz_dist_near': u_dist - r_ulvz,
                     'ulvz_dist_far': u_dist + r_ulvz,
                     'ulvz_lat': u_lat,
                     'ulvz_lon': u_lon,
                     'event_lat': e_lat,
                     'event_lon': e_lon
                     }
    json.dump(args_stations, open(out_dir / "args.json", 'w'))

    ############
    # plot map #
    ############
    plt.figure(dpi=200)
    # map
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='ivory', lake_color='lightblue')
    m.drawmapboundary(fill_color='lightblue')

    # earthquake
    m.scatter(e_lon, e_lat, c='red', s=50, marker='*')

    # ulvz
    m.scatter(u_lon, u_lat, c='orange', s=200)
    m.scatter(u_lon_far, u_lat_far, c='blue', s=10, marker='+')
    m.scatter(u_lon_near, u_lat_near, c='blue', s=10, marker='+')

    # stations
    s_data = np.loadtxt(out_dir / f'STATIONS_INCIDENT_FLUID', dtype=str)
    s_lat = s_data[:, 2].astype(float)
    s_lon = s_data[:, 3].astype(float)
    m.scatter(s_lon, s_lat, c='green', s=1, marker='v', linewidth=1)
    plt.savefig(out_dir / 'map_incident.png', bbox_inches='tight', pad_inches=0.0)

    ######################
    # wave extrapolation #
    ######################
    if args['array']['wave_extrapolation']:
        n_dist = int(np.ceil(dist_half_range / dist_delta)) + 1
        grid_dist_WE = np.linspace(0, dist_half_range, n_dist)
        grid_azim_WE = np.radians(np.linspace(0, 360, 2 * args_mesh['nu_to_use'] + 1)[:-1] * 1.)
        to_station_file(grid_depth_solid, 'solid', grid_dist_WE, grid_azim_WE, u_lat, u_lon, "EXTRAPOLATION")
        to_station_file(grid_depth_fluid, 'fluid', grid_dist_WE, grid_azim_WE, u_lat, u_lon, "EXTRAPOLATION")
        np.savetxt(out_dir / 'grid_dist_extrapolation.txt', grid_dist_WE)
        np.savetxt(out_dir / 'grid_azim_extrapolation.txt', grid_azim_WE)
        fout = open(out_dir / 'STATIONS_ARRAY_EXTRAPOLATION', 'w')
        fin = open(out_dir / 'STATIONS_ARRAY')
        fout.write(fin.read())
        fin = open(out_dir / 'STATIONS_EXTRAPOLATION_SOLID')
        fout.write(fin.read())
        fin = open(out_dir / 'STATIONS_EXTRAPOLATION_FLUID')
        fout.write(fin.read())
        fout.close()
        # plot
        plt.figure(dpi=200)
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=0, urcrnrlon=360, resolution='c')
        m.drawcoastlines()
        m.fillcontinents(color='ivory', lake_color='lightblue')
        m.drawmapboundary(fill_color='lightblue')
        m.scatter(e_lon, e_lat, c='red', s=50, marker='*')
        m.scatter(u_lon, u_lat, c='orange', s=200)
        m.scatter(u_lon_far, u_lat_far, c='blue', s=10, marker='+')
        m.scatter(u_lon_near, u_lat_near, c='blue', s=10, marker='+')
        s_data = np.loadtxt(out_dir / f'STATIONS_EXTRAPOLATION_FLUID', dtype=str)
        s_lat = s_data[:, 2].astype(float)
        s_lon = s_data[:, 3].astype(float)
        m.scatter(s_lon, s_lat, c='green', s=1, marker='v', linewidth=1)
        plt.savefig(out_dir / 'map_extrapolation.png', bbox_inches='tight', pad_inches=0.0)
