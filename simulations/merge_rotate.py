import json

from seismogram import merge_nc, rotate_nc

if __name__ == "__main__":
    merge_nc(f'./@@_solve_prem/output/stations/ARRAY')
    merge_nc(f'./@@_solve_1d/output/stations/ARRAY')
    merge_nc(f'./@@_solve_2d/output/stations/ARRAY')
    args_stations = json.load(open(f'./@@_stations/args.json'))
    rotate_nc(f'./@@_solve_3d/output/stations/ARRAY',
              f'./@@_stations/STATIONS_ARRAY',
              e_lat=args_stations['event_lat'],
              e_lon=args_stations['event_lon'],
              u_lat=args_stations['ulvz_lat'],
              u_lon=args_stations['ulvz_lon'])
