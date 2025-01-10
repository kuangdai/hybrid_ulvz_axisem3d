import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

from utils import replace_in_file

if __name__ == "__main__":
    run_name = sys.argv[1]
    args = json.load(open(f'inputs/{run_name}/args.json'))
    out_dir = Path(f'outputs/{run_name}/@@_exodus')
    out_dir.mkdir(parents=True, exist_ok=True)

    #################
    # generate mesh #
    #################
    # bm file
    r_cmb = 3480.
    ulvz_height = args['ulvz']['height']
    box_height = args['box']['height']
    # max velocity reduction
    layers = args['ulvz']['layers']
    dvs = np.min(np.array([layer.split(' ')[1] for layer in layers]).astype(float))
    vs_coefficients_ulvz = np.array([6.9254, 1.4672, -2.0834, 0.9783]) * (1 + dvs)
    bm_file = out_dir / f"mesh.bm"
    replace_in_file('templates/prem_iso_smooth_ulvz.bm',
                    {'__BOX_BOT__': r_cmb - box_height,
                     '__ULVZ_TOP__': r_cmb + ulvz_height,
                     '__BOX_TOP__': r_cmb + ulvz_height + box_height,
                     '__VS__': ' '.join(vs_coefficients_ulvz.astype(str))},
                    dest=bm_file)

    # mesh
    T = args['mesh']['period']
    epw = args['mesh']['epw']
    mesh_file = out_dir / f"mesh.e"
    cmd = f'python -m salvus_mesh_lite.interface AxiSEM ' \
          f'--basic.model {bm_file} --basic.period {T:.1f} ' \
          f'--output_filename {mesh_file} --overwrite_file ' \
          f'--advanced.elements_per_wavelength {epw}'
    os.system(cmd)

    ##########################
    # change VS back to PREM #
    ##########################
    # open
    f = h5py.File(mesh_file, 'r+')

    # read connect and coords
    connect = f['connect1'][:, :] - 1
    s = f['coordx'][:]
    z = f['coordy'][:]
    r = np.sqrt(np.power(s, 2) + np.power(z, 2))
    t = np.arccos(z / np.maximum(r, 1e-12))

    # compute center
    r_cen = np.ndarray((len(connect),))
    t_cen = np.ndarray((len(connect),))
    for iq in np.arange(len(connect)):
        r_cen[iq] = (r[connect[iq, 0]] + r[connect[iq, 1]] + r[connect[iq, 2]] + r[connect[iq, 3]]) / 4
        t_cen[iq] = (t[connect[iq, 0]] + t[connect[iq, 1]] + t[connect[iq, 2]] + t[connect[iq, 3]]) / 4

    # read vs
    vs = np.array([f['vals_elem_var17eb1'][0][:],
                   f['vals_elem_var18eb1'][0][:],
                   f['vals_elem_var19eb1'][0][:],
                   f['vals_elem_var20eb1'][0][:]])

    # change vs and find centers of ULVZ elements
    t_cen_ulvz = []
    n_upper = 0
    n_lower = 0
    for iq in np.arange(len(connect)):
        if r_cmb * 1e3 < r_cen[iq] < r_cmb * 1e3 + ulvz_height * 1e3:
            t_cen_ulvz.append(t_cen[iq])
            vs[:, iq] /= (1 + dvs)
        if (r_cmb - box_height) * 1e3 < r_cen[iq] < r_cmb * 1e3:
            n_lower += 1
        if (r_cmb + ulvz_height) * 1e3 < r_cen[iq] < (r_cmb + ulvz_height + box_height) * 1e3:
            n_upper += 1
    t_cen_ulvz = np.array(t_cen_ulvz)

    # copy to file
    f['vals_elem_var17eb1'][0, :] = vs[0, :]
    f['vals_elem_var18eb1'][0, :] = vs[1, :]
    f['vals_elem_var19eb1'][0, :] = vs[2, :]
    f['vals_elem_var20eb1'][0, :] = vs[3, :]
    NEX_surface = len(f['elem_ss1'])
    f.close()

    #########################
    # ULVZ and box geometry #
    #########################
    # half-length of ULVZ in degrees
    half_length_ulvz_deg = args['ulvz']['radius_deg']

    # determine NEX on CMB
    theta_tol = np.pi / NEX_surface / 10.
    NEX = len(np.unique(np.round(t_cen_ulvz / theta_tol).astype(int)))
    assert len(t_cen_ulvz) % NEX == 0, f'NEX={NEX}, NE_ULVZ={len(t_cen_ulvz)}'
    assert n_lower == NEX and n_upper == NEX, f'Only one layer is needed for injection boundary.'
    n_elem_layer_ulvz = len(t_cen_ulvz) // NEX

    # box range
    NEX_U = int(np.ceil(NEX / 180 * half_length_ulvz_deg))
    rb_min = (r_cmb - box_height) * 1e3
    rb_max = (r_cmb + ulvz_height + box_height) * 1e3
    tb_max = np.pi / NEX * (NEX_U + 1)

    # ulvz range
    ru_min = r_cmb * 1e3
    ru_max = (r_cmb + ulvz_height) * 1e3
    tu_max = np.pi / NEX * NEX_U

    # find quads
    quad_bound = []
    quad_ulvz = []
    for iq in np.arange(len(connect)):
        if rb_min < r_cen[iq] < rb_max and t_cen[iq] < tb_max:
            if ru_min < r_cen[iq] < ru_max and t_cen[iq] < tu_max:
                quad_ulvz.append(iq)
            else:
                quad_bound.append(iq)
    print('# quad in ulvz = %d' % len(quad_ulvz))
    print('# quad on boundary = %d' % len(quad_bound))
    print('# NEX_U = %d' % NEX_U)
    print('# elem layers in ULVZ = %d' % n_elem_layer_ulvz)

    print("quads on boundary:")
    print(quad_bound)
    print()

    # Equation (4)
    # 6. is reduced Vs on CMB
    nr_needed = 2 * np.pi * r_cmb * np.sin(tb_max) / (6. * T) * epw * 4
    nu_needed = nr_needed / 2
    nu_to_use = int(np.ceil(nu_needed / 10) * 10)
    print('nc needed = %.18f' % nu_needed)
    print('nc to use = %d\n' % nu_to_use)

    print('box r min = %d' % rb_min)
    print('box r max = %d' % rb_max)
    print('box t max = %.18f' % tb_max)

    print('ulvz t max = %.18f' % (np.pi / NEX * NEX_U))
    print('ulvz t max degree = %.18f' % (180 / NEX * NEX_U))

    r0 = 1e10
    r1 = 0
    t1 = 0
    for iq in quad_bound:
        r0 = min(r0, r[connect[iq, 0]], r[connect[iq, 1]], r[connect[iq, 2]], r[connect[iq, 3]])
        r1 = max(r1, r[connect[iq, 0]], r[connect[iq, 1]], r[connect[iq, 2]], r[connect[iq, 3]])
        t1 = max(t1, t[connect[iq, 0]], t[connect[iq, 1]], t[connect[iq, 2]], t[connect[iq, 3]])

    print('\nCheck their values in mesh:')
    print('box r min', r0)
    print('box r max', r1)
    print('box t max', t1)

    # write derived info for later use
    args_mesh = {'NEX': NEX,
                 'NEX_U': NEX_U,
                 'n_elem_layer_ulvz': n_elem_layer_ulvz,
                 'nu_to_use': nu_to_use,
                 'rb_min': rb_min,
                 'rb_max': rb_max,
                 'tb_max': tb_max}
    json.dump(args_mesh, open(out_dir / "args.json", 'w'))
    # quads on boundary for later use
    np.savetxt(out_dir / 'quad_bound.txt', np.array(quad_bound), fmt='%d')

    #############
    # demo mesh #
    #############
    os.system('cp %s %s' % (mesh_file, str(mesh_file)[:-2] + '__demo.e'))
    f = h5py.File(str(mesh_file)[:-2] + '__demo.e', 'r+')
    # read fluid
    fluid = np.array(f['vals_elem_var25eb1'][0][:])
    # change fluid
    # ulvz
    for iq in quad_ulvz:
        fluid[iq] = -0.5
    # bound
    for iq in quad_bound:
        fluid[iq] = 0.5
    # copy to file
    f['vals_elem_var25eb1'][0, :] = fluid[:]
    f.close()
