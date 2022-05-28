import sys
import h5py
import numpy as np
import lib_topology as es


def package(flow_filenames, hdf5_filename):
    with h5py.File(hdf5_filename, 'w') as h5file:
        flows = h5file.require_group('flows')
        plaqs = h5file.require_group('plaqs')
        byensemble = h5file.require_group('byensemble')
        for fname in flow_filenames:
            print('Processing', fname, end='\r')
            N, L, beta, flow_data = es.topo_load_raw_data(fname)

            flow_ds = flows.create_dataset(name=f'{N}_{L}_{beta}',
                                           data=flow_data)

            plaq_fname = fname + '_plaq'
            plaq_data = np.genfromtxt(plaq_fname)
            plaq_ds = plaqs.create_dataset(name=f'{N}_{L}_{beta}',
                                           data=plaq_data)

            for ds in plaq_ds, flow_ds:
                ds.attrs['beta'] = beta
                ds.attrs['N'] = N
                ds.attrs['L'] = L

            ensemble_group = byensemble.require_group(f'{N}_{L}_{beta}')
            ensemble_group['flows'] = flow_ds
            ensemble_group['plaqs'] = plaq_ds


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('flow_filenames', nargs='+')
    parser.add_argument('--hdf5_filename', default='datapackage.h5')
    args = parser.parse_args()

    package(args.flow_filenames, args.hdf5_filename)


if __name__ == '__main__':
    main()
