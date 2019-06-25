config = {
    'outcars_folder': 'Fe12Sn4',
    'magnetic_elements': ['Fe'],
    'spin_values': [1.0] * 12,
    'cutoff_radius': 4.3,
    'delta_tolerance': 1.,
    'bounds': [(-100, -0.1, -0.1, -0.1, -0.1, -0.1), (90, 0.1, 0.1, 0.1, 0.1, 0.1)],
    'x0': [-9.29128027e+01, 4.48035228e-02, 5.40010960e-02, 2.98253280e-02, 2.79007547e-02, -4.46021585e-03]
}

config['read_structures_log'] = config['outcars_folder'] + '_outcars.log'
config['data_load_file'] = config['outcars_folder'] + '.pkl'
