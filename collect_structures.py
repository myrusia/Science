import glob
from pymatgen.io.vasp import Poscar, Outcar
from config import config
from pymatgen import Structure
import numpy as np
import os
import pickle
import logging


def data_from_pickle():
    with open(config['data_load_file'], 'rb') as f:
        different_structures = pickle.load(f)
    return different_structures


def initLogging():
    logLevel = logging.INFO
    logging.basicConfig(level=logLevel)
    log_file = config['read_structures_log']
    fileHandler = logging.FileHandler(log_file)
    logging.root.addHandler(fileHandler)


def is_unique(array, value, tolerance):
    assert tolerance > 0
    if not len(array):
        return True
    else:
        return np.all(np.abs(np.array(array) - value) > tolerance)


def data_from_outcars_in_folder(folder):
    initLogging()
    if not os.path.exists(config['data_load_file']):
        different_energies = []
        different_structures = []
        poscar = Poscar.from_file(folder+'/POSCAR')
        logging.info('Assigned Spin values Energy Unique energy FilePath')
        for outcar_file in sorted(glob.glob(folder+'/OUTCAR*')):
            structure = poscar.structure.copy()
            magnetic_element_index = [n for n, e in enumerate(structure) if e.species_string in config['magnetic_elements']]
            magnetic_substructure = Structure.from_sites([e for e in structure if e.species_string in config['magnetic_elements']])
            outcar = Outcar(outcar_file)
            try:
                magnetization = np.array([entry['tot'] for entry in outcar.magnetization])[magnetic_element_index]
            except IndexError:
                logging.info(f'Bad OUTCAR: {outcar_file}')
                continue
            spins = np.zeros_like(magnetization)
            for i, m in enumerate(magnetization):
                spins[i] = np.sign(m) * config['spin_values'][i]
                #if abs(m) > config['spin_threshold']:
                #    spins[i] = config['spin_value'] * np.sign(m)
            magnetic_substructure.add_spin_by_site(spins)
            magnetic_substructure.vasp_energy = outcar.final_energy
            magnetic_substructure.name = outcar_file.split('/')[-1]
            isUnique = False
            if is_unique(different_energies, magnetic_substructure.vasp_energy, 1e-15):
                different_energies.append(magnetic_substructure.vasp_energy)
                isUnique = True
                different_structures.append(magnetic_substructure)
            #energy_checked = np.round(magnetic_substructure.vasp_energy, decimals=3)
            #if not np.in1d(energy_checked, different_energies).any():
            #    different_energies.append(energy_checked)
            #    isUnique = True
            #    different_structures.append(magnetic_substructure)
            logging.info((spins, magnetic_substructure.vasp_energy, isUnique, outcar_file))
        logging.info(f'Data reading complete, saving to {config["data_load_file"]}')
        with open(config['data_load_file'], 'wb') as f:
            pickle.dump(different_structures, f)
        logging.info('Data load complete. Rerun script')
        import sys
        sys.exit()
    else:
        with open(config['data_load_file'], 'rb') as f:
            different_structures = pickle.load(f)
            for s in different_structures:
                logging.info(s.vasp_energy)
    return different_structures


if __name__ == '__main__':

    structures = data_from_outcars_in_folder(config['outcars_folder'])
    print(len(structures))
