from config import config
from collect_structures import data_from_pickle
from pymatgen.analysis.energy_models import IsingModel
import numpy as np


class MultiIsingModel:
    """
         A little complicated Ising model, with j series.

         Args:
             j (list): The interaction parameter. E = J * spin1 * spin2 + ...
             radius (float): max_radius for the interaction.
         """

    def __init__(self, j, max_radius, structure, magnetic_element='Fe'):
        self.j = j[1:]
        self.max_radius = max_radius
        self.magnetic_lattice = structure.copy()
        """
        self.magnetic_lattice.remove_species(
            [element for element in self.magnetic_lattice.species if element.value != magnetic_element])
        self.magnetic_lattice.add_spin_by_site(spins)
        """
        self.distances = np.unique(
            np.round([dist for site in self.magnetic_lattice.get_all_neighbors(r=self.max_radius) for nn, dist in site],
                     decimals=3))

    def recommended_j(self):
        print(self.distances)
        print(self.j[0] * self.distances[0] / self.distances)

    def get_energy(self):
        all_nn = self.magnetic_lattice.get_all_neighbors(r=self.max_radius)
        energy = 0
        for i, nn in enumerate(all_nn):
            s1 = getattr(self.magnetic_lattice[i].specie, "spin", 0)
            for site, dist in nn:
                idx = int(np.where(np.isclose(self.distances, np.round(dist, decimals=3)))[0])
                energy -= 0.5 * self.j[idx] * s1 * getattr(site.specie, "spin", 0)
        return energy


def objective_function_Ising(x, structures):
    """
    :param x: list of E0 and J
    :param structures: vasp calculated structures
    :return: error = (property_predicted - property_target) ** 2 / delta_tolerance
    """

    model = IsingModel(x[1], config['cutoff_radius'])
    error = sum([(model.get_energy(st) + x[0] - st.vasp_energy) ** 2 for st in structures]) / config['delta_tolerance']
    return error


def objective_function_multi(x, structures):
    error = 0
    for st in structures:
        model = MultiIsingModel(x, config['cutoff_radius'], st)
        error += (model.get_energy() + x[0] - st.vasp_energy) ** 2 / config['delta_tolerance']
    return error


class DeMulti:
    def __init__(self, bounds):
        self.structures = data_from_pickle()[:]
        with open('Fe12Sn4/dofit', 'r') as f:
            dofit = f.read().split()
        structures = [s for s in self.structures if s.name in dofit]
        self.structures = structures
        self.bnds = bounds
        self.fm = self.structures[0].copy()
        self.fm.add_spin_by_site(config['spin_values'])

    def objective_function_IsingMulti_de(self, x):
        error = 0
        model_fm = MultiIsingModel(x, config['cutoff_radius'], self.fm)
        #E0 = config['FM_energy'] - model_fm.get_energy()
        E0 = x[0]
        for st in self.structures:
            model = MultiIsingModel(x, config['cutoff_radius'], st)
            #model.recommended_j()
            #raise  RuntimeError
            error += (model.get_energy() + E0 - st.vasp_energy) ** 2 / config['delta_tolerance']
        return error

    def run(self):
        from scipy.optimize import differential_evolution, least_squares
        from scipy.optimize.optimize import fmin_bfgs
        #res = fmin_bfgs(self.objective_function_IsingMulti_de, np.array([0.02, 0.02, 0.000, 0.000, 0.000]))
        #x0 = np.array(np.array([0.02, 0.02, 0.000]))
        x0 = np.array(config['x0'])
        x0 = np.array([-90, 1, -1, 1, -1, 1, -1])
        #x0 = np.array(np.array([0.02, 0.02, 0.000, 0.000, 0.000]))
        res = least_squares(self.objective_function_IsingMulti_de, x0, jac='3-point', bounds=self.bnds, ftol=1e-12,
                            xtol=1e-12, gtol=1e-12, max_nfev=100000, verbose=2)
        #res = differential_evolution(self.objective_function_IsingMulti_de, self.bnds, strategy='randtobest1bin', popsize=1000, disp=True, workers=6)
        print(res)

    def print_energy(self, x):
        model_fm = MultiIsingModel(x, config['cutoff_radius'], self.fm)
        E0 = config['FM_energy'] - model_fm.get_energy()
        print(E0)
        for st in self.structures:
            model = MultiIsingModel(x, config['cutoff_radius'], st)
            E0 = st.vasp_energy - model.get_energy()
            print(E0)


if __name__ == '__main__':

    """
    structures = data_from_pickle()[:]
    error = objective_function_Ising([config['E0'], config['guess_J']], structures)
    print(error)

    probe_j = [-91.9098722] + [-0.02372642, -0.02264993, -0.01941252, -0.01638308, -0.01473944, -0.01357223,
 -0.012308, -0.01214998]
    error = objective_function_multi(probe_j, structures)
    print(error)

    from scipy.optimize import minimize
    import numpy as np

    bnds = ((-92, -88), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001), (-0.1, -0.0001))
    res = minimize(objective_function_multi, np.array(probe_j), args=structures, bounds=bnds, options={'disp': True})
    print(res)

    error = objective_function_multi(res.x, structures)
    print(error)
    """

    """
    from scipy.optimize import minimize
    import numpy as np
    bnds = ((-92, -88), (-2.0, 0.0))
    res = minimize(objective_function_Ising, np.array([config['E0'], config['guess_J']]), args=structures, bounds=bnds, options={'disp': True})
    print(res)
    error = objective_function_Ising([res.x[0], res.x[1]], structures)
    print(error)
    model = IsingModel(res.x[1], config['cutoff_radius'])
    for st in data_from_pickle()[:]:
        print(st.vasp_energy, model.get_energy(st) + res.x[0], model.get_energy(st) + res.x[0] - st.vasp_energy)
    import matplotlib.pyplot as plt
    plt.plot([st.vasp_energy for st in data_from_pickle()[:]], [model.get_energy(st) + res.x[0] for st in data_from_pickle()[:]], 'o')
    plt.xlabel('Vasp energy')
    plt.ylabel('Ising model energy')
    plt.savefig('Fe12_Ising model.png')
    """


    bounds_all = config['bounds']
    bounds_all = [(-100, -1, -1, -1, -1, -1, -1), (-90, 1, 1, 1, 1, 1, 1)]

    de_all = DeMulti(bounds_all)

    de_all.run()





