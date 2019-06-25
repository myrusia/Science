from ase.build import bulk
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from pymatgen.analysis.energy_models import IsingModel
from time import time
from copy import deepcopy
import pickle
import os


class MCModel:
    def __init__(self, config):
        self.config = config
        self.mode = self.config['mode']
        self.savefile = self.config['savefile']
        if self.mode == 'calculate':
            self.structure = self.make_simulation_grid(self.config['dim'])
            self.temperature = self.config['temperature'] * 8.6173303e-05  # Boltzmann constant in eV/K
            self.init_weight = np.array([x for x in self.config['initial_spin_weight']]) / sum(self.config['initial_spin_weight'])
            self.spins = None
            self.initialize_spins()
            self.history = [self.structure.copy()]
            self.size = len(self.structure)
            self.model = IsingModel(j=self.config['J'], max_radius=self.config['max_radius'])
            self.step = 0
            self.maxsteps = self.config['maxsteps']
            self.energy = self.model.get_energy(self.structure)
            self.get_magnetization()
            self.start_time = None
            self.saveevery = self.config['saveevery']
            self.force_stop = False
        elif self.mode == 'plot':
            with open(self.savefile, 'rb') as f:
                self.history = pickle.load(f)

    def make_simulation_grid(self, dim):
        atoms = bulk("Fe")
        structure = AseAtomsAdaptor().get_structure(atoms)
        structure.make_supercell(dim)
        return structure

    def initialize_spins(self):
        self.spins = np.random.choice([self.config['spin'], -self.config['spin']], len(self.structure), p=self.init_weight)
        self.structure.add_spin_by_site(self.spins)

    def make_step(self):
        self.step += 1
        test_step = deepcopy(self.structure.copy())
        idx = np.random.randint(self.size)
        self.spins[idx] *= -1
        test_step.add_spin_by_site(self.spins)
        energy_after_step = self.model.get_energy(test_step)
        change = energy_after_step - self.energy
        p = np.random.rand()
        print(np.exp((-1.0 * change)/self.temperature), p)
        if change < 0:
            self.structure = deepcopy(test_step.copy())
            self.energy = energy_after_step
            self.structure.add_spin_by_site(self.spins)
        elif np.exp((-1.0 * change)/self.temperature) > p:
            self.structure = deepcopy(test_step.copy())
            self.energy = energy_after_step
            self.structure.add_spin_by_site(self.spins)
        else:
            self.spins[idx] *= -1
        self.history.append(self.structure.copy())
        self.get_magnetization()
        if self.step % self.saveevery == 0:
            with open(self.savefile, 'wb') as f:
                pickle.dump(self.history, f)
            if os.path.exists('stop'):
                self.force_stop = True

    def get_magnetization(self):
        m = sum([getattr(self.structure[idx].specie, "spin") for idx in range(self.size)]) / self.size
        print(f'Step {self.step:05d}: magnetization {m:.8f} energy {self.energy:0.6f}')

    def run(self):
        self.start_time = time()
        while self.step < self.maxsteps:
            self.make_step()
            if self.force_stop:
                break
        runtime = time() - self.start_time
        print(f'Run time {runtime:.1f} seconds, {runtime/self.maxsteps:.3f} second per step')

    def plot(self):
        pass


if __name__ == '__main__':
    config_1 = {
        'dim': (11, 11, 1),
        'J': -0.00532566,
        'max_radius': 2.5,
        'spin': 2.2,
        'initial_spin_weight': (1.0, 0.0),
        'maxsteps': 500,
        'saveevery': 500,
        'savefile': 'config1.pkl',
        'temperature': 61.8,
        'mode': 'calculate'  # calculate, plot
    }
    config_2 = {
        'dim': (21, 21, 3),
        'J': -0.00532566,
        'max_radius': 10.0,
        'spin': 2.2,
        'initial_spin_weight': (0.5, 0.5),
        'maxsteps': 13250,
        'saveevery': 250,
        'savefile': 'config2.pkl',
        'temperature': 0.00000001,
        'mode': 'plot'  # calculate, plot
    }
    model = MCModel(config_1)
    if model.mode == 'calculate':
        model.run()
    if model.mode == 'plot':
        model.plot()

