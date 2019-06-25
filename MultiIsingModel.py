from ase.build import sort
from ase.io.vasp import write_vasp, read_vasp
from pymatgen import Structure
from pymatgen.analysis.energy_models import EnergyModel, IsingModel
import matplotlib.pyplot as plt
import numpy as np


class MultiIsingModel:
    """
     A little complicated Ising model, with j series.

     Args:
         j (list): The interaction parameter. E = J * spin1 * spin2 + ...
         radius (float): max_radius for the interaction.
     """
    def __init__(self, j, max_radius, structure, spins, magnetic_element='Fe'):
        self.j = j
        self.max_radius = max_radius
        self.magnetic_lattice = structure.copy()
        self.magnetic_lattice.remove_species([element for element in self.magnetic_lattice.species if element.value !=magnetic_element])
        self.magnetic_lattice.add_spin_by_site(spins)
        self.distances = np.unique(np.round([dist for site in self.magnetic_lattice.get_all_neighbors(r=self.max_radius) for nn, dist in site], decimals=6))

    def recommended_j(self):
        print(self.distances)
        print(self.j[0] * self.distances[0] / self.distances)

    def get_energy(self):
        all_nn = self.magnetic_lattice.get_all_neighbors(r=self.max_radius)
        energy = 0
        for i, nn in enumerate(all_nn):
            s1 = getattr(self.magnetic_lattice[i].specie, "spin", 0)
            for site, dist in nn:
                idx = int(np.where(np.isclose(self.distances, dist))[0])
                energy -= self.j[idx] * s1 * getattr(site.specie, "spin", 0)
        return energy

#atoms = read_vasp('POSCAR2')
#write_vasp('POSCAR', sort(atoms), vasp5=True, direct=True, label=atoms.get_chemical_formula())


if __name__ == '__main__':
    """
    J = [0.00240319, 0.00224361, 0.00146472, 0.00125072, 0.00122649]
    model_uuu = MultiIsingModel(J, 5.0, Structure.from_file('2FEPOSCAR'), [3, 3])
    model_uuu.recommended_j()
    model_uud = MultiIsingModel(J, 5.0, Structure.from_file('2FEPOSCAR'), [3, -3])
    uuu = model_uuu.get_energy()
    uud = model_uud.get_energy()
    print(uuu - uud)
    """
    J = -0.03489
    uu = Structure.from_file('2FEPOSCAR')
    ud = uu.copy()
    uu.add_spin_by_site([3.291, 3.291])
    ud.add_spin_by_site([3.285, -3.285])
    model = IsingModel(J, 5.0)
    print(model.get_energy(uu) - model.get_energy(ud))
    uuu = Structure.from_file('POSCAR3')
    uud = uu.copy()
    uuu.add_spin_by_site([3, 3, 3])
    uud.add_spin_by_site([3, 3,-3])
    model = IsingModel(J, 5.0)
    print(model.get_energy(uuu) - model.get_energy(uud))