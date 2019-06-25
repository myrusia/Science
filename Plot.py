import matplotlib.pyplot as plt
from config import config
from collect_structures import data_from_pickle
from objective_function import MultiIsingModel
import matplotlib.patches as mpatches
from pymatgen import Structure


def plot_energy_prediction():
    x = [-9.28649456e+01,  4.65876504e-02,  5.58262486e-02,  3.54432955e-02, 3.03638793e-02, -5.52253644e-03,  5.00219072e-03]
    with open(f'{config["outcars_folder"]}/dofit', 'r') as f:
        dofit = f.read().split()
    all_structures = data_from_pickle()[:]
    fm = Structure.from_file(f'{config["outcars_folder"]}/fm.vasp')
    model = MultiIsingModel(x, config['cutoff_radius'], fm)
    fm_model_energy = model.get_energy() + x[0]
    vasp_energies = []
    model_energies = []
    colors = []
    print('Structure Error')
    for struct in all_structures:
        model = MultiIsingModel(x, config['cutoff_radius'], struct)
        vasp_energies.append(struct.vasp_energy)
        model_energies.append(model.get_energy() + x[0])
        if struct.name in dofit:
            colors.append('b')
        else:
            colors.append('r')
        print(f'{struct.name} {vasp_energies[-1] - model_energies[-1]}')
    emin = min(vasp_energies + model_energies)
    emax = max(vasp_energies + model_energies)
    plt.figure(figsize=(20, 20), dpi=600)
    fig, ax = plt.subplots()
    plt.plot([emin, emax], [emin, emax], 'k')
    patch1 = mpatches.Patch(color='#0000ff', label='Used in fit')
    patch2 = mpatches.Patch(color='#ff0000', label='Predicted')
    patch3 = mpatches.Patch(color='#000000', label='Ferromagnetic')
    plt.scatter(vasp_energies, model_energies, color=colors, edgecolors=['k'] * len(vasp_energies), alpha=0.66)
    plt.plot(-95.10789834, fm_model_energy, 'ok', label='FM')
    all_handles = (patch1, patch2, patch3)
    leg = ax.legend(handles=all_handles)
    plt.xlabel('vasp energy')
    plt.ylabel('model energy')
    plt.savefig(f'{config["outcars_folder"]}.png')


if __name__ == '__main__':
    plot_energy_prediction()
