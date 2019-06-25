import numpy as np
from pymatgen import Structure
from scipy.constants import e


def write_vampire(structure, J, max_radius, S):
    SCALEFACTOR = 100
    distances = np.unique(
        np.round([dist for site in structure.get_all_neighbors(r=max_radius) for nn, dist in site],
                 decimals=4))
    print(distances)
    output_head = '# Vampire input structure file. Unit cell size\n'
    output_head += '{} {} {}\n'.format(*structure.lattice.abc)
    output_head += '# Unit cell lattice vectors:\n'
    for i, v in enumerate(structure.lattice.abc):
        output_head += '{} {} {}\n'.format(*structure.lattice.matrix[i, :] / v)
    output_head += '# Atoms\n'
    output_head += '{} {} \n'.format(len(structure), len(structure.composition))
    for i, site in enumerate(structure):
        output_head += '{} {} {} {}\n'.format(i, *site.frac_coords)
    output_head += '# Interactions\n'
    output_tail = ''
    IID = 0
    """
    for i, site in enumerate(structure):
        for j, site in enumerate(structure):
            if i != j:
                dist = structure[i].distance(structure[j])
                idx = np.where(np.isclose(distances, np.round(dist, decimals=4)))[0]
                if idx:
                    print(idx, dist)
                    output_tail += '{} {} {} {} {} {} {}\n'.format(IID, i, j, 0, 0, 0, *J[idx] * e * SCALEFACTOR)
                    IID += 1
    """
    for n, site in enumerate(structure):
        unique_d = []
        interactions = structure.get_neighbors(site, max_radius, include_index=True)
        ijlist = []
        for remote_site, dist, origin_index in interactions:
            supercell_index = (remote_site.coords - structure[origin_index].coords) / structure.lattice.abc
            if all([(np.isclose(abs(a), 1.0, atol=0.01) or np.isclose(a, 0.0, atol=0.01)) for a in supercell_index]):
                idx = np.where(np.isclose(distances, np.round(dist, decimals=4)))[0]
                if idx not in unique_d:
                    unique_d.append(int(idx))
                if (n, origin_index) not in ijlist:

                    output_tail += '{} {} {} {} {} {} {}\n'.format(IID, n, origin_index, int(supercell_index[0]), int(supercell_index[1]), int(supercell_index[2]), S * S *J[idx] * e * SCALEFACTOR)
                    ijlist.append((n, origin_index))
                    IID += 1
                    output_tail += '{} {} {} {} {} {} {}\n'.format(IID, origin_index, n, -int(supercell_index[0]), -int(supercell_index[1]), -int(supercell_index[2]), S * S *J[idx] * e * SCALEFACTOR)
                    ijlist.append((origin_index, n))
                    IID += 1
    print(unique_d)
    output_zhong = f'{IID} 0\n'
    #print(output)
    with open('Fe12.ucf', 'w') as f:
        f.write(output_head + output_zhong + output_tail)


def get_all_interactions(structure, J, max_radius):
    distances = np.unique(
        np.round([dist for site in structure.get_all_neighbors(r=max_radius) for nn, dist in site],
                 decimals=4))
    for site in structure:
        interactions = structure.get_neighbors(site, max_radius, include_index=True)
        for remote_site, dist, origin_index in interactions:
            supercell_index = (remote_site.coords - structure[origin_index].coords) / structure.lattice.abc
            if all([(np.isclose(abs(a), 1.0) or np.isclose(a, 0.0)) for a in supercell_index]):
                idx = np.where(np.isclose(distances, np.round(dist, decimals=4)))[0]
                print(*supercell_index, np.round(dist, decimals=4), J[idx])


            #print((remote_site.coords - structure[origin_index].coords) / structure.lattice.abc, dist, origin_index)
        break


if __name__ == '__main__':
    structure = Structure.from_str("""Fe3Sn
1.0
        5.5219998360         0.0000000000         0.0000000000
       -2.7610005548         4.7821917700         0.0000000000
        0.0000000000         0.0000000000         8.6879997253
   Fe
   12
Direct
     0.149999996         0.300000009         0.125000000
     0.149999996         0.300000009         0.625000000
     0.849999988         0.699999988         0.375000000
     0.849999988         0.699999988         0.875000000
     0.699999983         0.850000042         0.125000000
     0.699999983         0.850000042         0.625000000
     0.300000004         0.150000004         0.375000000
     0.300000004         0.150000004         0.875000000
     0.149999996         0.850000042         0.125000000
     0.149999996         0.850000042         0.625000000
     0.849999970         0.150000004         0.375000000
     0.849999970         0.150000004         0.875000000""", fmt='poscar')
    J = np.array([1.79794272e-04,  2.96712523e-03,  7.18270723e-03,
        7.59394248e-04,  3.04796980e-04,  0.00000000e+00,  2.46742713e-04,
        6.25025769e-04])
    write_vampire(structure, J, 5.0, 3.0)


    #get_all_interactions(structure, J, 5.0)
