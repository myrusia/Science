from config import config
import os
import numpy as np
import glob


def make_fit_list(N, rewrite=False):
    outcars = [x.split(f'{config["outcars_folder"]}/')[1] for x in sorted(glob.glob(f'{config["outcars_folder"]}/OUTCAR-*'))]
    assert len(outcars) > N
    selected = np.random.choice(outcars, N, replace=False)
    output = ''
    for outcar in sorted(selected):
        output += outcar + '\n'
    if rewrite or not os.path.exists(f'{config["outcars_folder"]}/dofit'):
        with open(f'{config["outcars_folder"]}/dofit', 'w') as f:
            f.write(output)


if __name__ == '__main__':
    make_fit_list(80)
