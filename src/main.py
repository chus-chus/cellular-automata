""" Cellular automata: competition dynamics of mutator and stable cell populations

Jesús Antoñanzas
Maria Elzaurdi
Jordan Ortiz """

import argparse
import numpy as np

from sparseSimulation import sparse_simulation


def main():
    parser = argparse.ArgumentParser('Launch simulation')
    popTypeHelp = """ Initial state of the cell population. 'sparse' only has a few stable cells and unstable cells, 
     while 'dense' has the board full of stable cells with a small cluster of already established mutator cells. """
    parser.add_argument('stateCase', cases=['sparse', 'dense'], help=popTypeHelp)
    popDensityHelp = """ If 'stateCase' = 'dense', it indicates the proportion of mutated cells. If 'stateCase' = 'sparse',
     it indicates... """
    parser.add_argument('den', type=float, default=0.01)
    parser.add_argument('epochs', default=500)
    parser.add_argument('randomSeed', default=888)
    # diffusion rate

    args = parser.parse_args()

    rng = np.random.default_rng(args['randomSeed'])

    if args['stateCase'] == 'sparse':
        sparse_simulation(args, rng)
    else:
        raise NotImplementedError

    return


if __name__ == "__main__":
    main()