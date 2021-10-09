""" Cellular automata: competition dynamics of mutator and stable cell populations

Jesús Antoñanzas
Maria Elzaurdi
Jordan Ortiz """

import argparse
import numpy as np

from sparseSimulation import sparse_simulation
from denseSimulation import dense_simulation
from utils import str2bool


def main():
    parser = argparse.ArgumentParser('Launch simulation', fromfile_prefix_chars='@')

    stateCaseHelp = """ Initial state of the cell population. 'sparse' only has a few stable cells and unstable cells, 
     while 'dense' has the board full of stable cells with a small cluster of already established mutator cells. """
    parser.add_argument('--stateCase', choices=['sparse', 'dense'], default='sparse', help=stateCaseHelp)

    totalPopDensityHelp = 'Total initial population density.'
    parser.add_argument('--totalPopDensity', type=float, default=0.2, help=totalPopDensityHelp)
    popDensityHelp = """ Population density of stable cells given the total pop. density, between 0 and 1. ' \
                     Mutator pop. will be (1 - arg)."""
    parser.add_argument('--stablePopDensity', type=float, default=0.5, help=popDensityHelp)

    parser.add_argument('--epochs', type=int, default=150, help='Number of generations to simulate.')
    parser.add_argument('--worldSize', type=int, default=50, help='Width and height of the grid.')

    parser.add_argument('--mutatorRR', type=float, default=2., help='Replication rate of mutator cells.')
    parser.add_argument('--stableRR', type=float, default=1., help='Replication rate of stable cells.')

    parser.add_argument('--diffusionRate', type=float, default=.01, help='Diffusion rate.')
    parser.add_argument('--damageProb', type=float, default=.02, help='Genetic hit rate.')

    parser.add_argument('--deathProb', type=float, default=.05, help='Probability of a damaged cell dying.')
    parser.add_argument('--mutationProb', type=float, default=.1, help='Probability of a damaged cell mutating.')

    sRProbHelp = 'Probability of a damaged stable cell repairing.'
    parser.add_argument('--stableRepairProb', type=float, default=.99, help=sRProbHelp)
    mRProbHelp = 'Probability of a damaged stable cell repairing.'
    parser.add_argument('--mutatorRepairProb', type=float, default=.1, help=mRProbHelp)

    parser.add_argument('--arrestProb', type=float, default=.3, help='Probability of cycle arrest during repair.')

    animationHelp = 'Create an animation depicting the evolution of the population?'
    parser.add_argument('--createAnimation', type=str2bool, default='True', help=animationHelp)
    parser.add_argument('--randomSeed', type=int, default=888, help='Random state for the simulation.')

    args = parser.parse_args()

    rng = np.random.default_rng(args.randomSeed)

    if args.stateCase == 'sparse':
        sparse_simulation(args, rng)
    else:
        dense_simulation(args, rng)

    return


if __name__ == "__main__":
    main()
