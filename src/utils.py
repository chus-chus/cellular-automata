""" Utility variables, methods and classes """

# celltypes
import argparse
from copy import deepcopy
from typing import List

import matplotlib.pyplot as plt
import numpy as np

cellTypes = {'stable': 0, 'mutator': 1}

# cellstates
cellStates = {'healthy': 0, 'repairing': 1, 'mutated': 2, 'damaged': 3}

# grid colors
__WHITE = [255, 255, 255]
__BLACK = [0, 0, 0]
__RED = [230, 80, 50]
__DARK_RED = [150, 30, 10]
__GREEN = [120, 230, 60]
__DARK_GREEN = [40, 100, 10]
__BLUE = [100, 190, 210]
# __DARK_BLUE = [0, 0, 255]

colors = {cellTypes['stable']: {cellStates['healthy']: __GREEN, cellStates['repairing']: __DARK_GREEN,
                                cellStates['mutated']: __DARK_GREEN, cellStates['damaged']: __BLUE},
          cellTypes['mutator']: {cellStates['healthy']: __RED, cellStates['repairing']: __DARK_RED,
                                 cellStates['mutated']: __DARK_RED, cellStates['damaged']: __BLUE},
          'empty': __WHITE}


class Cell(object):
    """ Represents a basic cell, of stable or mutant tyes and with healthy, repairing or mutated states. """

    def __init__(self,
                 cellType: str = 'stable',
                 cellState: str = 'healthy',
                 replicationRate: float = 1,
                 repairProb: float = 0.95):

        self.type = cellTypes[cellType]
        self.state = cellStates[cellState]
        self.replicationRate = replicationRate
        self.repairProb = repairProb


def init_world(world: List[List[Cell]], totalDensity: float, stableDensity: float,
               args, rng: np.random.default_rng) -> None:
    """ Initialise grid with stable or mutator cells. There are density * len(world) cells placed, half of them
     being stable and half of them being mutator. """

    rows = rng.choice(len(world), size=round(len(world) * len(world) * totalDensity), replace=True)
    cols = rng.choice(len(world), size=round(len(world) * len(world) * totalDensity), replace=True)

    for i in range(len(rows)):
        if rng.uniform() < stableDensity:
            world[rows[i]][cols[i]] = Cell('stable', replicationRate=args.stableRR,
                                           repairProb=args.stableRepairProb)
        else:
            world[rows[i]][cols[i]] = Cell('mutator', replicationRate=args.mutatorRR,
                                           repairProb=args.mutatorRepairProb)

    return


def precompute_moore_domain(world):
    mooreDomain = deepcopy(world)
    for i in range(len(mooreDomain)):
        for j in range(len(mooreDomain)):
            domain = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                      [i,     j - 1],             [i,     j + 1],
                      [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]

            # parse invalid indexes (edges)
            for t in range(len(domain)):
                if domain[t][0] >= len(world) or domain[t][0] < 0:
                    domain[t][0] = None
                if domain[t][1] >= len(world) or domain[t][1] < 0:
                    domain[t][1] = None

            # delete invalid indexes (edges)
            mooreDomain[i][j] = [index for index in domain if index[0] is not None and index[1] is not None]

    return mooreDomain


def choose_moore_domain(row, col, mooreDomain, n, rng) -> np.ndarray:
    """ Choose n indexes w/o replacement in the moore domain centered around i, j. Does not consider grid edges.

    :return: a list of indexes as lists. """

    return rng.choice(mooreDomain[row][col], size=n, replace=False)


def color_world(world: List[List[Cell]]):
    # translate the world into a np.ndarray with the specified color values
    coloredWorld = [[None for _ in range(len(world))] for _ in range(len(world))]

    for i in range(len(world)):
        for j in range(len(world)):
            cell = world[i][j]
            if cell is None:
                coloredWorld[i][j] = colors['empty']
            else:
                coloredWorld[i][j] = colors[cell.type][cell.state]

    return coloredWorld


def clear_plots():
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def gen_barplot_two_lines(x, y1, y2, line1Legend, line2Legend, xlabel, ylabel, plotTitle,
                          filepath, filename, dpi=300):
    """ linechart with two lines.  """

    plt.plot(x, y1, label=line1Legend)
    plt.plot(x, y2, label=line2Legend)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(plotTitle)
    plt.legend()
    plt.savefig(f'{filepath}/{filename}.png', dpi=dpi)


def gen_save_plots(epochs, stats, path):
    clear_plots()

    x = list(range(epochs))

    # total population numbers
    gen_barplot_two_lines(x=x, y1=stats['stable']['total'],
                          y2=stats['mutator']['total'],
                          xlabel='Epoch', ylabel='Number of cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Total population evolution',
                          filepath=path, filename='total_population_evolution')

    clear_plots()
    # healthy cells
    gen_barplot_two_lines(x=x, y1=stats['stable']['healthy'],
                          y2=stats['mutator']['healthy'],
                          xlabel='Epoch', ylabel='Number of healthy cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Healthy population evolution',
                          filepath=path, filename='healthy_population_evolution')
    clear_plots()
    # mutated cells
    gen_barplot_two_lines(x=x, y1=stats['stable']['mutated'],
                          y2=stats['mutator']['mutated'],
                          xlabel='Epoch', ylabel='Number of mutated cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Mutated population evolution',
                          filepath=path, filename='mutated_population_evolution')
    # damaged cells
    clear_plots()
    gen_barplot_two_lines(x=x, y1=stats['stable']['damaged'],
                          y2=stats['mutator']['damaged'],
                          xlabel='Epoch', ylabel='Number of damaged cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Damaged population evolution',
                          filepath=path, filename='damaged_population_evolution')


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def parse_all_params(args):
    params = str(vars(args)).split(',')  # ["{'a': 1", ..., "'b': 2}",]
    params = [p[:p.find(':')] + '\n' + p[p.find(':'):] for p in params]
    params = ['--' + param for param in params]  # ["--{'a': \n 1", ..., "--'b': \n 2}"]
    params = '\n'.join(params)  # "--{'a': \n 1 \n ... \n --'b': \n 2}"
    params = params.replace('{', '').replace('}', '')  # "--'a': \n 1 \n ... \n --'b': \n 2"
    params = params.replace("'", '').replace(':', '')  # "--a \n 1 \n ... \n --b \n 2"
    return params.replace(' ', '')