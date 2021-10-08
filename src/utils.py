""" Utility variables, methods and classes """

# celltypes
import pathlib
from typing import List

import numpy as np
import matplotlib.pyplot as plt

cellTypes = {'stable': 0, 'mutator': 1}

# cellstates
cellStates = {'healthy': 0, 'repairing': 1, 'mutated': 2, 'damaged': 3}

# grid colors
__WHITE = [255, 255, 255]
__BLACK = [0, 0, 0]
__RED = [255, 0, 0]
__DARK_RED = [150, 0, 0]
__GREEN = [0, 255, 0]
__DARK_GREEN = [0, 150, 0]
__BLUE = [0, 0, 150]
__DARK_BLUE = [0, 0, 255]

colors = {cellTypes['stable']: {cellStates['healthy']: __GREEN, cellStates['repairing']: __GREEN,
                                cellStates['mutated']: __DARK_GREEN, cellStates['damaged']: __BLUE},
          cellTypes['mutator']: {cellStates['healthy']: __RED, cellStates['repairing']: __RED,
                                 cellStates['mutated']: __DARK_RED, cellStates['damaged']: __DARK_BLUE},
          'empty': __WHITE}


class Cell(object):
    """ Represents a basic cell, of stable or mutant tyes and with healthy, repairing or mutated states. """

    def __init__(self, cellType: str = 'stable',
                 cellState: str = 'healthy',
                 replicationRate: float = 1,
                 damageProb: float = 0.005,
                 repairProb: float = 0.05,
                 deathProb: float = 0.05,
                 mutationProb: float = 0.1):

        self.type = cellTypes[cellType]
        self.state = cellStates[cellState]
        self.replicationRate = replicationRate
        self.damageProb = damageProb
        self.repairProb = repairProb
        self.deathProb = deathProb
        self.mutationProb = mutationProb


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

    # total population numberes
    gen_barplot_two_lines(x=x, y1=stats['stable']['total'],
                          y2=stats['mutator']['total'],
                          xlabel='Epoch', ylabel='Number of cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Pop. evolution. Params go here',
                          filepath=path, filename='total_population_evolution')

    clear_plots()
    # healthy cells
    gen_barplot_two_lines(x=x, y1=stats['stable']['healthy'],
                          y2=stats['mutator']['healthy'],
                          xlabel='Epoch', ylabel='Number of healthy cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Healhty pop. evolution. Params go here',
                          filepath=path, filename='healthy_population_evolution')
    clear_plots()
    # mutated cells
    gen_barplot_two_lines(x=x, y1=stats['stable']['mutated'],
                          y2=stats['mutator']['mutated'],
                          xlabel='Epoch', ylabel='Number of mutated cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Mutated pop. evolution. Params go here',
                          filepath=path, filename='mutated_population_evolution')
    # damaged cells
    clear_plots()
    gen_barplot_two_lines(x=x, y1=stats['stable']['damaged'],
                          y2=stats['mutator']['damaged'],
                          xlabel='Epoch', ylabel='Number of damaged cells',
                          line1Legend='Stable cells', line2Legend='Mutator cells',
                          plotTitle='Damaged pop. evolution. Params go here',
                          filepath=path, filename='damaged_population_evolution')