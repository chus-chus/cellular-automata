""" Utility variables, methods and classes """

# celltypes
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
                 repairRate: float = 0.5):

        self.type = cellTypes[cellType]
        self.state = cellStates[cellState]
        self.replicationRate = replicationRate
        self.repairRate = repairRate


def plot_world(world: List[List[Cell]]):
    # translate the world into a np.ndarray with the specified color values
    coloredWorld = [[None for _ in range(len(world))] for _ in range(len(world))]

    for i in range(len(world)):
        for j in range(len(world)):
            cell = world[i][j]
            if cell is None:
                coloredWorld[i][j] = colors['empty']
            else:
                coloredWorld[i][j] = colors[cell.type][cell.state]

    plt.imshow(coloredWorld)
    plt.show()
