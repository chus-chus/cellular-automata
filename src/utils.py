""" Utility variables, methods and classes """

# celltypes
cellTypes = {'stable': 0, 'mutator': 1}

# cellstates
cellStates = {'healthy': 0, 'repairing': 1, 'mutant': 2}

# grid colors
__WHITE = [255, 255, 255]
__BLACK = [0, 0, 0]
__RED = [255, 0, 0]
__DARK_RED = [150, 0, 0]
__GREEN = [0, 255, 0]
__DARK_GREEN = [0, 150, 0]

colors = {'stable': {'healthy': __GREEN, 'repairing': __GREEN, 'mutant': __DARK_GREEN},
          'mutant': {'healthy': __RED, 'repairing': __RED, 'mutant': __DARK_RED},
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

