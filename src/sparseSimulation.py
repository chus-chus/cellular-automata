import numpy as np
from typing import List

from utils import Cell, cellStates, plot_world


def init_world(world: List[List[Cell]], density: float, rng: np.random.default_rng) -> None:
    """ Initialise grid with stable or mutator cells. There are density * len(world) cells placed, half of them
     being stable and half of them being mutator. """
    rows = rng.choice(len(world), size=round(len(world) * density), replace=False)
    cols = rng.choice(len(world), size=round(len(world) * density), replace=False)

    for i in range(0, len(cols) // 2):
        world[rows[i]][cols[i]] = Cell('stable')

    for i in range(len(cols) // 2, len(cols)):
        world[rows[i]][cols[i]] = Cell('mutator')

    return


def compute_statistics(world):
    return 0


def choose_moore_domain(i: int, j: int, n: int, worldSize: int, rng: np.random.default_rng) -> np.ndarray:
    """ Choose n indexes w/o replacement in the moore domain centered around i, j. Does not consider grid edges.

    :return: a list of indexes as lists. """

    mooreDomain = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                   [i,     j - 1],             [i,     j + 1],
                   [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]

    # parse invalid indexes (edges)
    for i in range(len(mooreDomain)):
        if mooreDomain[i][0] >= worldSize or mooreDomain[i][0] < 0:
            mooreDomain[i][0] = None
        if mooreDomain[i][1] >= worldSize or mooreDomain[i][1] < 0:
            mooreDomain[i][1] = None

    # delete invalid indexes (edges)
    mooreDomain = [index for index in mooreDomain if index[0] is not None and index[1] is not None]

    return rng.choice(mooreDomain, size=n, replace=False)


def autoreplicate(world: List[List[Cell]], row: int, col: int, rng: np.random.default_rng) -> None:
    """ Autoreplication step. Two random cells are chosen in the Moore domain of the cell, and only
     if one of those is of the same type and the other one is empty will the cell replicate, the new one
     being of the same type. """

    centerCell = world[row][col]

    # replication only happens for mutated or healthy cells
    if centerCell.state == cellStates['damaged']:
        return

    targetCellIndexes = choose_moore_domain(row, col, 2, len(world), rng)

    # check what is in each of the picked grid spots and replicate accordingly
    nullFirstCell = True
    for i, cellIndex in enumerate(targetCellIndexes):
        try:
            targetCell = world[cellIndex[0]][cellIndex[1]]
        except:
            a = 1
        if i == 0 and targetCell is not None:
            nullFirstCell = False
        elif i == 1:
            if nullFirstCell and targetCell is None:
                # this cell is not empty and neither is the previous one, so replication cannot happen
                return
            elif nullFirstCell and targetCell is not None:
                if targetCell.type == centerCell.type and targetCell.state == centerCell.state:
                    # replication!
                    world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = centerCell
            elif not nullFirstCell and targetCell is not None:
                # both cells are full, no replication
                return
            elif not nullFirstCell and targetCell is None:
                firstCell = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
                if firstCell.type == centerCell.type and firstCell.state == centerCell.state:
                    # replication!
                    world[cellIndex[0]][cellIndex[1]] = centerCell
    return


def forward_generation(world: List[List[Cell]], statistics, rng: np.random.default_rng) -> None:
    """ An evolution epoch. That is, :math:`L^2` random grid spots are chosen. If it's empty, do nothing.
     If a cell is in the grid spot, perform the following steps:

       - Autoreplication
       - Mutation
       - Genetic damage
       - Degradation
       - Diffusion

     , according to the specified probabilities. Returns a statistical description of the current
     world. """

    L = len(world)
    rows = rng.choice(L, L * L, replace=True)
    cols = rng.choice(L, L * L, replace=True)
    indexes = list(zip(rows, cols))
    for index in indexes:
        row, col = index[0], index[1]
        if world[row][col] is None:
            # empty grid spot!
            pass
        else:
            # autoreplication
            autoreplicate(world, row, col, rng)
            # mutation
            # damage
            # degradation
            # diffusion

    # save statistics
    statistics = compute_statistics(world)

    return


def sparse_simulation(args, rng):
    """ In this simulation, cells can only be in healthy, mutated or damaged states. """

    # worldSize = args['worldSize']
    worldSize = 300
    density = 1
    rrMutant = 1.5
    rrStable = 1.2
    epochs = 50

    world = [[None for _ in range(worldSize)] for _ in range(worldSize)]

    init_world(world, density, rng)

    statistics = None
    for i in range(epochs):
        forward_generation(world, statistics, rng)
        if i % 10 == 0:
            plot_world(world)


def main():
    sparse_simulation(None, np.random.default_rng(888))


if __name__ == "__main__":
    main()
