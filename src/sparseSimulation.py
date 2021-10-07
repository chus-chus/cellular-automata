import os
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import Cell, cellStates, color_world, cellTypes

SIMULATIONSTATS = {'stable': [], 'mutator': []}


def init_world(world: List[List[Cell]], density: float, stableRR: float, mutatorRR: float,
               rng: np.random.default_rng) -> None:
    """ Initialise grid with stable or mutator cells. There are density * len(world) cells placed, half of them
     being stable and half of them being mutator. """

    rows = rng.choice(len(world), size=round(len(world) * len(world) * density), replace=True)
    cols = rng.choice(len(world), size=round(len(world) * len(world) * density), replace=True)

    for i in range(0, len(cols) // 2):
        world[rows[i]][cols[i]] = Cell('stable', replicationRate=stableRR)

    for i in range(len(cols) // 2, len(cols)):
        world[rows[i]][cols[i]] = Cell('mutator', replicationRate=mutatorRR)

    return


def update_statistics(world):
    """ Updates or creates a dict with statistics of the current world """

    stableCells = 0
    mutatorCells = 0
    for i in range(len(world)):
        for j in range(len(world)):
            if world[i][j] is not None:
                if world[i][j].type == cellTypes['stable']:
                    stableCells += 1
                else:
                    mutatorCells += 1

    SIMULATIONSTATS['stable'].append(stableCells)
    SIMULATIONSTATS['mutator'].append(mutatorCells)

    return


def choose_moore_domain(i: int, j: int, n: int, worldSize: int, rng: np.random.default_rng) -> np.ndarray:
    # todo precompute moore domain for all cells
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

    # try to replicate according to the cell replication rate
    replicationTimes = (centerCell.replicationRate // 1) + float(rng.uniform() < centerCell.replicationRate % 1)
    for _ in range(int(replicationTimes)):
        targetCellIndexes = choose_moore_domain(row, col, 2, len(world), rng)
        # check what is in each of the picked grid spots and replicate accordingly
        nullFirstCell = True
        for i, cellIndex in enumerate(targetCellIndexes):
            targetCell = world[cellIndex[0]][cellIndex[1]]
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


def mutate(world, row, col, rng):
    centerCell = world[row][col]

    # mutation only happens for healthy cells
    if centerCell.state == cellStates['damaged'] or centerCell.state == cellStates['mutated']:
        return

    targetCellIndexes = choose_moore_domain(row, col, 2, len(world), rng)
    # check what is in each of the picked grid spots and replicate accordingly
    nullFirstCell = True
    for i, cellIndex in enumerate(targetCellIndexes):
        targetCell = world[cellIndex[0]][cellIndex[1]]
        if i == 0 and targetCell is not None:
            nullFirstCell = False
        elif i == 1:
            if nullFirstCell and targetCell is None:
                # this cell is not empty and neither is the previous one, so replication cannot happen
                return
            elif nullFirstCell and targetCell is not None:
                if targetCell.type == centerCell.type and targetCell.state == cellStates['mutated']:
                    # mutation with some probability!
                    if rng.uniform() < centerCell.mutationProb:
                        world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = targetCell
            elif not nullFirstCell and targetCell is not None:
                # both cells are full, no mutation
                return
            elif not nullFirstCell and targetCell is None:
                firstCell = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
                if firstCell.type == centerCell.type and firstCell.state == 'mutated':
                    # mutation with some probability!
                    if rng.uniform() < centerCell.mutationProb:
                        world[cellIndex[0]][cellIndex[1]] = firstCell
    return


def damage(world, row, col, rng):
    centerCell = world[row][col]

    # only healthy cells get genetic damage
    if centerCell.state == cellStates['damaged'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < centerCell.damageProb:
        world[row][col].state = cellStates['damaged']

    return


def degradation_or_repair(world, row, col, rng):
    """ In degradation, either the damaged cell dies or it gets repaired. """

    centerCell = world[row][col]

    # only damaged cells get degradated or repaired
    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < centerCell.repairProb:
        world[row][col].state = cellStates['healthy']
    else:
        # dies
        world[row][col] = None

    return


def diffusion(world, row, col, rng, diffusionRate):
    """ Diffusion (swap) stage """

    centerCell = world[row][col]
    targetCellIndexes = choose_moore_domain(row, col, 1, len(world), rng)

    if targetCellIndexes[0] is not None:
        if rng.uniform() < diffusionRate:
            world[row][col] = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
            world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = centerCell


def forward_generation(world: List[List[Cell]], diffusionRate: float, rng: np.random.default_rng) -> dict:

    """ An evolution epoch. That is, :math:`L^2` random grid spots are chosen. If it's empty, do nothing.
     If a cell is in the grid spot, perform the following steps:

       - Autoreplication
       - Mutation
       - Genetic damage
       - Degradation or Repair
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
            continue
        else:
            autoreplicate(world, row, col, rng)
            # mutate(world, row, col, rng)
            damage(world, row, col, rng)
            degradation_or_repair(world, row, col, rng)
            diffusion(world, row, col, rng, diffusionRate)

    # save statistics
    update_statistics(world)

    return


def gen_save_plots(epochs, dirName):
    # clear all plots
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    # population plot
    epochs = list(range(epochs))
    plt.plot(epochs, SIMULATIONSTATS['stable'], label='Stable cells')
    plt.plot(epochs, SIMULATIONSTATS['mutator'], label='Mutator cells')
    plt.xlabel('Epoch')
    plt.ylabel('Number of cells')
    plt.title('Pop. evolution. Params go here')
    plt.legend()
    plt.savefig(f'{dirName}/population_evolution.png', dpi=300)


def sparse_simulation(args, rng):
    """ In this simulation, cells can only be in healthy, mutated or damaged states. """

    # worldSize = args['worldSize']
    worldSize = 50
    density = 0.1
    rrMutant = 1.5
    rrStable = 1.1
    epochs = 100
    diffRate = 0.05
    animate = True

    dirName = f'experiment_{datetime.now().strftime("%d%m%Y_%H%M%S")}'
    os.mkdir(dirName)

    world = [[None for _ in range(worldSize)] for _ in range(worldSize)]

    init_world(world, density, rrStable, rrMutant, rng)

    if animate:
        fig = plt.figure()
        data = np.zeros((worldSize, worldSize, 3))
        im = plt.imshow(data)

        def init():
            im.set_data(np.zeros((worldSize, worldSize, 3)))
            return im

        def animate_frame(_):
            forward_generation(world, diffRate, rng)
            im.set_data(color_world(world))
            return im

        anim = animation.FuncAnimation(fig, animate_frame, init_func=init, frames=epochs)
        anim.save(f'{dirName}/system_evolution.gif', fps=5)
    else:
        for i in range(epochs):
            forward_generation(world, diffRate, rng)

    gen_save_plots(epochs, dirName)


def main():
    sparse_simulation(None, np.random.default_rng(0))


if __name__ == "__main__":
    main()
