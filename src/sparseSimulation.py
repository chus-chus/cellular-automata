import os
from copy import copy
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import Cell, cellStates, color_world, cellTypes, gen_save_plots

import pathlib

SIMULATIONSTATS = {'stable': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []},
                   'mutator': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []}}


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
    stableMutated = 0
    stableDamaged = 0
    stableHealthy = 0

    mutatorCells = 0
    mutatorMutated = 0
    mutatorDamaged = 0
    mutatorHealthy = 0

    for i in range(len(world)):
        for j in range(len(world)):
            if world[i][j] is not None:
                if world[i][j].type == cellTypes['stable']:
                    stableCells += 1
                    if world[i][j].state == cellStates['mutated']:
                        stableMutated += 1
                    elif world[i][j].state == cellStates['damaged']:
                        stableDamaged += 1
                    elif world[i][j].state == cellStates['healthy']:
                        stableHealthy += 1
                else:
                    mutatorCells += 1
                    if world[i][j].state == cellStates['mutated']:
                        mutatorMutated += 1
                    elif world[i][j].state == cellStates['damaged']:
                        mutatorDamaged += 1
                    elif world[i][j].state == cellStates['healthy']:
                        mutatorHealthy += 1

    SIMULATIONSTATS['stable']['total'].append(stableCells)
    SIMULATIONSTATS['stable']['healthy'].append(stableHealthy)
    SIMULATIONSTATS['stable']['damaged'].append(stableDamaged)
    SIMULATIONSTATS['stable']['mutated'].append(stableMutated)

    SIMULATIONSTATS['mutator']['total'].append(mutatorCells)
    SIMULATIONSTATS['mutator']['healthy'].append(mutatorHealthy)
    SIMULATIONSTATS['mutator']['damaged'].append(mutatorDamaged)
    SIMULATIONSTATS['mutator']['mutated'].append(mutatorMutated)

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

    centerCell = copy(world[row][col])

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
    centerCell = copy(world[row][col])

    # mutation only happens for –--healthy--- damaged cells
    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        return

    targetCellIndexes = choose_moore_domain(row, col, 2, len(world), rng)
    # check what is in each of the picked grid spots and replicate accordingly
    nullFirstCell = True
    mutated = False
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
                        mutated = True
            elif not nullFirstCell and targetCell is not None:
                # both cells are full, no mutation
                return
            elif not nullFirstCell and targetCell is None:
                firstCell = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
                if firstCell.type == centerCell.type and firstCell.state == 'mutated':
                    # mutation with some probability!
                    if rng.uniform() < centerCell.mutationProb:
                        world[cellIndex[0]][cellIndex[1]] = firstCell
                        mutated = True
    # if mutation by moore domain has not happened, a small chance of mutation exists to jumpstart
    if not mutated:
        if rng.uniform() < 0.01 * centerCell.mutationProb:
            world[row][col].state = cellStates['mutated']
    return


def damage(world, row, col, damageProb, rng):
    centerCell = copy(world[row][col])
    # only healthy cells get genetic damage
    if centerCell.state == cellStates['damaged'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < damageProb:
        centerCell.state = cellStates['damaged']
        world[row][col] = centerCell


def degradation_or_repair(world, row, col, rng):
    """ In degradation, either the damaged cell dies or it gets repaired. """

    centerCell = copy(world[row][col])

    # only damaged cells get degradated or repaired
    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < centerCell.repairProb:
        world[row][col].state = cellStates['healthy']
    elif rng.uniform() < centerCell.deathProb:
        # dies
        world[row][col] = None
    return


def diffusion(world, row, col, rng, diffusionRate):
    """ Diffusion (swap) stage """

    centerCell = copy(world[row][col])
    targetCellIndexes = choose_moore_domain(row, col, 1, len(world), rng)

    if rng.uniform() < diffusionRate:
        if targetCellIndexes[0] is not None:
            world[row][col] = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
        else:
            world[row][col] = None
        world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = centerCell


def forward_generation(world: List[List[Cell]], diffusionRate: float,
                       damageProb: float, rng: np.random.default_rng) -> dict:

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
            damage(world, row, col, damageProb, rng)
            mutate(world, row, col, rng)
            degradation_or_repair(world, row, col, rng)
            diffusion(world, row, col, rng, diffusionRate)

    update_statistics(world)

    return


def sparse_simulation(args, rng):
    """ In this simulation, cells can only be in healthy, mutated or damaged states. """

    # worldSize = args['worldSize']
    worldSize = 50
    density = 0.1
    rrMutant = 3
    rrStable = 1.1
    epochs = 150
    diffRate = 0.05
    damageProb = 0.02
    animate = True

    dirName = f'experiment_{datetime.now().strftime("%d%m%Y_%H%M%S")}'
    currentPath = pathlib.Path(__file__).parent.resolve()

    os.mkdir(currentPath / dirName)

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
            forward_generation(world, diffRate, damageProb, rng)
            im.set_data(color_world(world))
            return im

        anim = animation.FuncAnimation(fig, animate_frame, init_func=init, frames=epochs)
        anim.save(f'{dirName}/system_evolution.gif', fps=5)
    else:
        for i in range(epochs):
            forward_generation(world, diffRate, damageProb, rng)

    gen_save_plots(epochs, SIMULATIONSTATS, currentPath / dirName)


def main():
    sparse_simulation(None, np.random.default_rng(0))


if __name__ == "__main__":
    main()
