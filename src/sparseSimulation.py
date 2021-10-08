import os
from copy import copy, deepcopy
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import Cell, cellStates, color_world, cellTypes, gen_save_plots, init_world

import pathlib

SIMULATIONSTATS = {'stable': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []},
                   'mutator': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []}}


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


def autoreplicate(world: List[List[Cell]], row: int, col: int, mooreDomain, rng: np.random.default_rng) -> None:
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
        targetCellIndexes = choose_moore_domain(row, col, mooreDomain, 2, rng)
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
                        world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = deepcopy(centerCell)
                elif not nullFirstCell and targetCell is not None:
                    # both cells are full, no replication
                    return
                elif not nullFirstCell and targetCell is None:
                    firstCell = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
                    if firstCell.type == centerCell.type and firstCell.state == centerCell.state:
                        # replication!
                        world[cellIndex[0]][cellIndex[1]] = deepcopy(centerCell)
    return


def mutate(world, row, col, mutationProb, mooreDomain, rng):
    centerCell = world[row][col]

    # mutation only happens for –--healthy--- damaged cells
    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        return

    targetCellIndexes = choose_moore_domain(row, col, mooreDomain, 2, rng)
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
                    if rng.uniform() < mutationProb:
                        world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = deepcopy(targetCell)
                        mutated = True
            elif not nullFirstCell and targetCell is not None:
                # both cells are full, no mutation
                return
            elif not nullFirstCell and targetCell is None:
                firstCell = world[targetCellIndexes[0][0]][targetCellIndexes[0][1]]
                if firstCell.type == centerCell.type and firstCell.state == 'mutated':
                    # mutation with some probability!
                    if rng.uniform() < mutationProb:
                        world[cellIndex[0]][cellIndex[1]] = deepcopy(firstCell)
                        mutated = True
    # if mutation by moore domain has not happened, a small chance of mutation exists to jumpstart
    if not mutated:
        if rng.uniform() < 0.01 * mutationProb:
            world[row][col].state = cellStates['mutated']
    return


def damage(world, row, col, damageProb, rng):
    centerCell = world[row][col]
    # only healthy cells get genetic damage
    if centerCell.state == cellStates['damaged'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < damageProb:
        centerCell.state = cellStates['damaged']
        world[row][col] = centerCell


def degradation_or_repair(world, row, col, deathProb, rng):
    """ In degradation, either the damaged cell dies or it gets repaired. """

    centerCell = world[row][col]

    # only damaged cells get degradated or repaired
    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        return

    if rng.uniform() < centerCell.repairProb:
        world[row][col].state = cellStates['healthy']
    elif rng.uniform() < deathProb:
        # dies
        world[row][col] = None
    return


def diffusion(world, row, col, rng, mooreDomain, diffusionRate):
    """ Diffusion (swap) stage """

    centerCell = copy(world[row][col])
    targetCellIndexes = choose_moore_domain(row, col, mooreDomain, 1, rng)

    if rng.uniform() < diffusionRate:
        if targetCellIndexes[0] is not None:
            world[row][col] = copy(world[targetCellIndexes[0][0]][targetCellIndexes[0][1]])
        else:
            world[row][col] = None
        world[targetCellIndexes[0][0]][targetCellIndexes[0][1]] = centerCell


def forward_generation(world: List[List[Cell]], diffusionRate: float, damageProb: float,
                       deathProb: float, mutationProb: float, mooreDomain, rng: np.random.default_rng) -> dict:

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
            autoreplicate(world, row, col, mooreDomain, rng)
            damage(world, row, col, damageProb, rng)
            mutate(world, row, col, mutationProb, mooreDomain, rng)
            degradation_or_repair(world, row, col, deathProb, rng)
            diffusion(world, row, col, rng, mooreDomain, diffusionRate)

    update_statistics(world)

    return


def sparse_simulation(args, rng):
    """ In this simulation, cells can only be in healthy, mutated or damaged states. """

    worldSize = args.worldSize
    totalPopDensity = args.totalPopDensity
    stableDensity = args.stablePopDensity
    epochs = args.epochs
    diffRate = args.diffusionRate
    damageProb = args.damageProb
    deathProb = args.deathProb
    mutationProb = args.mutationProb
    animate = args.createAnimation

    expName = f'experiment_{datetime.now().strftime("%d%m%Y_%H%M%S")}'
    currentPath = pathlib.Path(__file__).parent.resolve()

    figurePath = currentPath.parent.resolve() / 'experiment_results/'
    if not figurePath.exists():
        os.mkdir(figurePath)

    os.mkdir(figurePath / expName)

    world = [[None for _ in range(worldSize)] for _ in range(worldSize)]

    mooreDomain = precompute_moore_domain(world)
    init_world(world, totalPopDensity, stableDensity, args, rng)

    if animate:
        fig = plt.figure()
        data = np.zeros((worldSize, worldSize, 3))
        im = plt.imshow(data)

        def init():
            im.set_data(np.zeros((worldSize, worldSize, 3)))
            return im

        def animate_frame(_):
            forward_generation(world, diffRate, damageProb, deathProb, mutationProb, mooreDomain, rng)
            im.set_data(color_world(world))
            return im

        anim = animation.FuncAnimation(fig, animate_frame, init_func=init, frames=epochs)
        anim.save(f'{figurePath}/{expName}/system_evolution.gif', fps=min(min(10, epochs / 20), 24))
    else:
        for i in range(epochs):
            forward_generation(world, diffRate, damageProb, deathProb, mutationProb, mooreDomain, rng)

    gen_save_plots(epochs, SIMULATIONSTATS, figurePath / expName)


def main():
    sparse_simulation(None, np.random.default_rng(0))


if __name__ == "__main__":
    main()
