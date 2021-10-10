import os
from copy import copy, deepcopy
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import Cell, cellStates, color_world, cellTypes, gen_save_plots, init_world, choose_moore_domain, \
    precompute_moore_domain, parse_all_params

import pathlib

SIMULATIONSTATS = {'stable': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []},
                   'mutator': {'total': [], 'healthy': [], 'damaged': [], 'mutated': []},
                   'stableVictories': 0, 'mutatorVictories': 0}


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


def replicate(world: List[List], row: int, col: int, mooreDomain, rng: np.random.default_rng):
    centerCell = world[row][col]

    if centerCell.state == cellStates['healthy'] or centerCell.state == cellStates['mutated']:
        # replicate according to the cell replication rate
        replicationTimes = int(centerCell.replicationRate // 1 + int(rng.uniform() < centerCell.replicationRate % 1))
        replicasIndexes = choose_moore_domain(row, col, mooreDomain, replicationTimes, rng)
        for cellIndex in replicasIndexes:
            world[cellIndex[0]][cellIndex[1]] = deepcopy(centerCell)


def apply_rules(world: List[List], row: int, col: int, damageProb: float, mutationProb: float,
                arrestProb: float, deathProb: float, rng: np.random.default_rng):
    centerCell = world[row][col]

    if centerCell.state == cellStates['healthy']:
        if rng.uniform() < damageProb:
            if rng.uniform() < centerCell.repairProb:
                world[row][col].state = cellStates['repairing']
            else:
                if rng.uniform() < mutationProb:
                    # mutation occurred
                    world[row][col].state = cellStates['mutated']
                else:
                    # dies
                    world[row][col] = None

    if centerCell.state == cellStates['repairing'] and rng.uniform() < arrestProb:
        # repairs itself
        world[row][col].state = cellStates['healthy']

    if centerCell.state == cellStates['mutated'] and rng.uniform() < deathProb:
        # dies
        world[row][col] = None

    return

def forward_generation(world: List[List], damageProb: float, deathProb: float,
                       mutationProb: float, arrestProb: float, mooreDomain, save_stats: bool, rng: np.random.default_rng) -> dict:

    """ An evolution epoch. That is, :math:`L^2` random grid spots are chosen. If it's empty, do nothing.
        If a cell is in the grid spot, it replicates and the model's rules are applied in order to change
        the cell's state, according to the specified probabilities.

        Returns a statistical description of the current world. """

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
            replicate(world, row, col, mooreDomain, rng)
            apply_rules(world, row, col, damageProb, mutationProb, arrestProb, deathProb, rng)

    if save_stats: update_statistics(world)

    return


def case_a_simulation(args, rng):

    worldSize = args.worldSize
    totalPopDensity = args.totalPopDensity
    stableDensity = args.stablePopDensity
    epochs = args.epochs
    damageProb = args.damageProb
    deathProb = args.deathProb
    mutationProb = args.mutationProb
    arrestProb = args.arrestProb
    animate = args.createAnimation
    iterations = args.iterations

    expName = f'experiment_{datetime.now().strftime("%d%m%Y_%H%M%S")}'
    currentPath = pathlib.Path(__file__).parent.resolve()

    figurePath = currentPath.parent.resolve() / 'experiment_results/'
    if not figurePath.exists():
        os.mkdir(figurePath)

    os.mkdir(figurePath / expName)

    # save parameters of the experiment
    with open(figurePath / expName / "params.txt", "w") as file:
        file.write(parse_all_params(args))

    world = [[None for _ in range(worldSize)] for _ in range(worldSize)]

    mooreDomain = precompute_moore_domain(world)

    if iterations == 1:
        save_stats = True
        init_world(world, totalPopDensity, stableDensity, args, rng)

        if animate:
            fig = plt.figure()
            data = np.zeros((worldSize, worldSize, 3))
            im = plt.imshow(data)

            plt.tick_params(axis='both',  # changes apply to both
                            which='both',  # both major and minor ticks are affected
                            bottom=False,  # ticks along the bottom edge are off
                            top=False,  # ticks along the top edge are off
                            left=False,
                            labelleft=False,
                            labelbottom=False)

            def init():
                im.set_data(np.zeros((worldSize, worldSize, 3)))
                return im

            def animate_frame(_):
                forward_generation(world, damageProb, deathProb, mutationProb, arrestProb, mooreDomain, save_stats, rng)
                im.set_data(color_world(world))
                return im

            anim = animation.FuncAnimation(fig, animate_frame, init_func=init, frames=epochs)
            anim.save(f'{figurePath}/{expName}/system_evolution.gif', fps=min(max(10, epochs / 20), 24))
        else:
            for i in range(epochs):
                forward_generation(world, damageProb, deathProb, mutationProb, arrestProb, mooreDomain, save_stats, rng)

        gen_save_plots(epochs, SIMULATIONSTATS, figurePath / expName)

    elif iterations > 1:
        save_stats = False
        for i in range(iterations):
            new_seed = rng.integers(10000000)
            new_rng = np.random.default_rng(new_seed)
            init_world(world, totalPopDensity, stableDensity, args, new_rng)

            for _ in range(epochs):
                forward_generation(world, damageProb, deathProb, mutationProb, arrestProb, mooreDomain, save_stats, rng)

            stableCells = 0
            mutatorCells = 0
            for i in range(len(world)):
                for j in range(len(world)):
                    if world[i][j] is not None:
                        if world[i][j].type == cellTypes['stable']:
                            stableCells += 1
                        else:
                            mutatorCells += 1

            if mutatorCells > stableCells:
                SIMULATIONSTATS['mutatorVictories'] += 1
            else:
                SIMULATIONSTATS['stableVictories'] += 1

        print(f"iterations: {iterations}\n"
              f"stable cell victories: {SIMULATIONSTATS['stableVictories']}"
              f"mutator cell victories: {SIMULATIONSTATS['mutatorVictories']}")


def main():
    case_a_simulation(None, np.random.default_rng(0))


if __name__ == "__main__":
    main()
