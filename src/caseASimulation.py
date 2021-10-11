import os
from copy import deepcopy
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import cellStates, color_world, cellTypes, gen_save_plots, init_world, choose_moore_domain, \
    precompute_moore_domain, parse_all_params, update_statistics, update_mean_stats, simulation_skeleton

import pathlib


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


def case_a_forward_function(world: List[List], damageProb: float, deathProb: float,
                            mutationProb: float, arrestProb: float, mooreDomain, prevStats: dict,
                            rng: np.random.default_rng) -> dict:
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

    stats = update_statistics(world, prevStats)

    return stats


def case_a_simulation(args, rng):
    ff_args = {'damageProb': args.damageProb, 'deathProb': args.deathProb,
               'mutationProb': args.mutationProb, 'arrestProb': args.arrestProb}

    simulation_skeleton(args, rng, case_a_forward_function, args.worldSize, args.iterations, args.totalPopDensity,
                        args.stablePopDensity, args.createAnimation, args.epochs, ff_args)

