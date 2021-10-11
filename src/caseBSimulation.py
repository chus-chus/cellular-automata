import os
from copy import copy, deepcopy
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from typing import List

from utils import Cell, cellStates, color_world, cellTypes, gen_save_plots, init_world, choose_moore_domain, \
    precompute_moore_domain, parse_all_params, update_statistics, simulation_skeleton

import pathlib


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


def case_b_forward_function(world: List[List[Cell]], diffusionRate: float, damageProb: float,
                            deathProb: float, mutationProb: float, mooreDomain, prevStats:dict,
                            rng: np.random.default_rng) -> dict:

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

    stats = update_statistics(world, prevStats)

    return stats


def case_b_simulation(args, rng):
    """ In this simulation, cells can only be in healthy, mutated or damaged states. """

    ff_args = {'damageProb': args.damageProb, 'deathProb': args.deathProb, 'mutationProb': args.mutationProb,
               'diffusionRate': args.diffusionRate}

    simulation_skeleton(args, rng, case_b_forward_function, args.worldSize, args.iterations, args.totalPopDensity,
                        args.stablePopDensity, args.createAnimation, args.epochs, ff_args)
