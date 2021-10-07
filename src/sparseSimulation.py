import numpy as np
from utils import Cell


def init_world(world, density, rng):
    rows = rng.choice(world.size, round(world.size * density))
    cols = rng.choice(world.size, round(world.size * density))
    for i in range(len(cols)):
        world[rows[i]][cols[i]].cellState = CELLSTATES
    world[rows[i]][cols[i]] for
    return world

def sparse_simulation(args, rng):
    # worldSize = args['worldSize']
    worldSize = 300
    density = 0.01
    rrMutant = 1.5
    rrStable = 1.2

    world = [[Cell for _ in range(worldSize)] for _ in range(worldSize)]



