PXL_SIZE = 4
FPS = 30

ITERATIONS = 15

WHITE = [255] * 3
BLACK = [0] * 3
RED = [255, 0, 0]
DARK_RED = [150, 0, 0]
GREEN = [0, 255, 0]
DARK_GREEN = [0, 150, 0]

# celltypes
STABLE = 1
MUTATOR = 2
NONE = 0

# cellstates
EMPTY = 0
HEALTHY = 1
REPAIRING = 2
MUTANT = 3
DEAD = 4


class CellState(object):
    def __init__(self):
        self.empty = 0
        self.healthy = 1
        self.repairing = 2
        self.mutant = 4


class Cell(object):
    def __init__(self, _type, r):
        self.type = _type
        self.state = CELL_STATES['healthy']

        if self.type == CELL_STATES['healthy']:
            self.r = r_s
            self.e = e_s
        else:
            self.r = r_m
            self.e = e_m

    def apply_rules(self):
        if self.state == HEALTHY:
            if chance(u):
                if chance(self.e):
                    self.state = REPAIRING
                else:
                    if chance(alpha):
                        self.state = MUTANT
                    else:
                        self.state = DEAD

        if self.state == REPAIRING and chance(beta):
            self.state = HEALTHY

        if self.state == MUTANT and chance(a):
            self.state = DEAD