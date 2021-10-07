import numpy as np
import pygame
import time



# EMPTY = 0
# HEALTHY = 1
# REPAIRING = 2
# MUTANT = 3
# DEAD = 4
#
# STABLE = 0
# MUTATOR = 1

COLORS = {CELL_TYPES['stable']: {CELL_STATES['healthy']: GREEN,
                                 CELL_STATES['repairing']: GREEN,
                                 CELL_STATES['mutant']: DARK_GREEN,
                                 }}


def color(cell):
    if cell == EMPTY:
        return WHITE
    else:
        if cell.type == STABLE:
            if cell.state == HEALTHY or cell.state == REPAIRING:
                return GREEN
            elif cell.state == MUTANT:
                return DARK_GREEN
            else:
                return BLACK

        if cell.type == MUTATOR:
            if cell.state == HEALTHY or cell.state == REPAIRING:
                return RED
            elif cell.state == MUTANT:
                return DARK_RED
            else:
                return BLACK


#   PARAMETERS

# GRID SIZE
L = 200
# REPLICATION RATES
r_s = 1
r_m = 1.3
# REPAIR PROB
e_s = 0.99
e_m = 0.7
# GENETIC ALTERATION PROB
u = 0.2
# CELL CYCLE ARREST
beta = 0.3
# VIABLE MUTANT PROB
alpha = 0.05
# APOPTOSIS PROB
a = 0.3


class CellState(object):
    def __init__(self):
        self.empty = 0
        self.healthy = 1
        self.repairing = 2
        self.mutant = 4


class Cell(object):
    def __init__(self, _type):
        self.type = _type
        self.state = HEALTHY

        if self.type == STABLE:
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


pos = [[-1, -1], [-1, 0], [-1, 1],
       [0, -1], [0, 1],
       [1, -1], [1, 0], [1, 1]]


def apply_rules(grid, rng):
    # new_grid = [[EMPTY for i in range(L)] for j in range(L)]
    for i in range(L):
        for j in range(L):
            grid[i][j].apply_rules()
            # replicate if not dead
            if grid[i][j].state != DEAD:
                # times the cell i,j replicates
                replicas = int(grid[i][j].r // 1) + np.random.rand(1) < grid[i][j].r % 1
                # replicate in random spots in its Moore domain
                for n in rng.choice(8, replicas):
                    grid[(i + pos[n][0]) % L][(j + pos[n][1]) % L] = Cell(grid[i][j].type)

    return grid



def put_pxl(position, window, color=WHITE):
    pygame.draw.rect(window, color, (position[1], position[0], PXL_SIZE, PXL_SIZE))


def print_window(window, grid):
    window.fill(BLACK)
    for i in range(L):
        for j in range(L):
            put_pxl([i * PXL_SIZE, j * PXL_SIZE], color(grid[i][j]))

    pygame.display.update()


def main():
    randomSeed = 888
    rng = np.random.default_rng(randomSeed)

    # RANDOM GRID
    # grid = [[np.random.choice(3) for i in range(L)] for j in range(L)]

    # HALF GRID
    # grid = [[STABLE_CELL if i < L/2 else MUTANT_CELL for i in range(L)] for j in range(L)]

    # CENTER CANCER
    grid = [[Cell(STABLE) for _ in range(L)] for _ in range(L)]

    for i in range(4):
        for j in range(4):
            grid[L // 2 - 2 + i][L // 2 - 2 + j] = Cell(MUTATOR)

    clock = pygame.time.Clock()
    pygame.init()
    window = pygame.display.set_mode((L * PXL_SIZE, L * PXL_SIZE))
    print_window(window, grid)

    counter = 0
    end = False

    while not end:
        clock.tick(30)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                end = True

        counter += 1
        grid = apply_rules(grid)
        print_window(window, grid)

        # if counter >= ITERATIONS:
        #    end = True

    # pygame.image.save(window, "screenshot_after.jpeg")
    # print("Screenshot saved")

    pygame.quit()


if __name__ == "__main__":
    main()
