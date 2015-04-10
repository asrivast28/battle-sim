import battlesim
import numpy
from random import randint, seed
from time import sleep

seed(0)

H, W = 600, 600

ff = battlesim.FloorField(H, W)
for y in xrange(0, W):
    ff.setSoldier(0, y, battlesim.Soldier(0, battlesim.Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(1, y, battlesim.Soldier(0, battlesim.Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 2, y, battlesim.Soldier(1, battlesim.Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 1, y, battlesim.Soldier(1, battlesim.Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
#ff.printGrid()

ff.setTarget(0, H, W / 2)
ff.setTarget(1, 0, W / 2)

ff.initializeNeighborhood()

army_grid = numpy.full((H, W), -1, dtype = numpy.uint8)

import gizeh
import moviepy.editor as mpy

i = 0

def make_frame(t):
    global i
    if i % 2 == 0:
        ff.move()
    else:
        ff.kill()
    i += 1
    surface = gizeh.Surface(H, W)
    image = surface.get_npimage()
    ff.armyGrid(army_grid)
    it = numpy.nditer(army_grid, op_flags = ['readwrite'], flags = ['multi_index'])
    while not it.finished:
        #import pdb;pdb.set_trace()
        if it[0] == 0:
            it[0] = -1
            image[it.multi_index[0]][it.multi_index[1]] = (255, 0, 0)
        elif it[0] == 1:
            it[0] = -1
            image[it.multi_index[0]][it.multi_index[1]] = (0, 0, 255)
        it.iternext()
    return image

clip = mpy.VideoClip(make_frame, duration=600)
clip.write_gif("circle.gif", fps=2, opt="OptimizePlus", fuzz=10)
