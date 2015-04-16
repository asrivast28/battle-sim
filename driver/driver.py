import battlesim
from battlesim import Soldier, FloorField

import numpy
from random import randint, seed
from time import sleep

import gizeh
import moviepy.editor as mpy


seed(0)

H, W = 200, 200

accessibility = numpy.full((H, W), 255, dtype = numpy.uint8)
for x in xrange(90, 110):
    for y in xrange(90, 110):
        accessibility[x][y] = 0

ff = FloorField(accessibility)
#ff = FloorField(H, W)

for y in xrange(0, W):
    ff.setSoldier(0, y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(1, y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 2, y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 1, y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
#ff.printGrid()

ff.setTarget(0, H, W / 2)
ff.setTarget(1, 0, W / 2)

ff.initializeNeighborhood()

army = numpy.full((H, W), -1, dtype = numpy.uint8)

obstructions = []
it = numpy.nditer(accessibility, op_flags = ['readonly'], flags = ['multi_index'])
while not it.finished:
    if it[0] != 255:
        color = list(it[0] / 255.0 for x in xrange(3))
        obstructions.append(gizeh.square(1, xy = (it.multi_index[1], it.multi_index[0]), fill = color))
    it.iternext()
obstructions = gizeh.Group(obstructions)

i = 0

def make_frame(t):
    global i
    if i % 2 == 0:
        ff.move()
    else:
        ff.kill()
    i += 1
    field = gizeh.Surface(H, W, bg_color = (1, 1, 1))
    obstructions.draw(field)
    ff.armyGrid(army)
    soldiers = []
    it = numpy.nditer(army, op_flags = ['readwrite'], flags = ['multi_index'])
    while not it.finished:
        if it[0] == 0:
            it[0] = -1
            soldiers.append(gizeh.square(1, xy = (it.multi_index[1], it.multi_index[0]), fill = (0, 1, 0)))
        elif it[0] == 1:
            it[0] = -1
            soldiers.append(gizeh.square(1, xy = (it.multi_index[1], it.multi_index[0]), fill = (0, 0, 1)))
        it.iternext()
    gizeh.Group(soldiers).draw(field)
    return field.get_npimage()

clip = mpy.VideoClip(make_frame, duration=100)
clip.write_videofile("battle.mp4", fps=5)
#clip.write_gif("circle.gif", fps=5, opt="OptimizePlus", fuzz=10)
