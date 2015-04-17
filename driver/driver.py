import battlesim
from battlesim import Soldier, BattleField

import numpy
from random import randint, seed
from time import sleep

import gizeh
import moviepy.editor as mpy


seed(0)

H, W = 80, 80

accessibility = numpy.full((H, W), 255, dtype = numpy.uint8)
killed = numpy.full((H, W), False, dtype = numpy.bool)
#for x in xrange(90, 110):
    #for y in xrange(90, 110):
        #accessibility[x][y] = 0


#ff.printGrid()

ff = BattleField(accessibility)
#ff = BattleField(H, W)

for y in xrange(0, W):
    ff.setSoldier(0, y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(1, y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 2, y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))
    ff.setSoldier(H - 1, y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255)))

ff.setTarget(0, H, W / 2)
ff.setTarget(1, 0, W / 2)

ff.initializeNeighborhood()


class PixelFill(object):
    """
    Metadata for filling the image with pixels.
    """
    def __init__(self):
        self.scale = 1
        self.size = 1
        self.mapping = {0 : {}, 1 : {}}

        self.mapping[0][Soldier.LEADER] = (0, 1, 0)
        self.mapping[1][Soldier.LEADER] = (0, 0, 1)

        self.mapping[0][Soldier.SWORDSMAN] = (0, 1, 0)
        self.mapping[1][Soldier.SWORDSMAN] = (0, 0, 1)

        self.mapping[0][Soldier.ARCHER] = (0, 1, 0)
        self.mapping[1][Soldier.ARCHER] = (0, 0, 1)

        self.dead = (1, 0, 0)


class IconFill(PixelFill):
    """
    Metadata for filling the image with icons.
    """
    def __init__(self):
        self.scale = 10
        self.size = 20
        self.mapping = {0: {}, 1: {}}

        for t in ('LEADER', 'SWORDSMAN', 'ARCHER'):
            for a in (0, 1):
                self.mapping[a][getattr(Soldier, t)] = self.getImage('%s_%d.png'%(t.lower(), a)).scale(10)

        self.dead = self.getImage('dead.png').scale(10)

    def getImage(self, name):
        import itertools, png
        image = png.Reader(filename = name)
        pngdata = image.asDirect()[2]
        image_2d = numpy.vstack(itertools.imap(numpy.uint8, pngdata))
        image_3d = numpy.reshape(image_2d, (image.height, image.width, image.planes))
        return gizeh.ImagePattern(image_3d)



class FrameBuilder(object):
    """
    Class that builds the actual frames in the video.
    """
    def __init__(self, fill):
        self.fill = fill
        self.field = self.createField(accessibility)
        self.soldiers = numpy.full((H, W), -1, dtype = numpy.uint8)
        self.killed = numpy.full((H, W), False, dtype = numpy.bool)
        self.i = 0

    def createField(self, accessibility):
        """
        Creates an image of the field with obstructions.
        """
        obstructions = []
        it = numpy.nditer(accessibility, op_flags = ['readonly'], flags = ['multi_index'])
        while not it.finished:
            if it[0] != 255:
                color = list(it[0] / 255.0 for x in xrange(3))
                obstructions.append(gizeh.square(self.fill.size, xy = (it.multi_index[1], it.multi_index[0]), fill = color))
            it.iternext()
        obstructions = gizeh.Group(obstructions)
        field = gizeh.Surface(H * self.fill.scale, W * self.fill.scale, bg_color = (1, 1, 1))
        obstructions.draw(field)
        return field.get_npimage()

    def move(self):
        # create new frame
        self.frame = gizeh.Surface.from_image(self.field)
        # move soldiers and obtain their new positions
        ff.move(self.soldiers)
        # draw soldiers, in new position, on the field
        soldiers = []
        it = numpy.nditer(self.soldiers, op_flags = ['readwrite'], flags = ['multi_index'])
        while not it.finished:
            if it[0] != 255:
                xy = [it.multi_index[a] * self.fill.scale for a in (1, 0)]
                fill = self.fill.mapping[int(it[0])][Soldier.SWORDSMAN]
                soldiers.append(gizeh.square(self.fill.size, xy = xy, fill = fill))
                it[0] = -1
            it.iternext()
        gizeh.Group(soldiers).draw(self.frame)

    def kill(self):
        # kill soldiers and find out which soldiers are killed
        ff.kill(self.killed)
        # draw killed soldiers
        killed = []
        it = numpy.nditer(self.killed, op_flags = ['readwrite'], flags = ['multi_index'])
        while not it.finished:
            if it[0]:
                xy = [it.multi_index[a] * self.fill.scale for a in (1, 0)]
                fill = self.fill.dead
                killed.append(gizeh.square(self.fill.size, xy = xy, fill = fill))
                it[0] = False
            it.iternext()
        gizeh.Group(killed).draw(self.frame)

    def __call__(self, t):
        if self.i % 2 == 0:
            self.move()
        else:
            self.kill()
        self.i += 1
        # create a new surface for each frame
        return self.frame.get_npimage()


#clip = mpy.VideoClip(FrameBuilder(IconFill()), duration = 80)
clip = mpy.VideoClip(FrameBuilder(PixelFill()), duration=50)
#clip.write_videofile("battle.mp4", fps=2)
clip.write_gif("circle.gif", fps=2, opt="OptimizePlus", fuzz=10)
