import battlesim
from battlesim import Soldier, FloorField

import numpy
from random import randint, seed
from time import sleep

import gizeh
import moviepy.editor as mpy


seed(0)

H, W = 20, 20 

accessibility = numpy.full((H, W), 255, dtype = numpy.uint8)
#for x in xrange(90, 110):
    #for y in xrange(90, 110):
        #accessibility[x][y] = 0


#ff.printGrid()

ff = FloorField(accessibility)
#ff = FloorField(H, W)

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

        self.mapping[0][Soldier.SWORDSMAN] = (0, 1, 0)
        self.mapping[1][Soldier.SWORDSMAN] = (0, 0, 1)

        self.mapping[0][Soldier.ARCHER] = (0, 1, 0)
        self.mapping[1][Soldier.ARCHER] = (0, 0, 1)

class IconFill(PixelFill):
    """
    Metadata for filling the image with icons.
    """
    def __init__(self):
        self.scale = 10
        self.size = 20
        self.mapping = {0: {}, 1: {}}

        self.mapping[0][Soldier.SWORDSMAN] = self.getImage('swordsman_0.png').scale(10)
        self.mapping[1][Soldier.SWORDSMAN] = self.getImage('swordsman_1.png').scale(10)

        self.mapping[0][Soldier.ARCHER] = self.getImage('archer_0.png').scale(10)
        self.mapping[1][Soldier.ARCHER] = self.getImage('archer_1.png').scale(10)

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

    def advance(self):
        """
        Advance to the next time step.
        This does move and kill alternatively.
        """
        if self.i == 0:
            pass
        elif self.i % 2 == 0:
            ff.move()
        else:
            ff.kill()
        self.i += 1

    def __call__(self, t):
        self.advance()
        # create a new surface for each frame
        field = gizeh.Surface.from_image(self.field) 
        ff.getSoldiers(self.soldiers)
        soldiers = []
        it = numpy.nditer(self.soldiers, op_flags = ['readwrite'], flags = ['multi_index'])
        while not it.finished:
            if it[0] != 255:
                xy = [it.multi_index[a] * self.fill.scale for a in (1, 0)]
                fill = self.fill.mapping[int(it[0])][Soldier.SWORDSMAN]
                soldiers.append(gizeh.square(self.fill.size, xy = xy, fill = fill))
                it[0] = -1
            it.iternext()
        gizeh.Group(soldiers).draw(field)
        return field.get_npimage()


clip = mpy.VideoClip(FrameBuilder(IconFill()), duration = 50)
#clip = mpy.VideoClip(FrameBuilder(IconFill()), duration=50)
#clip.write_videofile("battle.mp4", fps=2)
clip.write_gif("circle.gif", fps=5, opt="OptimizePlus", fuzz=10)
