import battlesim
from battlesim import Soldier 

import numpy
from random import randint, seed

import gizeh
import moviepy.editor as mpy


seed(0)

H, W = 40, 80


class PixelFill(object):
    """
    Metadata for filling the image with pixels.
    """
    def __init__(self):
        self.scale = 1
        self.size = 1
        self.mapping = {0 : {}, 1 : {}}

        self.win = {0 : (0, 1, 0), 1 : (0, 0, 1)}

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
        super(IconFill, self).__init__()
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


class BattleField(battlesim.BattleField):
    """
    Light wrapper around C++ side BattleField class.
    """
    @staticmethod
    def createField():
        accessibility = numpy.full((H, W), 255, dtype = numpy.uint8)
        #for x in xrange(90, 110):
            #for y in xrange(90, 110):
                #accessibility[x][y] = 0
        return accessibility

    def __setSoldiers(self):
        soldiers = []
        # army 0
        soldiers.extend((0 * W + y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255))) for y in xrange(0, W))
        soldiers.extend((1 * W + y, Soldier(0, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255))) for y in xrange(0, W))
        # army 1
        soldiers.extend(((H - 2) * W + y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255))) for y in xrange(0, W))
        soldiers.extend(((H - 1) * W + y, Soldier(1, Soldier.SWORDSMAN, randint(0, 255), randint(0, 255))) for y in xrange(0, W))

        self.setSoldiers(soldiers)

    def __setTarget(self):
        self.setFlag(0, H, W / 2)
        self.setFlag(1, 0, W / 2)

        self.setTarget(self.ANNIHILATE_ENEMY)

    def __init__(self, accessibility):
        super(BattleField, self).__init__(accessibility)
        self.__setTarget()
        self.__setSoldiers()
        count = self.soldierCount(0) + self.soldierCount(1)
        # list for storing current position of soldiers
        self.__soldiers = battlesim.SoldierPositionVector(count)
        # list for storing position of killed soldiers
        self.__killed = battlesim.KilledPositionVector()

    def move(self):
        super(BattleField, self).move(self.__soldiers)
        count = self.soldierCount(0) + self.soldierCount(1)
        it = iter(self.__soldiers)
        for x in xrange(count):
            yield it.next()

    def kill(self):
        count = super(BattleField, self).kill(self.__killed)
        it = iter(self.__killed)
        for x in xrange(count):
            yield it.next()


class FrameBuilder(object):
    """
    Class that builds the actual frames in the video.
    """
    def __createFieldFrame(self, accessibility):
        """
        Creates an empty frame of the field with obstructions.
        """
        obstructions = []
        it = numpy.nditer(accessibility, op_flags = ['readonly'], flags = ['multi_index'])
        while not it.finished:
            if it[0] != 255:
                color = list(it[0] / 255.0 for x in xrange(3))
                obstructions.append(gizeh.square(self.__fill.size, xy = (it.multi_index[1], it.multi_index[0]), fill = color))
            it.iternext()
        obstructions = gizeh.Group(obstructions)
        field = gizeh.Surface(W * self.__fill.scale, H * self.__fill.scale, bg_color = (1, 1, 1))
        obstructions.draw(field)
        return field.get_npimage()

    def __init__(self, fill):
        # fill specifications
        self.__fill = fill

        # setup battle
        accessibility = BattleField.createField()
        self.__bf = BattleField(accessibility)
        # iteration counter
        self.__iter = 0

        # create frame for empty field 
        self.__field = self.__createFieldFrame(accessibility)

    def move(self):
        # create new frame from the field snapshot
        self.__frame = gizeh.Surface.from_image(self.__field)
        # move soldiers and draw them in their new positions
        soldiers = []
        for pos, army in self.__bf.move():
            xy = [i * self.__fill.scale for i in (pos % W, pos / W)]
            fill = self.__fill.mapping[army][Soldier.SWORDSMAN]
            soldiers.append(gizeh.square(self.__fill.size, xy = xy, fill = fill))
        gizeh.Group(soldiers).draw(self.__frame)

    def kill(self):
        # kill soldiers and indicate killed soldiers on the previous frame
        killed = []
        for pos in self.__bf.kill():
            xy = [i * self.__fill.scale for i in (pos % W, pos / W)]
            fill = self.__fill.dead
            killed.append(gizeh.square(self.__fill.size, xy = xy, fill = fill))
        gizeh.Group(killed).draw(self.__frame)

    def __call__(self, t):
        if self.__bf.status() != BattleField.ONGOING:
            # declare that the war is finished
            gizeh.text('FIN!', 'Amiri', H, xy = ((W / 2) * self.__fill.scale, (H / 2) * self.__fill.scale)).draw(self.__frame)
        elif self.__iter % 2 == 0:
            # start off by moving
            self.move()
        else:
            # then killing
            self.kill()
        self.__iter += 1
        # return the image of the current frame 
        return self.__frame.get_npimage()

clip = mpy.VideoClip(FrameBuilder(IconFill()), duration = 80)
#clip = mpy.VideoClip(FrameBuilder(PixelFill()), duration=50)
clip.write_videofile("battle.mp4", fps=2)
#clip.write_gif("circle.gif", fps=2, opt="OptimizePlus", fuzz=10)
