import battlesim
import numpy
from random import normalvariate, seed


class Soldier(battlesim.Soldier):
    randomParameters = {}

    @classmethod
    def configure(cls, config):
        configMap = battlesim.SoldierTypeMap()
        for section in ('KillRadius', 'DynamicField'):
            options = config.options(section)
            for option in options:
                configMap[getattr(cls, option.upper())] = int(config.get(section, option))
            getattr(cls, 'set%sMap'%(section))(configMap)
            configMap.clear()

        for army in xrange(2):
            section = str(army)
            options = config.options(section)
            thisDict = {}
            for kind in ('ARCHER', 'SWORDSMAN', 'LEADER'): 
                thisParams = [float(config.get(section, kind.title() + param)) for param in ('Mu', 'Sigma')]
                # assert six sigma rule
                assert (thisParams[0] + 6 * thisParams[1] <= 1.0) and (thisParams[0] - 6 * thisParams[1] >= 0.0)
                thisDict[getattr(cls, kind)] = thisParams
            cls.randomParameters[army] = thisDict

    @classmethod
    def normalRandom(cls, army, kind, maximum):
        return int(round(normalvariate(*cls.randomParameters[army][kind]) * maximum))

    def __init__(self, army, kind):
        skill = self.normalRandom(army, kind, 255)
        aggression = self.normalRandom(army, kind, 255)
        super(Soldier, self).__init__(army, kind, skill, aggression)


class BattleField(battlesim.BattleField):
    """
    Light wrapper around C++ side BattleField class.
    """
    @classmethod
    def configure(cls, config):
        for param in ('ExtendedNeighborhoodSize', 'DynamicFieldDecayFactor'):
            getattr(cls, 'set%s'%(param))(int(config.get('Common', param)))
        for param in ('FollowPreviousProbability', ):
            getattr(cls, 'set%s'%(param))(float(config.get('Common', param)))

    def __setSoldiers(self, H, W):
        soldiers = []
        # army 0
        soldiers.extend((0 * W + y, Soldier(0, Soldier.ARCHER)) for y in xrange(0, W))
        soldiers.extend((1 * W + y, Soldier(0, Soldier.SWORDSMAN)) for y in xrange(0, W))
        soldiers.extend((4 * W + y, Soldier(0, Soldier.LEADER)) for y in xrange(0, W, 10))
        # army 1
        soldiers.extend(((H - 5) * W + y, Soldier(1, Soldier.LEADER)) for y in xrange(0, W, 10))
        soldiers.extend(((H - 2) * W + y, Soldier(1, Soldier.SWORDSMAN)) for y in xrange(0, W))
        soldiers.extend(((H - 1) * W + y, Soldier(1, Soldier.ARCHER)) for y in xrange(0, W))

        self.setSoldiers(soldiers)

    def __setTarget(self, H, W):
        self.setFlag(0, H, W / 2)
        self.setFlag(1, 0, W / 2)

        self.setTarget(self.ANNIHILATE_ENEMY)

    def __init__(self, accessibility, soldiers=None, initTarget=False):
        """
        Initializes the battlefield with given accessibility matrix.
        Also sets up soldiers and targets.
        """
        super(BattleField, self).__init__(accessibility)
        H, W = accessibility.shape
        if soldiers:
            self.setSoldiers(soldiers)
        else:
            self.__setSoldiers(H, W)
        if initTarget:
            self.__setTarget(H, W)
        count = self.soldierCount(0) + self.soldierCount(1)
        # list for storing current position of soldiers
        self.__soldiers = battlesim.SoldierPositionVector(count)
        # list for storing position of killed soldiers
        self.__killed = battlesim.KilledPositionVector()

    def move(self):
        """
        Returns an iterator to soldier positions after the move.
        """
        super(BattleField, self).move(self.__soldiers)
        count = self.soldierCount(0) + self.soldierCount(1)
        it = iter(self.__soldiers)
        for x in xrange(count):
            pos, info = it.next()
            army = (info & (1 << 7)) >> 7
            kind = info ^ (army << 7)
            yield pos, army, kind

    def kill(self):
        """
        Returns an iterator to killed soldiers after the kill.
        """
        count = super(BattleField, self).kill(self.__killed)
        it = iter(self.__killed)
        for x in xrange(count):
            yield it.next()
