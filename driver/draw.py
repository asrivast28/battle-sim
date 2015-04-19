from pygame import *

#from battlesimwrap import *
import battlesim

import numpy
import ConfigParser


skill = 100
aggression = 100


def battlefield_from_arr(image3d):
    access = numpy.full((W, H), fill_value=255, dtype=numpy.uint8)
    field = battlesim.BattleField(access)
    # set soldiers
    soldiers = [] # battlesim.SoldierVector()
    print (W,H)
    for x in range(0, W):
        for y in range(0, H):
            if image3d[x,y,0] == 255 and image3d[x,y,1] == 255 and image3d[x,y,2] == 255:
                # white
                pass
            elif image3d[x,y,1] == 255:
                # green army!
                soldiers.append((x * W + y, battlesim.Soldier(0, battlesim.Soldier.ARCHER, skill, aggression)))
            elif image3d[x,y,2] == 255:
                # blue army
                soldiers.append((x * W + y, battlesim.Soldier(1, battlesim.Soldier.ARCHER, skill, aggression)))
    field.setSoldiers(soldiers)
    # set target (annihilation)
    # TODO: support adding user targets
    field.setFlag(0, H, W / 2)
    field.setFlag(1, 0, W / 2)
    field.setTarget(field.ANNIHILATE_ENEMY)
    return field

# read config
# soldier configuration
config.read('soldier.config')
#Soldier.configure(config)

# battlefield configuration
config.read('battlefield.config')
#BattleField.configure(config)

# init screen
init()
display.init()
W = 40
H = 40
screen = display.set_mode((W, H))

display.update(screen.fill(Color('white'), Rect((0, 0, W, H))))

draw_color = 'black'
# draw buttons

q = False
# two modes:
play = False
edit = True
firstit = True
i = 0
while not q:

    if play:
        arr = surfarray.array3d(screen)
        if firstit:
            firstit = False
            #arr = surfarray.array2d(screen)
            bf = battlefield_from_arr(arr)
            #print arr
            i = 0
            pass
        soldier_count = bf.soldierCount(0) + bf.soldierCount(1)
        moved_soldiers = battlesim.SoldierPositionVector(soldier_count)
        if i % 2 == 0:
            bf.move(moved_soldiers)
            soldier_count = bf.soldierCount(0) + bf.soldierCount(1)
            screen.fill(Color('white'), Rect(0, 0, W, H))
            arr = surfarray.pixels3d(screen)
            it = iter(moved_soldiers)
            for s in xrange(soldier_count):
                pos, info = it.next()
                army = (info & (1 << 7)) >> 7
                kind = info ^ (army << 7)
                # draw!
                if army == 0:
                    arr[pos[0], pos[1]] = numpy.array([0, 255, 0])
                else:
                    arr[pos[0], pos[1]] = numpy.array([0, 0, 255])
        else:
            bf.kill()
        i = i + 1

        display.update()

        for e in event.get():
            if e.type == KEYUP:
                print e.key
                if e.key == K_SPACE:
                    edit = True
                    firstit = True
                    play = False
            if e.type in (MOUSEBUTTONDOWN, MOUSEMOTION):
                if mouse.get_pressed() == (1, 0, 0):
                    screen.fill(Color(draw_color), Rect(mouse.get_pos(), (20, 20)))
                    display.update()
            if e.type == QUIT:
                q = True

            # get matrix and make into battlefield
        # create new frame

    if edit:
        for e in event.get():
            if e.type == KEYUP:
                print e.key
                if e.key == K_r:
                    draw_color = 'red'
                elif e.key == K_g:
                    draw_color = 'green'
                elif e.key == K_b:
                    draw_color = 'blue'
                elif e.key == K_s:
                    draw_color = 'black'
                elif e.key == K_SPACE:
                    edit = False
                    firstit = True
                    play = True
            if e.type in (MOUSEBUTTONDOWN, MOUSEMOTION):
                if mouse.get_pressed() == (1, 0, 0):
                    screen.fill(Color(draw_color), Rect(mouse.get_pos(), (20, 20)))
                    display.update()
            if e.type == QUIT:
                q=True
quit()
