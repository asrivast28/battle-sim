from pygame import *

from battlesimwrap import *
#import battlesim

import numpy
import ConfigParser


blue_skill = 255
green_skill = 0
#blue_aggression = 255



def battlefield_from_arr(image3d):
    access = numpy.full((W, H), fill_value=255, dtype=numpy.uint8)
    access_rot = numpy.copy(access)
    #field = battlesim.BattleField(access)
    # set soldiers
    soldiers = [] # battlesim.SoldierVector()
    #print (W,H)
    for x in range(0, W):
        for y in range(0, H):
            if image3d[x,y,0] == 255 and image3d[x,y,1] == 255 and image3d[x,y,2] == 255:
                # white
                pass
            elif image3d[x,y,1] == 255:
                # green army!
                soldiers.append((y * W + x, battlesim.Soldier(0, battlesim.Soldier.ARCHER, green_skill, 100)))
            elif image3d[x,y,2] == 255:
                # blue army
                soldiers.append((y * W + x, battlesim.Soldier(1, battlesim.Soldier.ARCHER, blue_skill, blue_skill)))
            else:
                access[y,x] = 0
                access_rot[x,y] = 0
    #field.setSoldiers(soldiers)
    field = BattleField(access, soldiers)
    # set target (annihilation)
    # TODO: support adding user targets
    field.setFlag(0, H-1, W-1)
    field.setFlag(1, H-1, W-1)
    field.setTarget(field.ANNIHILATE_ENEMY)
    return (field, access_rot)

# read config
config = ConfigParser.ConfigParser()
# soldier configuration
config.read('soldier.config')
Soldier.configure(config)

# battlefield configuration
config.read('battlefield.config')
BattleField.configure(config)

# init screen
init()
display.init()
W = 400
H = 400
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
            i = 0
            (bf,access) = battlefield_from_arr(arr)
            access_arr = numpy.zeros((W, H, 3),dtype=numpy.uint8)
            for x in range(0, W):
                for y in range(0, H):
                    access_arr[x,y] = numpy.array([access[x,y],access[x,y],access[x,y]])
            #print arr
            pass
        arr = surfarray.pixels3d(screen)

        if i % 2 == 0:
            numpy.copyto(arr, access_arr)
            #screen.fill(Color('white'), Rect(0, 0, W, H))
            for p, army, kind in bf.move():
                # draw!
                pos = (p % W, p / W)
                #print arr[pos[0], pos[1]]
                if army == 0:
                    arr[pos[0], pos[1]] = numpy.array([0, 255, 0])
                else:
                    arr[pos[0], pos[1]] = numpy.array([0, 0, 255])
        else:
            for p in bf.kill():
                pos = (p % W, p / W)
                arr[pos[0], pos[1]] = numpy.array([255,0,0])
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
