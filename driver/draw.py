from pygame import *

from battlesimwrap import *
#import battlesim

import numpy
import ConfigParser
import time

blue_skill = 255
green_skill = 150
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
screen = display.set_mode((2*W, 2*H), RESIZABLE)


draw_color = 'black'
draw_type = 'field'

green_target = (H/2, W/2)
blue_target = (H/2, W/2)


# draw buttons

q = False
# two modes:
play = False
edit = True
firstit = True
i = 0
# initialize font; must be called after 'pygame.init()' to avoid 'Font not Initialized' error
myfont = font.SysFont("monospace", 15)

# render text
label = myfont.render("BattleSim 0.1", 1, (255, 0, 0))
#label_ttl = 10
label_active = False
label_starttime = time.time()
pfollow = 1.0

def draw_target(screen, pos, col):
    size = 4
    maxx = screen.get_size()[0]
    maxy = screen.get_size()[1]
    draw.line(screen, Color(col), [max(0,pos[0]-size), max(0,pos[1]-size)], [min(maxx,pos[0]+size),min(maxy,pos[1]+size)], 2)
    draw.line(screen, Color(col), [max(0,pos[0]-size), min(maxy,pos[1]+size)], [min(maxx,pos[0]+size),max(0,pos[1]-size)], 2)


def msg(text):
    global label
    global label_active
    global label_starttime
    label = myfont.render(text, 1, (255, 0, 0))
    label_starttime = time.time()
    label_active = True
    #label_ttl = 10


def pos_transform(pos):
    target = (W, H)
    relto = screen.get_size()
    result = [0, 0]
    for i in [0, 1]:
        result[i] = pos[i]*target[i] / relto[i]
    return (result[0], result[1])



#screen.blit(label, (100, 100))
backbuf = Surface((W,H))
fieldbuf = Surface((W,H))
fieldbuf.fill(Color('white'), Rect((0, 0, W, H)))
#display.update()
while not q:

    if play:
        arr = surfarray.array3d(fieldbuf)
        if firstit:
            firstit = False
            i = 0
            (bf, access) = battlefield_from_arr(arr)
            access_arr = numpy.zeros((W, H, 3), dtype=numpy.uint8)
            for x in range(0, W):
                for y in range(0, H):
                    access_arr[x,y] = numpy.array([access[x,y],access[x,y],access[x,y]])
            #print arr
            pass
        arr = surfarray.pixels3d(fieldbuf)
        bf.setFlag(0, green_target[1], green_target[0])
        bf.setFlag(1, blue_target[1], blue_target[0])

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
                arr[pos[0], pos[1]] = numpy.array([255, 0, 0])
        i = i + 1

        del arr

    # process key events
    for e in event.get():
        if e.type == KEYUP:
            print e.key
            if e.key == K_g and edit:
                msg("Draw Green army")
                draw_color = 'green'
                draw_type = 'field'
            elif e.key == K_b and edit:
                msg("Draw Blue army")
                draw_color = 'blue'
                draw_type = 'field'
            elif e.key == K_s and edit:
                msg("Draw static field")
                draw_color = 'black'
                draw_type = 'field'
            elif e.key == K_h:
                msg("Set Green Target")
                draw_type = 'target'
                draw_color = 'green'
            elif e.key == K_n:
                msg("Set Blue Target")
                draw_type = 'target'
                draw_color = 'blue'
            elif e.key == K_p and play:
                if pfollow >= 1.0:
                    pfollow = 0.0
                else:
                    pfollow += 0.2
                msg("Set Follow = %f.1"%pfollow)
                bf.setFollowPreviousProbability(pfollow)
            elif e.key == K_SPACE:
                if edit:
                    msg("Play (press SPACE to pause)")
                    edit = False
                    firstit = True
                    play = True
                else:
                    edit = True
                    firstit = True
                    play = False
                    msg("Edit Mode (press SPACE to start)")
        if e.type in (MOUSEBUTTONDOWN, MOUSEMOTION):
            if mouse.get_pressed() == (1, 0, 0):
                if draw_type == 'field' and edit:
                    fieldbuf.fill(Color(draw_color), Rect(pos_transform(mouse.get_pos()), (20, 20)))
                if draw_type == 'target':
                    # set target!
                    if draw_color == 'blue':
                        blue_target = pos_transform(mouse.get_pos())
                    if draw_color == 'green':
                        green_target = pos_transform(mouse.get_pos())
        elif e.type == VIDEORESIZE:
            screen = display.set_mode(e.dict['size'],RESIZABLE)
        if e.type == QUIT:
            q = True
    # draw
    backbuf.blit(fieldbuf, (0,0))
    draw_target(backbuf, green_target, 'lightgreen')
    draw_target(backbuf, blue_target, 'lightblue')
    transform.smoothscale(backbuf, screen.get_size(), screen)
    #screen.blit(backbuf, (0, 0))
    # draw targets!
    # draw msg?
    if label_active:
        if time.time() - label_starttime > 5:
            label_active = False
        else:
            screen.blit(label, (0, 0))
    display.update()
quit()
