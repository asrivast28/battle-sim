from pygame import *

from battlesimwrap import *

import numpy



def battlefield_from_arr(image3d):

    


# init screen
init()
display.init()
W = 20
H = 20
screen = display.set_mode((W, H))

display.update(screen.fill(Color('white'), Rect((0, 0, W, H))))

draw_color = 'black'
# draw buttons

q = False
# two modes:
play = False
edit = True
firstit = True
while not q:

    if play:
        if firstit:
            firstit = False
            #arr = surfarray.array2d(screen)
            arr = surfarray.array3d(screen)
            bf = battlefield_from_arr(arr)
            print arr
            pass
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
