#!/usr/bin/env python
################################################################################
##                                                                            ##
##  hpRNA_ms2_draw.py                                                 MODULE  ##
##  ------------------------------------------------------------------------  ##
##                                                                            ##
##  v. 08/10/2014                                                             ##
##                                                                            ##
##  James Geraets, University of York                                         ##
##  jg923@york.ac.uk                                                          ##
##                                                                            ##
##  Module for generating Hamiltonian path images in 2d for bacteriophage     ##
##  ms2.                                                                      ##
##                                                                            ##
##  (C) University of York 2014                                               ##
##                                                                            ##
################################################################################

### MODULE IMPORTS

from math import sqrt, atan, sin, cos, pi
import matplotlib
import cairo
import os.path

### CONSTANTS

r3b2 = sqrt(3.0)/2.0
yoffset = 2 * r3b2 / 3
br3 = 1.0/sqrt(3.0)

### MODIFIERS

# various drawing parameters
c_rad = 0.1
scaling = 200
offset = 250, 250
linew = 10
partlinew = 6
text_size = 40
text_offsetx = 20
text_offsety = -10
gap_greater = 1

# dashes for the forced unoccupied edges
dashed_no = True
dashes_list = [10,10]

# principal colour yes
pr, pg, pb = 0, 0.7, 0

# dashed colour no
rr, rg, rb = 0.7, 0, 0

# inferred colour yes
ir, ig, ib = 0, 0, 0

# background
background = False
br, bg, bb = 1, 1, 1

# include homodimer/heterodimer colours
include_dimers = True

# labels for vertices and PS positions
include_vertices = False
include_labels = False

# scaffold for path
include_scaffold = True
sr, sg, sb = 0.75, 0.75, 0.75

### GEOMETRY SPECIFICATION

move_C5 = [('a', 't', 'N', 'u', 'd'),
           ('L', '6', 'R', 'O', 'M'),
           ('p', 'I', 'K', 's', 'r'),
           ('b', 'f', 'j', 'm', 'q'),
           ('c', 'w', 'y', 'g', 'e'),
           ('v', 'P', 'S', 'W', 'x'),
           ('Q', '2', 'Y', 'V', 'T'),
           ('G', '5', '3', '7', 'J'),
           ('l', 'F', 'H', 'o', 'n'),
           ('h', 'A', 'C', 'k', 'i'),
           ('z', 'X', 'U', '0', 'B'),
           ('D', '1', 'Z', '4', 'E')]

# assign numbers to each 12 vertices, for ease when working in 3d
vertex_names = {('a', 't', 'N', 'u', 'd'): 1,
                ('b', 'f', 'j', 'm', 'q'): 2,
                ('c', 'w', 'y', 'g', 'e'): 3,
                ('v', 'P', 'S', 'W', 'x'): 4,
                ('L', '6', 'R', 'O', 'M'): 5,
                ('p', 'I', 'K', 's', 'r'): 6,
                ('h', 'A', 'C', 'k', 'i'): 7,
                ('z', 'X', 'U', '0', 'B'): 8,
                ('Q', '2', 'Y', 'V', 'T'): 9,
                ('G', '5', '3', '7', 'J'): 10,
                ('l', 'F', 'H', 'o', 'n'): 11,
                ('D', '1', 'Z', '4', 'E'): 12}

move_DS = [('a', 'b'),
           ('c', 'd'),
           ('e', 'f'),
           ('g', 'h'),
           ('i', 'j'),
           ('k', 'l'),
           ('m', 'n'),
           ('o', 'p'),
           ('q', 'r'),
           ('s', 't'),
           ('u', 'v'),
           ('w', 'x'),
           ('y', 'z'),
           ('A', 'B'),
           ('C', 'D'),
           ('E', 'F'),
           ('G', 'H'),
           ('I', 'J'),
           ('K', 'L'),
           ('M', 'N'),
           ('O', 'P'),
           ('Q', 'R'),
           ('S', 'T'),
           ('U', 'V'),
           ('W', 'X'),
           ('Y', 'Z'),
           ('0', '1'),
           ('2', '3'),
           ('4', '5'),
           ('6', '7')]

# find the vertex number for each protein easily
protein_vertex = dict(((k, v) for (ks, v) in vertex_names.items() for k in ks))

# each face is composed of three heterodimers
face_name_list = [
                  ('a', 'f', 'c'),
                  ('d', 'w', 'v'),
                  ('u', 'P', 'M'),
                  ('N', 'L', 's'),
                  ('t', 'r', 'b'),
                  ('g', 'A', 'z'),
                  ('W', 'U', 'T'),
                  ('R', '2', '7'),
                  ('I', 'G', 'o'),
                  ('m', 'l', 'i'),
                  ('D', 'k', 'F'),
                  ('E', 'H', '5'),
                  ('4' ,'3', 'Y'),
                  ('Z', 'V', '0'),
                  ('1', 'B', 'C'),
                  ('n', 'q', 'p'),
                  ('J', 'K', '6'),
                  ('Q', 'O', 'S'),
                  ('X', 'x', 'y'),
                  ('h', 'e', 'j')]

# map each face to its position in the 2d projection
face_2d_mapping = [
                  (0.0, 0.0, False),
                  (1.0, 0.0, False),
                  (2.0, 0.0, False),
                  (3.0, 0.0, False),
                  (4.0, 0.0, False),
                  (0.5, r3b2, False),
                  (1.5, r3b2, False),
                  (2.5, r3b2, False),
                  (3.5, r3b2, False),
                  (4.5, r3b2, False),
                  (4.5, (r3b2 + br3), True),
                  (3.5, (r3b2 + br3), True),
                  (2.5, (r3b2 + br3), True),
                  (1.5, (r3b2 + br3), True),
                  (0.5, (r3b2 + br3), True),
                  (4.0, br3, True),
                  (3.0, br3, True),
                  (2.0, br3, True),
                  (1.0, br3, True),
                  (0.0, br3, True)
                  ]

### FUNCTION DEFINITIONS

def clockwise(pos, move):
    '''
    Cycle through the moves clockwise.
    
    '''
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) + 1) % len(cyc)]

def anticlockwise(pos, move):
    '''
    Cycle through the moves anticlockwise.
    
    '''
    for cyc in move:
        if pos in cyc:
            return cyc[(cyc.index(pos) - 1) % len(cyc)]
            
def m1(pos):
    '''
    Dimer switch is move 1.
    
    '''
    return clockwise(pos, move_DS)
    
def m2(pos):
    '''
    Clockwise around five-fold is move 2.
    
    '''
    return clockwise(pos, move_C5)

def m3(pos):
    '''
    Anticlockwise around five-fold is move 3.
    
    '''
    return anticlockwise(pos, move_C5)


class Val_Creator(object):
    '''
    classdocs
    '''


    def __init__(self, ratio, rotation=0):
        '''
        Constructor
        '''
        lam = 1 / sqrt((ratio*ratio) + (2*ratio) + 4)
        theta = atan(sqrt(3)/(ratio + 1))
        extend = lam*sin(pi / 3) / sin((2 * pi / 3) - theta)
        
        self.vals = dict()
        
        self.vals['vertex'] = (0,
                        yoffset)
        
        self.vals['stem_pos'] = (-lam * sin((pi / 6) - theta),
                        yoffset - (lam * cos((pi / 6) - theta)))
        
        self.vals['+3 -1'] = (extend * sin(pi / 6),
                        yoffset - (extend * cos(pi / 6)))
        
        self.vals['+1 -3'] = (- extend * sin(pi / 6),
                        yoffset - (extend * cos(pi / 6)))
        
        self.vals['+2 -2'] = (0.25,
                        yoffset - (r3b2 / 2))

        self.vals['centre'] = (0,0)

        self.vals['homodimer 1'] = ((1.0 / 3.0), 0)

        self.vals['homodimer 2'] = ((1.0 / 3.0) * sin(pi / 6.0), (1.0 / 3.0) * cos(pi / 6.0))
    
        self.vals['homodimer mid'] = ((self.vals['homodimer 1'][0] + self.vals['homodimer 2'][0]) / 2.0, (self.vals['homodimer 1'][1] + self.vals['homodimer 2'][1]) / 2.0)

        self.vals['homodimer split'] = (sum([self.vals['homodimer 1'][0], self.vals['centre'][0], self.vals['centre'][0]]) / 3.0, sum([self.vals['homodimer 1'][1], self.vals['centre'][1], self.vals['centre'][1]]) / 3.0)

        self.vals['hetdimer 1'] = ((1.0 / 3.0) * sin(pi / 6.0), (1.0 / 3.0) * cos(pi / 6.0))

        self.vals['hetdimer 2'] = (-(1.0 / 3.0) * sin(pi / 6.0), (1.0 / 3.0) * cos(pi / 6.0))

        self.vals['hetdimer split 1'] = (sum([self.vals['hetdimer 1'][0], self.vals['centre'][0], self.vals['centre'][0]]) / 3.0, sum([self.vals['hetdimer 1'][1], self.vals['centre'][1], self.vals['centre'][1]]) / 3.0)

        self.vals['hetdimer split 2'] = (sum([self.vals['hetdimer 2'][0], self.vals['vertex'][0], self.vals['vertex'][0]]) / 3.0, sum([self.vals['hetdimer 2'][1], self.vals['vertex'][1], self.vals['vertex'][1]]) / 3.0)
      
        for i in range(rotation * 2):
            self.rotate_grid()
            
    def rotate_grid(self):
        for i in self.vals:
            j, face = self.vals[i]
            m = (0.5 * j) - (r3b2 * face)
            n = (r3b2 * j) + (0.5 * face)
            self.vals[i] = (m, n)

    def flip_grid(self):
        for i in range(3):
            self.rotate_grid()

    def translate(self, xt, yt):
        valslist = self.vals.items()
        for i, (x,y) in valslist:
            self.vals[i] = (x + xt, y + yt)

def give_face_name(v):
    for f in face_name_list:
        if v in f:
            return f

def draw_scaffold():
    '''
    Draw the baseline of the MS2 architecture to lay the path description above.
    
    '''
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 1400, 800)
    if background:
        bg = cairo.Context(surface)
        bg.set_source_rgb(1,1,1)
        bg.paint()
    
    draw_dict = dict()
    for ind, face in enumerate(face_name_list):
        f0, f1, f2 = face
        g0 = Val_Creator(1.95, rotation=0)
        g1 = Val_Creator(1.95, rotation=2)
        g2 = Val_Creator(1.95, rotation=1)
        if not face_2d_mapping[ind][2]:
            g0.flip_grid()
            g1.flip_grid()
            g2.flip_grid()
        g0.translate(face_2d_mapping[ind][0], face_2d_mapping[ind][1])
        g1.translate(face_2d_mapping[ind][0], face_2d_mapping[ind][1])
        g2.translate(face_2d_mapping[ind][0], face_2d_mapping[ind][1])
        g0.position = f0
        g1.position = f1
        g2.position = f2
        g0.face = face
        g1.face = face
        g2.face = face
        draw_dict[f0] = g0
        draw_dict[f1] = g1
        draw_dict[f2] = g2

    for face in face_name_list:
        f0, f1, f2 = face
        vl = cairo.Context(surface)
        v0x, v0y = draw_dict[f0].vals['vertex']
        v1x, v1y = draw_dict[f1].vals['vertex']
        v2x, v2y = draw_dict[f2].vals['vertex']
        vl.set_source_rgb(0.3, 0.3, 0.3)
        vl.move_to(offset[0] + scaling*v0x, offset[1] + scaling*v0y)
        vl.line_to(offset[0] + scaling*v1x, offset[1] + scaling*v1y)
        vl.line_to(offset[0] + scaling*v2x, offset[1] + scaling*v2y)
        vl.line_to(offset[0] + scaling*v0x, offset[1] + scaling*v0y)
        vl.stroke()


    if include_dimers:
        # include the colourful hetero and homodimers
        for e1, e2 in move_DS:
            testc = cairo.Context(surface)
            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer split 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer 1'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['vertex'][0], offset[1] + scaling*draw_dict[e1].vals['vertex'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer split 2'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer split 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer split 1'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(0.8, 1, 0.8)
            testc.fill()
            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['centre'][0], offset[1] + scaling*draw_dict[e1].vals['centre'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer split 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer split 2'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer split 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['hetdimer 2'][0], offset[1] + scaling*draw_dict[e1].vals['hetdimer 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['centre'][0], offset[1] + scaling*draw_dict[e1].vals['centre'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(0.8, 0.8, 1)
            testc.fill()

            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer split 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer 1'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['vertex'][0], offset[1] + scaling*draw_dict[e2].vals['vertex'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer split 2'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer split 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer split 1'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(0.8, 1, 0.8)
            testc.fill()
            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['centre'][0], offset[1] + scaling*draw_dict[e2].vals['centre'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer split 1'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer split 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer split 2'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer split 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['hetdimer 2'][0], offset[1] + scaling*draw_dict[e2].vals['hetdimer 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['centre'][0], offset[1] + scaling*draw_dict[e2].vals['centre'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(0.8, 0.8, 1)
            testc.fill()

            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['homodimer 2'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['centre'][0], offset[1] + scaling*draw_dict[e1].vals['centre'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['homodimer split'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer split'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['homodimer mid'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer mid'][1])
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['homodimer 2'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer 2'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(1, 0.8, 0.8)
            testc.fill()

            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['homodimer 1'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['homodimer split'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer split'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e1].vals['homodimer mid'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer mid'][1])
            testc.move_to(offset[0] + scaling*draw_dict[e1].vals['homodimer 1'][0], offset[1] + scaling*draw_dict[e1].vals['homodimer 1'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(1, 0.8, 0.8)
            testc.fill()

            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['homodimer 2'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer 2'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['centre'][0], offset[1] + scaling*draw_dict[e2].vals['centre'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['homodimer split'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer split'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['homodimer mid'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer mid'][1])
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['homodimer 2'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer 2'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(1, 0.8, 0.8)
            testc.fill()

            testc.set_source_rgb(0, 0, 0)
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['homodimer 1'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer 1'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['homodimer split'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer split'][1])
            testc.line_to(offset[0] + scaling*draw_dict[e2].vals['homodimer mid'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer mid'][1])
            testc.move_to(offset[0] + scaling*draw_dict[e2].vals['homodimer 1'][0], offset[1] + scaling*draw_dict[e2].vals['homodimer 1'][1])
            testc.stroke_preserve()
            testc.set_source_rgb(1, 0.8, 0.8)
            testc.fill()

        testc.set_source_rgb(0, 0, 0)
        testc.move_to(offset[0] + scaling*draw_dict['L'].vals['homodimer 2'][0], offset[1] + scaling*draw_dict['L'].vals['homodimer 2'][1])
        testc.line_to(offset[0] + scaling*draw_dict['L'].vals['centre'][0], offset[1] + scaling*draw_dict['L'].vals['centre'][1])
        testc.line_to(offset[0] + scaling*draw_dict['L'].vals['homodimer 1'][0], offset[1] + scaling*draw_dict['L'].vals['homodimer 1'][1])
        testc.move_to(offset[0] + scaling*draw_dict['L'].vals['homodimer 2'][0], offset[1] + scaling*draw_dict['L'].vals['homodimer 2'][1])
        testc.stroke_preserve()
        testc.set_source_rgb(0.2, 0.8, 0.8)
        testc.fill()

        testc.set_source_rgb(0, 0, 0)
        testc.move_to(offset[0] + scaling*draw_dict['K'].vals['homodimer 2'][0], offset[1] + scaling*draw_dict['K'].vals['homodimer 2'][1])
        testc.line_to(offset[0] + scaling*draw_dict['K'].vals['centre'][0], offset[1] + scaling*draw_dict['K'].vals['centre'][1])
        testc.line_to(offset[0] + scaling*draw_dict['K'].vals['homodimer 1'][0], offset[1] + scaling*draw_dict['K'].vals['homodimer 1'][1])
        testc.move_to(offset[0] + scaling*draw_dict['K'].vals['homodimer 2'][0], offset[1] + scaling*draw_dict['K'].vals['homodimer 2'][1])
        testc.stroke_preserve()
        testc.set_source_rgb(0.2, 0.8, 0.8)
        testc.fill()

        testc.set_source_rgb(0.3, 0.3, 0.3)
        coordmp = (offset[0] + (scaling*draw_dict['L'].vals['vertex'][0] + scaling*draw_dict['K'].vals['vertex'][0])/2.0, offset[1] +(scaling*draw_dict['L'].vals['vertex'][1] + scaling*draw_dict['K'].vals['vertex'][1])/2.0 )
        testc.set_font_size(30)
        testc.move_to(coordmp[0] - 12, coordmp[1] - 10)
        testc.show_text('M')
        testc.move_to(coordmp[0] - 9, coordmp[1] + 35)
        testc.show_text('P')
        testc.fill()


    if include_vertices:
        # label the vertices from 1 to 12
        for i,j in draw_dict.items():
            cc = cairo.Context(surface)
            cr = cairo.Context(surface)
            cc.set_source_rgb(1, 1, 1)
            cr.set_source_rgb(0,0,0)
            x, y = j.vals['vertex']
            #cc.move_to(offset[0] + (scaling*x), offset[1] + (scaling*(y + c_rad)))
            cc.arc(offset[0] + (scaling*x), offset[1] + (scaling*y), scaling*c_rad*1.5, 0, 2 * pi)
            cr.arc(offset[0] + (scaling*x), offset[1] + (scaling*y), scaling*c_rad*1.5, 0, 2 * pi)
            cc.fill()
            cc.stroke()
            cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            cr.set_font_size(text_size)
            (ex, ey, ewidth, eheight, edx, edy) = cr.text_extents(str(protein_vertex[i[0]])) 
            cr.move_to(offset[0] + scaling*x - ((ewidth / 2.0) + 2.5), offset[1] + scaling*y + (eheight / 2.0))
            cr.show_text(str(protein_vertex[i[0]]))
            cr.stroke()

    if include_scaffold:
        for e1, e2 in move_DS:
            st1x, st1y = draw_dict[e1].vals['stem_pos']
            p1x, p1y = draw_dict[e1].vals['+2 -2']
            p2x, p2y = draw_dict[e2].vals['+2 -2']
            st2x, st2y = draw_dict[e2].vals['stem_pos']
            sc = cairo.Context(surface)
            sc.set_line_width(linew)
            sc.set_source_rgb(sr, sg, sb)
            sc.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
            sc.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
            sc.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
            sc.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
            sc.stroke()
            sc.arc(offset[0] + scaling*st1x, offset[1] + scaling*st1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p1x, offset[1] + scaling*p1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p2x, offset[1] + scaling*p2y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*st2x, offset[1] + scaling*st2y, linew/2.0, 0, 2 * pi)
            sc.fill()

        for e1, e2 in move_DS:
            e0 = m1(e1)
            st1x, st1y = draw_dict[e0].vals['stem_pos']
            p1x, p1y = draw_dict[e0].vals['+3 -1']
            p2x, p2y = draw_dict[e0].vals['+1 -3']
            st2x, st2y = draw_dict[e0].vals['stem_pos']
            sc = cairo.Context(surface)
            sc.set_line_width(linew)
            sc.set_source_rgb(sr, sg, sb)
            sc.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
            sc.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
            sc.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
            sc.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
            sc.stroke()
            sc.arc(offset[0] + scaling*st1x, offset[1] + scaling*st1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p1x, offset[1] + scaling*p1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p2x, offset[1] + scaling*p2y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*st2x, offset[1] + scaling*st2y, linew/2.0, 0, 2 * pi)
            sc.fill()
            
        for e1, e2 in move_DS:
            e0 = m1(e2)
            st1x, st1y = draw_dict[e0].vals['stem_pos']
            p1x, p1y = draw_dict[e0].vals['+3 -1']
            p2x, p2y = draw_dict[e0].vals['+1 -3']
            st2x, st2y = draw_dict[e0].vals['stem_pos']
            sc = cairo.Context(surface)
            sc.set_line_width(linew)
            sc.set_source_rgb(sr, sg, sb)
            sc.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
            sc.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
            sc.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
            sc.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
            sc.stroke()
            sc.arc(offset[0] + scaling*st1x, offset[1] + scaling*st1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p1x, offset[1] + scaling*p1y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*p2x, offset[1] + scaling*p2y, linew/2.0, 0, 2 * pi)
            sc.fill()
            sc.arc(offset[0] + scaling*st2x, offset[1] + scaling*st2y, linew/2.0, 0, 2 * pi)
            sc.fill()

    #if not os.path.exists('draw'):
    #    os.makedirs('draw')
    #surface.write_to_png(os.path.join('draw', 'scaffold.png'))
    
    return surface, draw_dict

def hami_draw(constrain_occ, constrain_unocc, draw, name, movestr, protstr):
    surface, draw_dict = draw_scaffold()
    
    for (e1, e2), ev in constrain_occ:
        if m1(e1) == e2:
            st1x, st1y = draw_dict[e1].vals['stem_pos']
            p1x, p1y = draw_dict[e1].vals['+2 -2']
            p2x, p2y = draw_dict[e2].vals['+2 -2']
            st2x, st2y = draw_dict[e2].vals['stem_pos']
        elif m2(e1) == e2:
            st1x, st1y = draw_dict[e1].vals['stem_pos']
            p1x, p1y = draw_dict[e1].vals['+3 -1']
            p2x, p2y = draw_dict[e2].vals['+1 -3']
            st2x, st2y = draw_dict[e2].vals['stem_pos']
        elif m3(e1) == e2:
            st1x, st1y = draw_dict[e2].vals['stem_pos']
            p1x, p1y = draw_dict[e2].vals['+3 -1']
            p2x, p2y = draw_dict[e1].vals['+1 -3']
            st2x, st2y = draw_dict[e1].vals['stem_pos']
        br = cairo.Context(surface)
        br.set_line_width(linew)
        br.set_source_rgb(pr, pg, pb)
        br.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
        br.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
        br.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
        br.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
        br.stroke()
        br.arc(offset[0] + scaling*st1x, offset[1] + scaling*st1y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*p1x, offset[1] + scaling*p1y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*p2x, offset[1] + scaling*p2y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*st2x, offset[1] + scaling*st2y, linew/2.0, 0, 2 * pi)
        br.fill()


    for (e1, e2), ev in constrain_unocc:
        if dashed_no:
            if m1(e1) == e2:
                st1x, st1y = draw_dict[e1].vals['stem_pos']
                p1x, p1y = draw_dict[e1].vals['+2 -2']
                p2x, p2y = draw_dict[e2].vals['+2 -2']
                st2x, st2y = draw_dict[e2].vals['stem_pos']
            elif m2(e1) == e2:
                st1x, st1y = draw_dict[e1].vals['stem_pos']
                p1x, p1y = draw_dict[e1].vals['+3 -1']
                p2x, p2y = draw_dict[e2].vals['+1 -3']
                st2x, st2y = draw_dict[e2].vals['stem_pos']
            elif m3(e1) == e2:
                st1x, st1y = draw_dict[e2].vals['stem_pos']
                p1x, p1y = draw_dict[e2].vals['+3 -1']
                p2x, p2y = draw_dict[e1].vals['+1 -3']
                st2x, st2y = draw_dict[e1].vals['stem_pos']
            br = cairo.Context(surface)
            br.set_line_width(linew)
            br.set_dash(dashes_list)
            br.set_source_rgb(rr, rg, rb)
            br.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
            br.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
            br.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
            br.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
            br.stroke()
            
    for (e1, e2) in draw:
        if m1(e1) == e2:
            st1x, st1y = draw_dict[e1].vals['stem_pos']
            p1x, p1y = draw_dict[e1].vals['+2 -2']
            p2x, p2y = draw_dict[e2].vals['+2 -2']
            st2x, st2y = draw_dict[e2].vals['stem_pos']
        elif m2(e1) == e2:
            st1x, st1y = draw_dict[e1].vals['stem_pos']
            p1x, p1y = draw_dict[e1].vals['+3 -1']
            p2x, p2y = draw_dict[e2].vals['+1 -3']
            st2x, st2y = draw_dict[e2].vals['stem_pos']
        elif m3(e1) == e2:
            st1x, st1y = draw_dict[e2].vals['stem_pos']
            p1x, p1y = draw_dict[e2].vals['+3 -1']
            p2x, p2y = draw_dict[e1].vals['+1 -3']
            st2x, st2y = draw_dict[e1].vals['stem_pos']
        br = cairo.Context(surface)
        br.set_line_width(linew)
        br.set_source_rgb(ir, ig, ib)
        br.move_to(offset[0] + scaling*st1x, offset[1] + scaling*st1y)
        br.line_to(offset[0] + scaling*p1x, offset[1] + scaling*p1y)
        br.move_to(offset[0] + scaling*p2x, offset[1] + scaling*p2y)
        br.line_to(offset[0] + scaling*st2x, offset[1] + scaling*st2y)
        br.stroke()
        br.arc(offset[0] + scaling*st1x, offset[1] + scaling*st1y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*p1x, offset[1] + scaling*p1y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*p2x, offset[1] + scaling*p2y, linew/2.0, 0, 2 * pi)
        br.fill()
        br.arc(offset[0] + scaling*st2x, offset[1] + scaling*st2y, linew/2.0, 0, 2 * pi)
        br.fill()

                
    if include_labels:
        cc = cairo.Context(surface)
        cr = cairo.Context(surface)
        cc.set_source_rgb(1, 1, 1)
        cr.set_source_rgb(0,0,0)
        cr.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(text_size/2.0)
        for v in move_C5:
            for p in v:
                x, y = draw_dict[p].vals['stem_pos']
                cc.arc(offset[0] + (scaling*x), offset[1] + (scaling*y), scaling*c_rad*0.7, 0, 2 * pi)
                cr.arc(offset[0] + (scaling*x), offset[1] + (scaling*y), scaling*c_rad*0.7, 0, 2 * pi)
                cc.fill()
                cc.stroke()
                (ex, ey, ewidth, eheight, edx, edy) = cr.text_extents(p) 
                cr.move_to(offset[0] + scaling*x - ((ewidth / 2.0) + 2.5), offset[1] + scaling*y + (eheight / 2.0))
                cr.show_text(p)
                cr.stroke()
                

    tx = cairo.Context(surface)
    tx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    tx.set_font_size(text_size/2.0)
    tx.move_to(55, 50)
    tx.show_text(movestr)
    tx.move_to(50, 75)
    tx.show_text(protstr)
    tx.stroke()
   
    surface.write_to_png(name)


## MAIN - not usually to be run

#if __name__ == '__main__':
#    longdict = dict()
#    shortdict = dict()
#    inf_longdict = dict()
#    inf_shortdict = dict()    
#    
#    hami_draw(longdict, shortdict, inf_longdict, inf_shortdict, 'blank_path')

### END OF MODULE