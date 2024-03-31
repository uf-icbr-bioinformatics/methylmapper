# (c) 2016, A. Riva, ariva@ufl.edu

import sys

class ColorMap():
    drawer = None
    cmap = {}
    colors = {}

    def __init__(self, d):
        self.cmap = {}
        self.drawer = d

    def allocate(self, name, r, g, b):
        idx = int(self.drawer.colorAllocate(r, g, b))
        # print "Alloc: {},{},{} => {}".format(r, g, b, idx)
        self.cmap[name] = idx
        self.colors[idx] = (r, g, b)
        return idx

    def decode(self, name):
        #print "{}, {}".format(name, self.cmap)
        #raw_input()
        if name in self.cmap:
            return self.cmap[name]
        else:
            return None

class StandardColorMap(ColorMap):
    """Some standard colors, taken from CSS 2."""
    stdcolors = ['255 255 255		white',
                 '  0   0   0		black',
                 '255   0   0		red',
                 '  0 255   0           green',
                 '  0   0 255           blue',
                 '128   0   0           darkRed',
                 '  0 128   0           darkGreen',
                 '  0   0 128           darkBlue',
                 '255 255   0           yellow',
                 '255   0 255           fuchsia',
                 '  0 255 255           aqua',
                 '128 128   0           olive',
                 '128   0 128           purple',
                 '  0 128 128           teal',
                 '128 128 128           grey',
                 '192 192 192           silver',
                 '255 165   0           orange',
                 '255 255 204           lightyellow']

    def __init__(self, d, debug=False):
        self.cmap = {}
        self.drawer = d
        for sc in self.stdcolors:
            sc = sc.split()
            if debug: print("({}, {}, {})".format(int(sc[0]), int(sc[1]), int(sc[2])))
            v = self.allocate(sc[3], int(sc[0]), int(sc[1]), int(sc[2]))
            if debug: print(v)

class RangeColorMap(ColorMap):
    ncolors = 0
    minc = None
    maxc = 0

    def allocate(self, n, red1, green1, blue1, red2, green2, blue2):
        self.ncolors += n
        crange = self.drawer.colorRangeAllocate(n, red1, green1, blue1, red2, green2, blue2)
        if self.minc == None:
            self.minc = crange[0]
        self.maxc = crange[1]
        return crange

    def decode(self, v):
        """Decode a value `v' in the range [0, 1] to a color index in this colormap."""
        idx = int(round(self.minc + v * (self.maxc - self.minc)))
        #sys.stderr.write("{} >>> {}\n".format(v, idx))
        return idx

