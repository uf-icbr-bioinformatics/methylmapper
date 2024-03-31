import colormaps

## Drawing primitives

class Drawer():
    cr = None                   # Creator (connection to gdcreate)
    cm = None                   # Default colormap

    def __init__(self, creator):
        self.cr = creator
        self.cm = colormaps.ColorMap(self)

    def setColormap(self, colormap):
        self.cm = colormap

    def color(self, c):
        #print "Decoding color {}".format(c)
        if isinstance(c, int):
            return c
        else:
            return self.cm.decode(c)

    def terminate(self):
        self.cr.terminate()

    def createImage(self, width, height):
        return self.cr.sendCommand("CR", width, height)

    def colorAllocate(self, red, green, blue):
        return self.cr.sendCommand("CA", red, green, blue)

    def colorRangeAllocate(self, n, red1, green1, blue1, red2, green2, blue2):
        """Allocate `n' different colors ranging from (red1, green1, blue1) to (red2, green2, blue2). Returns the index of the
first and last allocated colors."""
        result = self.cr.sendCommand("C*", n, red1, green1, blue1, red2, green2, blue2)
        crange = result.split("-")
        return (int(crange[0]), int(crange[1]))

    def drawPixel(self, x, y, color):
        return self.cr.sendCommand("PI", x, y, self.color(color))

    def drawManyPixels(self, pixels, color):
        words = ["P*", len(pixels), color]
        for p in pixels:
            words.append(p[0])
            words.append(p[1])
        return self.cr.sendCommand(*words)

    def drawRectangle(self, x1, y1, x2, y2, color):
        return self.cr.sendCommand("RE", x1, y1, x2, y2, self.color(color))

    def drawFilledRectangle(self, x1, y1, x2, y2, color):
        return self.cr.sendCommand("RF", x1, y1, x2, y2, self.color(color))

    def drawPolygon(self, coordinates, color):
        return self.cr.sendCommandList(["PO", len(coordinates) / 2] + coordinates + [self.color(color)])

    def drawFilledPolygon(self, coordinates, color):
        args = ["PF", len(coordinates) / 2] + coordinates + [self.color(color)]
        return self.cr.sendCommandList(args)

    def drawDot(self, x, y, color, size=1):
        return self.cr.sendCommand("DO", x, y, size, self.color(color))

    def drawManyDots(self, dots, color, size=1):
        words = ["D*", len(dots), size, color]
        for p in dots:
            words.append(p[0])
            words.append(p[1])
        return self.cr.sendCommand(*words)

    def drawLine(self, x1, y1, x2, y2, color):
        # print "Drawing line: {}".format((x1, y1, x2, y2, color))
        return self.cr.sendCommand("LI", x1, y1, x2, y2, self.color(color))

    def drawArrow(self, x1, y1, x2, y2, size, color):
        # print "Drawing line: {}".format((x1, y1, x2, y2, color))
        return self.cr.sendCommand("AW", x1, y1, x2, y2, size, self.color(color))

    def setThickness(self, th):
        return self.cr.sendCommand("TH", th)

    def setStyle(self, colors):
        ncolors = len(colors)
        colors = [ self.color(c) for c in colors ]
        return self.cr.sendCommand("SS", ncolors, *colors)

    def setViewport(self, bx1, by1, bx2, by2, vx1, vy1, vx2, vy2):
        return self.cr.sendCommand("VI", bx1, by1, bx2, by2, vx1, vy1, vx2, vy2)

    def cancelViewport(self):
        return self.cr.sendCommand("VO")

    def setFont(self, font):
        """Sets the current font to `font'. If `font' is a number, it is interpreted as one of the
five built-in GD fonts. Otherwise it should be the name of a TrueType font."""
        # print "setting font to {}".format(font)
        if isinstance(font, str):
            return self.cr.sendCommand("SF", '0', font)
        else:
            return self.cr.sendCommand("SF", str(font))

    def drawString(self, string, x, y, color, anchor=1):
        # print "GD printing `{}' at {},{}, anchor={}".format(string, x, y, anchor)
        return self.cr.sendCommand("ST", x, y, self.color(color), anchor, string)

    def drawStringFT(self, string, x, y, color, anchor=1, pointsize=10.0, angle=0.0):
        # print "FT printing `{}' at {},{}, anchor={}".format(string, x, y, anchor)
        return self.cr.sendCommand("S*", x, y, self.color(color), anchor, pointsize, angle, string)

    def saveImage(self, filename):
        return self.cr.sendCommand("SA", filename)

    def close(self):
        return self.cr.sendCommand("ZZ")

    # Special purpose
    def drawGene(self, start, end, y, strand, color, smallboxes, largeboxes):
        words = ['GE', str(start), str(end), str(y), str(strand), str(color)]
        words.append(str(len(smallboxes)))
        for b in smallboxes:
            words.append(str(b[0]))
            words.append(str(b[1]))
        words.append(str(len(largeboxes)))
        for b in largeboxes:
            words.append(str(b[0]))
            words.append(str(b[1]))
        return self.cr.sendCommand(*words)
    
## High-level objects

class Point():
    x = None
    y = None

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return "P({}, {})".format(self.x, self.y)

###
###   bp1---------+        +---------vp2
###   |           |        |           |
###   +---------bp2        vp1---------+
###

class Viewport():
    vp1 = None                  # Bottom left (virtual)
    vp2 = None                  # Top right (virtual)
    bp1 = None                  # Top left (real)
    bp2 = None                  # Bottom right (real)
    vd = None                   # Dimensions of viewport (virtual)
    bd = None                   # Dimensions of viewport (real)

    def __init__(self, bp1, bp2, vp1, vp2):
        self.vp1 = vp1
        self.vp2 = vp2
        self.bp1 = bp1
        self.bp2 = bp2
        self.vd = Point(vp2.x - vp1.x, vp2.y - vp1.y)
        self.bd = Point(bp2.x - bp1.x, bp2.y - bp1.y)

    def virtualToReal(self, x, y):
        """Convert coordinates (x, y) in virtual space into a point in real (bitmap) space."""
        return Point(int(self.bp1.x + (1.0 * (x - self.vp1.x) / self.vd.x) * self.bd.x),
                     int(self.bp1.y + (1.0 * (self.vp2.y - y) / self.vd.y) * self.bd.y))

    def realToVirtual(self, x, y):
        """Convert coordinates (x, y) in bitmap space into a point in virtual space."""
        # print "RtoV {}: {}".format(y, self.vp1.y + (1.0 * (self.bp2.y - y) / self.bd.y) * self.vd.y)
        return Point(self.vp1.x + (1.0 * (x - self.bp1.x) / self.bd.x) * self.vd.x,
                     self.vp1.y + (1.0 * (self.bp2.y - y) / self.bd.y) * self.vd.y)

###
###   p1---------+
###   |          |
###   +---------p2
###

class plotPanel():
    p1 = None                   # Top left of rectangle
    p2 = None                   # Bottom right of rectangle
    viewport = None
    color = None                # Color of rectangle (inside)
    border = None               # Color for rectangle and tick lines
    gridx = None
    gridy = None
    gridstyle = None
    tickx = None
    ticky = None
    xformat = None
    xtemplate = "{:.2f}"
    yformat = None
    ytemplate = "{:.2f}"
    xlabels = None
    ylabels = None

    def __init__(self, x1, y1, x2, y2, color):
        self.p1 = Point(x1, y1)
        self.p2 = Point(x2, y2)
        self.border = color
#        self.xformat = utils.defaultFormat
#        self.yformat = utils.defaultFormat

    def setgrid(self, nx, ny, style):
        self.gridx = nx
        self.gridy = ny
        self.gridstyle = style

    def setViewport(self, vx1, vy1, vx2, vy2):
        self.viewport = Viewport(self.p1, self.p2, Point(vx1, vy1), Point(vx2, vy2))

    def viewportOn(self, d):
        vp1 = self.viewport.vp1
        vp2 = self.viewport.vp2
        return d.setViewport(self.p1.x, self.p1.y, self.p2.x, self.p2.y, 
                             vp1.x, vp1.y, vp2.x, vp2.y)

    def viewportOff(self, d):
        self.viewport = None
        return d.cancelViewport()

    def drawGrid(self, d):
        d.setStyle(self.gridstyle)
        if self.gridx:
            w = (self.p2.x - self.p1.x) / self.gridx
            for i in range(1, self.gridx):
                px = self.p1.x + (w * i)
                d.drawLine(px, self.p1.y+1, px, self.p2.y-1, -1)

        if self.gridy:
            h = (self.p2.y - self.p1.y) / self.gridy
            for i in range(1, self.gridy):
                py = self.p1.y + (h * i)
                d.drawLine(self.p1.x+1, py, self.p2.x-1, py, -1)
        
    def drawTicks(self, d):
        d.setFont(2)
        if self.tickx:
            w = 1.0 * (self.p2.x - self.p1.x) / self.gridx
            labelfunc = self.xformat
            for i in range(0, self.gridx+1):
                px = self.p1.x + (w * i)
                d.drawLine(px, self.p2.y+1, px, self.p2.y+5, self.border)
                if self.xlabels:
                    label = self.xlabels[i]
                else:
                    vpx = self.viewport.realToVirtual(px, 0)
                    label = labelfunc(self, vpx.x, True)
                d.drawString(label, px, self.p2.y+6, 1, 2)
        if self.ticky:
            h = 1.0 * (self.p2.y - self.p1.y) / self.gridy
            # print "h = {}-{}/{} = {}".format(self.p2.y, self.p1.y, self.gridy, h)
            labelfunc = self.yformat
            for i in range(0, self.gridy+1):
                py = self.p1.y + (h * i)
                # print (i, py)
                d.drawLine(self.p1.x-5, py, self.p1.x-1, py, self.border)
                if self.ylabels:
                    label = self.ylabels[i]
                else:
                    vpy = self.viewport.realToVirtual(0, py)
                    label = labelfunc(self, vpy.y, False)
                d.drawString(label, self.p1.x-6, py, 1, 6)
            
    def draw(self, d):
        """Draw this plotPanel on drawer `d'."""
        if self.color:
            d.drawFilledRectangle(self.p1.x, self.p1.y, self.p2.x, self.p2.y, self.color)
        d.drawRectangle(self.p1.x-1, self.p1.y-1, self.p2.x+1, self.p2.y+1, self.border)
        if self.gridstyle:
            self.drawGrid(d)
        if self.tickx or self.ticky:
            self.drawTicks(d)

