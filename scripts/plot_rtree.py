import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import collections as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys, os

class BBox:
    def __init__(self, lowleft, upright, level):

        (x_low, y_low) = lowleft
        (x_top, y_top) = upright

        dx = x_top - x_low
        dy = y_top - y_low

        self.coords = [ \
                (x_low,    y_low), \
                (x_low+dx, y_low), \
                (x_low+dx, y_low+dy), \
                (x_low,    y_low+dy) ]

        self.dx      = dx
        self.dy      = dy
        self.center  = (x_low+0.5*dx, y_low+0.5*dy)
        self.lowleft = lowleft
        self.upright = upright
        self.level   = level



class RTree:

    def __init__(self, file_name):
        ''' Initialize the RTree from a text file
        '''
        with open(file_name, 'r') as f:
            lines = f.readlines()

        self.__io_clear_comments( lines )
        self.__io_read_tree( lines )


    def __io_clear_comments(self, lines):
        ''' Clears all comments from the input string
        '''
        for i in range(len(lines)-1, -1, -1):
            line = lines[i].replace(' ','')
            if line[0] == '#':
                lines.pop( i )

    def __io_read_tree(self, lines):
        ''' Read the entire tree data from the input string
        '''
        bboxes = []
        i_line = 0
        max_level = 0

        while True:
            line = lines[i_line].split(",")

            (n, level, index, parent_id) = [int(l) for l in line]

            if level > max_level:
                max_level = level

            for i in range(1,n+1):
                line = lines[i_line+i].split(",")
                (x_low, y_low, x_top, y_top) = [float(l) for l in line]
                bboxes.append( BBox((x_low, y_low), (x_top, y_top), level) )

            i_line += n + 1

            if i_line >= len(lines)-1:
                break

        self.bboxes = bboxes
        self.height = max_level + 1


    def plot_bboxes(self, ax):
        ''' Plot all bounding boxes of the R-Tree
        '''
        # Gather bbox coordinate array
        levels, coords = [], []
        for b in self.bboxes:
            coords.append( b.coords )
            levels.append( b.level )
        coords = np.array(coords)
        levels = np.array(levels).astype(int)

        # Create bbox colors
        colors = plt.cm.Set1( np.linspace(0.0,1.0,self.height) )
        color_map = np.zeros( levels.shape ).astype(int)
        for i, c in enumerate(np.unique(levels)):
            color_map[levels == c] = i
        bbox_colors = colors[color_map]

        alphas = 1. - levels / float(self.height-1.)
        bbox_colors[:,-1] = 0.3 + 0.5 * alphas

        linewidths = 0.3 + 0.8 * (1. - alphas)

        bbox_collection = mc.PolyCollection(coords, edgecolor=['k'],
                                            facecolors=bbox_colors,
                                            lw=linewidths)

        ax.add_collection( bbox_collection )

        xlim = (np.min(coords[:,:,0]), np.max(coords[:,:,0]))
        ylim = (np.min(coords[:,:,1]), np.max(coords[:,:,1]))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)




def main(argv):
    if len(argv) < 2:
        print("plot_rtree.py <rtree.txt> <flags>")
        sys.exit(1)

    file_name = sys.argv[1]

    rtree = RTree( file_name )


    # Create the plot
    fig, ax = plt.subplots(1,1,dpi=200)

    rtree.plot_bboxes( ax )

    plt.show()



if __name__ == '__main__': main(sys.argv)
