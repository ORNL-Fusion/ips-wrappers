"""
plt_XY_curves.py

Contains classes and functions to plot one figure containing one or more x,y curves.

Class XY_curve(x, y, **kwargs) contains the data for one curve.  Keyword arguments allow one to
set curve label, curve color, or curve linestyle.

Class XY_Curves_Fig(curve_list, title = '', xlabel = '', ylabel = '', **kwargs) produces a one
page figure from a list of instances of XY_curve.  Keyword arguments allow one to set the
figure size.

Function plot_XY_Curves_Fig(fig) takes an instance of XY_Curves_Fig and plots it.  Currently it
is set to write a page in a pdf file.

Functions open_file_XY_Curves_Fig(filename) and close_file_XY_Curves_Fig(filename) open and close
the file into which the figures are written.  They should be called before the first call to
plot_XY_Curves_Fig and after the last call to plot_XY_Curves_Fig respectively

This script requires external modules:
    matplotlib
    numpy

change log:
 1/5/2011
 version 1.0

 1/25/2014
 Starting with matplotlib 1.3.1 a warning is issued if more than 20 figures are open at
 once.  For use with the monitor component this would happen when the number of pdf
 pages exceeds 20, which is all the time.  The warning can be eliminated by raising the
 rcParam "max_open_warning".  However I found that just closing each figure after it is
 saved in plot_XY_Curves_Fig also does it.  So that's what I did.

 11/21/2021 (DBB)
 Modified XY_Curves_Fig to accept matplotlib xlim and ylim as keyword args

 11/2/2022 (DBB)
 Modified XY_Curves_Fig to accept figsize and aspect_ratio as keyword args

"""

from matplotlib import use
# PTB use('TkAgg')
#use('MacOSX')
use('pdf')

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_mgr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import numpy as np
import math
#import string

__all__ = ['XY_Curves_Fig', 'XY_curve', 'plot_XY_Curves_Fig', 'open_file_XY_Curves_Fig',\
           'close_file_XY_Curves_Fig', 'plot_summary', 'plot_index' ]

debug = False
sci_threshold = 1.e4
figure_count = 0
page_count = 0

#_________________________________________________________________________________________________
def scale_curve_list(curve_list):

# Scale x, y values if max(abs(x)), max(abs(y)) is greater than sci_threshold or less than
# 1/sci_threshold.  Returns the scale exponents so can use scientific notation in the x,y
# axis labels

    if debug:
        print('curve_list = ', curve_list)
    fig_xmax = 0.
    fig_ymax = 0.
    for curve in curve_list:
        max_abs_x = max(list(map(abs,curve.x)))
        if max_abs_x >= fig_xmax:
            fig_xmax = max_abs_x

        max_abs_y = max(list(map(abs,curve.y)))
        if max_abs_y >= fig_ymax:
            fig_ymax = max_abs_y

    xscale = yscale = 1

    if fig_xmax >= sci_threshold:
        xscale = int(math.log10(fig_xmax))
    if fig_xmax <= 1.0/sci_threshold and fig_xmax > 0.:
        xscale = int(math.log10(fig_xmax)-1.)
    if xscale != 1:
        scale = pow(10.0, xscale)
        for curve in curve_list:
            curve.x = [x/scale for x in curve.x]

    if fig_ymax >= sci_threshold:
        yscale = int(math.log10(fig_ymax))
    if fig_ymax <= 1.0/sci_threshold and fig_ymax > 0.:
        yscale = int(math.log10(fig_ymax)-1.)
    if yscale != 1:
        scale = pow(10.0, yscale)
        for curve in curve_list:
            curve.y = [y/scale for y in curve.y]

    return (xscale, yscale)

#_________________________________________________________________________________________________

def plot_summary(summary, **kwargs):

# summary is a dictionary containing the file name and global_attributes
# global_attributes is a list of length 2 lists -> [[att name, att value],...]

    entries_per_line = 4
    max_lines_per_column = 20
    column_x = [0.05, 0.25, 0.55, 0.75]

    file_name = summary['file_name']
    global_attributes = summary['global_attributes']

    figsize = (9.0, 6.0)
    if 'figsize' in kwargs:
        figsize = kwargs['figsize']
        print('figsize = ', figsize)

    fig = plt.figure(figsize=figsize)
    ax = plt.axes([0,0,1,1])

    fig.suptitle( file_name, fontsize = 18, fontweight = 'bold')

    lines = ''
    for i in range(len(global_attributes)):
        lines = lines + global_attributes[i][0] + global_attributes[i][1] + '\n'
    plt.annotate(lines, (0.1, 0.9), xycoords = 'figure fraction', fontsize = 12,\
            fontweight = 'bold', verticalalignment = 'top')

    plt.draw()
    #plt.show()

    plt.savefig(plot_file, format = 'pdf' )

def plot_index(index, entries_per_line = 2, **kwargs):

# index is a list of entries which are themselves length 2 lists -> [plot #, title]

    max_lines_per_column = 25
    entries_per_page = entries_per_line * max_lines_per_column
    column_x = [0.05, 0.5]
    y_top_row = 0.9

    figsize = (9.0, 6.0)
    if 'figsize' in kwargs:
        figsize = kwargs['figsize']
        print('figsize = ', figsize)

    fig = plt.figure(figsize=figsize)
    ax = plt.axes([0,0,1,1])

    fig.suptitle( 'Index', fontsize = 18, fontweight = 'bold')

    lines = ''
    line_count = 1
    column_count = 0
    page_line_count = 0
    for i in range(len(index)):
        lines = lines + str.ljust(str(index[i][0]), 3, ' ') + ' ' +\
                index[i][1] + '\n'
        line_count = line_count + 1
        page_line_count = page_line_count + 1

        if i == len(index) - 1: # write column if this is last entry
            plt.annotate(lines, (column_x[column_count], y_top_row), xycoords = 'figure fraction',\
                          verticalalignment = 'top')
            plt.draw()
            plt.savefig(plot_file, format = 'pdf' )
            break

        if line_count > max_lines_per_column:
            plt.annotate(lines, (column_x[column_count], y_top_row), xycoords = 'figure fraction',\
                          verticalalignment = 'top')
            lines = ''
            line_count = 1
            column_count = column_count + 1

        if page_line_count >= entries_per_page :
            plt.draw()
            plt.savefig(plot_file, format = 'pdf' )

            line_count = 1
            column_count = 0
            page_line_count = 0
            fig = plt.figure(figsize=figsize)
            ax = plt.axes([0,0,1,1])
            fig.suptitle( 'Index', fontsize = 18, fontweight = 'bold')

#_________________________________________________________________________________________________

class XY_Curves_Fig:
    def __init__(self, curve_list, title = '', xlabel = '', ylabel = '', **kwargs):

        if debug:
            print('XYCurves_Fig: curve_list = ', curve_list)
# Check the arguments
        if isinstance(curve_list, XY_curve): # It's legal to pass in a single XY_curve instance
            self.curve_list = [curve_list]
        elif type(curve_list) == list:
            for curve in curve_list: # check that everything is an instance of XY_curve
                if not  isinstance(curve, XY_curve):
                    print('curve type not = instance')
                    print('plt_XY_Curves: arg curve_list must be a list of instances of XY_curve')
                    raise Exception('plt_XY_Curves: arg curve_list must be a list of instances of\
                    XY_curve')
                self.curve_list = curve_list
        else:
            print('curve_list type = ', type(curve_list))
            print('plt_XY_Curves: arg curve_list must be a list of instances of XY_curve')
            raise Exception('plt_XY_Curves: arg curve_list must be a list of instances of XY_curve')

        self.title = str(title)
        self.xlabel = str(xlabel)
        self.ylabel = str(ylabel)
        self.kwargs = kwargs
# Arguments are OK

# Create figure
        global figure_count
        figure_count = figure_count + 1
        self.fig_number = figure_count

        self.figsize = (9.0, 6.0)
        if 'figsize' in kwargs:
            self.figsize = kwargs['figsize']
            print('figsize = ', self.figsize)

        fig = plt.figure(figsize=self.figsize)

        if 'aspect_ratio' in kwargs:
            self.aspect_ratio = kwargs['aspect_ratio']
            print('aspect_ratio = ', self.aspect_ratio)
            plt.axes((0.1, 0.1, 0.75, 0.75), aspect = self.aspect_ratio)
        else:
           plt.axes((0.1, 0.1, 0.75, 0.75))

#        plt.axes((0.15, 0.1, 0.3, 0.75))

        if 'xlim' in kwargs:
            self.xlim = kwargs['xlim']
            print('xlim = ', self.xlim)
            plt.xlim(self.xlim)

        if 'ylim' in kwargs:
            self.ylim = kwargs['ylim']
            print('ylim = ', self.ylim)
            plt.ylim(self.ylim)

        # Note to DBB: Annotating figure number below is confusing because it doesn't coincide with
        # page number.  Change this to add page number instead.
        #str_fig_number = str(figure_count)
        #plt.annotate(str_fig_number, (0.9, 0.05), xycoords = 'figure fraction')

        # Check if the maximum of data in the x or y of any curve in the curve_list is outside
        # the threshold to use scientific notation in the axis labels.  If so divide the values
        # by the appropriate scale factor and change the x and/or y axis labels to reflect this.

        scaleX, scaleY = scale_curve_list(self.curve_list)
        if scaleX != 1:
            power_string = '  ' + r'$' + r'\times' + '10^{' + str(scaleX) + r'})$'
            self.xlabel = xlabel + power_string
        if scaleY != 1:
            power_string = '  (' + r'$' + r'\times' + '10^{' + str(scaleY) + r'})$'
            self.ylabel = ylabel + power_string

        if debug:
            print('scaleX = ', scaleX)
            print('scaleX = ', scaleX)
            print('self.xlabel = ', self.xlabel)
            print('self.ylabel = ', self.ylabel)

        return
#_________________________________________________________________________________________________

class XY_curve:

    count = 0

    def __init__(self, x, y, **kwargs):

        if (type(x) == list or type(x) == np.ndarray or type(x) == np.ma.core.MaskedArray) and\
                (type(y) == list or type(y) == np.ndarray or type(y) == np.ma.core.MaskedArray):
            self.x = x
            self.y = y
            self.kwargs = kwargs
        else:
            print('type(x) = ', type(x))
            print('type(y) = ', type(y))
            print('XY_curve: required args x and y must be lists or numpy.ndarrays')
            raise Exception('XY_curve: required args x and y must be lists or numpy.ndarrays')
        if len(x) != len(y):
            print('XY_curve: required args x and y must be the same length')
            raise Exception('XY_curve: required args x and y must be the same length')

        if 'label' in kwargs:
            self.label = kwargs['label']

        if 'color' in kwargs:
            self.color = kwargs['color']

        if 'linestyle' in kwargs:
            self.linestyle = kwargs['linestyle']

        self.__class__.count = self.__class__.count +1

    def __getattr__(self, attrname):

        if attrname in ['label', 'color','linestyle'] :
            return None
        else:
            raise AttributeError(attrname)

    def setLabel(self, label):
        self.label = label
        self.kwargs['label'] = label

    def setColor(self, color):
        self.color = color
        self.kwargs['color'] = color

    def setLinestyle(self, linestyle):
        self.linestyle = linestyle
        self.kwargs['linestyle'] = linestyle

#_________________________________________________________________________________________________

def plot_XY_Curves_Fig(fig):

    if not isinstance(fig, XY_Curves_Fig):
        print('plot_XY_Curves_FIG: arg fig must be an instances of XY_Curves_Fig')
        raise Exception('plot_XY_Curves_FIG: arg fig must be an instances of XY_Curves_Fig')

    fig_has_legend = False
    for curve in fig.curve_list:
        if debug:
            print('plot_XY_Curves_Fig: curve = ', curve)

        plt.plot(curve.x, curve.y, **curve.kwargs)
        if 'label' in curve.kwargs:
            fig_has_legend = True

    plt.title(fig.title)
    plt.xlabel(fig.xlabel)
    plt.ylabel(fig.ylabel)

    # print 'dealing with legend'

    # deal with legend format only if curves on this fig have labels
    if fig_has_legend:
        prop = font_mgr.FontProperties(size = 'small')
        leg = plt.legend(shadow=True,loc = 'upper left', bbox_to_anchor = (1.01, 1.0), \
              borderaxespad=0., prop = prop)

    plt.draw()
    #plt.show()

    plt.savefig(plot_file, format = 'pdf' )
    plt.close()

#_________________________________________________________________________________________________

def open_file_XY_Curves_Fig(filename):

    global plot_file

    if (type(filename) != str):
        print('open_file_XY_Curves_Fig: filename must be of type str')
        raise Exception('open_file_XY_Curves_Fig: filename must be of type str')

    plot_file = PdfPages(filename)

#_________________________________________________________________________________________________

def close_file_XY_Curves_Fig():

    global plot_file

    plot_file.close()

#_________________________________________________________________________________________________


if __name__ == '__main__':

    import math

    open_file_XY_Curves_Fig('plot_output.pdf')

    x = [0.1*i for i in range(11)]

    curve_list = []
    for i in range(1,6):
        lbl = 'curve(' + str(i) + ')'
        y = [-2.0e4 + 1.e5*i*r for r in x]
        new_curve = XY_curve(x, y, label = lbl)
        curve_list.append(new_curve)

    curve_list[0].setLinestyle('dashed')

    title = 'Some real good stuff'
    xlabel = 'x(cm)'
    ylabel = 'volts'
    plot1 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)

    plot_XY_Curves_Fig(plot1)

    title = 'Even better stuff'
    plot2 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
    plot_XY_Curves_Fig(plot2)

    title = 'Try to set xlim'
    plot3 = XY_Curves_Fig(curve_list, title, xlabel, ylabel, xlim = [0.2, 0.8])
#    plot3 = XY_Curves_Fig(curve_list, title, xlabel, ylabel)
    plot_XY_Curves_Fig(plot3)

    title = 'Try to set aspect_ratio'
    plot4 = XY_Curves_Fig(curve_list, title, xlabel, ylabel, aspect_ratio = 'equal')
    plot_XY_Curves_Fig(plot4)

    x = []
    y = []
    vx = []
    vy = []
    vscale = []
    for i in range(16):
        t = 0.1*3.1415926*i
        x.append(math.cos(t))
        y.append(math.sin(t))
        vx.append(math.cos(t))
        vy.append(math.sin(t))
        vscale.append(0.1*t)

    title = 'Parametric Plot - different fig size'
    xlabel = 'x(cm)'
    ylabel = 'y(cm)'
    curve_list = [XY_curve(x, y)]
    plot_parametric2 = XY_Curves_Fig(curve_list, title, xlabel, ylabel, aspect_ratio = 'equal')
    plot_XY_Curves_Fig(plot_parametric2)


    title = 'Parametric Plot'
    xlabel = 'x(cm)'
    ylabel = 'y(cm)'
    curve_list = [XY_curve(x, y)]
    plot_parametric = XY_Curves_Fig(curve_list, title, xlabel, ylabel,figsize = (8., 8.),\
                      aspect_ratio = 'equal')

    for i in range(16):
        plt.arrow(x[i], y[i], vscale[i]*vx[i], vscale[i]*vy[i], shape='full', head_width = 0.02)
    plot_XY_Curves_Fig(plot_parametric)


    title = 'Parametric Plot - different fig size'
    xlabel = 'x(cm)'
    ylabel = 'y(cm)'
    curve_list = [XY_curve(x, y)]
    plot_parametric2 = XY_Curves_Fig(curve_list, title, xlabel, ylabel, aspect_ratio = 'equal')

    plot_XY_Curves_Fig(plot_parametric2)

    global_attributes = [['Global_label = ', 'Global_label'], ['RunID = ', 'RunID'],\
                      ['tokamak_id = ', 'tokamak_id'], ['shot_number = ', str(2)] ]

#_________________________________________________________________________________________

    print('adding index')
    index = [ [1,'fig 1'], [2,'fig 2'], [3,'fig 3'], [4,'fig 4'] ]
    summary ={'global_attributes': global_attributes, 'index': index, 'file_name': 'no file'}

    plot_summary(summary)

    close_file_XY_Curves_Fig()
