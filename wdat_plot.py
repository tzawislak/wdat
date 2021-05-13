import numpy as np
import struct
import cmath, math

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# this class unpacks .wtxt files
# as a parameter in constructor data prefix should be passed
class Wdat:
    types = {'real' : 1, 'complex' : 2, 'vector' : 3}
    N = [0,0,0]  # dimensions
    D = [0,0,0]  # lattice spacing
    consts = {} # dictionary of constants
    var = {}   # dictionary of variables
    datadim = 3
    cycles =1   # number of cycles (measurements)
    t0=-1   # time value for the first cycle
    dt=1
    prefix=""
    ######################################################################
    #####   Prints all class elements
    ######################################################################
    def __str__(self):
        text="prefix:\t{0}\nN={1}\nD={2}\ndatadim:\t{3}\n".format(self.prefix, self.N, self.D, self.datadim)
        text+="cycles:\t{0}\nt0:\t{1}\ndt:\t{2}\n### constants ###\n".format(self.cycles, self.t0, self.dt)
        for key, value in self.consts.items():  text += "{0}\t{1}\n".format(key, value)
        text += "### variables ###\n"
        for key, value in self.var.items():  text += "{0}\t\t{1}\n".format(key, value)
        return text
    ######################################################################
    ##### Reads .wtxt file
    ##### @param filepath - path to .wtxt file
    ######################################################################
    def __init__(self, filepath : str):
        print( "Reading: ", filepath+".wtxt")
        self.prefix = filepath

        #
        # READING .wtxt file from file
        #
        data = open( filepath + ".wtxt", "r")
        # split into lines and remove comments
        lines = data.read().splitlines()
        for i in range(0, len(lines)):
            lines[i] = lines[i].split("#", 1)[0]
            lines[i] = lines[i].split(" ")
            # remove empty elements
            while "" in lines[i] : lines[i].remove("")
        while [] in lines : lines.remove([])

        metadata = []
        for i in lines:
            for j in i:
                metadata.append(str(j))

        #
        # Unpack variable types
        #
        self.readN(metadata)
        self.readD(metadata)
        self.readConsts(metadata)
        self.readVars(metadata)
        self.readOthers(metadata)
    # done

    ######################################################################
    ##### INITIALIZERS
    ##### 5 functions below fill in class elements with content from .wtxt
    ######################################################################
    def readN(self, lines : list):
        self.N[0] = int(lines[ lines.index('NX') + 1 ])
        self.N[1] = int(lines[ lines.index('NY') + 1 ])
        self.N[2] = int(lines[ lines.index('NZ') + 1 ])
    def readD(self, lines : list):
        self.D[0] = int(lines[ lines.index('DX') + 1 ])
        self.D[1] = int(lines[ lines.index('DY') + 1 ])
        self.D[2] = int(lines[ lines.index('DZ') + 1 ])
    def readConsts(self, lines : list):
        for i in range(len(lines)):
            if lines[i] == "const":
                self.consts[ lines[i+1] ] = float( lines[i+2] )
    def readVars(self, lines : list):
        for i in range(len(lines)):
            if lines[i] == "var":
                self.var[ lines[i+1] ] = str( lines[i+2] )
    def readOthers(self, lines : list):
        self.datadim = int(lines[ lines.index('datadim') + 1 ])
        self.cycles = int(lines[ lines.index('cycles') + 1 ])
        self.t0 = float(lines[ lines.index('t0') + 1 ])
        self.dt = float(lines[ lines.index('dt') + 1 ])
    ######################################################################
    ##### Reads data from prefix_variable.wdat file
    ##### @param variable - variable to be read
    ##### @param cycle - number of cycle, if not passed
    #####                last cycle will be read
    ##### WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
    ##### For the time being 2D data can be plotted
    ##### returns:
    #####   2D: 3D numpy array
    #####       dimension of outermost list is equal self.datadim
    ######################################################################
    def read( self, variable : str, cycle=-1 ):
        if not self.var[variable] : return
        # SET BLOCKLENGHT
        typelen = self.types[self.var[variable]]
        blocklength = self.N[0]
        if self.datadim == 2 :
            blocklength *= self.N[1]
        if self.datadim == 3 :
            blocklength *= self.N[1]*self.N[2]
        # SET BLOCKSIZE
        blocksize = blocklength * 8 * typelen # number of bytes

        if cycle==-1 : cycle = self.cycles-1

        # open prefix_variable.wdat file
        datafile = open( self.prefix + "_" + variable + ".wdat", "rb")
        datafile.seek( cycle*blocksize ) # moves to the part which is about to be read
        data = datafile.read( blocksize )
        #
        # 2D data support only!!!
        #

        if self.datadim == 2 :
            ar = np.zeros( (typelen, self.N[0], self.N[1]) )
            sliced = [data[i:i+typelen*8] for i in range(0, len(data), typelen*8)]
            for x in range(0, self.N[0]):
                for y in range(0, self.N[1]):
                    payload = sliced[ x*self.N[1] + y  ]
                    for t in range(0, typelen):
                        [number] = struct.unpack('d', [payload[i:i+8] for i in range(0, len(payload), 8)][t])
                        ar[t][y][x] = number
            return ar
    ######################################################################
    ##### Generally, complex variables are written in Re Im format.
    ##### This function converts data[][][] to Abs Arg format
    ##### @param data - 3D numpy array read for
    #####               complex data
    ##### returns 3D numpy array
    #####         [0] - absolute value
    #####         [1] - argument
    ######################################################################
    def ReIm2AbsArg(self, data):
        if self.datadim == 2 :
            ar = np.zeros( (2, self.N[0], self.N[1]) )
            for x in range(0, self.N[0]):
                for y in range(0, self.N[1]):
                    c = complex(data[0][y][x], data[1][y][x])
                    ar[0][y][x] = abs(c)
                    ar[1][y][x] = cmath.phase(c)
            return ar

    ######################################################################
    ##### Draws 2D data
    ##### @param data - [][] list (double nested list) of data
    ##### @param figpath - name of figure to be saved,
    #####                  if not specified, the plot will not be saved
    ######################################################################
    def Draw2D(self, data, figpath=""):
        ##### X Y setup #####
        xlist = np.linspace(-1*self.N[0]/2. * self.D[0], self.N[0]/2. * self.D[0], self.N[0])
        ylist = np.linspace(-1*self.N[1]/2. * self.D[1], self.N[1]/2. * self.D[1], self.N[1])
        X, Y = np.meshgrid(xlist, ylist)


        fig,ax=plt.subplots(1,1, figsize=(12, 10 ))
        #ax = plt.axes(projection ='3d')

        ##### COLOR SCALE CONFIG #####
        cp = ax.contourf(Y, X, data, levels=100, cmap=plt.get_cmap("jet") )
        # assign contour and tics in legend to fig
        print( "MAX:\t", np.max(data), "\nMIN:\t", np.min(data) )
        min_scale = np.min(data)-1e-6
        max_scale = np.max(data)-1e-6
        N_of_labels = 5
        prec = 2
        # Add a colorbar to a plot
        fig.colorbar(cp, ticks=np.round(np.arange(min_scale, max_scale+(max_scale-min_scale)/(N_of_labels-1), (max_scale-min_scale)/(N_of_labels-1)), prec ))

        # change scale's ticks size
        ax.figure.axes[1].tick_params(labelsize=22)
        # set plot title and padding,
        #title = r'$\Delta_{\varphi}\quad [\pi]$'
        title = r'$|\Delta|\quad [\epsilon_{F}]$'
        ax.set_title(title, fontsize=28, pad=30)
        plt.setp(ax.get_yticklabels(), fontsize=24)
        plt.setp(ax.get_xticklabels(), fontsize=24)

        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_facecolor('w')
        plt.show()
        fig.patch.set_facecolor('xkcd:white')
        if figpath != "": fig.savefig( figpath )
        fig.clf()
