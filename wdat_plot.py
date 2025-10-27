import numpy as np
import struct
import cmath, math



class Wdat:
    """
    A class that parses the .wdat files and prepares them for plotting.
    Pass the data prefix as the constructor argument
    """



    types = {'real' : 1, 'complex' : 2, 'vector' : 3}
    N = [0.,0.,0.]  # dimensions
    D = [0.,0.,0.]  # lattice spacing
    L = [0.,0.,0.]
    V = 1
    bx = [0.,0.,0.]
    by = [0.,0.,0.]
    bz = [0.,0.,0.]
    X,Y,Z=0,0,0
    consts = {} # dictionary of constants
    var = {}   # dictionary of variables
    datadim = 3
    cycles =1   # number of cycles (measurements)
    t0=-1   # time value for the first cycle
    dt=1
    prefix=""
    lenscale=r' [a_{\text{ho}}]$'

    ######################################################################
    def __str__(self):
        text="prefix:\t{0}\nN={1}\nD={2}\ndatadim:\t{3}\n".format(self.prefix, self.N, self.D, self.datadim)
        text+="cycles:\t{0}\nt0:\t{1}\ndt:\t{2}\n### constants ###\n".format(self.cycles, self.t0, self.dt)
        for key, value in self.consts.items():  
            text += "{0}\t\t{1}\t\t{2}\n".format(key, value[0], value[1])
        text += "### variables ###\n"
        for key, value in self.var.items():  
            text += "{0}\t\t{1}\t\t{2}\n".format(key, value[0], value[1])
        return text

    ######################################################################
    def __init__(self, filepath : str, ahoscale=False, _mute=False):
        self.types = {'real' : 1, 'complex' : 2, 'vector' : 3}
        self.N = [0.,0.,0.]  # dimensions
        self.D = [0.,0.,0.]  # lattice spacing
        self.L = [0.,0.,0.]  # lengths
        self.V = 1
        self.bx = [0.,0.,0.]
        self.by = [0.,0.,0.]
        self.bz = [0.,0.,0.]
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.consts = {} # dictionary of constants
        self.var = {}   # dictionary of variables
        self.datadim = 3
        self.cycles =1   # number of cycles (measurements)
        self.t0=-1   # time value for the first cycle
        self.dt=1
        self.prefix=""
        self.lenscale=r' [a_{\text{ho}}]$'
        if not _mute:
            print( "Reading: ", filepath)
        self.prefix = filepath[:-10]
        #
        # READING .wtxt file from file
        #
        data = open( filepath, "r")
        # split into lines and remove comments
        lines = data.read().splitlines()
        data.close()
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
        try:
            self.readn(metadata)
        except:
            self.readN(metadata)
        try:
            self.readd(metadata)
        except:
            self.readD(metadata)
        self.readConsts(metadata)
        self.readVars(metadata)
        self.readOthers(metadata)
        
        
        if not ahoscale:
            self.D[0] = self.D[0] * self.consts['aho'][0]
            self.D[1] = self.D[1] * self.consts['aho'][0]
            self.D[2] = self.D[2] * self.consts['aho'][0]
            self.lenscale  = r' [\mu m]$'
        ddx = 0# -self.D[0]/2.0
        ddy = 0#-self.D[1]/2.0    
        ddz = 0#-self.D[2]/2.0
        # origin not shifted
        #self.bx=[ -0.5*(self.N[0]-1)*self.D[0] + ddx, 0.5*(self.N[0]-1)*self.D[0] + ddx, self.D[0] -1e-6 ]
        #self.by=[ -0.5*(self.N[1]-1)*self.D[1] + ddy, 0.5*(self.N[1]-1)*self.D[1] + ddy, self.D[1] -1e-6 ]
        #self.bz=[ -0.5*(self.N[2]-1)*self.D[2] + ddz, 0.5*(self.N[2]-1)*self.D[2] + ddz, self.D[2] -1e-6 ]
        self.bx=[ -0.5*self.N[0]*self.D[0], (0.5*self.N[0]-1)*self.D[0], self.D[0] -1e-6 ]
        self.by=[ -0.5*self.N[1]*self.D[1], (0.5*self.N[1]-1)*self.D[1], self.D[1] -1e-6 ]
        self.bz=[ -0.5*self.N[2]*self.D[2], (0.5*self.N[2]-1)*self.D[2], self.D[2] -1e-6 ]
        self.X = np.arange( -0.5*self.N[0]*self.D[0] + 0.5*self.D[0], (0.5*self.N[0])*self.D[0] + 0.5*self.D[0], self.D[0]  )
        self.Y = np.arange( -0.5*self.N[1]*self.D[1] + 0.5*self.D[1], (0.5*self.N[1])*self.D[1] + 0.5*self.D[1], self.D[1]  )
        self.Z = np.arange( -0.5*self.N[2]*self.D[2] + 0.5*self.D[2], (0.5*self.N[2])*self.D[2] + 0.5*self.D[2], self.D[2]  )

        self.L[0] = self.N[0]*self.D[0]
        self.L[1] = self.N[1]*self.D[1]
        self.L[2] = self.N[2]*self.D[2]
        self.V = self.L[0]*self.L[1]*self.L[2]
    # done

    ################
    # INITIALIZERS #
    ################
    ######################################################################
    def readN(self, lines : list):
        self.N[0] = int(lines[ lines.index('NX') + 1 ])
        self.N[1] = int(lines[ lines.index('NY') + 1 ])
        self.N[2] = int(lines[ lines.index('NZ') + 1 ])


    def readD(self, lines : list):
        self.D[0] = int(lines[ lines.index('DX') + 1 ])
        self.D[1] = int(lines[ lines.index('DY') + 1 ])
        self.D[2] = int(lines[ lines.index('DZ') + 1 ])

    def readn(self, lines : list):
        self.N[0] = int(lines[ lines.index('nx') + 1 ])
        self.N[1] = int(lines[ lines.index('ny') + 1 ])
        self.N[2] = int(lines[ lines.index('nz') + 1 ])


    def readd(self, lines : list):
        self.D[0] = float(lines[ lines.index('dx') + 1 ])
        self.D[1] = float(lines[ lines.index('dy') + 1 ])
        self.D[2] = float(lines[ lines.index('dz') + 1 ])


    def readConsts(self, lines : list):
        for i in range(len(lines)):
            if lines[i] == "const":
                self.consts[ lines[i+1] ] = (float( lines[i+2] ), lines[i+3])

    def readVars(self, lines : list):
        for i in range(len(lines)):
            if lines[i] == "var":
                self.var[ lines[i+1] ] = (str( lines[i+2] ), lines[i+3])

    def readOthers(self, lines : list):
        self.datadim = int(lines[ lines.index('datadim') + 1 ])
        self.cycles = int(lines[ lines.index('cycles') + 1 ])
        self.t0 = float(lines[ lines.index('t0') + 1 ])
        self.dt = float(lines[ lines.index('dt') + 1 ])


    ######################################################################
    def read( self, variable : str, cycle=-1 ):
        """ Reads a binary file and formats it accoring to the data type

        """
        
        if ('_final' in variable):
            typelen = self.types[self.var[variable.replace('_final', '')][0]]
        elif ('init' in variable):
            typelen = self.types[self.var[variable.replace('_init', '')][0]]
        else:
            typelen = self.types[self.var[variable][0]]
            
        blocklength = self.N[0]
        if self.datadim == 2 :
            blocklength *= self.N[1]
        if self.datadim == 3 :
            blocklength *= self.N[1]*self.N[2]
        # SET BLOCKSIZE
        blocksize = blocklength * 8 # number of bytes

        if cycle==-1 : cycle = self.cycles-1


        # open prefix_variable.wdat file
        datafile = open( self.prefix + "_" + variable + ".wdat", "rb")
        datafile.seek( cycle*blocksize ) # moves to the part which is about to be read
        data = datafile.read( blocksize )
        datafile.close()
        
        if self.datadim == 1 :
            ar = np.fromfile( self.prefix + "_" + variable + ".wdat",
                             count = typelen*self.N[0]*self.N[1]*self.N[2],
                             offset= typelen*cycle*blocksize          )
            ar = ar.reshape( (self.N[2], self.N[1], self.N[0], typelen) )
            
            if typelen==2:
                return self.ReIm3AbsArg(ar)[:,:,0,0]    
            if typelen==1:
                return np.transpose(ar, (3,2,1,0))[:,:,0,0]


        if self.datadim == 2 :
            ar = np.zeros( (typelen, self.N[1], self.N[0]) )
            sliced = [data[i:i+typelen*8] for i in range(0, len(data), typelen*8)]
            for x in range(0, self.N[0]):
                for y in range(0, self.N[1]):
                    payload = sliced[ x*self.N[1] + y  ]
                    for t in range(0, typelen):
                        [number] = struct.unpack('d', [payload[i:i+8] for i in range(0, len(payload), 8)][t])
                        ar[t][y][x] = number
            return ar

        if self.datadim == 3 :
            ar = np.fromfile( self.prefix + "_" + variable + ".wdat",
                             count = typelen*self.N[0]*self.N[1]*self.N[2],
                             offset= typelen*cycle*blocksize          )
            
            ar = ar.reshape( (self.N[2], self.N[1], self.N[0], typelen) )
            
            
            
            if typelen==2:
                return self.ReIm3AbsArg(ar)
            
            if typelen==1:
                return np.transpose(ar, (3,2,1,0))
            
        
    ######################################################################
    def ReIm2AbsArg(self, data):
        """ converts 2D arrays of complex numbers saved in [Re,Im] format to [Abs,Arg] """
        ar = np.zeros( (2, self.N[1], self.N[0]) )
        for x in range(0, self.N[0]):
            for y in range(0, self.N[1]):
                c = complex(data[0][y][x], data[1][y][x])
                ar[0][y][x] = abs(c)
                ar[1][y][x] = cmath.phase(c)
        return ar

    ######################################################################
    def ReIm3AbsArg(self, data):
        """ converts 3D arrays of complex numbers saved in [Re,Im] format to [Abs,Arg] """
        cmplx = 1j*data[:,:,:,1]; cmplx += data[:,:,:,0]
        ar = np.array([ np.abs(cmplx[:,:,:]), np.angle(cmplx[:,:,:]) ])
        ar =  np.transpose(ar, (0,3,2,1))
        return ar
    ######## Pass 2D array only ##########################################
    def CalcXYGrad(self, data):
        """ returns a 2D gradient of data """
        gx = np.arccos(np.cos( (data[1:,1:] - data[1:,:-1])*np.pi ))
        gy = np.arccos(np.cos( (data[1:,1:] - data[:-1,1:])*np.pi ))
        
        return ((gx))**2 + ((gx))**2
    ######## Pass 2D array only ##########################################
    def ReturnXYGrad(self, data):
        """ returns a 2D gradient of data """
        gx = np.arcsin(np.sin( (data[1:,1:] - data[1:,:-1])*np.pi ))
        gy = np.arcsin(np.sin( (data[1:,1:] - data[:-1,1:])*np.pi ))
        
        return gx, gy
        

    ######################################################################
    def Integrate2D( self, data, _type="X" ):
        # pass 'iX' if you want to integrate over the first index but make sure it's the 
        # x-th direction
        if _type[0]=="i":
            data = data.transpose()
            _type = _type[1:]
            
        if _type=="X":
            NX, dx = self.N[0], self.D[0]
        elif _type=="Y":
            NX, dx = self.N[1], self.D[1]
        elif _type=="Z":
            NX, dx = self.N[2], self.D[2]
        
        ar = np.zeros(len(data))
        for n1 in range(0, len(data)):
            nsum=0
            for x in range(0, NX):
                if _type=="X":
                    nsum += data[n1][x] * dx
                elif _type=="Y":
                    nsum += data[n1][x] * dx
                elif _type=="Z":
                    nsum += data[n1][x] * dx
            ar[n1] = nsum
        return ar
        
    
    ######################################################################
    def Integrate3D( self, data, _type="XY" ):
        """Integrates a 3D array over one direction

        Args:
            data  - a 3D data array
            _type - specifies the integration direction:
                        XY - over Z
                        XZ - over Y
                        YZ - over X
        Returns:
            A 2D array integrated over the 3rd direction

        Prints the sum of all elements
        """
        if _type=="XY":
            N1, N2, N3 = self.N[0], self.N[1], self.N[2]
        elif _type=="XZ":
            N1, N2, N3 = self.N[0], self.N[2], self.N[1]
        elif _type=="YZ":
            N1, N2, N3 = self.N[1], self.N[2], self.N[0]
            
        ar = np.zeros( (N2, N1) )
        integral=0
        for n1 in range(0, N1):
            for n2 in range(0, N2):
                xysum=0
                for n3 in range(0, N3):
                    if _type=="XY":
                        xysum += data[n1][n2][n3] * self.D[2]
                    elif _type=="XZ":
                        xysum += data[n1][n3][n2] * self.D[1]
                    elif _type=="YZ":
                        xysum += data[n3][n1][n2] * self.D[0]
                    
                ar[n2][n1] = xysum
                if _type=="XY":
                    integral += xysum * self.D[0] * self.D[1]
                elif _type=="XZ":
                    integral += xysum * self.D[0] * self.D[2] 
                elif _type=="YZ":
                    integral += xysum * self.D[2] * self.D[1]
        #check:
        print("Data integral: ", integral)# * (self.consts['aho'][0])**3 )
        return ar 
    
    ######################################################################
    def NumericalVolume(self):
        # Returns the volume of the simulated space
        return self.N[0]*self.D[0] * self.N[1]*self.D[1] *self.N[2]*self.D[2]
    
    ######################################################################
    def Cut3D(self, values, _x, _y, _z):
        """ Removes marginal array rows from the 3D array 
        
        Args:
            values      - the 3D array to trim
            _x, _y, _z  - number of rows to remove from the beginning
                          and the end of the X,Y,Z row respectively

        Returns:
            values[-_x:_x, -_y:_y, -_z:_z]

        """
        nnx = len(np.arange(-1*_x, _x, self.bx[2]))
        nny = len(np.arange(-1*_y, _y, self.by[2]))
        nnz = len(np.arange(-1*_z, _z, self.bz[2]))
        if (nnx % 2) == 1:
            nnx += 1
        if (nny % 2) == 1:
            nny += 1
        if (nnz % 2) == 1:
            nnz += 1
        print(nnx, nny, nnz)
        ix = [ int((self.N[0]-nnx)/2), int((self.N[0]+nnx)/2)  ]
        iy = [ int((self.N[1]-nny)/2), int((self.N[1]+nny)/2)  ]
        iz = [ int((self.N[2]-nnz)/2), int((self.N[2]+nnz)/2)  ]

        return values[ ix[0]:ix[1], iy[0]:iy[1], iz[0]:iz[1]  ]
        
    

   
    #######################################################################
    def Section(self, data, _type="XY"):
        # returns a 2D central section from a 3D array 
        if _type=="XY":
            return data[:,:,int(self.N[2]/2.)].transpose()
        elif _type=="XZ":
            return data[:,int(self.N[1]/2.),:].transpose()
        elif _type=="YZ":
            return data[int(self.N[0]/2.),:,:].transpose()
        
        


    ######################################################################
    """ The following ReturnSomething() functions:
    
    Args:
        _dim - dimensions of the data to apply the units to

    Return:
        c - the scale units in which the code operates
        l - string with unit
        m - Colormap for drawing
    
    Note:
        If variable type is complex, then lists are returned for amplitude and phase 
    """
    def ReturnDensityScale(self, _dim):
        c = [self.consts['aho'][0]**(-3)]
        if _dim==1:
            l = [r"$\large n\ [\mu m^{-1}]$"]
        elif _dim==2:
            l = [r"$\large n(x,y)\ [\mu m^{-2}]$"]
        else:
            l = [r"$\large n(\textbf{r})\ [\mu m^{-3}]$"]
        m = ["Bluyl"]
        return c, l, m
    
    def ReturnPsiScale0(self, _dim):
        c = [(self.consts['aho'][0])**(-3), 1/np.pi]
        if _dim==1:
            l = [r"$\large |\Psi(\textbf{r})|^{2}\ [\mu m^{-1}]$", r"$\large \varphi(\textbf{r})\ [\pi]$"]
        elif _dim==2:
            l = [r"$\large |\Psi(\textbf{x,y})|^{2}\ [\mu m^{-2}]$", r"$\large \varphi(\textbf{r})\ [\pi]$"]
        else:
            l = [r"$\large |\Psi(\textbf{r})|^{2}\ [\mu m^{-3}]$", r"$\large \varphi(\textbf{r})\ [\pi]$"]
            
        m = ["speed", "IceFire"]
        return c, l, m
    
    def ReturnPsiScale(self, _dim, _n):
        c = [(self.consts['aho'][0])**(-3), 1/np.pi]
        if _dim==1:
            l = [r"$\large |\Psi_{NN}(\textbf{r})|^{2}\ [\mu m^{-1}]$".replace("NN",str(_n)), r"$\large \varphi_{NN}(\textbf{r})\ [\pi]$".replace("NN",str(_n))]
        elif _dim==2:
            l = [r"$\large |\Psi_{NN}(\textbf{x,y})|^{2}\ [\mu m^{-2}]$".replace("NN",str(_n)), r"$\large \varphi_{NN}(\textbf{r})\ [\pi]$".replace("NN",str(_n))]
        else:
            l = [r"$\large |\Psi_{NN}(\textbf{r})|^{2}\ [\mu m^{-3}]$".replace("NN",str(_n)), r"$\large \varphi_{NN}(\textbf{r})\ [\pi]$".replace("NN",str(_n))]
            
        m = ["speed", "IceFire"]
        return c, l, m
    
    def ReturnVextScale(self):
        c = [self.consts['Eho'][0]/self.consts['kB'][0], 1/np.pi]
        l = [r"$\large V_{ext}(\textbf{r})\ [nK]$", r"$\large \varphi_{V_{ext}}(\textbf{r})\ [\pi]$"]
        m = ['Aggrnyl_r', "twilight"]
        return c, l, m  
    
    def ReturnScale(self, _type, _dim=2):
        switch={
           'density':  self.ReturnDensityScale(_dim),
           'psi':      self.ReturnPsiScale0(_dim),
           'psi_final':self.ReturnPsiScale0(_dim),
           'psi_init':self.ReturnPsiScale0(_dim),
           '1psi':      self.ReturnPsiScale(_dim, 1),
           '1psi_final':self.ReturnPsiScale(_dim, 1),
           '1psi_init':self.ReturnPsiScale(_dim, 1),
           '2psi':      self.ReturnPsiScale(_dim, 2),
           '2psi_final':self.ReturnPsiScale(_dim, 2),
           '2psi_init':self.ReturnPsiScale(_dim, 2),
           'vext':     self.ReturnVextScale(),
           '1vext':     self.ReturnVextScale(),
           '2vext':     self.ReturnVextScale(),
        }
        return switch.get(_type,'Choose observable')




