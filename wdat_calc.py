import numpy as np
import struct
import cmath, math
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.special as sc
from scipy.interpolate import interp1d

from scipy import optimize
import pandas as pd



######################################################################
def CalcCoM(self, _type="psi", _y='Y', _d=1 ):
    """Calculates Center-of-mass in the XY or XZ plane
    
    Args:
        _type  - type of the data to be read
        _y     - choose between 'Y' and 'Z' as the second direction to be determined 
        _d     - file tag number 
        
    """ 
    filename = self.prefix
    filename = filename[filename.rfind('/')+1:]
    tag = '_COM'
    if _d != 1: 
        tag = tag +str(_d)
    f = open( filename + tag + ".txt", "w")
    print("Time Xcm "+_y+"cm", file=f)

    for it in range( self.cycles-1 ): 
        values = self.read(_type, it)
        if "psi" in _type: # Convert |psi| to its square
            values = values[0]**2
            
        vsc, unit, colorsc = self.ReturnScale(_type, 2)
        # scale 
        values = values*vsc[0]
        vsc, unit, colorsc = vsc[0], unit[0], colorsc[0]        
        pos = np.unravel_index(values.argmax(), values.shape)
        print(values.shape)        
        print(pos)
        pos = np.unravel_index(values.argmin(), values.shape)
        print(pos)    
        if _y == 'Y':
            vxy = values.sum(axis=2)*self.D[2]
            vx = vxy.sum(axis=1)*self.D[1]
            vy = vxy.sum(axis=0)*self.D[0]

            xs = np.mgrid[ self.bx[0]:self.bx[1]:self.bx[2] ]
            ys = np.mgrid[ self.by[0]:self.by[1]:self.by[2] ]


            # now calculate the center of mass:
            xcm = (vx*(xs)**(_d)).sum() / vx.sum()
            ycm = (vy*(ys)**(_d)).sum() / vy.sum()
        
        elif _y == 'Z':
            vxz = values.sum(axis=1)*self.D[1]
            vx = vxz.sum(axis=1)*self.D[2]
            vz = vxz.sum(axis=0)*self.D[0]

            xs = np.mgrid[ self.bx[0]:self.bx[1]:self.bx[2] ]
            zs = np.mgrid[ self.bz[0]:self.bz[1]:self.bz[2] ]

            # now calculate the center of mass:
            xcm = (vx*(xs)**(_d)).sum() / vx.sum()
            ycm = (vz*(zs)**(_d)).sum() / vz.sum()
        
        time = self.t0 + it*self.dt*self.consts['itmod'][0]
        timetxt = "time: "+format(time, '04.3f')+"ms"
        print(time, xcm, ycm, file=f)
        
    f.close()





######################################################################
def AddVariableToWtxt(self, _name, _value, _unit, _replace=False):
    """Appends a new variable to the .wtxt file
    
    Args:
        _name    - name of the variable
        _value   - its value 
        _unit    - its unit
        _replace - to be appended (False) or replaced (True)
        
    """ 
    test = Wdat(self.prefix+"_info.wtxt", _mute=True)
    if not _name in test.consts:
        print("Writing {n} to {p} ...".format(n=_name, p=self.prefix+"_info.wtxt"))
        with open(self.prefix+"_info.wtxt", "a") as myfile:
            myfile.write("const           {n}                      {v}            {u}\n".format(n=_name, v=_value, u=_unit) )

    




def CountDroplets(self, _type="psi_final"): 
    """Counts the number of peaks in the density profile

    Args:
        _type  - type of the data to be read

    """ 
    filename = self.prefix[self.prefix.rfind('/')+1:]
    vsc, unit, colorsc = self.ReturnScale(_type, 1)
    # scale 
    vsc, unit, colorsc = vsc[0], unit[0], colorsc[0]

    xX = np.mgrid[ (self.bx[0]):(self.bx[1]):(self.bx[2]) ]
    
    values = self.read(_type, 0)
    if "psi" in _type: # Convert |psi| to its square
        values = values[0]**2
    values = values*vsc
    # integrate over z
    vx = values.sum(axis=2).sum(axis=1)*self.D[1]*self.D[2]


    # shrink the array to  a meaningful size
    maxima_indices = []
    n = len(vx)
    for i in range(1, n - 1):
        if vx[i] > vx[i - 1] and vx[i] > vx[i + 1]:
            maxima_indices.append(i)


    print("N droplets: ", len(maxima_indices) )
    self.AddVariableToWtxt("Ndrpl", len(maxima_indices), "1")



#################################################################################################3
def CalcfsTube(self, _type="psi_final", _it=0, _int=[], _InterpolationFactor=4):
    """Calculates Leggett's upper bound on the superfluid fraction in a tube
    
    Args:
        _type  - type of the data to be read
        _it    - cycle number (for intermediate wavefunctions from *_psi.wdat)
        _int   - x-axis range limits 
        _InterpolationFactor - Interpolates the density profile with (by default) 4 times more points for a better precision
        
    Returns:
        fs     - superfluid fraction
    """ 
    filename = self.prefix[self.prefix.rfind('/')+1:]
    psi = self.read(_type, _it)
    
    if "psi" in _type: # Convert |psi| to its square
        psi = psi[0]**2
        
    vsc, unit, colorsc = self.ReturnScale(_type, 1)
    # scale 
    psi = psi*vsc[0]
    vsc, unit, colorsc = vsc[0], unit[0], colorsc[0]    

    # integrate over zy
    rhox = psi.sum(axis=2).sum(axis=1)*self.D[2]*self.D[1] / self.consts['npart'][0]
    inv_rhox = np.reciprocal(rhox) 

    if _int != [] :
        # interpolate the density
        x_raw = np.linspace(0, (self.N[0])*self.D[0], self.N[0] , endpoint=False)  # 16 points over one period of cosine

        x_int = np.linspace(x_raw.min(), x_raw.max(), (self.N[0])*_InterpolationFactor, endpoint=False)  # 64 points within the same range
        interpolator = interp1d(x_raw, rhox, kind='cubic')    # Linear interpolation
        rhox_int = interpolator(x_int)
        # limit the range
        rhox_int = rhox_int[_int[0]:_int[1]] * self.consts['Ndrpl'][0]

        inv_rhox = np.reciprocal(rhox_int) 
        
        L = abs(_int[0]-_int[1])* self.D[0] / _InterpolationFactor**(0.5)
        fs = inv_rhox.sum() * self.D[0] / (L)**2 # has to be squared, see Legget and change primed variables to not primed
        fs = 1/fs

       
    else:
        fs = inv_rhox.sum() * self.D[0] / (2*np.pi*self.consts['ringRho'][0])**2 # has to be squared, see Legget and change primed variables to not primed
        fs = 1/fs
        print("fs = ", fs)
        self.AddVariableToWtxt("Lfs", fs, "1")
    
    return fs



from wdat_plot import Wdat
Wdat.CalcCoM = CalcCoM
Wdat.AddVariableToWtxt = AddVariableToWtxt
Wdat.CountDroplets = CountDroplets
Wdat.CalcfsTube = CalcfsTube




