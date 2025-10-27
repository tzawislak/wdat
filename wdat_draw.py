import numpy as np
import struct
import cmath, math
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
from template import *

#######################################################################
#######################################################################
                        # DRAWING FUNCTIONS #
#######################################################################
#######################################################################

def ReadPsi(self, _type="psi", _cycle=0, _dim=3):
    """Reads a wavefunction (complex binary file)

    Args:
        _type  - type of the data 
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read
        _dime  - dimensionality of the data

    Returns:
        values - values of the wavefunction in physical units
        vsc    - scale units
        unit   - units (strings to be put as labels)
        colorsc- colormaps
    """
    values = self.read(_type, _cycle)
    if ("psi" in _type) :  # Convert |psi| to its square
        values[0] = values[0]**2
    
    vsc, unit, colorsc = self.ReturnScale(_type, _dim)
    # scale 
    values[0] = values[0]*vsc[0]
    values[1] = values[1]*vsc[1]

    return  values, vsc, unit, colorsc


def Draw1Dphase_1c(self, _type="psi", _cycle=0, _zmin=-0.111, _zmax=0.111, _z0=False, _logz=False ):   
    """Draws 1D density and phase for one-component BECs
    
    Args:
        _type  - type of the data
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin  - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax  - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0    - does nothing 
        _logz  - log scale along y

    Returns:
        fig  - the plotly figure
    """
    values1 = self.read(_type, _cycle)

    if ("psi" in _type) :  # Convert |psi| to its square
        values1[0] = values1[0]**2

    # set The scale   
    dim=1
    vsc1, unit, colorsc1 = self.ReturnScale(_type, dim)
    # scale 
    v1 = values1[0].flatten()*vsc1[0]
    p1 = values1[1].flatten()*vsc1[1]

    fig = make_subplots(rows=2, cols=1, horizontal_spacing=0.1, column_widths=[1], shared_xaxes=True)
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v1, mode="lines", line=dict(color='#0043c9', width=4, ), name="density"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p1, mode="lines", line=dict(color='#ca5050', width=4, ), name="phase  "), row=2, col=1  )


    fig.update_layout(title="", showlegend=True, plot_bgcolor='white', )
    fig['layout']['xaxis2'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text="Density") #, range=(_zmin*0.995, _zmax*1.005))
    fig['layout']['yaxis2'].update(title_text="Phase", range=(-1.01, 1.01))


    return fig


def Draw1Dphase_2c(self, _mode="12", _type="psi", _cycle=0, _zmin=-0.111, _zmax=0.111, _z0=False, _logz=False ):   
    """Draws 1D density and phase for two-component BECs
    
    Args:
        _mode  - "12": the two component density will be plotted
                 "ds": the spin and density will be plotted
        _type  - type of the data
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin  - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax  - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0    - does nothing 
        _logz  - log scale along y
        
    Returns:
        fig  - the plotly figure
    """
    val1, vsc1, unit1, colorsc1 =  ReadPsi(self, "1"+_type, _cycle, 1)
    val2, vsc2, unit2, colorsc2 =  ReadPsi(self, "2"+_type, _cycle, 1)

    _1,_2="1","2"
    if _mode=="ds":
        val1 = 0.5*(val1 + val2)
        val2 = 0.5*(val1 - val2)
        _1,_2="d","s"

    vxy1 = self.Section( val1[0], _type="XY")
    pxy1 = self.Section( val1[1], _type="XY")
    v1 = vxy1[int(len(vxy1)/2)]
    p1 = pxy1[int(len(pxy1)/2)]

    vxy2 = self.Section( val2[0], _type="XY")
    pxy2 = self.Section( val2[1], _type="XY")
    v2 = vxy2[int(len(vxy2)/2)]
    p2 = pxy2[int(len(pxy2)/2)]

    if _zmin==-0.111:
       _zmin = np.min(np.concatenate((v1, v2)))
    if _zmax==0.111:
       _zmax = np.max(np.concatenate((v1, v2)))


    ttl = "Sections at surfaces 0"

    fig = make_subplots(rows=2, cols=1, horizontal_spacing=0.1, column_widths=[1], shared_xaxes=True)
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v1, mode="lines", line=dict(color='#0043c9', width=4, ), name="density 1"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v2, mode="lines", line=dict(color='#3971e1', width=4, dash='dot'), name="density 2"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=(v2+v1)/2., mode="lines", line=dict(color='#71c8f4', width=3, dash='dot'), name="n/2"), row=1, col=1  )

    
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p1, mode="lines", line=dict(color='#ca5050', width=4, ), name="phase 1"), row=2, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p2, mode="lines", line=dict(color='#9b1a45', width=4, dash='dot'), name="phase 2"), row=2, col=1  )


    fig.update_layout(title=ttl, showlegend=True, plot_bgcolor='white', )
    fig['layout']['xaxis2'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text="Density", range=(_zmin*0.995, _zmax*1.005))
    fig['layout']['yaxis2'].update(title_text="Phase", range=(-1.01, 1.01))


    return fig

############################################################################################
# 1D data only
def Draw1Dphase_2c_1D(self, _type="psi", _cycle=0, _zmin=0, _zmax=1, _z0=False, _logz=False ):  
    """Draws 1D density and phase for two-component BECs
    
    Args:
        _type   - type of the data
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin  - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax  - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0    - does nothing 
        _logz  - log scale along y
        
    Returns:
        fig  - the plotly figure
    """ 
    values1 = self.read("1"+_type, _cycle)
    values2 = self.read("2"+_type, _cycle)

    if ("psi" in _type) :  # Convert |psi| to its square
        values1[0] = values1[0]**2
        values2[0] = values2[0]**2

    # set The scale   
    dim=1
    vsc1, unit, colorsc1 = self.ReturnScale(_type, dim)
    # scale 
    v1 = values1[0].flatten()*vsc1[0]
    p1 = values1[1].flatten()*vsc1[1]

    v2 = values2[0].flatten()*vsc1[0]
    p2 = values2[1].flatten()*vsc1[1]

    

    fig = make_subplots(rows=2, cols=1, horizontal_spacing=0.1, column_widths=[1])
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v1, mode="lines", line=dict(color='#0043c9', width=4, ), name="density 1"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v2, mode="lines", line=dict(color='#3971e1', width=4, dash='dot'), name="density 2"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=(v2+v1)/2., mode="lines", line=dict(color='#71c8f4', width=3, dash='dot'), name="n/2"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p1, mode="lines", line=dict(color='#ca5050', width=4, ), name="phase 1"), row=2, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p2, mode="lines", line=dict(color='#9b1a45', width=4, dash='dot'), name="phase 2"), row=2, col=1  )


    fig.update_layout(title="", showlegend=True, plot_bgcolor='white', )
    fig['layout']['xaxis2'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text=unit[0]) #, range=(_zmin, _zmax))
    fig['layout']['yaxis2'].update(title_text=unit[1], range=(-1.01, 1.01))


    return fig


def Draw1Dphase(self, _type="psi", _cycle=0, _zmin=0, _zmax=1, _z0=False, _logz=False ): 
    """Draws 1D sections of density and phase for one-component BECs
    
    Args:
        _type   - type of the data
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin  - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax  - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0    - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                 the middle of the array 
        _logz  - log scale along y
        
    Returns:
        fig  - the plotly figure
        v3   - 1D section of the density at the trap center
        p3   - 1D section of the phase at the trap center
    """  
    values = self.read(_type, _cycle)
    if ("psi" in _type) :  # Convert |psi| to its square
        values[0] = values[0]**2

    # set The scale   
    dim=1
    
    vsc, unit, colorsc = self.ReturnScale(_type, dim)
    # scale 
    values[0] = values[0]*vsc[0]
    values[1] = values[1]*vsc[1]

    phases, pvsc, punit, pcolorsc = values[1], vsc[1], unit[1], colorsc[1]  
    values, vsc, unit, colorsc = values[0], vsc[0], unit[0], colorsc[0]  


    vxy = self.Section( values, _type="XY")
    vxz = self.Section( values, _type="XZ" )
    vzy = self.Section( values, _type="YZ" )

    pxy = self.Section( phases, _type="XY")
    pxz = self.Section( phases, _type="XZ" )
    pzy = self.Section( phases, _type="YZ" )
    if _z0==True:
        # find maximum value of YZ plane:
        flat = vzy.flatten()
        indices = np.argpartition(flat, -1)[-1:]
        indices = indices[np.argsort(-flat[indices])]
        (y0, z0) =  np.unravel_index(indices, vzy.shape) 
        v1 = vxz[:, z0[0]]
    else:
        v1 = vxz[:, int(len(vzy[0])/2)]
    v2 = vxy[:, int(len(vxy[0])/2)]
    v3 = vxy[int(len(vxy)/2)]
    p3 = pxy[int(len(pxy)/2)]
    ttl = "Sections at surfaces 0"



    fig = make_subplots(rows=2, cols=1, horizontal_spacing=0.1, column_widths=[1])
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v3, mode="markers+lines"), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=p3, mode="markers+lines"), row=2, col=1  )


    fig.update_layout(title=ttl, showlegend=False, plot_bgcolor='white', )
    fig['layout']['xaxis2'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text=unit, range=(_zmin, _zmax))
    fig['layout']['yaxis2'].update(title_text=punit, range=(-1.01, 1.01))


    return fig, v3, p3
        
#######################################################################
def Draw1Dof2D(self, _mod="sec", _type="psi", _cycle=0, _v=0, _zmin=0, _zmax=1, _z0=False, _logz=False ):   
    """Draws 1D sections of density and phase for one-component BECs
    
    Args:
        _mod   - "sec" - draw a 1D section, "int" - integrate over transverse directions
        _type  - type of the data
        _cycle - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _v     - 0 - draw density, 1 - draw the phase 
        _zmin  - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax  - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0    - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                 the middle of the array 
        _logz  - log scale along y
        
    Returns:
        fig        - the plotly figure
        [v1,v2,v3] - 1D sections (or integrated values) of the density along the three directions
    
    """  
    values = self.read(_type, _cycle)
    if ("psi" in _type) and _v==0:  # Convert |psi| to its square
        values[0] = values[0]**2

    # set The scale   
    dim=3
    if _mod=="int":
        dim=1
    vsc, unit, colorsc = self.ReturnScale(_type, dim)
    # scale 
    values[_v] = values[_v]*vsc[_v]
    values, vsc, unit, colorsc = values[_v], vsc[_v], unit[_v], colorsc[_v]  


    if _mod=="int":
        vxy = values.sum(axis=2)*self.D[2]
        vxz = values.sum(axis=1)*self.D[1]
        vzy = values.sum(axis=0)*self.D[0]
        v3 = vxy.sum(axis=1)*self.D[1]
        v2 = vxy.sum(axis=0)*self.D[0]
        v1 = vzy.sum(axis=0)*self.D[1]
        integral = v3.sum()*self.D[0]
        print("Norm: ", integral)
        ttl = "Integrated View"

    elif _mod=="sec":
        vxy = self.Section( values, _type="XY")
        vxz = self.Section( values, _type="XZ" )
        vzy = self.Section( values, _type="YZ" )
        if _z0==True:
            # find maximum value of YZ plane:
            flat = vzy.flatten()
            indices = np.argpartition(flat, -1)[-1:]
            indices = indices[np.argsort(-flat[indices])]
            (y0, z0) =  np.unravel_index(indices, vzy.shape) 
            v1 = vxz[:, z0[0]]
        else:
            v1 = vxz[:, int(len(vzy[0])/2)]
        v2 = vxy[:, int(len(vxy[0])/2)]
        v3 = vxy[int(len(vxy)/2)]
        ttl = "Sections at surfaces 0"



    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.1, column_widths=[0.2, 0.2, 0.6])
    fig.add_trace(  go.Scatter(y=np.linspace(self.bz[0], self.bz[1], self.N[2]), x=v1), row=1, col=1  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.by[0], self.by[1], self.N[1]), y=v2), row=1, col=2  )
    fig.add_trace(  go.Scatter(x=np.linspace(self.bx[0], self.bx[1], self.N[0]), y=v3, mode="markers"), row=1, col=3  )


    fig.update_layout(title=ttl, showlegend=False, plot_bgcolor='white', )
    fig['layout']['yaxis'].update(title_text=r'$\large Z\ '+ self.lenscale)
    fig['layout']['xaxis2'].update(title_text=r'$\large Y\ '+ self.lenscale)
    fig['layout']['xaxis3'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['xaxis'].update(title_text=unit)
    fig['layout']['yaxis2'].update(title_text=unit)
    fig['layout']['yaxis3'].update(title_text=unit)
    return fig, [v1,v2,v3]

#######################################################################
def Draw2Dof3D(self, values=[], _mod="sec", _type="psi", _cycle=-1, _v=0, _zmin=0, _zmax=1, _logz=False, _passData=False, soloXY=False ):   
    """Draws 2D sections of density and phase for one-component BECs
    
    Args:
        values    - pass your 3D data array, _passData must be set to True
        _mod      - "sec" - draw a 1D section, "int" - integrate over transverse directions
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _v        - 0 - draw density, 1 - draw the phase 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0       - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                    the middle of the array 
        _logz     - log scale along y
        _passData - if True, values passed will be drawn, if False, _type will be read and drawn 
        soloXY    - draw only XY plane 
        
    Returns:
        fig        - the plotly figure
    
    """  


    # set The scale   
    # _mod: int, sec, none
    dim=3
    if _mod=="int":
        dim=2


    vsc, unit, colorsc = self.ReturnScale(_type, dim)
    
    if not _passData:
        values = self.read(_type, _cycle)

        if "psi" in _type:            # Convert |psi| to its square
            values[0] = values[0]**2
        # scale 
        values[_v] = values[_v]*vsc[_v]
    # create reference array for conditional drawing
    reference = values[_v] 

    values, vsc, unit, colorsc = values[_v], vsc[_v], unit[_v], colorsc[_v]        
    if _mod=="int":
        vxy = values.sum(axis=2)*self.D[2]
        vxz = values.sum(axis=1)*self.D[1]
        vzy = values.sum(axis=0)*self.D[0]

        integral = vxy.sum(axis=1)*self.D[1]
        integral = integral.sum()*self.D[0]
        print("Norm: ", integral)

        vxz, vxy = vxz.transpose(), vxy.transpose()

        ttl = "Integrated plots of "+_type
        rvxy = reference.sum(axis=2)*self.D[2]
        rvxy = rvxy.transpose()

    elif _mod=="sec":
        vxy = self.Section( values, _type="XY")
        vxz = self.Section( values, _type="XZ" )
        vzy = self.Section( values, _type="YZ" )
        vzy = vzy.transpose()
        print("Max: ", np.max(vxy))
        print("Min: ", np.min(vxy))
        ttl = "2D Cuts of "+_type
        rvxy = self.Section( reference, _type="XY")
    elif _mod=="none":
        ttl = "data"
        vxy = values

    # remove all entries to the vxy matrix for which absolute value is tiny (useful for phase of the wf)
    #vxy = vxy * ( rvxy>1e-6 )
    #if not _mod=='none':
    #    vxy = np.where( ( rvxy>rvxy/1000 ), vxy, np.zeros((self.N[1], self.N[0])) )

    if _logz:
        vxy, vxz, vzy = np.log10(vxy), np.log10(vxz), np.log10(vzy)
    if _zmax == 0.111:
        _zmax = np.max(vxy)
    if _zmin == -0.111:
        _zmin = np.min(vxy)


    xmian = -1*self.bx[0]+self.bx[1]  -1*self.bz[0]+self.bz[1]
    ymian = -1*self.by[0]+self.by[1]  -1*self.bz[0]+self.bz[1]
    xx=[ (-1*self.bx[0]+self.bx[1])/xmian, (-1*self.bz[0]+self.bz[1])/xmian ]
    yy=[ (-1*self.bz[0]+self.bz[1])/ymian, (-1*self.by[0]+self.by[1])/ymian ]
    # IF only XY plot is desired 
    if soloXY:
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(
        go.Heatmap( z=vxy, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2], showscale=True,  zmin=_zmin, zmax=_zmax,  colorscale=colorsc ), # colorbar={"title": unit, "titleside": "top"},
        row=1, col=1  )
        fig['layout']['xaxis'].update(title_text=r'$\large X\ '+ self.lenscale)
        fig['layout']['yaxis'].update(title_text=r'$\large Y\  '+ self.lenscale)
        # add unit to the colorbar
        #fig.add_annotation(x=1.1, y=1.1,xref="paper",yref="paper",text=unit, showarrow=False,)
        fig.update_layout(height=500, width=540, title_text=ttl, margin=dict(l=80, r=80, t=60, b=60), paper_bgcolor="white",)
        return fig

    # Draw all XY, XZ and YZ
    fig = make_subplots(rows=2, cols=2, shared_yaxes=True, shared_xaxes=True,
                vertical_spacing=0.02, horizontal_spacing=0.02, column_widths=[xx[0], xx[1]], row_heights=[yy[0], yy[1]]) 
    fig.add_trace(
        go.Heatmap( z=vxy, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2], showscale=True,  zmin=_zmin, zmax=_zmax,  colorscale=colorsc ), # colorbar={"title": unit, "titleside": "top"},
        row=2, col=1  )
    
    fig.add_trace(
        go.Heatmap( z=vzy, x0=self.bz[0], dx=self.bz[2], y0=self.by[0], dy=self.by[2], showscale=False, zsmooth='best', zmin=_zmin, zmax=_zmax, colorscale=colorsc ),
        row=2, col=2  ) 

    fig.add_trace(
        go.Heatmap( z=vxz, x0=self.bx[0], dx=self.bx[2], y0=self.bz[0], dy=self.bz[2], showscale=False, zsmooth='best', zmin=_zmin, zmax=_zmax, colorscale=colorsc ),
        row=1, col=1  )

    fig['layout']['xaxis3'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['xaxis4'].update(title_text=r'$\large Z\ '+ self.lenscale)

    fig['layout']['yaxis3'].update(title_text=r'$\large Y\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text=r'$\large Z\  '+ self.lenscale)
    # add unit to the colorbar
    fig.add_annotation(x=1.1, y=1.05,xref="paper",yref="paper",text=unit, showarrow=False,)
    # add arrow indicating the magnetic field orientation
    angl = self.prefix.find('angl', 1)
    if angl  > 0:
        self.DrawMagneticFieldDirection(fig)

    if 'dipole_x' in self.consts:
        print("Dipole x component: ", self.consts['dipole_x'][0])
    fig.update_layout(height=700, width=700, title_text=ttl, margin=dict(l=80, r=120, t=100, b=100), paper_bgcolor="white",)

    return fig

#######################################################################
def DrawComplex(self, _mod="sec", _type="psi", _cycle=0, _zmin=0, _zmax=1, _logz=False, _dG=False, _zn=-1, _zx=1, _yn=-1, _yx=1, _xn=-1, _xx=1, _remove_init=False  ):   
    """Draws 2D sections of density and phase for one-component BECs
    
    Args:
        _mod      - "sec" - draw a 1D section, "int" - integrate over transverse directions
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _v        - 0 - draw density, 1 - draw the phase 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0       - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                    the middle of the array 
        _logz     - log scale along y
        _remove_init - removes initial state (hard coded!) 
        
    Returns:
        fig        - the plotly figure
    
    """ 
    if _remove_init:
        
        #name = self.prefix.replace("R_", '')
        name = self.prefix.replace("pin0.1_", '').replace("R_", '')
        p = Wdat(name + "_info.wtxt", _mute=False)
        init = p.read("psi_final", 0)
        init[0] = init[0]**2    
        
        
    # set The scale   
    dim=3
    if _mod=="int":
        dim=2
    vsc, unit, colorsc = self.ReturnScale(_type, dim)
    values = self.read(_type, _cycle)
    
    
    if "psi" in _type:
        # Convert |psi| to its square
        values[0] = values[0]**2

    if _remove_init:
         values[0] = values[0] - init[0]


    # scale 
    values[0] = values[0]*vsc[0]
    values[1] = values[1]*vsc[1] 
    N = 5

    flattened_array = values[0].flatten()
    indices = np.argpartition(flattened_array, -N)[-N:]
    # Step 3: Sort these indices by their corresponding values
    indices = indices[np.argsort(-flattened_array[indices])]
    # Step 4: Get the top N values
    top_N_values = flattened_array[indices]
    print("Top N values:", top_N_values)
    
    print("Max: ",np.max(values[0]) )

    if _mod=="int":
        mxy = values[0].sum(axis=2)*self.D[2]
        integral = mxy.sum(axis=1)*self.D[1]
        integral = integral.sum()*self.D[0]
        print("Norm: ", integral)
        mxy = mxy.transpose()
        axy = values[1].sum(axis=2)*self.D[2]
        axy = axy.transpose()
        ttl = "Integrated plots of "+_type

    elif _mod=="sec":
        mxy = self.Section( values[0], _type="XY")
        axy = self.Section( values[1], _type="XY")
        ttl = "2D Cuts of "+_type


    if _logz:
        mxy = np.log10(mxy)

    if _zmax == 0.111:
        _zmax = np.max(mxy)
    if _zmin == -0.111:
        _zmin = np.min(mxy)

    mxymax = np.max(mxy)
    # remove all entries to the vxy matrix for which absolute value is tiny (useful for phase of the wf)
    axy = np.where( ( mxy>mxymax/1000 ), axy, np.zeros(axy.shape) )

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.2 ) 

    fig.add_trace(
        go.Heatmap( z=mxy, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2], showscale=True, zsmooth='best', zmin=_zmin, zmax=_zmax,  colorscale=colorsc[0], colorbar={ "x": 0.4, "len": 0.95}),
        row=1, col=1  )

    fig.add_trace(
        go.Heatmap( z=axy, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2], zmax=_zx, zmin=_zn, showscale=True,colorscale=colorsc[1], colorbar={ "x": 1.0, "len": 0.95} ),
        row=1, col=2  ) 
    contour, shade = self.GetContourTrace(mxy)
    fig.add_trace( shade  , row=1, col=2  )
    fig.add_trace( contour, row=1, col=2  )

    if _dG:
        import plotly.figure_factory as ff
        axy = np.where( ( mxy>mxymax/100 ), axy, np.zeros(axy.shape) )
        mxy = mxy[:-1, :-1]
        qx,qy = np.meshgrid(np.arange(self.bx[0], self.bx[1]-self.bx[2], self.bx[2]), 
                            np.arange(self.by[0], self.by[1]-self.by[2], self.by[2]))
        
        gx, gy = self.ReturnXYGrad(axy)
        print(mxy.shape)
        gx = np.where( ( mxy>mxymax/100 ), gx, np.zeros(gx.shape) )
        gy = np.where( ( mxy>mxymax/1 ), gy, np.zeros(gy.shape) )
        
        # ring
        #gx, gy = gx[8:-8:8, 8:-8:8], gy[8:-8:8, 8:-8:8]
        #qx, qy = qx[1:-1, 1:-1], qy[1:-1, 1:-1]
        # tube
        #gx, gy = gx[8:-8:8, 4:-4:4], gy[8:-8:8, 4:-4:4] 
        #qx, qy = qx[1:-1:1, 4:-4:4], qy[1:-1:1, 4:-4:4]
        xx, yy = 2, 64
        print( (gx[::xx, ::yy]).shape , (qx[::xx, ::yy]).shape )
        
        fiq = ff.create_quiver(qx[::xx, ::yy], qy[::xx, ::yy], gx[::xx, ::yy], gy[::xx, ::yy], scale = 1.15, arrow_scale = 1.15, 
                               name = 'quiver', line = dict(width=1, color="red"))
        fig.add_traces(data = fiq.data)

    sx=[_xn, _xx]
    sy=[_yn, _yx]
    fig['layout']['xaxis'].update(title_text=r'$\large X\ '+ self.lenscale, range=sx)
    fig['layout']['xaxis2'].update(title_text=r'$\large X\ '+ self.lenscale, range=sx)
    fig['layout']['yaxis'].update(title_text=r'$\large Y\  '+ self.lenscale, range=sy)
    fig['layout']['yaxis2'].update(title_text=r'$\large Y\ '+ self.lenscale, range=sy)

    # add unit to the colorbar
    fig.add_annotation(x=0.44, y=1.1,xref="paper",yref="paper",text=unit[0], showarrow=False,)
    fig.add_annotation(x=1.08, y=1.1,xref="paper",yref="paper",text=unit[1], showarrow=False,)

    fig.update_layout(height=500, width=1000, title_text=ttl, margin=dict(l=80, r=120, t=100, b=100), paper_bgcolor="white",)

    return fig
######################################################################
def GetContourTrace(self, _d, _N=2):
    """Draws 2D contours
    
    Args:
        _d - density
        _N - number of contours
    Returns:
        trc   - plotly contour figure
        shade - shade to mask irrelevant values
    
    """ 
    _d = _d - np.mean(_d)
    dmax = np.max(_d)
    shade = go.Contour(z=_d, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2],
                     line_width=0, autocontour=False, opacity=0.4, showlegend=False, fillcolor="black",
                     contours=dict(
                        type="constraint", 
                        coloring='none',
                        operation=">",
                        value=0.1*dmax
                    ))
    trc =  go.Contour(z=_d, x0=self.bx[0], dx=self.bx[2], y0=self.by[0], dy=self.by[2],
                      line_width=1, line_color='white', autocontour=False, showlegend=False,
                      ncontours=_N, contours=dict(
                               coloring='none', start= 0.3*dmax, end  = 0.92*dmax,
                      ))
    return trc, shade

######################################################################
def DrawRT(self, _mod="sec", _type="psi", _cycle=0, _zmin=-0.111, _zmax=0.111, _logz=False, _dG=False, _zn=-1, _zx=1, _yn=-1, _yx=1, _xn=-1, _xx=1, _remove_init=False  ): 
    """Draws 2D sections of density and phase time evolution for one-component BECs
    
    Args:
        _mod      - "sec" - draw a 1D section, "int" - integrate over transverse directions
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0       - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                    the middle of the array 
        _logz     - log scale along y
        _remove_init - removes initial state (hard coded!) 
        
    Returns:
        fig        - the plotly figure
    
    """ 

    is1D=False
    if _mod == "1dsec":
        _mod = "sec"
        is1D = True
    
    for it in range(0, self.cycles, 1):
        print("\rCreate frame: ", it+1, "/", self.cycles)

        if is1D:
            f, st = self.Draw1Dphase(_type=_type, _cycle=it, _zmin=_zn, _zmax=_zx, _z0=False, _logz=_logz )
        else:
            f = self.DrawComplex( _mod=_mod, _type=_type, _cycle=it,
                             _zmin=_zmin, _zmax=_zmax, _logz=_logz, _dG=_dG, _zn=_zn, _zx=_zx, _yn=_yn, _yx=_yx, _xn=_xn, _xx=_xx, _remove_init=_remove_init )
        time = self.t0 + it*self.dt*self.consts['itmod'][0]
        timetxt = "time: "+format(time, '04.3f')+"ms"
        f.add_annotation(x=0.5, y=1.15,xref="paper",yref="paper",text=timetxt, showarrow=False)
        prfx = self.prefix[self.prefix.rfind('/')+1:]
        f.write_image(prfx + _type+ "_t"+format(time, '08.3f')+".png")#, width=800, height=800)

######################################################################
def DrawRT_2c(self, _mod="sec", _type="psi", _cycle=0, _zmin=-0.111, _zmax=0.111, _logz=False, _dG=False, _zn=-1, _zx=1, _yn=-1, _yx=1, _xn=-1, _xx=1, _remove_init=False  ): 
     """Draws 2D sections of density and phase time evolution for one-component BECs
    
    Args:
        _mod      - "sec" - draw a 1D section, "int" - integrate over transverse directions
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0       - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                    the middle of the array 
        _logz     - log scale along y
        _remove_init - removes initial state (hard coded!) 
                
    """ 
    is1D=False
    if _mod == "1dsec":
        _mod = "sec"
        is1D = True
    
    for it in range(0, self.cycles, 1):
        print("\rCreate frame: ", it+1, "/", self.cycles)

        if is1D:
            f = self.Draw1Dphase_2c(_type=_type, _cycle=it, _zmin=_zmin, _zmax=_zmax, _z0=False, _logz=_logz )
        else:
            f = self.DrawComplex( _mod=_mod, _type=_type, _cycle=it,
                             _zmin=_zmin, _zmax=_zmax, _logz=_logz, _dG=_dG, _zn=_zn, _zx=_zx, _yn=_yn, _yx=_yx, _xn=_xn, _xx=_xx, _remove_init=_remove_init )
        time = self.t0 + it*self.dt*self.consts['itmod'][0]
        timetxt = "time: "+format(time, '04.3f')+"ms"
        f.add_annotation(x=0.5, y=1.15,xref="paper",yref="paper",text=timetxt, showarrow=False)
        prfx = self.prefix[self.prefix.rfind('/')+1:]
        f.write_image(prfx + _type+ "_t"+format(time, '08.3f')+".png")#, width=800, height=800)

######################################################################
def Draw3D(self, vals, _x, _y, _z, figpath="", _isomin=0, _isomax=1000, _opc=0.4, _sc=20):
    """Draws 3D data
    
    Args:
        vals    - the 3D data to plot
        _x      - trim the vals from both sides along the X direction
        _y      - trim the vals from both sides along the Y direction
        _z      - trim the vals from both sides along the Z direction
        figpath - does nothing 
        _isomin - minimum value of the isosurfaces
        _isomax - maximum value of the isosurfaces
        _opc    - opacity
        _sc     - does nothing
        
    Returns:
        fig     - the plotly figure
    
    """ 
    X, Y, Z = np.mgrid[ -1*_x:_x:(self.bx[2]), -1*_y:_y:self.by[2], -1*_z:_z:self.bz[2]]
    print(len(X), len(vals), "   ", len(Y[0]), len(vals[0]), "   ", len(Z[0][0]), len(vals[0][0]) )
    fig = go.Figure(data=go.Volume(  x=X.flatten(), y=Y.flatten(), z=Z.flatten(), value=vals.flatten(),
        isomin=_isomin,
        isomax=_isomax,
        opacity=_opc,     # needs to be small to see through all surfaces
        surface_count=90, # needs to be a large number for good volume rendering
        ))
    fig.update_layout(
        plot_bgcolor='white',
        title="",
        scene = dict(  xaxis_title=r'$\large X\ '+ self.lenscale,
                       yaxis_title=r'$\large Y\ '+ self.lenscale,
                       zaxis_title=r'$\large Z\ '+ self.lenscale),
                       width=900,
                       margin=dict(r=20, b=10, l=10, t=100))
    return fig

######################################################################
def DrawStatus(self, isreal=False):
    """Draws status data from *status.out
    
    Args:
        isreal - defines the data structure from either real or imaginary time evolution
        
    Returns:
        fig     - the plotly figure
    
    """ 
    stat_file = self.prefix + '_status.out'

    stat = pd.read_csv(stat_file, delimiter =r"\s+")
    xax = 'Iter.'
    
    if len(stat.axes[1]) > 4: 
        sttls = ("Norm", "Chem. pot.", "Energy", "E convergence", "Lz/N/hbar")
        if isreal:
            xax = "Time"
            sttls = ("Norm", "Chem. pot.", "Energy", "dE", "Px/N/hbar", "cos")
        fig = make_subplots(rows=3, cols=2, subplot_titles=sttls)
        #stat[xax]=stat[xax][1:]
        #stat['Norm']=stat['Norm'][1:]
        #stat['dE']=stat['dE'][1:]
        #stat['Chem']=stat['Chem'][1:]
        #stat['Energy']=stat['Energy'][1:]
        #stat['Lz']=stat['Lz'][1:]
        #stat['Px>'] = stat['<Px>'] #/ self.consts['hbar'][0]
        fig.add_trace(go.Scatter(x=stat[xax], y=100*(stat['Norm']-stat['Norm'][0])/stat['Norm'][0], name=r'$\Delta N/N_{0}\ \%$'), 1, 1) 
        fig.add_trace(go.Scatter(x=stat[xax], y=100*(stat['mu']-stat['mu'][0])/stat['mu'][0], name=r'$\Delta \mu/\mu_{0}\ \%$'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=100*(stat['Energy']-stat['Energy'][0])/stat['Energy'][0], name=r'$\Delta E/E_{0}\ \%$'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=100*(stat['Ediff']/stat['Energy'][0]), name=r'$(E(t_{i})-E(t_{i-1}))/E_{0}\ \%$'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Px'], name=r'$ <Px>\ [\hbar]$'), 3, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['<cos>'], name=r'$ <cos>\ [?]$'), 3, 2)

        #fig['layout']['yaxis3'].update(type="log")
        fig.update_layout(legend=dict( yanchor="top",  y=0.3, xanchor="right",x=0.99))
    else:
        sttls = ("Norm", "Chem. pot.", "Energy", "E convergence")
        if isreal:
            xax = 'Time'
            sttls = ("Norm", "Chem. pot.", "Energy", "dE")
        fig = make_subplots(rows=2, cols=2, subplot_titles=sttls)

        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm'], name='Norm'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu'], name='Chem. pot.'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy'], name='Energy'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff'], name='dE'), 2, 2)
        fig['layout']['yaxis3'].update(type="log")

    return fig
######################################################################
def DrawStatus_2c(self, _obs="", isreal=False, mode="Simple"):
    """Draws status data from *status.out for 2-component BECs
    
    Args:
        _obs   - additional observable name to be plotted
        isreal - defines the data structure from either real or imaginary time evolution
        mode   - "Simple" - normal 2-component gas, "Rabi" - Rabi coupled gas
    Returns:
        fig     - the plotly figure
    
    """ 
    stat_file = self.prefix + '_status.out'

    stat = pd.read_csv(stat_file, delimiter=r"\s+")

    if isreal:
        # REAL TIME EVOLUTION STATUS

        sttls = ("Norm", "Chem. pot.", "Energy", "dE", "iter. time") #, r"$\langle\cos\rangle$")
        if _obs != "":
            sttls += (_obs,)

        fig = make_subplots(rows=3, cols=2, subplot_titles=sttls)

        xax = 'Time'
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm1'], mode="lines", name='1', legend='legend1'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend1'), 1, 1)
        
        if mode == "Simple":    
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu1'], mode="lines", name='1', legend='legend2'), 1, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend2'), 1, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy1'], mode="lines", name='1', legend='legend3'), 2, 1)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend3'), 2, 1)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff1'], mode="lines", name='1', legend='legend4'), 2, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='1', legend='legend4'), 2, 2)
        if mode == "Rabi":
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm1']  +stat['Norm2'], mode="lines", name='sum', legend='legend1'), 1, 1)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu1']    +stat['mu2'],     mode="lines", name='mu', legend='legend2'), 1, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy1']+stat['Energy2'], mode="lines", name='E', legend='legend3'), 2, 1)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff1'] +stat['Ediff2'] , mode="lines", name='dE', legend='legend4'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['ItTime'], mode="lines", name='IT', legend='legend5'), 3, 1)
       
        if _obs=="<cos_d>":
            fig.add_trace(go.Scatter(x=stat[xax], y=stat["<cos_d>"], mode="markers", name="<cos_d>", legend='legend6'), 3, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat["<cos_s>"], mode="lines", name="<cos_s>", legend='legend6'), 3, 2)
        elif _obs=="<cos_1>":
            fig.add_trace(go.Scatter(x=stat[xax], y=stat["<cos_1>"], mode="markers", name="<cos_1>", legend='legend6'), 3, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat["<cos_2>"], mode="lines", name="<cos_2>", legend='legend6'), 3, 2)
        else:
            try:
                fig.add_trace(go.Scatter(x=stat[xax], y=stat[_obs+"1"], mode="lines", name=_obs+"1", legend='legend6'), 3, 2)
                fig.add_trace(go.Scatter(x=stat[xax], y=stat[_obs+"2"], mode="lines", name=_obs+"2", legend='legend6'), 3, 2)
            except:
                fig.add_trace(go.Scatter(x=stat[xax], y=stat[_obs], mode="lines", name=_obs, legend='legend6'), 3, 2)

     
        fig['layout']['yaxis4'].update(type="log")
        fig.update_layout({
            'legend1': { 'x': 0.4, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend2': { 'x': 0.94, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend3': { 'x': 0.4, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend4': { 'x': 0.94, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend5': { 'x': 0.4, 'y': 0.3, 'xanchor': 'left', 'yanchor': 'top', },
        })

        if _obs != "":
            fig.update_layout({'legend6': { 'x': 0.94, 'y': 0.2, 'xanchor': 'left', 'yanchor': 'top', },})

    else:
         # IMAGINARY TIME EVOLUTION STATUS

        sttls = ("Norm", "Chem. pot.", "Energy", "E convergence", "time per. it.")
        if _obs != "":
            sttls += (_obs,)

        fig = make_subplots(rows=3, cols=2, subplot_titles=sttls)
        xax="Iter."
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm1'], mode="lines", name='1', legend='legend1'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend1'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu1'], mode="lines", name='1', legend='legend2'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend2'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy1'], mode="lines", name='1', legend='legend3'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend3'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff1'], mode="lines", name='1', legend='legend4'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff2'], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend4'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Time'], mode="lines", name='IT', legend='legend5'), 3, 1)

        if _obs != "":
            fig.add_trace(go.Scatter(x=stat[xax], y=stat[_obs+"1"], mode="lines", name='1', legend='legend6'), 3, 2)
            fig.add_trace(go.Scatter(x=stat[xax], y=stat[_obs+"2"], mode="markers", marker=dict(symbol='circle-open', size=6), name='2', legend='legend6'), 3, 2)

        fig['layout']['yaxis4'].update(type="log")
        fig.update_layout({
            'legend1': { 'x': 0.4, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend2': { 'x': 0.94, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend3': { 'x': 0.4, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend4': { 'x': 0.94, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend5': { 'x': 0.4, 'y': 0.2, 'xanchor': 'left', 'yanchor': 'top', },
        })

        if _obs != "":
            fig.update_layout({'legend6': { 'x': 0.94, 'y': 0.2, 'xanchor': 'left', 'yanchor': 'top', },})

    return fig

######################################################################
def DrawStatus_1c(self, isreal=False):
    """Draws status data from *status.out
    
    Args:
        isreal - defines the data structure from either real or imaginary time evolution
        
    Returns:
        fig     - the plotly figure
    
    """ 
    stat_file = self.prefix + '_status.out'

    stat = pd.read_csv(stat_file, delimiter=r"\s+")

    if isreal:
        # REAL TIME EVOLUTION STATUS

        sttls = ("Norm", "Chem. pot.", "Energy", "dE", "iter. time", r"$\langle\cos\rangle$")
        fig = make_subplots(rows=3, cols=2, subplot_titles=sttls)

        xax = 'Time'
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm'], mode="lines", name='1', legend='legend1'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu'], mode="lines", name='1', legend='legend2'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy'], mode="lines", name='1', legend='legend3'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff'], mode="lines", name='1', legend='legend4'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['ItTime'], mode="lines", name='IT', legend='legend5'), 3, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['<cos>'], mode="lines", name='<cos>', legend='legend6'), 3, 2)

        fig['layout']['yaxis4'].update(type="log")
        fig.update_layout({
            'legend1': { 'x': 0.4, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend2': { 'x': 0.94, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend3': { 'x': 0.4, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend4': { 'x': 0.94, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend5': { 'x': 0.4, 'y': 0.3, 'xanchor': 'left', 'yanchor': 'top', },
            'legend6': { 'x': 0.94, 'y': 0.3, 'xanchor': 'left', 'yanchor': 'top', },

        })
    else:
         # IMAGINARY TIME EVOLUTION STATUS
        sttls = ("Norm", "Chem. pot.", "Energy", "E convergence", "time per. it.")
        fig = make_subplots(rows=3, cols=2, subplot_titles=sttls)
        xax="Time"
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Norm'], mode="lines", name='1', legend='legend1'), 1, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['mu'], mode="lines", name='1', legend='legend2'), 1, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Energy'], mode="lines", name='1', legend='legend3'), 2, 1)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['Ediff'], mode="lines", name='1', legend='legend4'), 2, 2)
        fig.add_trace(go.Scatter(x=stat[xax], y=stat['ItTime'], mode="lines", name='IT', legend='legend5'), 3, 1)

        fig['layout']['yaxis4'].update(type="log")
        fig.update_layout({
            'legend1': { 'x': 0.4, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend2': { 'x': 0.94, 'y': 1.1, 'xanchor': 'left', 'yanchor': 'top', },
            'legend3': { 'x': 0.4, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend4': { 'x': 0.94, 'y': 0.6, 'xanchor': 'left', 'yanchor': 'top', },
            'legend5': { 'x': 0.4, 'y': 0.2, 'xanchor': 'left', 'yanchor': 'top', },
        })

    return fig
        
######################################################################
def DrawMagneticFieldDirection(self, figs):
    """Draws an arrow with the direction of the magnetic field for dipolar BECs
    
    Args:
        figs - and existing plotly figure
    
    """ 
    angl = float(self.prefix[self.prefix.find('angl', 1)+4:self.prefix.find('angl', 1)+6])
    print( angl )
    figs.add_shape(type="circle", xref="paper", yref="paper", fillcolor="PaleTurquoise", line_color="LightSeaGreen",
        x0=0.02, y0=0.02, x1=0.12, y1=0.12, line_width=2    )
    figs.add_shape(type="circle", xref="paper", yref="paper", fillcolor="PaleTurquoise", line_color="LightSeaGreen",
        x0=0.02, y0=0.98, x1=0.12, y1=0.88, line_width=2    )
    figs.add_shape(type="circle", xref="paper", yref="paper", fillcolor="PaleTurquoise", line_color="LightSeaGreen",
        x0=0.98, y0=0.02, x1=0.88, y1=0.12, line_width=2    )
    x=np.sin(angl/360*2*np.pi) * 0.05
    z=np.cos(angl/360*2*np.pi) * 0.05
    print(x, z)
    figs.add_shape(type="line",xref="paper", yref="paper", x0= 0.07+x, y0= 0.93+z, x1=0.07, y1=0.93,
        line=dict(color="LightSeaGreen",width=2))
    figs.add_shape(type="line",xref="paper", yref="paper", x0= 0.07+x, y0= 0.07, x1=0.07, y1=0.07,
        line=dict(color="LightSeaGreen",width=2))
    figs.add_shape(type="line",xref="paper", yref="paper", x0= 0.93+z, y0= 0.07, x1=0.93, y1=0.07,
        line=dict(color="LightSeaGreen",width=2))

    

    
#######################################################################
def Draw1D(self, _type="psi", _cycle=0, _v=0, _zmin=0, _zmax=1, _logz=False ):  
    """Draws 1D density of a 1D one-component BECs
    
    Args:
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _v        - 0 - draw density, 1 - draw the phase 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _logz     - log scale along y
        
    Returns:
        fig        - the plotly figure
    
    """ 
    values = self.read(_type, _cycle)
    if ("psi" in _type) and _v==0:  # Convert |psi| to its square
        values[0] = values[0]**2
    print(values[0])
    dim=1
    vsc, unit, colorsc = self.ReturnScale(_type, dim)
    # scale 
    values[_v] = values[_v]*vsc[_v]
    values, vsc, unit, colorsc = values[_v], vsc[_v], unit[_v], colorsc[_v]  

    fig = make_subplots(rows=1, cols=1, horizontal_spacing=0.1)
    fig.add_trace(  go.Scatter(y=np.linspace(self.bx[0], self.bx[1], self.N[0]), x=values), row=1, col=1  )

    fig.update_layout(title="", showlegend=False, plot_bgcolor='white' )
    fig['layout']['xaxis'].update(title_text=r'$\large X\ '+ self.lenscale)
    fig['layout']['yaxis'].update(title_text=unit)
    return fig




#######################################################################
def Draw1DVelocity(self, _rm_vs=0, _rm_vn=0, _type="psi", _cycle=0, _v=0, _zmin=0, _zmax=1, _logz=False ):   
    """Draws a 1D velocity field 
    
    Args:
        _rm_vs    - value of the v_s to be removed
        _rm_vn    - value of the v_n to be removed
        _type     - type of the data
        _cycle    - if the *psi.wdat contains more than one wavefunction, specifies the index to be read 
        _v        - 0 - draw density, 1 - draw the phase 
        _zmin     - sets minimum value of y axis ploted ( default value indicates lowest from the data set will be choosen)
        _zmax     - sets maximum value of y axis ploted ( default value indicates highest from the data set will be choosen)
        _z0       - find the index along the Z axis, at which the densit is the highest. Otherwise asumes this index is in 
                    the middle of the array 
        _logz     - log scale along y
        _remove_init - removes initial state (hard coded!) 
        
    Returns:
        fig        - the plotly figure
    
    """ 
    fig, v3, p3 = self.Draw1Dphase(_type=_type, _cycle=_cycle, _zmin=0, _zmax=1, _z0=False, _logz=_logz )
    hbar, mass, Eho = self.consts['hbar'][0], self.consts['mass'][0], self.consts['Eho'][0]
    x = np.arange( -self.N[0]*self.D[0]/2, self.N[0]*self.D[0]/2, self.D[0])
    vncr = hbar/(2*mass * (96/(2*np.pi))) 

    v_moving = _rm_vn/2.5696508420412703 * vncr
    vs       = _rm_vs * 2*np.pi*hbar/(96*mass)
    ddp3     = np.gradient(p3, self.D[0]) * hbar/mass *2*np.pi - v_moving - vs
    

    fig = go.Figure()
    fig.add_trace( go.Scatter(x=x, y=ddp3 /vncr, name="v_s = 0", line=dict(width=2, color="red")))

    fig.update_xaxes(title=r"$x\ [\mu\text{m}]$", range=(-10, 10))

    fig.update_layout(showlegend=True,width=800, height=300, template=sc_template, title=r"", )  # xgap and ygap for spacing
    fig.update_yaxes(title=r"$\frac{\hbar}{m}\nabla\phi(x)\ [v_n^{\text{cr}}]$")
    return fig, ddp3 /vncr

from wdat_plot import Wdat

Wdat.ReadPsi = ReadPsi
Wdat.Draw1Dphase_1c = Draw1Dphase_1c
Wdat.Draw1Dphase_2c = Draw1Dphase_2c
Wdat.Draw1Dphase_2c_1D = Draw1Dphase_2c_1D
Wdat.DrawMagneticFieldDirection = DrawMagneticFieldDirection
Wdat.DrawStatus = DrawStatus
Wdat.DrawStatus_1c = DrawStatus_1c
Wdat.DrawStatus_2c = DrawStatus_2c
Wdat.Draw1DVelocity = Draw1DVelocity
Wdat.Draw3D = Draw3D
Wdat.DrawRT = DrawRT
Wdat.DrawRT_2c = DrawRT_2c
Wdat.GetContourTrace = GetContourTrace
Wdat.DrawComplex = DrawComplex
Wdat.Draw2Dof3D = Draw2Dof3D
Wdat.Draw1Dof2D = Draw1Dof2D
Wdat.Draw1D = Draw1D
Wdat.Draw1Dphase = Draw1Dphase
