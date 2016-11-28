# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 19:44:52 2016

@author: mark
"""
#this is a bunch of stuff needed for 3d plotting...
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
#end, plot support libraries.

import numpy as np

# This tool is intended to plot a beam and a few nulls,
#arrayed across a 2-d circular grid 
#maybe later with "space taper"
#applying a gradient search in direction determined by SVD of the 
#"Jacobian" that is, the partial derivatives of nulls w.r.t. each of the 
#AIP phasors.




   
#-----------Supporting class-------------------------------------------

                         
                         
class beam(object):
    def __init__(self):
        # here's the aperture AIP locations, in m.
        # Here are some coordinates to start with, for a 7x7 array.
        # Note the center element is at Xgrid=0, Ygrid=0.
        # Later, I should pass these in, instead.
    
        X = np.arange(-10, 10.1, 0.5);
        Y = np.arange(-10, 10.1, 0.5);
        X, Y = np.meshgrid(X, Y);

        self.Ygrid=Y;  self.Xgrid=X;
        
#        self.Ygrid =np.matrix([[3.,3.,3.,3.,3.,3.,3.],
#                         [2.,2.,2.,2.,2.,2.,2.],
#                         [1.,1.,1.,1.,1.,1.,1.],
#                         [0.,0.,0.,0.,0.,0.,0.],
#                         [-1.,-1.,-1.,-1.,-1.,-1.,-1.],
#                         [-2.,-2.,-2.,-2.,-2.,-2.,-2.],
#                         [-3.,-3.,-3.,-3.,-3.,-3.,-3.]]);
#        self.Xgrid= np.matrix([[-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.],
#                         [-3.,-2.,-1.,0.,1.,2.,3.]]);
        
        self.elements = self.Ygrid.size;  #number of elements
        
        self.wavelength = 1.;  #meters
               
        #This initial condition sends a beam out the boresight.  Remember, the
        #zeros are the phase angles of the unit power vectors at each AIP.
        self.BeamDelays= np.zeros(self.elements);
      
    def calcGain(self,az,el):
        yVals=self.AzRotate(az); #get new y-values at this Az.
        
        #Now get the geometric delays from this off-axis sidelobe direction.
        SideLobeDelays = np.matrix(self.BeamDelays);
        
        c= 3.e8;  #m/sec, speed of light
        
        #how much addtl delay must we add for every meter in -x
        delayRamp = np.sin(el)/c;  #something like this
    
        #Now go through the AIPs, adding geometric delay proportionate to the 
        #y-coord to the value stored in the AIP to get a net delay projected 
        #along this sidelobe's vector.
        SideLobeDelays += delayRamp*yVals;
        #N.B. "+" because El is +x rotation: +y is farther so, real delay > 0.
       
        #now compute the total beam gain. This is a complex sum across all the 
        #AIP phasors.
        w=2.*np.pi*c/self.wavelength;
        SideLobeComplexDelays=np.exp(1.j*w*np.matrix(SideLobeDelays,dtype=complex));
               
        complexGain = np.sum(SideLobeComplexDelays);
        gain = np.absolute(complexGain);
        phase = np.angle(complexGain);
        #Now we have the net delay per AIP, so sum over the Array
        #Calculate the 
        return(gain);
        

    
    def AzRotate(self,az):
        # provides the y coordinates of the array AIPs, in
        # a frame rotated to align -y with the Az direction
        
        # The rotation is az(radians), pos. rotation = CCW.
        # Sense of the rotation is +z, with x right, y up & z out of the page.
        
        #First, precalc the sin and cos:
        cosAz=np.cos(az);
        sinAz=np.sin(az);
            
        #I need a view or something to index through with one iterator. 
        #   Below,  .A1 method returns a collapsed Array view of Matrix.
        #   Later, ndarray.flat collapses an array
        Yarray=self.Ygrid.flat;  #Danger: no copy. ObjRef. (was A1)
        Xarray=self.Xgrid.flat; 
        
        #initialize space for the results for AIP y Values in Az frame:
        AIPyValsInAzFrame = np.ndarray([self.elements]);
        
        for i in range(0, self.elements) :
            AIPyValsInAzFrame[i] = -Xarray[i]*sinAz + Yarray[i]*cosAz;
            
        return(AIPyValsInAzFrame);    
       
        
    def IncrementAIPDelays(self,el,AIPYs):    
        
        c= 3e8;  #m/sec, speed of light
        
        #how much addtl delay must we add for every meter in -x
        delayRamp = np.sin(el)/c;  #something like this
    
        #Now go through the AIPs, adding delay proportionate to the y-coord
        for i in range(0,self.BeamDelays.size):
            self.BeamDelays[i] -= delayRamp*AIPYs[i];
            #N.B. "-" 'because'cause El is +x rotation: +y is farther so correction < 0.


    def SetAIPDelays(self,el,AIPYs):
        #Caution: you better have performed AzRotate first!
        #We are going to load the BeamDelays, as needed to create
        #the squint angle el, in the az frame defined by AzGrid.
        #(el is the "squint" angle away from aperture boresight)
        self.BeamDelays[:]=0.;
        self.IncrementAIPDelays(el,AIPYs);
        
    def aim(self,az,el):  #Loads the Delays for a commanded beam direction.
        #Transform all the array AIP grid locations to the Az frame 
        #Python Q: is self already local in AzRotate?
        self.BeamAzAngle=az;  # make these part of the class state.
        self.BeamElAngle=el;
        #RFF will this syntax work?
        AIPyValsInAzFrame=self.AzRotate(az);  #x-form to az frame.
        
        # Now, in this new frame, where we can...
        # apply the needed phase delays to all AIP to squint through el.
        self.SetAIPDelays(el,AIPyValsInAzFrame);
    

          

#------------------------------------------------------------------


def beamFormExample(az=0.5,el=15/57.3,Nulls=np.matrix([[0],[0]])):
       #compute the ramp delay in sec/meter, that will generate El
       B = beam();  # TODO I would like to pass in the AIP grid someday.
       
       B.aim(az,el);  #This calculates the delays needed to aim the beam.
    
       #Now compute and plot the beam power for all Az-El
    
       #Here, chosing the domain for the ArrayFactor calculation.
       elBins = 80; azBins=90;
       ElRange=np.linspace(-40.,41.,elBins)/57.3; 
       AzRange=np.linspace(-90.,91.,azBins)/57.3;
       ArrayFactor = np.zeros([azBins,elBins]);  #answer will go here.
       
       for i in range(0,AzRange.size):
           for j in range(0, ElRange.size):
       
               ArrayFactor[i,j]=B.calcGain(AzRange[i],ElRange[j]);
       
       fig = plt.figure();
       Z = np.absolute(ArrayFactor);      
       plt.plot(ArrayFactor);        
            
       fig = plt.figure()
       ax = fig.gca(projection='3d')
       X, Y = np.meshgrid(ElRange,AzRange)
       maxGain = np.max(Z);
      
       surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)
       ax.set_zlim(-1.01, maxGain);
       
       ax.zaxis.set_major_locator(LinearLocator(10))
       ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
       
       fig.colorbar(surf, shrink=0.5, aspect=5)
       
       plt.show()

    
       # for i in range (0,Nulls.shape[0])  #iterate through list of Nulls 
       i=0;  # leave iteration for later, hardcode to just one Null for now.
       
       #Next, we will add a null.
       NullAz=Nulls[0,i];  NullEl=Nulls[1,i];  #pull the values from the input.
       
       #For a first step, we will want to know the power radiated along the 
       #null direction, as a consequence of the beam forming already completed.
    
       
       
       #Transform the ApertureGrids to the Null Az coord frame.
       #Hey, are we going to just want a long vector of AIP coords?
       #(That is, I have not yet found any use for the 2-d grid...)
       #[XnullGrid,YnullGrid]=B.AzRotate(ApertureXgrid,ApertureYgrid,NullAz);
       
       # Seteup memory for the increments we want to create nulls, and the
       # time delays coming from the actual null El angle.
       #NullDelays=np.matrix(BeamDelays);  # this syntax does a deep copy. 
       #GeomDelays=np.matrix(BeamDelays);
       
       #Now, to each of the AIP array elements, add the perceived delay added
       #by the geometry pertaining to the Null's Elevation angle.
       #Remember, this fn will write to data @ NullDelays...
       #HEY  do I need RampSubtract, instead?  I'm not changing the delay in 
       #AIP memory, but instead adding the physical space delay that will happen
       #from geometry when we look into the array at this El angle.
       #GeomDelays=RampAdd(NullEl,Xgrid,Ygrid) + NullDelays;
       
       #Now that we have a null, compute the beam power for it.  It should match
       #one of the points in the plot made previously.  Now, however, we can fix it
       #by storing some changes.  Later we will solve mulitple nulls, but for now, 
       #just this one.
       #Say we have calculated the power: that is just the vector sum of all the 
       #GeomDelays (which function I haven't yet defined) and we want to kill it.
       
       #This next function calculates whether it will help to ootch the phase
       #(delay) of a particular AIP, and how much. The direction of change of 
       #these unit vectors is perpendicular to the vectors' direction, basically a
       #90deg rotation. Like the vector, the delay sensitivity is also 2-d, a 
       #dX/dPhase and dy/dPhase pair, and the effectivity of adjusting the phase
       #of an AIP will be the dot product of this sensitivity and the error vector.
       #DelaySensitivity=calcDelaySensitivity(GeomDelays);
       
    
    # here's the Main.  Takes Az and El of the commanded beam direction, and 
    # arrays of the null direction commands.
    # N.B. This is not a class method.
    # N.B. fn definition syntax that puts defaults into the parameter list.
    
ArrayFactor= beamFormExample();
        
        
