# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 19:44:52 2016

@author: mark
"""
#As of 10:40AM, 25Nov'16, I pushed this into github.

import numpy as np
# This tool is intended to plot a beam and a few nulls,
#arrayed across a 2-d circular grid 
#maybe later with "space taper"
#applying a gradient search in direction determined by SVD of the 
#"Jacobian" that is, the partial derivatives of nulls w.r.t. each of the 
#AIP phasors.



def beamFormExample(az=0.,el=30/57.3,Nulls=np.matrix([[0],[0]])):
   #compute the ramp delay in sec/meter, that will generate El
   B = beam();  # TODO I would like to pass in the AIP grid someday.
   
   B.aim(az,el);  #This calculates the delays needed to aim the beam.
   
  
   #Now that we know the delays, compute and plot the beam power for all Az-El
   #Also note a basic metric, the max ratio of all sidelobes to main beam pwr.
   #TODO: do that.
  
   #To calculate the ArrayFactor (an array itself, of  expressing the gain of 
   #of the antenna array, at different angles) we need to define the domain
  
   #Somewhat arbitrarily, we'll make ArrayFactor the same shape.
   ArrayFactor = np.zeros(B.shape); 
   ArrayFactorAz = np.zeros(B.shape);
   ArrayFactorEl = np.zeros(B.shape);
   ArrayFactorArray=ArrayFactor.reshape(B.elements); #Danger: no copy. ObjRef.
   for i in range(0,B.elements):
       AzAFA[i]=B.AzAngle(i);
       
  
  

   # for i in range (0,Nulls.shape[0])  #iterate through list of Nulls 
   i=0;  # leave iteration for later, hardcode to just one Null for now.
   
   #Next, we will add a null.
   #For a first step, we will want to know the power radiated along the 
   #null direction, as a consequence of the beam forming already completed.
   #(We could get lucky and be in a null already!  ...but, it's not likely.)
   
   NullAz=Nulls[0,i];  NullEl=Nulls[1,i];  #pull the values from the input.
   
   #Transform the ApertureGrids to the Null Az coord frame.
   #Hey, are we going to just want a long vector of AIP coords?
   #(That is, I have not yet found any use for the 2-d grid...)
   [XnullGrid,YnullGrid]=AzRotate(ApertureXgrid,ApertureYgrid,NullAz);
   
   # Seteup memory for the increments we want to create nulls, and the
   # time delays coming from the actual null El angle.
   NullDelays=np.matrix(BeamDelays);  # this syntax does a deep copy. 
   GeomDelays=np.matrix(BeamDelays);
   
   #Now, to each of the AIP array elements, add the perceived delay added
   #by the geometry pertaining to the Null's Elevation angle.
   #Remember, this fn will write to data @ NullDelays...
   #HEY  do I need RampSubtract, instead?  I'm not changing the delay in 
   #AIP memory, but instead adding the physical space delay that will happen
   #from geometry when we look into the array at this El angle.
   GeomDelays=RampAdd(NullEl,Xgrid,Ygrid) + NullDelays;
   
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
   DelaySensitivity=calcDelaySensitivity(GeomDelays);
   
   
#-----------Supporting class-------------------------------------------

                         
wavelength = 2.5;  #meters, wavelength we're working with


                         
class beam(object):
    def __init__(self):
        # here's the aperture AIP locations, in m.
        # Here are some coordinates to start with, for a 7x7 array.
        # Note the center element is at Xgrid=0, Ygrid=0.
        # Later, I should pass these in, instead.
        self.Ygrid =np.matrix([[-3.,-3.,-3.,-3.,-3.,-3.,-3.],
                         [-2.,-2.,-2.,-2.,-2.,-2.,-2.],
                         [-1.,-1.,-1.,-1.,-1.,-1.,-1.],
                         [0.,0.,0.,0.,0.,0.,0.],
                         [1.,1.,1.,1.,1.,1.,1.],
                         [2.,2.,2.,2.,2.,2.,2.],
                         [3.,3.,3.,3.,3.,3.,3.]]);
        self.Xgrid= np.matrix([[-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.],
                         [-3.,-2.,-1.,0.,1.,2.,3.]]);
        
        self.gridshape = self.Ygrid.shape; 
        self.elements = self.Ygrid.size;  #number of elements
        
        #These next variables are locations transformed into Az coords.
        self.AzYgrid = np.zeros(self.gridshape);
        self.AzYarray =AzYgrid.reshape(elements);
        
        self.AzXgrid = np.zeros(self.gridshape); 
        self.AzXarray= AzXgrid.reshape(elements);
               
        #This initial condition sends a beam out the boresight.  Remember, the
        #zeros are the phase angles of the unit power vectors at each AIP.
        self.BeamDelays= np.zeros(self.gridshape);
        
    def AzRotate(self,az):
        # provides the x and y positions of the array AIPs, in
        # coordinate frame rotated to align with the Az direction
        # First, allocate space for the answer: that's a class data element
        # so, as it's done, we can immediately load it.
        cosAz=np.cos(az);
        sinAz=np.sin(az);
        # Sense of this rotation is CCW, with x right & y up, on the page
        # So the rotation is +z, out of the page, towards you.
        # Python question: can I assume the "self." eg on self.Xgrid? A:no.
        
        #I need a view or something to index through with one iterator. In 
        #the below, matrix.reshape does not remove a dimension, but .A1 method 
        #returns a collapsed Array view. (view means objref: not a copy)
        Xarray=self.Xgrid.A1;       #Danger: no copy. ObjRef.
        Yarray=self.Ygrid.A1;
        for i in range(0, self.elements) :
            self.AzYgrid.A1[i] = -Xarray[i]*sinAz + Yarray[i]*cosAz;
            self.AzXgrid.A1[i] = Xarray[i]*cosAz + Yarray[i]*sinAz;
   
    def RampAdd(el):
        #Caution: you better have performed AzRotate first!
        #We are going to load the BeamDelays, as needed to create
        #the squint angle el, in the az frame defined by AzGrid.
        #(el is the "squint" angle away from aperture boresight)
        
        c= 3e8;  #m/sec, speed of light
        
        #how much addtl delay must we add for every meter in -x
        delayRamp = sin(el)/c;  #something like this
    
        #First, make a temp linear array pointer into BeamDelays:
        BeamDelayArray=self.BeamDelays.reshape(self.elements)
        #Now go through the AIPs, adding delay proportionate to the x-coord
        for i in range(0,self.elements):
            #N.B. BeamDelayArray == BeamDelays so we are loading "both."
            BeamDelayArray[i]= -delayRamp*AzYgrid[i];
            #N.B. "-" because El is +x rotation: +y is farther so delay < 0.

        
    def aim(self,az,el):  #Loads the Delays for a commanded beam direction.
        #Transform all the array AIP grid locations to the Az frame 
        #Python Q: is self already local in AzRotate?
        self.BeamAzAngle=az;  # make these part of the class state.
        self.BeamElAngle=el;
        self.AzRotate(az);  #x-form to az frame, where we can...
        # apply the needed phase delays to all AIP to squint through el.
        # (This function will set the BeamDelays)
        self.RampAdd(el);
        
    def AzAngle(self,whichElement):
        Xarray=self.Xgrid.reshape(self.elements); #Danger: no copy. ObjRef.
        Yarray=self.Ygrid.reshape(self.elements);
        angle = atan2(Yarray[i],Xarray[i]);
        return(angle);
                   

#------------------------------------------------------------------


# here's the Main.  Takes Az and El of the commanded beam direction, and 
# arrays of the null direction commands.
# N.B. This is not a class method.
# N.B. fn definition syntax that puts defaults into the parameter list.

#beamFormExample();
    
    
