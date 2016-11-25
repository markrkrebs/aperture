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



                         
wavelength = 2.5;  #meters, wavelength we're working with
                         
class beam(object):
    def __init__(self):
        # here's the aperture AIP locatoins, in m.
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
        self.AzYgrid = np.zeros(self.gridshape);
        self.AzXgrid = np.zeros(self.gridshape);  
               
        #This initial condition sends a beam out the boresight.  Remember, the
        #zeros are the phase angles of the unit power vectors at each AIP.
        self.BeamDelays= np.zeros(self.gridshape);
        self.elements = self.Ygrid.size;  #number of elements
        
    def AzRotate(self,az):
        # provides the x and y positions of the array AIPs, in
        # coordinate frame rotated to align with the Az direction
        # First, allocate space for the answer: that's a class data element
        # so, as it's done, we can immediately load it.
        cosAz=np.cos(az);
        sinAz=np.sin(az);
        # Sense of this rotation is CCW, with x right & y up, on the page
        # So the rotation is +z, out of the page, towards you.
        # Python question: can I assume the "self." eg on self.Xgrid?
        Xarray=Xgrid.reshape(self.elements); #Danger: no copy. ObjRef.
        Yarray=Ygrid.reshape(self.elements);
        for i in range(0, self.elements) :
            self.AzYgrid[i] = -Xarray[i]*sinAz + Yarray[i]*cosAz;
            self.AzXgrid[i] = Xarray[i]*cosAz + Yarray[i]*sinAz;
   
    def RampAdd(el):
        #Caution: you better have performed AzRotate first!
        #We are going to load the BeamDelays, as needed to create
        #the squint angle el, in the az frame defined by AzGrid.
        #(el is the "squint" angle away from aperture boresight)
        
        c= 3e8;  #m/sec, speed of light
        
       #how much addtl delay must we add for every meter in -x
       delayRamp = sin(el)/c;  #something like this
    
        #First, make a temp linear array pointer into BeamDelays
        BeamDelayArray=self.BeamDelays.reshape(self.elements)
        #Now go through the AIPs, adding delay proportionate to the x-coord
        for i in range(0,self.elements):
            BeamDelayArray[i]= delayRamp*AzYgrid[i];

        
    def aim(self,az,el):  #Loads the Delays for a commanded beam direction.
        #Transform all the array AIP grid locations to the Az frame 
        #Python Q: is self already local in AzRotate?
        self.BeamAzAngle=az;  # make these part of the class state.
        self.BeamElAngle=el;
        self.AzRotate(az);  #x-form to az frame, where we can...
        # apply the needed phase delays to all AIP to squint through el.
        # (This function will set the BeamDelays)
        self.RampAdd(el);
                   

#------------------------------------------------------------------
# here's the Main.  Takes Az and El of the commanded beam direction, and 
# arrays of the null direction commands.
# N.B. This is not a class method.
# N.B. fn definition syntax that puts defaults into the parameter list.

def beamformExample(az=0.,el=30/57.3,Nulls=np.matrix([[0],[0]])):
   #compute the ramp delay in sec/meter, that will generate El
   B = beam();  # TODO I would like to pass in the AIP grid someday.
   
   B.aim(az,el);  #This calculates the delays needed to aim the beam
   
  
   #Now that we know the delays, compute and plot the beam power for all Az-El
   #Also note a basic metric, the max ratio of all sidelobes to main beam pwr.
   #TODO: do that.
  
  

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
   
   
   

    
    
#==============================================================================
#    Here are comments motivating this program.
#    
#    For n beams…  This first note is just reminder that we will do this 
#     several times, completely independently, and later sum them at the AIP. 
#     
# Define the beam as a phase ramp, applied across the aperture in a particular 
# azimuthal direction.  I’m tempted to call these angles Az and Squint, since the
# second one is deviation from mechanical boresight. Call Az, El however, as 
# that is a little more familiar, though this is kind of the complement of 
# canonical elevation.  Imagine starting in beam co-ords with 100% vertically 
# aligned, zero-angle phasors.
# 
# Transform to antenna coordinates by applying the phase ramp, in the Az 
# direction. (Can I be clearer about this? Given an AIP coordinate, rotate 
# through Az halfway to beam frame. In this intermediate frame, the AIP x-coord 
# is perpendicular to LOS, and y-coordinate needs phase delay proportionate to 
# El command and the value of y.)
# 
# 
# Next transform to null LOS coordinates, where the newly assigned beam delay 
# adds to another geometric delay instantiated by the null LOS. Here we hope for 
# a random looking set, which aggregates to zero. Whatever is left is sidelobe, 
# which must be nulled.  The value of that phase sum is negated and “applied” 
# to only those AIPs whose phasor, viewed in the null frame, is substantially 
# perpendicular to the error.  This will demand the least total modifications 
# to the ideal beam.
# 
# 
# Repeat this process for all the other nulls.
# 
# 
# Consider an iterative procedure, where the main beam is viewed again and 
# moderately corrected?  Every error from zero phase will have been applied in 
# pursuit of achieving a null so “repairing” the main beam would be both 
# trivial, and catastrophic to the null formation.  However, it may be 
# worthwhile to sweep through the nulls again making another round of even 
# smaller adjustments to converge (?) the AIP phasors to values mutually 
# satisfactory to all the null beams.
# 
# 
# I think all the null beam settings can be computed at once with a pseudo-
# inverse.  It’s trivial to compute d(phasorSum)/d(AIP phase) for each AIP, 
# which represents a vector of weights. Assembled, these make a matrix as wide 
# as element count and tall as # of nulls. In [A]x=b notation, where b is the 
# command vector, solve for x using pseudoinverse where b is -1 times the 
# residual phase in each null.  The command vector could be chosen to achieve 
# sufficient reduction (eg -40dB, vs the main beam) without requiring perfect 
# nulling.  Ideally the inverse would allow inequalities like <= in this 
# calculation. Hopefully in the presumably rare case where a deep null is found 
# in a desired direction, it will not be really bad to make the tiny corrections 
# required to “lift” it up to the specified attenuation.
#==============================================================================
