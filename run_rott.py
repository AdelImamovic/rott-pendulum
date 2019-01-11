# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:42:07 2016

@author: adeli
"""




from utils_rott import *
import sys


FIPATH='./Figures/'


if __name__=='__main__': 
    """This script demonstrates the use of the classes associated with the Rott pendulum.
    For class and function documentation see the ./rott_documentation.html file
    """
    #set up default pendulum and show the parameters 
    s1=RottSetup()
    s1.print_params()      

    #initialize pendulum and show the parameters
    pp1=RottProperties(s1)    
    pp1.print_prop()
    
    #prepare the trajectory in phasespace (set of points (alpha,gamma,alphadot,gammadot))
    t1=RottTrajectoryInPhasespace(s1,pp1)
    ischeme=0
    
    #integrate the pendulum with a default scheme   
    t1.integrate(int_scheme=0)
   
    #set up a new pendulum with changed initial conditions
    s2=RottSetup(g0=90.)
    pp2=RottProperties(s2)
    t2=RottTrajectoryInPhasespace(s2,pp2)    
    t2.integrate()
    
    #plot and save the output   
    t1.plot(savefig=True,figname='t1')
    #t2.plot(savefig=True,figname='t2')
  
    #overloading of +operator in trajectory_in_phasespace class results 
    #t1+t2
    
    #animate the pendulum and save the output as test.mp4
    videoname='animation_test'
    t2.animate_rott(videoname)

    #t1.energetics(plotme=False)
    
    
