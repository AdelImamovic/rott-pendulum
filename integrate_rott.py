# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:42:07 2016

@author: adeli
"""




from utils_rott import *
import sys

workingdir='/home/adeli/Documents/Vorlesungen/tutorials/Numerische_Methoden_Umweltphysik/mat2py/rott/rottpy'
sys.path.append(workingdir)


if __name__=='__main__': 
    """This script demonstrates the use of the classes associated with the Rott pendulum.
    For class and function documentation see the ./rott_documentation.html file
    """
    #set up default pendulum and show the parameters 
    s1=setup(adot0=5)
    s1.print_params()      

    #initialize pendulum and show the parameters
    pp1=pendulum_properties(s1)    
    pp1.print_prop()
    
    #prepare the trajectory in phasespace (set of points (alpha,gamma,alphadot,gammadot))
    t1=trajectory_in_phasespace(s1,pp1)
    ischeme=0
    
    #integrate the pendulum with a default scheme   
    t1.integrate(int_scheme=0)
   
    #set up a new pendulum with changed initial conditions
    s2=setup(a0=60)
    pp2=pendulum_properties(s2)
    t2=trajectory_in_phasespace(s2,pp2)    
    t2.integrate()
    
    #plot and save the output   
    t1.plot(savefig=True,figname='t1')
    t2.plot(savefig=True,figname='t2')
  
    #overloading of +operator in trajectory_in_phasespace class results 
    #t1+t2
    
    #animate the pendulum and save the output as test.mp4
    videoname='animation_test'
    t1.animate_pendulum(videoname)

    #t1.energetics(plotme=False)
    
    
