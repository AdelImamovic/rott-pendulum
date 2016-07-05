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
    #check for input size
    #write help function

    #set initial conditions, integration properties, 
    s1=setup()
    s1.print_params()      

    #initialize pendulum geometry / properties
    pp1=pendulum_properties(s1)    
    pp1.print_prop()
    
    #prepare objects for calculations:
    t1=trajectory_in_phasespace(s1,pp1)
    ischeme=0
    
    #available integratos
    
    t1.integrate(int_scheme=0)
    
    s2=setup(adot0=s1.adot0+1e-4)
    pp2=pendulum_properties(s2)
    t2=trajectory_in_phasespace(s2,pp2)    
    t2.integrate()
    #t2.integrate(s1,pp1)
    
    
    t1.plot(savefig=True,figname='t1')
    t2.plot(savefig=True,figname='t2')
  
    #t1+t2 #overload operator+
    
    
    #t1.energetics(plotme=False)
    
    
