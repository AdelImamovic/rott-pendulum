# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 14:42:31 2016

@author: adeli
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
import re
pi=np.pi
# maker of figures

class setup:
    _geomparams={
    'M':'documentation',
    'R':'documentation'}
    
    _initvalues={
    'a0':'documentation',
    'g0':'documentation',
    'adot':'documentation',
    'gdot0':'documentation'}
    
    _integrationparams={
    'ischeme':'documentation',
    'time':'documentation',
    'delt':'documentation'}
    
    def_time=60.
    def_M=1.
    def_R=2.0638
    def_a0=0.
    def_g0=0.
    def_adot0=15.
    def_gdot0=0.
    def_ischeme=1.
    def_delt=0.005

    def __init__(self,time=def_time,M=1.,R=2.0638,a0=0.,g0=15.,adot0=0,gdot0=0,ischeme=6,delt=0.005):
        self.time=time        
        self.M=M
        self.R=R
        self.a0=np.deg2rad(a0)
        self.g0=np.deg2rad(g0)
        self.adot0=adot0
        self.gdot0=gdot0
        self.initvec=(self.a0,self.g0,self.adot0,self.gdot0)
        self.ischeme=ischeme
        self.delt=delt
        
        
    def print_params(self):
        print 'M,R'
        print (self.M,self.R)
        print 'a0,g0,adot0,gdot0'
        print (self.a0,self.g0,self.adot0,self.gdot0)
        print 'ischeme,time,delt'
        print (self.ischeme,self.time,self.delt)
        
        
        
class pendulum_properties:
    def __init__(self,setup,ini_prop=6,ndim=1):
        self.PendulumType=ini_prop
        self.ndim=1
        if ini_prop==6:
            self.M=setup.M
            if setup.M<0.:
                print('WARNING: you chose M<0 (invalid). M=1.')
                self.M=1.
            self.R=setup.R
            self.g=9.81
            
            self.rho=1.
            self.rhoL=1.
            self.rhoR=1.
            self.b0=0.2
            self.hL=self.b0*1.5
            self.hR=self.b0*1.5 #!Redundant
            
            self.hR=np.sqrt((self.rhoL/self.rhoR)*(self.hL/self.M)**2)
            if self.R>0.:
                self.find6_h()
                #nasty: self.hL,self.hR=find6_h(self.R,self.M,self.rhoL,self.rho,g=g,self.b0,self.hL,self.hR)
            #self.aux=design_rott6(self.rhoL,self.rho,self.rhoR,g=g,self.b0,self.hL,self.hR) #physical auxillarly properties
            #self.mc,self.mt,self.s,self.a,self.b,self.bL,self.c,self.Ic,self.It=self.aux
            self.beta0=np.pi/2.
            self.design_rott6()
            
            
            #unused parameters?
            self.d=0
            self.mL=0.
            self.mR=0.            
            self.t=0.          
            
            #pendulum parameters 
            self.lambda2=self.g*self.mt*self.a/self.It
            self.omega2=self.g*self.mc*self.c/self.Ic
            self.eta=self.mc*self.b*self.c/self.It
            self.zeta=self.mc*self.b*self.c/self.Ic
            
            self.R=np.sqrt(self.omega2/self.lambda2)
            self.M=np.sqrt(self.a*self.mt/(self.c*self.mc))
            self.Q=self.c*self.mc*self.b/np.sqrt(self.Ic*self.It)
            
            if self.ndim==1:
                self.tau2=4*(pi**2)/self.omega2
                self.tau=np.sqrt(self.tau2)
                self.g=self.tau2*self.g
                self.lambda2=self.tau2*self.lambda2
                self.omega2=self.tau2*self.omega2
                
           
        else:
            print 'Pendulum Type not defined -- > Assert'
            assert 0
            
    def bl2o2ez(self):
        return (self.beta0,self.lambda2,self.omega2,self.eta,self.zeta)
         
    def print_prop(self):
        print 20*'-'
        print 'moment of inertia: It, Ic, mt, mc'
        print (self.It,self.Ic,self.mt,self.mc)
        print 20*'-'
        print 'Amplitude-ratio of Rott oscillations: a*mt, c*mc'
        print (self.a*self.mt,self.c*self.mc)
        print 20*'-'
        print 'pendulum parameters: lambda^2, omega^2, eta, zeta'
        print (self.lambda2,self.omega2,self.eta,self.zeta)
        print 20*'-'  
        print 'dimensionless parameters: R, M, Q'
        print (self.R,self.M,self.Q)
        print 20*'-'   
        print 'oscillations periods [s],ratio'
        print (2*pi/np.sqrt(self.lambda2),2*pi/np.sqrt(self.omega2),np.sqrt(self.omega2/self.lambda2))
        #print 20*'-' 
        #print 'oscillations periods [1],ratio'
        #print (2*pi/np.sqrt(self.lambda2),2*pi/np.sqrt(self.omega2),np.sqrt(self.omega2/self.lambda2))
        #print 20*'-' 
     
     #this function designs a physical pendulum (pendulum type 5) from geometrical 
     #and material parameters. The mass of the pendulum is concentrated in a
     #line with densities rhoL, rho and rhoR [kg/m], respectively.      
    def design_rott6(self):
        self.mc=self.rhoR*self.hR
        self.mt=self.rhoL*self.hL + self.rho*(2*self.b0) + self.mc     
        #compute center of gravity and geometry
        self.s=(-self.rhoL*self.hL*self.b0+self.rhoR*self.hR*self.b0)/self.mt
        self.b =self.b0-self.s
        self.bL=self.b0+self.s
        self.a=self.rhoL*self.hL*(self.hL/2) / self.mt
        self.c=self.hR/2
        #compute momenta of inertia
        self.It= (self.rhoL*(self.hL)**3)/3+self.rhoL*self.hL*(self.bL**2)+ (self.rho*(self.bL**3))/3 + (self.rho*(self.b**3))/3. + self.rhoR*self.hR*(self.b)**2
        self.Ic= (self.rhoR*(self.hR)**3)/3.
     
    #this function designs a physical pendulum (pendulum type 5) from geometrical 
    #and material parameters such that it has a pre-defined ratio R between the periods of
    #the large and small pendulas 
    def find6_h(self):
        pass
#        self.hrat=self.hR/self.hL
#        self.hL=self.fzero(self.funct,self.hL0)
#        if self.M>0:
#            self.HR=np.sqrt((self.rhoL/self.rhoR)*(self.hL/self.M)**2)
#        else:
#            self.hhR=self.hrat*self.hL
            

        
        
        
  
        
        
        
        
     
     
def rott_dynamics(state,tr,beta0,lambda2,omega2,eta,zeta): #state in phasespace, timerange pend_properties
    a,g,a1,g1=state
    #beta0,lambda2,omega2,eta,zeta=bl2o2ez    
#    a0=svec.a #zeroth derivative of alpha
#    g0=svec.g
#    a1=svec.adot #first derivative of alpha
#    g1=svec.gdot
    
    #aux variables
    C=np.cos(g-a-beta0)
    S=-np.sin(g-a-beta0)
    sa=np.sin(a)
    sg=np.sin(g)
    
    #second derivatives of alpha and gamma
    a2=(-lambda2*sa+eta*C*omega2*sg-eta*(zeta*C*S*a1**2+S*g1**2))/(1-eta*zeta*C**2)
    g2=(-omega2*sg+zeta*C*lambda2*sa+zeta*(eta*C*S*g1**2+S*a1**2))/(1-eta*zeta*C**2)
    
    
    return (a1,g1,a2,g2)
    




#type for point in phase space 
class state:
    def __init__(self,a,g,adot,gdot):
        self.a=a
        self.g=g
        self.adot=adot
        self.gdot=gdot
        self.vec=(a,g,adot,gdot)
    
        
        
#class containing all the output after calculations     
class trajectory_in_phasespace:
    def __init__(self,s0,pp):
        self.timeinsts=np.arange(0.,s0.time,s0.delt)
        self.setup=s0
        self.time=s0.time
        self.pp=pp
        #self.init_state=state_of_pendulum(s0.a0,s0.g0,s0.adot0,s0.gdot0)
        #self.p_states=[self.init_state]    
        self.a=np.zeros(len(self.timeinsts))
        self.g=np.copy(self.a)
        self.adot=np.copy(self.a)
        self.gdot=np.copy(self.a)

    def __setitem__(self,i,newstate):
        #if not (t in self.timeinsts):
        #    print 'time not in timeinsts --> interpolation to grid'
            
        #i=np.where(self.timeinsts==t)[0][0]
            
        (self.a)[i]=newstate.a
        (self.g)[i]=newstate.g
        (self.adot)[i]=newstate.adot
        (self.gdot)[i]=newstate.gdot

    def __getitem__(self,i):
        #if not (t in self.timeinsts):
        #    print 'time not in timeinsts --> interpolation to grid'
        #    print t
        #    print self.timeinsts
        #i=np.where(self.timeinsts==t)[0][0]
        return state(self.a[i],self.g[i],self.adot[i],self.gdot[i])
        
    def __add__(self,traj2):
        figname='addition'
        self.plot(figname=figname,savefig=False,holdagain=True)
        traj2.plot(figname=figname,savefig=True,holdagain=False)
        
    def integrate(self,int_scheme=0):
        initstate=(self.setup).initvec
        params=(self.pp).bl2o2ez()
        output=odeint(rott_dynamics,initstate,self.timeinsts,args=params)
        output=np.transpose(output)        
        self.a,self.g,self.adot,self.gdot=output
        
    def energetics(self, plotme=False):
        a=self.a
        g=self.g
        adot=self.adot
        gdot=self.gdot
        
        #pendprop
        beta0=(self.pp).beta0
        It=(self.pp).It
        Ic=(self.pp).Ic
        mc=(self.pp).mc
        b=(self.pp).b
        c=(self.pp).c
        gpp=(self.pp).g
        app=(self.pp).a
        mt=(self.pp).mt     
        C=np.cos(g-a-beta0)
        
        #energies
        self.ekin=0.5*It*adot**2+0.5*Ic*gdot**2+mc*b*c*C*adot*gdot
        self.epot=-gpp*(mt*app*np.cos(a)+mc*c*np.cos(g))      
        self.etot=self.ekin+self.epot
        #measure for dissipation
        self.dE=max(self.etot[0]-np.min(self.etot),np.max(self.etot)-self.etot[0])
        
        if plotme:
            plt.plot(self.timeinsts,self.etot)
            plt.savefig('energetics.pdf')
            print 'the numerical energy dissipation is '+str(self.dE)
        #self.ekin=0.5*setup.It*self.adot**2+0.5*setup
            
    def plot(self,figname,var='adotgdot',savefig=True,holdagain=False):
        if not (var=='adotgdot'):
            print 'no valid variable to plot'
        plt.plot(self.timeinsts,self.adot,linewidth=0.5)
        plt.hold(True)
        plt.plot(self.timeinsts,self.gdot,linewidth=0.5)
        
        if not holdagain:
            #plt.legend()
            plt.savefig(figname+str(self.time)+'.pdf',dpi=300)
            
        plt.hold(holdagain)
        
    def animate(self):#,framerate=20*(self.setup).dt):
        #ntimesteps=len(self.a)
        for alp,gam in zip(self.a,self.g):
            self.snapshot(alp=alp,gam=gam)
            
        #check if trajectory exists
        #compile the pendulum into snapshots
        #
    
    def er(angle):
        """Returns the vector that points to point on unit circle."""
        return np.array([np.cos(angle),np.sin(angle)])
    
    def snapshot(alp=0,gam=0,timeinst=10):#self,time):
        M=1
        R=1
        #alp=np.deg2rad(alp)
        #gam=np.deg2rad(gam)
        Lcolor='green'
        Icolor='red' 
        Llwidth=3
        Ilwidth=3
        pivot=np.zeros(2)
        
        #geometric helper
        ea=er(alp)
        eg=er(gam)
        ea_p=er(alp-np.pi/2)
        eg_p=er(gam-np.pi/2)
        
        l1=l2=l3=l4=1. #for testing                
        
        #five pivot points
        p1=pivot        
        p2=p1+ea*l1
        p3=p1-ea*l2
        p4=p2+eg_p*l3
        p5=p3+ea_p*l4
        allpivots=[p1,p2]#,p3,p4,p5]
        
        fig=plt.figure()        
        plt.hold(True)
        
        plt.axes().set_aspect('equal')
        #L shape:
        #p2->p3
        t=zip(p2,p3)
        plt.plot(t[0],t[1],color=Lcolor,linewidth=Llwidth)
        #p3->p5
        t=zip(p3,p5)
        plt.plot(t[0],t[1],color=Lcolor,linewidth=Llwidth)
        #plt.plot(,color=Lcolor)
        #I shape:
        #p2->p4
        t=zip(p2,p4)
        plt.plot(t[0],t[1],color=Icolor,linewidth=Ilwidth)
        
        #splt.plot(,color=Icolor)
        for piv in allpivots:
            circle=plt.Circle(piv,l1*0.05,color='black')
            fig.gca().add_artist(circle)
        
        #cosmetics        
        fact=2.5
        plt.xlim([pivot[0]-fact*l1,pivot[1]+fact*l1])
        plt.ylim([pivot[0]-fact*l3,pivot[1]+fact*l3])
        plt.xticks([],[])
        plt.yticks([],[])
        plt.grid()        
        plt.show()
        
       
        #return mpllibobject
        #creates a snapshot in $Format of the current pendulum state
           
   
#0      
def dummy_integrator(oldstate): #dummy integrator
    newstate=oldstate #do the integration here
    return newstate
    
def ode2cntr(old):
    pass

def ode3abs(old):
    pass

def ode3rk(old):
    pass

def ode6rk(old):
    pass

def ode8rk_hairer(old):
    pass

def ode8rk_scalar(old):
    pass

    
        
