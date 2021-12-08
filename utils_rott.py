# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 14:42:31 2016

@author: Adel Imamovic (adeli@ethz.ch)
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
from matplotlib import animation
# importing movie py libraries
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
 
pi=np.pi


class RottSetup:
    """The RottSetup class contains information on the initial conditions,
    numerical integration paramters, and pendulum geometries. 
    It is the Front-End for students of the ETH lecture (Christoph Schär)
    "Numerische Methoden der Umweltphysik"
    
    Usage:
    -----
        s1=RottSetup();
    
    Parameters:
    ----------
        
    TODO:
    ----
        DONE
    """
    
    def_time=60.
    def_M=1.
    def_R=2.0638
    def_a0=0.
    def_g0=0.
    def_adot0=15.
    def_gdot0=0.
    def_ischeme=1.
    def_delt=0.005
    def __init__(self,a0=0.,g0=90.,adot0=0.,gdot0=0.,
                 time=60.,delt=0.005,ischeme=6,
                 M=1.,R=2.0638):
        """Initilializes the Rott Pendulum (L+I) with default arguments.
        Converts initial angles/angular speeds a0,g0,adot0,gdot0 to radians.
        
        Keyword Arguments:
        -----------------
        Initial Conditions
            a0 :    float
                initial deflection angle of L part in degree
            g0 :    float
                initial deflection angle of I part in degree per second
            adot0:  float
                initial angular speed of L part in degree
            gdot0:  float
                initial angular speed of I part in degree per second
        
        Integration Parameters:
            time:   float
                total integration time in seconds
            delt:   float
                time step size in seconds
            ischeme: int
                integration scheme index
                
        Pendulum Parameters:
            M: float
                geometric parameter
            R: float
                geometric parameter   
        """   
        #initial state
        self.a0=np.deg2rad(a0)
        self.g0=np.deg2rad(g0)
        self.adot0=np.deg2rad(adot0)
        self.gdot0=np.deg2rad(gdot0)
        self.initvec=(self.a0,self.g0,self.adot0,self.gdot0)        
        
        #integration parameters
        self.time=time  
        self.delt=delt  
        self.ischeme=ischeme
        
        #pendulum parameters
        self.M=M
        self.R=R
        

class RottProperties:
    """ Contains the geometric information of the Rott Pendulum.
    Calculates further quantities for the integration.
    
    Parameters:
    ----------
        pend_type:
        
        pend_type=1 Standard resonant Rott pendulum. In this case you need 
        to specify very little information about the pendulum.
        
        pend_type=2 A generalized Rott pendulum, which is composed by 
        three arbitrary mass lines
       
        pend_type=3 Fully general double pendulum. In this case one needs 
        to specify for both pendulum elements, the momenta of inertia, 
        the mass and the center of gravity.
        
        pend_type=4 A physical pendulum made of metal plates. Given 
        some initial geometrical information, the initialization solves 
        for a pendulum with a given frequency ratio R.
        
        pend_type=5 A physical pendulum made of metal plates 
        and rounded corners. Given some initial geometrical information, 
        the initialization solves for a pendulum with a given frequency ratio R.
        
        pend_type=6 almost like a Rott pendulum composed of three mass lines, 
        full generality in definition, but b is adjusted such that beta0=pi/2     
    
    TODO:
    ----
        So far only pend_prop==6 is implemented
    """
    
    #some hard coded constants
    pend_geoms=range(1,7)
   
    def __init__(self,setup,pend_type=6,make_ndim=1):
        """Initializes plot properties.

        Arguments:
        ---------
            setup:      RottSetup
                        Init properties for Pendulum
            
            pend_type:  int in 1-6
                        switch for selecting different pendulum geometries
                        
            make_ndim:  bool
                        shall we make the pendulum non-dimensional?                       
        """      
        
        if pend_type not in [6]:
            print("Pendulum geometry is not implemented.")
            raise NotImplemented
        self.setup=setup    
        self.pend_type=pend_type
        self.make_ndim=make_ndim
        #calcualte the properties once pendulum type is known
        self.calculate_properties()
        self.init_general_params()
         
    def print_prop(self):
        """Prints all parameters used for the integration to stdout.
        
        TODO:
            Nicer and more readable print layout ifneedbe
            Add decorators for pend_type specific properties (instead of if/elseif)


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
        if self.pend_type==6:
            print 20*'-'
            print 'geometry: hL, bL, b, hR, a, c, beta0[deg], g'
            print [self.hL,self.bL,self.b,self.hR,self.a,self.c,self.beta0*180./pi,self.g]								
            print 'density: rhoL, rho, rhoR'								
            print [self.rhoL,self.rho,self.rhoR]        
        """
    def calculate_properties(self):
        """Decides at runtime which calculate_properties method to call.
        
        Use of dispatch method        
        """
        calc_method_name='calc_pendtype_'+str(self.pend_type)
        calc_method=getattr(self,calc_method_name,lambda:'nothing')
        calc_method()   
     
    def calc_pendtype_6(self):
        """Calculates properties of pendulum of pend_type 6        
        Is called by calculate_properties if self.pend_type is 6.
        
        Calculates pendulum properties.        
        """
        setup=self.setup
        self.M=setup.M
        if setup.M<0.:
            print('WARNING: you chose M<0 (invalid). M=1.')
            raise RuntimeWarning
            self.M=1.
        self.R=setup.R
        self.g=9.81 # g is a property of the instance, because it can be rescaled
        
        self.rho=1.
        self.rhoL=1.
        self.rhoR=1.
        self.b0=0.2
        self.hL=self.b0*1.5
        self.hR=self.b0*1.5 #!Redundant
        
        self.hR=np.sqrt((self.rhoL/self.rhoR)*(self.hL/self.M)**2)
        if self.R>0.:
            self.find6_h()
        
        self.design_rott6()
        self.beta0=np.pi/2.      
        
        #unused parameters?
        self.d=0.
        self.mL=0.
        self.mR=0.            
        self.t=0.
    
    def init_general_params(self):
        #general pendulum parameters 
        self.lambda2=self.g*self.mt*self.a/self.It
        self.omega2=self.g*self.mc*self.c/self.Ic
        self.eta=self.mc*self.b*self.c/self.It
        self.zeta=self.mc*self.b*self.c/self.Ic
        
        self.R=np.sqrt(self.omega2/self.lambda2)
        self.M=np.sqrt(self.a*self.mt/(self.c*self.mc))
        self.Q=self.c*self.mc*self.b/np.sqrt(self.Ic*self.It)
        
        if self.make_ndim==1:
            self.tau2=4*(pi**2)/self.omega2
            self.tau=np.sqrt(self.tau2)
            self.g=self.tau2*self.g
            self.lambda2=self.tau2*self.lambda2
            self.omega2=self.tau2*self.omega2
            
            
    def calc_pendtype_5(self):
        #function to be called if pend_prop==5
        pass
    def calc_pendtype_4(self):
        pass
    def calc_pendtype_3(self):
        pass
    def calc_pendtype_2(self):
        pass
    def calc_pendtype_1(self):
        pass
     
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
        
        
    def bl2o2ez(self):
        return (self.beta0,self.lambda2,self.omega2,self.eta,self.zeta)
     
    def find6_h(self):
        pass
        """Designs parameters for a physical pendulum (pendulum type 6) 
        from geometrical and material parameters such that it has a pre-defined 
        ratio R between the periods of the large and small pendulas.
        
        POST:
        ----
            changes hR,hL values of instance only
        
        
        TODO:
        ----
            design_rott6 call is unfortunate design_rott6_var workaround qnd
            unclear what the root search is needed for
        
        """
        from scipy.optimize import newton
        
        def design_rott6_var(rhoL,rho,rhoR,g,b0,hL,hR):
            mc=rhoR*hR
            mt=rhoL*hL + rho*(2.*b0) + mc     
            #compute center of gravity and geometry
            s=(-rhoL*hL*b0+rhoR*hR*b0)/mt
            b =b0-s
            bL=b0+s
            a=rhoL*hL*(hL/2.) / mt
            c=hR/2.
            #compute momenta of inertia
            It= (rhoL*(hL)**3)/3.+rhoL*hL*(bL**2)+ (rho*(bL**3))/3. + (rho*(b**3))/3. + rhoR*hR*(b)**2
            Ic= (rhoR*(hR)**3)/3.
            return (mc,mt,s,a,b,bL,c,Ic,It)

        def aux_funct(hhL):
            """auxillary function."""
            if (self.M>0):
                hhR=np.sqrt((self.rhoL/self.rhoR)*(hhL/self.M)**2)
            else:
                hhR=hrat*hhL
            
            mc,mt,s,a,b,bL,c,Ic,It=design_rott6_var(self.rhoL,self.rho,self.rhoR,self.g,self.b0,hhL,hhR)
            lambda2=self.g*mt*a/It
            omega2=self.g*mc*c/Ic                
            return omega2/lambda2-self.R**2   
      
        hL0=(self.hL)*1.
        hR0=(self.hR)*1.
        hrat=hL0,hR0 #ratios
        #find zeros of aux_funct, starting from hL0 guess
        self.hL=newton(aux_funct,hL0)
        
        if self.M>0:
            self.hR=np.sqrt((self.rhoL/self.rhoR)*(self.hL/self.M)**2)
        else:
            self.hR=hrat*self.hL  
        
            

      
     
     
def RottDynamics(state,tr,beta0,lambda2,omega2,eta,zeta): #state in phasespace, timerange pend_properties
    """Casts the equations into a set of first order ODEs."""
    a,g,a1,g1=state
    C=np.cos(g-a-beta0)
    S=-np.sin(g-a-beta0)
    sa=np.sin(a)
    sg=np.sin(g)
    
    #second derivatives of alpha and gamma
    a2=(-lambda2*sa+eta*C*omega2*sg-eta*(zeta*C*S*a1**2+S*g1**2))/(1-eta*zeta*C**2)
    g2=(-omega2*sg+zeta*C*lambda2*sa+zeta*(eta*C*S*g1**2+S*a1**2))/(1-eta*zeta*C**2)
    
    return (a1,g1,a2,g2)




class RottState:
    """Contains the State of the Pendulum as Phase Space coordinate.
    TODO:
    ----
        Redundant?
    """
    def __init__(self,a,g,adot,gdot):
        self.a=a
        self.g=g
        self.adot=adot
        self.gdot=gdot
        self.vec=(a,g,adot,gdot)
      
    
class RottTrajectoryInPhasespace:
    """Contains all the information on the phase space trajectory of the Rott
    Pendulum.
    
    Parameters:
    ----------
        to be filled in!
    
    TODO:
    ----
        Seperate Integration schemes from it
        Make function RottDynamics a method of the class    
    
    """
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
        return self.a[i],self.g[i],self.adot[i],self.gdot[i]
        
    def __add__(self,traj2):
        figname='addition'
        self.plot(figname=figname,savefig=False,holdagain=True)
        traj2.plot(figname=figname,savefig=True,holdagain=False)
        
    def integrate(self,int_scheme=0):
        """Given the initial state and a scheme, this integrates the trajectory in phasespace.
        Parameters
        ----------
            int_scheme : number
            Numer of integrationscheme to be used
        """           
        initstate=(self.setup).initvec
        params=(self.pp).bl2o2ez()
        output=odeint(RottDynamics,initstate,self.timeinsts,args=params)
        output=np.transpose(output)
        #cast the output to the values a,g,adot,gdot
        self.a,self.g,self.adot,self.gdot=output
        #"diagnose" the energies
        self.energetics()

    def energetics(self):
        """Given the evolution in phasespace x(t)=(a,g,adot,gdot) calculates the energies.

        -ekin (kinetic energies)
        -epot (potential energies
        -etot=ekin+epot (the conserved sum of both, i.e. the total energy)
        are computed and made a property of the instance.
        """
        #plotme=False
        a=self.a
        g=self.g
        adot=self.adot
        gdot=self.gdot
        #alias for pendulum properties
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
        #energy calculations
        self.ekin=0.5*It*adot**2+0.5*Ic*gdot**2+mc*b*c*C*adot*gdot
        self.epot=-gpp*(mt*app*np.cos(a)+mc*c*np.cos(g))      
        self.etot=self.ekin+self.epot
        #evolution of energy loss (numerical dissipation)
        self.dE=max(self.etot[0]-np.min(self.etot),np.max(self.etot)-self.etot[0])
            
    def plot(self,figname,var='adotgdot',savefig=True,holdagain=False):
        """Function to plot the time series of the variables in var.

        TODO:
        ----
            so far only adotgdot can be plotted
            add: energies etc
        """

        if not (var=='adotgdot'):
            print('no valid variable to plot')
            raise UserWarning
        f,ax = plt.subplots()
        ax.plot(self.timeinsts,self.adot,linewidth=0.5)
        #plt.hold(True)
        ax.plot(self.timeinsts,self.gdot,linewidth=0.5)
        
        if not holdagain:
            #plt.legend()
            f.savefig(figname+str(self.time)+'.pdf',dpi=300)
            
        #plt.hold(holdagain)
        
    def animate_rott(self,videoname='basic_animation'):#,framerate=20*(self.setup).dt):
        """Creates an animation of the Rott pendulum
        and creates video 'videoname.mp4'.
       
        PRE:
        ----
            Trajectory must already be calculated.

        Parameters:
        ----------
            videoname: string
            videoname (without format suffices)

        TODO:
        ----
            Add different geometries of the Pendulum
            Add capability to modify framerate
    
        """
     
        #set up the figure,axes,line,text boxes for time and energy
        fig, ax = plt.subplots()
        line, = ax.plot([], [],'-', lw=3)
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        ax.set_aspect('equal')
        time_text = ax.text(0.02, 0.95, '')#, transform=ax.transAxes)
        energy_text = ax.text(0.02, 0.90, '')#, transform=ax.transAxes)
        
        #function to produce the axes background that does not change
        def init():
            line.set_data([],[])
            time_text.set_text('')
            energy_text.set_text('')
            return line, time_text, energy_text
        
        #function to add lines for every snapshot 
        def animate(i):
            ax.clear()
            ax.set_xlim([-1.5,1.5])
            ax.set_ylim([-1.5,1.5])
            tinst = 16*i*(self.setup).delt
            #time_text.set_text('time=%.1f'%(self.timeinsts)[4*i])
            #energy_text.set_text('energy=%.1f'%(self.etot)[4*i])
            p5,p3,p2,p4=self.pivot_coordinates(tinst)
            #retrieve position of pivot points
            line, = ax.plot([], [],'-', lw=3)
            line.set_data(zip(p5,p3,p2,p4))
            line.set_c('b')
            
            #line.set_data(zip(p2,p4))
            #line.set_c('r')
            ax.plot([p2[0],p4[0]],[p2[1],p4[1]])
            
            ax.set_axis_off()
            return mplfig_to_npimage(fig)

        #timing for the animation
        from time import time
        t0=time()
        animate(0)
        t1=time()
        
        #create and save the animation
        interval=4*1200*(self.setup).delt
        #myanim=animation.FuncAnimation(fig,animate,init_func=init,frames=1200,blit=True)         
        animation = VideoClip(animate, duration = 20)

        animation.ipython_display(fps = 20)
       
            
   
    def pivot_coordinates(self,timeinstance=10):#self,time):
        """Compute the coordinates of the four pivot points of the Rott Pendulum
        at the timeslice before timeinstance.
        p3______p2
        |        \
        |         \
        |          \
        |           p4
        p5			  
        Parameters
        ----------
        timeinstance : float
            desired time instance
        
        Returns
        ------
        coordinates : four coordinates in arrays
                p5,p3,p2,p4=((xp_i),(yp_i)) for i=1,2,3,4
        Raises
        ------
        """
        if timeinstance>self.time:
            print('timeindex given out of integration range')
            pass
        
        #get timeslice index idx befor the time instance
        idx=int(timeinstance//(self.setup).delt)
        alp=(self.a)[idx]
        gam=(self.g)[idx]
      
        pivot=np.zeros(2)

        def er(angle):
            """Returns the vector that points to unit circle at an 
            angle ''angle'' from the abscissa."""
            return np.array([np.cos(angle),np.sin(angle)])


        #geometric helper
        ea=er(alp)
        ea_p=er(alp-np.pi/2)
        eg_p=er(gam-np.pi/2)
        
        #geometry of pendulum
        l1=l2=l3=l4=1.                  
        #five pivot points
        p1=pivot        
        p2=p1+ea*l1
        p3=p1-ea*l2
        p4=p2+eg_p*l3
        p5=p3+ea_p*l4
                
        return p5,p3,p2,p4

    def snapshot(self):
        """Function creates a snapshot of the pendulum.
        
        
        TODO:
        ----
            Still a scelecton.
        
        """
        
        
        Lcolor='green'
        Icolor='red'
        Llwidth=3
        ILwidth=3

        fig=plt.figure() 
        ax = plt.axes(xlim=(-3,3), ylim=(-3, 3))
        plt.hold(True)
        
        plt.axes().set_aspect('equal')
        #L shape:
        #p2->p3
        
        
        t=zip(self.p2,p3)
        ax.plot(t[0],t[1],color=Lcolor,linewidth=Llwidth)
        #p3->p5
        t=zip(p3,p5)
        ax.plot(t[0],t[1],color=Lcolor,linewidth=Llwidth)

        #plt.plot(,color=Lcolor)
        #I shape:
        #p2->p4
        t=zip(p2,p4)
        ax.plot(t[0],t[1],color=Icolor,linewidth=Ilwidth)
        
        #splt.plot(,color=Icolor)
        for piv in allpivots:
            circle=plt.Circle(piv,l1*0.05,color='black')
            fig.gca().add_artist(circle)
        
        #cosmetics        
        #fact=2.5
        #plt.xlim([pivot[0]-fact*l1,pivot[1]+fact*l1])
        #plt.ylim([pivot[0]-fact*l3,pivot[1]+fact*l3])
        #plt.xticks([],[])
        #plt.yticks([],[])
        #plt.grid()        
        #plt.show()
        #return ax.get_lines()
        
       
        #return mpllibobject
        #creates a snapshot in $Format of the current pendulum state
           
   
#0      
class Integrator:
    """
    Encapsulation of different integrator functionalities.
    
    TODO:
        implementation of different integrators (@andreaarteaga)
        Maybe with dispatch methods?    
    """
    
    def __init__(self,ischeme):
        _integratorlist={0:self.dummy_integrator, 1:self.ode2cntr}

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

        
            
