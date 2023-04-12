#Single ramp system for the inlet of an engine of a vehicle going at supersonic speeds
from math import sin, cos, tan, pi, acos, atan
import numpy as np

from modopt.api import Problem

class SingleRampInlet(Problem):
    def initialize(self, ):
        self.problem_name = 'SingleRampInlet'
        
        # Single ramp inlet variable initialization
        self.options.declare('M1_Single_Ramp', default=2.0, types=float)
        self.options.declare('gamma', default=1.4, types=float)
        self.options.declare('rho_pen', default=5.0, types=float)

    def setup(self):
        self.add_design_variables('x',
                                  shape=(1, ),
                                  lower=None,
                                  upper=None,
                                  equals=None,
                                  vals=np.array([0.1]))

        self.add_objective('f')

    def setup_derivatives(self):
        # Declare objective gradient and it's shape
        self.declare_objective_gradient(wrt='x', vals=None)
        #self.declare_objective_hessian(of='x', wrt='x', vals=None)

    def compute_objective(self, dvs, obj):
        delta=dvs['x'][0]
        obj['f']=self.my_compute_objective(delta)
        
    def my_compute_objective(self, delta):
        M1 = self.options['M1_Single_Ramp']
        gamma = self.options['gamma']

        """Calculating oblique shock angle"""

        # These equations are provided by Professor Sanchez in his MAE 212 Notes
        print('delta',delta*180/pi)
        lamda = np.sqrt((M1**2-1)**2-3*(1+((gamma-1)/2)*M1**2)*(1+((gamma+1)/2)*M1**2)*tan(delta)**2)
        print('lamda',lamda)
        X = ((M1**2 - 1)**3 - 9*(1 + ((gamma - 1)/2)*M1**2)*(1+((gamma-1)/2)*M1**2+((gamma+1)/4)*M1**4)*(tan(delta)**2))/(lamda**3)
        print('X',X)
        beta = atan(((M1**2-1)+2*lamda*cos((1/3)*(4*pi+acos(X))))/(3*(1+(gamma-1)/2*M1**2)*tan(delta)))
        print('this is beta',beta*180/pi)

        """Mach number calcs for oblique shock"""
        
        # These equations are provided by Professor Sanchez in his MAE 212 Notes for oblique shocks
        M1n = M1*sin(beta)
        M2n = np.sqrt(((gamma-1)*M1n**2+2)/((2*gamma*M1n**2)-(gamma-1)))
        M2 = M2n/(sin(beta-delta))
        
        """Mach number calcs for normal shock"""
        M3=np.sqrt(((gamma-1)*M2**2+2)/((2*gamma*M2**2)-(gamma-1))) # MAE 212 Normal Shock
        #Note M3 has to be <1
        print('This is M2',M2)
        print('This is M3',M3)

        
        """Pressure Calcs"""

        P2P1 = (2*gamma*(M1**2)*((sin(beta))**2)-(gamma-1))/(gamma+1) # NASA oblique
        Po2Po1 = ((((gamma+1)*(M1**2)*(sin(beta))**2)/((gamma-1)*(M1**2)*((sin(beta))**2)+2))**(gamma/(gamma-1)))*(1/P2P1)**(1/(gamma-1)) # NASA oblique
        Po3Po2 = ((2*gamma*M2**2-gamma+1)/(gamma+1)*((2+(gamma-1)*M2**2)/((gamma+1)*M2**2))**gamma)**(-1/(gamma-1)) # MAE 212 Normal Shock
        
        Po1Po2 = 1/Po2Po1
        Po2Po3 = 1/Po3Po2

        print('Po2Po1', Po2Po1)
        print('Po3Po2', Po3Po2)


        """Applying Exterior method"""
        
        rho_pen = self.options['rho_pen']
        
        if ((delta*180/pi)>0):
            if((delta*180/pi) <90):
                phi=0
            else:
                phi=(delta*180/pi)**2
        else:
            phi=(delta*180/pi)**2

        #Minimizing Po1Po3
        return Po1Po2*Po2Po3+ rho_pen*phi

    def oblique_p0dp0u(self, Mu, Md, beta, alpha):
        p0dp0u = Mu

    def compute_objective_gradient(self, dvs, grad):
        delta = dvs['x'][0]
        F = self.my_compute_objective(delta)
        h = 1e-5

        xh = delta+h
        Fh = self.my_compute_objective(xh)
        dfdx_FOFD = (Fh-F)/h

        grad['x']= dfdx_FOFD
