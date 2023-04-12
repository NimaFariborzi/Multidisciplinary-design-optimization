#Double ramp system for the inlet of an engine of a vehicle going at supersonic speeds

from math import sin, cos, tan, pi, acos, atan
import numpy as np

from modopt.api import Problem


class DoubleRampInlet(Problem):
    def initialize(self, ):
        self.problem_name = 'DoubleRampInlet'

        # Double ramp inlet variable initialization
        self.options.declare('M1_Double_Ramp', default=2.0, types=float)
        self.options.declare('gamma', default=1.4, types=float)
        self.options.declare('rho_pen', default=5.0, types=float)

    def setup(self):
        self.add_design_variables('x',
                                  shape=(2, ),
                                  lower=None,
                                  upper=None,
                                  equals=None,
                                  vals=np.array([.05,.05])) #.01, .01

        self.add_objective('f')

    def setup_derivatives(self):
        # Declare objective gradient and its shape
        self.declare_objective_gradient(wrt='x', vals=None)

    def compute_objective(self, dvs, obj):
        delta1=dvs['x'][0]
        delta2=dvs['x'][1]
        obj['f']=self.my_compute_objective(delta1,delta2)
        
    def my_compute_objective(self, delta1,delta2):
        M1 = self.options['M1_Double_Ramp']
        gamma = self.options['gamma']
        print('delta1',delta1*180/pi)
        print('delta2',delta2*180/pi)
        print('delta1+delta2', (delta1+delta2)*180/pi)

        """Calculating first oblique shock angle"""
        
        # These equations are provided by Professor Sanchez in his MAE 212 Notes for oblique shocks
        lamda1 = np.sqrt((M1**2-1)**2-3*(1+((gamma-1)/2)*M1**2)*(1+((gamma+1)/2)*M1**2)*tan(delta1)**2)
        # print('lamda1',lamda1)
        X1 = ((M1**2-1)**3-9*(1+((gamma-1)/2)*M1**2)*(1+((gamma-1)/2)*M1**2+((gamma+1)/4)*M1**4)*(tan(delta1)**2))/(lamda1**3)
        # print('X1',X1)
        beta1 = atan(((M1**2-1)+2*lamda1*cos((1/3)*(4*pi+acos(X1))))/(3*(1+(gamma-1)/2*M1**2)*tan(delta1)))
        # print('this is beta1',beta1*180/pi)

        """Mach number calcs for first oblique shock"""

        # These equations are provided by Professor Sanchez in his MAE 212 Notes for oblique shocks
        M1n = M1*sin(beta1)
        M2n = np.sqrt(((gamma-1)*M1n**2+2)/((2*gamma*M1n**2)-(gamma-1)))
        M2 = M2n/(sin(beta1-delta1))
        print('this is M2',M2)

        """Calculating second oblique shock angle"""

        # These equations are provided by Professor Sanchez in his MAE 212 Notes for oblique shocks
        lamda2 = np.sqrt((M2**2-1)**2-3*(1+((gamma-1)/2)*M2**2)*(1+((gamma+1)/2)*M2**2)*tan(delta2)**2)
        # print('lamda2',lamda2)
        X2 = ((M2**2-1)**3-9*(1+((gamma-1)/2)*M2**2)*(1+((gamma-1)/2)*M2**2+((gamma+1)/4)*M2**4)*(tan(delta2)**2))/(lamda2**3)
        print('X2',X2)
        beta2 = atan(((M2**2-1)+2*lamda2*cos((1/3)*(4*pi+acos(X2))))/(3*(1+(gamma-1)/2*M2**2)*tan(delta2)))
      
        # print('this is beta2',beta2*180/pi)

        """Mach number calcs for second oblique shock"""

        # These equations are provided by Professor Sanchez in his MAE 212 Notes for oblique shocks
        M2nn = M2*sin(beta2)
        # M3n = ((gamma-1)*M2nn**2+2)/(2*gamma*M2nn**2+1-gamma)
        M3n = np.sqrt(((gamma-1)*M2nn**2+2)/((2*gamma*M2nn**2)-(gamma-1)))
        M3 = M3n/(sin(beta2-delta2))

        """Mach number calcs for normal shock"""
        # M4 = np.sqrt(((gamma-1)*M2**2+2)/(2*gamma*M2**2+1-gamma))
        M4 = np.sqrt(((gamma-1)*M3**2+2)/((2*gamma*M3**2)-(gamma-1))) # MAE 212 Normal Shock
        print('This is M3', M3)
        print('This is M4', M4)

        """Pressure calcs"""
                
        # Oblique Shock 1
        P2P1 = (2*gamma*(M1**2)*((sin(beta1))**2)-(gamma-1))/(gamma+1) # NASA oblique
        Po2Po1 = ((((gamma+1)*(M1**2)*(sin(beta1))**2)/((gamma-1)*(M1**2)*((sin(beta1))**2)+2))**(gamma/(gamma-1)))*(1/P2P1)**(1/(gamma-1)) # NASA oblique
        
        Po1Po2 = 1/Po2Po1

        print('Po2Po1', Po2Po1)

        # Oblique Shock 2
        P3P2 = (2*gamma*(M2**2)*((sin(beta2))**2)-(gamma-1))/(gamma+1) # NASA oblique
        Po3Po2 = ((((gamma+1)*(M2**2)*(sin(beta2))**2)/((gamma-1)*(M2**2)*((sin(beta2))**2)+2))**(gamma/(gamma-1)))*(1/P3P2)**(1/(gamma-1)) # NASA oblique
        
        Po2Po3 = 1/Po3Po2

        print('Po3Po2', Po3Po2)

        # Normal Shock
        Po4Po3 = ((2*gamma*M3**2-gamma+1)/(gamma+1)*((2+(gamma-1)*M3**2)/((gamma+1)*M3**2))**gamma)**(-1/(gamma-1)) # MAE 212 Normal Shock
        Po3Po4 = 1/Po4Po3
        print('Po4Po3', Po4Po3)

        """Exterior point penalty method"""

        rho_pen = self.options['rho_pen']
        
        # phi 1: delta 1 >= 0
        if ((delta1*180/pi)>0):
            if((delta1*180/pi)<90):
                phi1=0
            else:
                phi1=(delta1*180/pi)**2
        else:
            phi1=(delta1*180/pi)**2

        # phi 2: delta 2 >= 0
        if ((delta2*180/pi)>0):
            if((delta2*180/pi)<90):
                phi2 = 0
            else:
                phi2 = (delta2*180/pi)**2
        else:
            phi2=(delta2*180/pi)**2

        #phi 3: delta 1 + delta 2 >= 0
        if(((delta1+delta2)*180/pi)>0):
            if(((delta1+delta2)*180/pi)<90):
                phi3 = 0
            else:
                phi3 = ((delta1+delta2)*180/pi)**2
        else:
            phi3 = ((delta1+delta2)*180/pi)**2

        #Minimizing Po1Po4
        return (Po1Po2*Po2Po3*Po3Po4)+(rho_pen*phi1)+(rho_pen*phi2)+(rho_pen*phi3)


    def compute_objective_gradient(self, dvs, grad):
        delta1 = dvs['x'][0]
        delta2 = dvs['x'][1]
        x = dvs['x']

        F = self.my_compute_objective(delta1, delta2)
        h = 1e-5
        e = np.zeros(len(x))
        dfdx_FOFD = np.zeros(len(x))

        for i in range(len(x)):
            e[i] = 1
            xh = x+h*e
            xh1 = xh[0]
            xh2 = xh[1]
            Fh = self.my_compute_objective(xh1, xh2)
            dfdx_FOFD[i] = (Fh-F)/h
            e[i] = 0

        grad['x']= dfdx_FOFD