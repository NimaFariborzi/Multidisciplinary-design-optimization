# double ramp system for the inlet of an engine of a vehicle going at supersonic speeds
from math import sin, cos, tan, pi, acos, atan
import numpy as np

def my_compute_objective(delta1, delta2):
        M1 = 2
        gamma = 1.4
        print('delta1',delta1*180/np.pi)
        print('delta2',delta2*180/np.pi)

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

        return Po4Po3*Po3Po2*Po2Po1