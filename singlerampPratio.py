#Single ramp system for the inlet of an engine of a vehicle going at supersonic speeds
from math import sin, cos, tan, pi, acos, atan
import numpy as np

def my_compute_objective(delta):
        M1 = 2
        gamma = 1.4

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

        print('Po2Po1', Po2Po1)
        print('Po3Po2', Po3Po2)

        return Po3Po2*Po2Po1