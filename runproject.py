#Run file for final project
import numpy as np
import matplotlib.pyplot as plot
import seaborn as sns
from singlerampPratio import my_compute_objective as PratioSing
from doublerampPratio import my_compute_objective as PratioDoub

from singleramp import SingleRampInlet as SRI
from doubleramp import DoubleRampInlet as DRI
from modopt.optimization_algorithms import QuasiNewton

# Solver Definitions
tol = 1E-8
max_itr = 100
gamma = 1.4 # normal atmospheric air is 1.4, (add .0 to the end if whole number)
rho_pen = 5.0 # penalty term for exterior point method of inequality constraint, (add .0 to the end if whole number)

# Single Ramp Inlet inputs
M1_Single_Ramp = 2.0 # make sure this is a float (add .0 to the end if whole number)

# Double Ramp Inlet inputs
M1_Double_Ramp = 2.0 # make sure this is a float (add .0 to the end if whole number)

"""Solving SingleRampInlet"""

prob = SRI(M1_Single_Ramp=M1_Single_Ramp, gamma = gamma)
optimizerSRI = QuasiNewton(prob, opt_tol=tol, max_itr=max_itr)
optimizerSRI.setup()
optimizerSRI.solve()
optimizerSRI.print_results(summary_table=True)

"""Plotting Single Ramp Inlet"""

# Plot of gradient norm vs iterations
sns.set()

plot.semilogy(optimizerSRI.outputs['itr'], optimizerSRI.outputs['opt'])
plot.xlabel('iterations')
plot.ylabel('gradient norm')
plot.title('Single Ramp Gradient Norm vs Iterations')
plot.show()

# #plot of change in stagnation pressure ratio p01/p04
# sns.set()

# plot.plot(optimizerSRI.outputs['x'],optimizerSRI.outputs['obj'])
# plot.xlabel('x')
# plot.ylabel('stagnation pressure')
# plot.title('change in stagnation pressure vs x')
# plot.show()

#plot of change in stagnation pressure ratio p03/p01
sns.set()
x_sing = optimizerSRI.outputs['x']
y_sing = np.zeros(len(x_sing))
for i in range(len(x_sing)):
    delta_sing = x_sing[i]
    y_sing[i] = PratioSing(delta_sing)

x_sing = x_sing*180/np.pi

plot.plot(x_sing, y_sing)
plot.xlabel('delta (degrees)')
plot.ylabel('P03/P01')
plot.title('Single Ramp Stagnation Pressure Ratio P03/P01 vs Ramp Angle delta')
plot.show()

#plot of objective change vs iterations

plot.plot(optimizerSRI.outputs['itr'], y_sing)
plot.xlabel('iterations')
plot.ylabel('P03/P01')
plot.title('Single Ramp Stagnation Pressure Ratio P03/P01 vs iterations')
plot.show()

"""Solving DoubleRampInlet"""

prob2=DRI(M1_Double_Ramp = M1_Double_Ramp, gamma = gamma)
optimizerDRI=QuasiNewton(prob2, opt_tol=tol, max_itr=max_itr)
optimizerDRI.setup()
optimizerDRI.solve()
optimizerDRI.print_results(summary_table=True)

# Plot of gradient norm vs iterations

sns.set()
plot.semilogy(optimizerDRI.outputs['itr'], optimizerDRI.outputs['opt'])
plot.xlabel('iterations')
plot.ylabel('gradient norm')
plot.title('Double Ramp Gradient Norm vs Iterations')
plot.show()


#plot of change in stagnation pressure ratio p03/p01

x_1 = optimizerDRI.outputs['x'][:,0]*180/np.pi
x_2 = optimizerDRI.outputs['x'][:,1]*180/np.pi
x_2 = x_2+x_1

y = 1/optimizerDRI.outputs['obj']

sns.set()
plot.plot(x_1, y, color = 'orange', label = 'Ramp 1')
plot.plot(x_2, y, color = 'blue', label = 'Ramp 2')
plot.legend()
plot.xlabel('delta (degrees)')
plot.ylabel('P04/P01')
plot.title('Double Ramp Stagnation Pressure Ratio P04/P01 vs Ramp Angle delta')
plot.show()

plot.plot(optimizerDRI.outputs['itr'], y)
plot.xlabel('iterations')
plot.ylabel('P04/P01')
plot.title('Double Ramp Stagnation Pressure Ratio P03/P01 vs iterations')
plot.show()

# x_tot = optimizerDRI.outputs['x']
# x_doub_1 = optimizerDRI.outputs['x'][:,0]
# x_doub_2 = optimizerDRI.outputs['x'][:,1]
# y_doub = np.zeros(len(x_doub_1))
# for i in range(len(x_doub_1)):
#     delta_doub_1 = x_doub_1[i]
#     delta_doub_2 = x_doub_2[i]
#     y_doub[i] = PratioDoub(delta1 = delta_doub_1, delta2 = delta_doub_2)

# x_doub_1 = x_doub_1*180/np.pi
# x_doub_2 = x_doub_2*180/np.pi

# sns.set()
# plot.plot(x_doub_1, y_doub, color = 'orange', label = 'Ramp 1')
# plot.plot(x_doub_2, y_doub, color = 'blue', label = 'Ramp 2')
# plot.legend()
# plot.xlabel('delta (degrees)')
# plot.ylabel('P04/P01')
# plot.title('Double Ramp Stagnation Pressure Ratio P04/P01 vs Ramp Angle delta')
# plot.show()

# """Plotting Contour Plot"""

# sns.set()
# plot.figure(figsize=(8, 5))

# N = 1000
# x1 = np.linspace(15, 20, N)#.reshape(1, -1)
# x2 = np.linspace(15, 20, N)#.reshape(-1, 1)

# x1 = x1*np.pi/180
# x2 = x2*np.pi/180
# obj = np.zeros(len(x1))

# for i in range(len(x1)):
#     delta_doub_1 = x1[i]
#     delta_doub_2 = x2[i]
#     obj[i] = PratioDoub(delta1 = delta_doub_1, delta2 = delta_doub_2)

# x1 = x1*180/np.pi
# x2 = x2*180/np.pi

# # x1, x2 = x1.flatten(), x2.flatten()

# cs = plot.contour(x1, x2, obj, 50)

# cbar = plot.colorbar(cs)

# # Plot Quadratic Penalty

# x_QP = optimizer_Penalty.outputs['x']

# x1_QP = x_QP[0:-1, 0]

# x2_QP = x_QP[0:-1, 1]

# plot.plot(x1_QP, x2_QP, color = 'red', marker = 'o', label = 'Quadratic Penalty Method')

# # Plot Method of Lagrange multipliers

# x_LM = optimizer_Lag.outputs['x']

# x1_LM = x_LM[0:-1, 0]

# x2_LM = x_LM[0:-1, 1]

# plot.plot(x1_LM, x2_LM, color = 'blue', marker = 'o', label = 'Method of Lagrange multipliers')

# plot.legend()
# plot.xlabel('x1')
# plot.ylabel('x2')
# plot.title('Contour Plot of Rosenbrock Problem with Equality Constraint Optimization Solvers')

# plot.show()