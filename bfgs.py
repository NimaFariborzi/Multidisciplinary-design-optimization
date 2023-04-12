import numpy as np
import time
from modopt.api import Optimizer

class BFGS2(Optimizer):
    def initialize(self):
        # Name your algorithm
        self.solver_name = 'BFGS'
        self.obj = self.problem._compute_objective
        self.grad = self.problem._compute_objective_gradient
        self.options.declare('max_itr', default=1000, types=int)
        self.options.declare('opt_tol', default=1e-5, types=float)
        # Specify format of outputs available from your optimizer after each iteration

        self.default_outputs_format = {
            'itr': int,
            'obj': float,
            # for arrays from each iteration, shapes need to be declared
            'xL': (float, (self.problem.nx, )),
            'opt': float,
            'time': float,
        }
        # Enable user to specify, as a list, which among the available outputs
        # need to be stored in memory and written to output files
        self.options.declare('outputs',
                             types=list,
                             default=['itr', 'obj', 'xL', 'opt', 'time'])
    def solve(self):
        nx = self.problem.nx
        x = self.problem.x.get_data()
        opt_tol = self.options['opt_tol']
        max_itr = self.options['max_itr']
        obj = self.obj
        grad = self.grad


        start_time = time.time()
        # Setting intial values for initial iterates
        x_k = x * 1.
        f_k = obj(x_k)
        g_k = grad(x_k)
        B_k= np.identity(len(x_k))
        #inital values for lamda
        lam=0
        con=x_k[0]+x_k[1]-1

        Jh=np.array([1,1])
        c_k= x_k[0] + x_k[1] - 1
        r= np.hstack([grad(x_k)+np.transpose(Jh)*lam,con])
        Jr1= np.vstack((B_k,Jh))
        Jr2= np.hstack((Jh,0)).reshape((3,1))
        Jr= np.block([Jr1, Jr2])      

        # Iteration counter
        itr = 0
        # Optimality
        opt = np.linalg.norm(g_k)

        # Initializing outputs
        self.update_outputs(itr=0,
                            xL=x_k,
                            obj=f_k,
                            opt=opt,
                            time=time.time() - start_time)
        while (opt > opt_tol and itr < max_itr):
            itr_start = time.time()
            itr += 1
            

            con=x_k[0]+x_k[1]-1


            p_k=np.transpose(-r)@np.linalg.inv(Jr)
            p_k=np.linalg.inv(Jr)@-r
            print('pk',p_k)
           
            itrwhile=1

            #Linesearch loop
            A=1
            rho=.1
            c1=1e-4
            a=A

            xnew=x_k+a*p_k[:2]
            LHS=obj(xnew)+3/2*np.transpose(xnew[0]+xnew[1]-1)*(xnew[0]+xnew[1]-1)+np.transpose(lam+a*p_k[2])*(xnew[0]+xnew[1]-1)
            RHS=(f_k+c1*a*np.transpose(g_k)@p_k[:2])

            while True:
                if (itrwhile>30):
                    break
                if (LHS<=RHS):
                    break
                itrwhile +=1
                a = a*rho
                xnew = x_k+a*p_k[:2]
                LHS=obj(xnew)+5/2*np.transpose(xnew[0]+xnew[1]-1)*(xnew[0]+xnew[1]-1)+np.transpose(lam+a*p_k[2])*(xnew[0]+xnew[1]-1) 
                
            lam += p_k[2]
            #print('this is lam', lam)
            y_k=grad(xnew)+np.transpose(Jh)*lam-(grad(x_k)+np.transpose(Jh)*lam)
            x_k=xnew
            g_k = grad(x_k)
            s_k = a*p_k[:2]

            #BFGS updating the B_K
            dk = s_k
            wk = y_k

            tol1 = 1e-14
            Bd = B_k.dot(dk)
            wTd = np.dot(wk, dk)
            sign = 1. if wTd >= 0. else -1.
            if abs(wTd) > tol1:
                B_k = B_k + np.outer(wk, wk) / wTd
            else:
                B_k = B_k + np.outer(wk, wk) / sign / tol1
            dTBd = np.dot(dk, Bd)
            sign = 1. if dTBd >= 0. else -1.
            if abs(dTBd) > tol1:
                B_k -= np.outer(Bd, Bd) / dTBd
            else:
                B_k -= np.outer(Bd, Bd) / sign / tol1

            r= np.hstack([grad(x_k)+np.transpose(Jh)*lam,con])
            opt = np.linalg.norm(r)

            Jr1= np.vstack((B_k,Jh))
            Jr2= np.hstack((Jh,0)).reshape((3,1))
            Jr= np.block([Jr1, Jr2])

            f_k = obj(x_k) 



            
        # Another way of doing BFGS 
            # first= np.outer(y_k,y_k)/(np.inner(y_k,s_k))
            # top= (np.dot(np.outer(np.dot(B_k,s_k),s_k),B_k))
            # bot=(np.dot(np.dot(s_k.T,B_k),s_k))
            # second=top/bot

            # B_k += first-second

            self.update_outputs(itr=itr,
                                xL=x_k,
                                obj=f_k,
                                opt=opt,
                                time=time.time() - start_time)
        # Run post-processing for the Optimizer() base class
        self.run_post_processing()
        end_time = time.time()
        self.total_time = end_time - start_time
