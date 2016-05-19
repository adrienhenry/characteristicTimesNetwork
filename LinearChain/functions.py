import numpy as np
from scipy import linalg
import scipy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib import rc
import time
######################################################################
##                         derivatives_MA                           ##
######################################################################
## Function computing the mass action derivative for the metabolite
## concentrations assuming the two extreme metabolites are constant.
## The function takes 3 arguments:
## - t: A lot of ode solver requires this variable
## - y: vector of metabolite concentrations
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,...,a_N-1,b_N-1]
def derivatives_MA(t,y,args):
    N = len(y)-2  
    a = [args[i] for i in range(0,2*(N+1)-1,2)]
    b = [args[i] for i in range(1,2*(N+1),2)]
    dy = np.zeros(N+2)
    v = np.zeros(N+1)
    for i in range(0,N+1):
        v[i] = a[i]*y[i] - y[i+1]*b[i]
    dy[0] = dy[N+1] = 0
    for i in range(1,N+1):
        dy[i] = v[i-1] - v[i]
    return dy    
######################################################################


######################################################################
##                          derivatives_MM                          ##
######################################################################
## Function computing the Michaelis-Menten derivative for the
## metabolite concentrations assuming the two extreme metabolites are
## constant.
## The function takes 3 arguments:
## - t: A lot of ode solver requires this variable
## - y: vector of metabolite concentrations
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,KmS_0,KmP_0,...,a_N-1,b_N-1,
## KmS_N+1,KmP_N+1]
def derivatives_MM(t,y,args):
    N = len(y)-2
    a = [args[i] for i in range(0,4*(N+1),4)]
    b = [args[i] for i in range(1,4*(N+1),4)]
    KmS = [args[i] for i in range(2,4*(N+1),4)]
    KmP = [args[i] for i in range(3,4*(N+1),4)]
    dy = np.zeros(N+2)
    v = np.zeros(N+1)
    for i in range(0,N+1):
        v[i] = (a[i]*y[i] - b[i]*y[i+1])/(1 + y[i]/KmS[i] + y[i+1]/KmP[i])
    dy[0] = dy[N+1] = 0
    for i in range(1,N+1):
        dy[i] = v[i-1] - v[i]
    return dy
######################################################################


######################################################################
##                          rates_MM                          ##
######################################################################
## Function computing the Michaelis-Menten derivative for the
## metabolite concentrations assuming the two extreme metabolites are
## constant.
## The function takes 3 arguments:
## - t: A lot of ode solver requires this variable
## - y: vector of metabolite concentrations
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,KmS_0,KmP_0,...,a_N-1,b_N-1,
## KmS_N+1,KmP_N+1]
def rates(y,args):
    N = len(y)-2
    a = [args[i] for i in range(0,4*(N+1),4)]
    b = [args[i] for i in range(1,4*(N+1),4)]
    KmS = [args[i] for i in range(2,4*(N+1),4)]
    KmP = [args[i] for i in range(3,4*(N+1),4)]
    dy = np.zeros(N+2)
    vp = np.zeros(N+1)
    vm = np.zeros(N+1)
    for i in range(0,N+1):
        vp[i] = a[i]/(1 + y[i]/KmS[i] + y[i+1]/KmP[i])
        vm[i] = b[i]/(1 + y[i]/KmS[i] + y[i+1]/KmP[i])
    return (vp,vm)
######################################################################



######################################################################
##                       evaluate_concST_MA                         ##
######################################################################
## Function that computes the steady states concentration in linear
## chain where the dynamic is given by mass action law.
## The generated function takes 3 arguments:
## - y: vector of metabolite concentrations, only the the first
## concentration and the vector length are use in the function.
## - args: the parameters for the function., the should be given in
# following order [a_0,B_0,...,a_N-1,b_N-1]
def evaluate_concST_MA(y,args,Flux):
    N = len(y)-2  
    a = [args[i] for i in range(0,2*(N+1)-1,2)]
    b = [args[i] for i in range(1,2*(N+1),2)]    
    conc_st = np.zeros(N+2)
    conc_st[0] = y[0]
    for i in range(1,N+2):
        conc_st[i] = a[i-1]/b[i-1]*conc_st[i-1] - Flux/b[i-1]
    return conc_st
######################################################################



######################################################################
##                       evaluate_concST_MM                         ##
######################################################################
## Function that computes the steady states concentration in linear
## chain where the dynamic is given by mass action law.
## The generated function takes 3 arguments:
## - y: vector of metabolite concentrations, only the the first
## concentration and the vector length are use in the function.
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,KmS_0,KmP_0,...,a_N-1,b_N-1,
## KmS_N+1,KmP_N+1]
def evaluate_concST_MM(y,args,Flux):
    N = len(y)-2    
    a = [args[i] for i in range(0,4*(N+1),4)]
    b = [args[i] for i in range(1,4*(N+1),4)]
    KmS = [args[i] for i in range(2,4*(N+1),4)]
    KmP = [args[i] for i in range(3,4*(N+1),4)]    
    conc_st = np.zeros(N+2)
    conc_st[0] = y[0]
    for i in range(1,N+2):
        conc_st[i] = (conc_st[i-1]*(a[i-1]-Flux/KmS[i-1]) - Flux)/(b[i-1] + Flux/KmP[i-1])
    return conc_st
######################################################################

######################################################################
##                         compute_jacobian                         ##
######################################################################
## Function that computes jacobian for the linear metabolite chain
## system. (mass action rate laws)
## The generated function takes 3 arguments:
## - y: vector of metabolite concentrations, only the length is used
## in the function.
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,...,a_N-1,b_N-1]
def compute_jacobian(size,args):
    N = size
    args = np.array(args)
    a = args[np.arange(0,(N+1)*2,2)]
    b = args[np.arange(1,(N+1)*2,2)]
    r_0 = np.zeros(N+2)
    r_last = np.zeros(N+2)
    jacob = [r_0]
    for i in range(1,N+1):
        r_i = np.zeros(N+2)
        r_i[[ i-1 , i , i+1 ]] = np.array([ a[i-1] , -(b[i-1] + a[i]), b[i] ])
        jacob.append(r_i)
    jacob.append(r_last)
    return np.array(jacob)
######################################################################


######################################################################
##                       compute_jacobian_MM                        ##
######################################################################
## Function that computes jacobian for the linear metabolite chain
## system. (Michaelis-Menten rate laws)
## The generated function takes 3 arguments:
## - y: vector of metabolite concentrations.
## - args: the parameters for the function., the should be given in
## following order [a_0,B_0,KmS_0,KmP_0,...,a_N-1,b_N-1,
## KmS_N+1,KmP_N+1]
def compute_jacobian_MM(y,args):
    N = len(y)-2
    a = [args[i] for i in range(0,4*(N+1)-1,4)]
    b = [args[i] for i in range(1,4*(N+1),4)]
    KmS = [args[i] for i in range(2,4*(N+1),4)]
    KmP = [args[i] for i in range(3,4*(N+1),4)]        
    r_0 = np.zeros(N+2)
    r_last = np.zeros(N+2)
    jacob = [r_0]
    for i in range(1,N+1):
        r_i = np.zeros(N+2)
        ma0 = a[i-1]*y[i-1] - b[i-1]*y[i]
        ma1 = a[i]*y[i] - b[i]*y[i+1]
        sat0 = (1+y[i-1]/KmS[i-1] + y[i]/KmP[i-1])
        sat1 = (1+y[i]/KmS[i] + y[i+1]/KmP[i])
        r_i[[ i-1 , i , i+1 ]] = np.array([ a[i-1]/sat0-ma0/(KmS[i-1]*sat0**2) ,-(b[i-1]/sat0 +ma0/(KmP[i-1]*sat0*sat0) + a[i]/sat1 - ma1/(KmS[i]*sat1*sat1)), b[i]/sat1 + ma1/(KmP[i]*sat1**2) ])
        jacob.append(r_i)
    jacob.append(r_last)
    return np.array(jacob)
######################################################################

######################################################################
##                             max_relax                            ##
######################################################################
def max_relax(param):
    A = param[0]
    B = param[1]
    return(1/(A+B - 2*np.sqrt(A*B)))
######################################################################

######################################################################
##                         crossover_size                           ##
######################################################################
def crossover_size(param):
    A = param[0]
    B = param[1]
    return(1 + 2*B*np.pi/(A-B))
######################################################################


######################################################################
##                           eval_transit                           ##
######################################################################
def eval_transit(size , param , norm_type=2):
    A = param[0]
    B = param[1]
    eigen_vec = []
    eigen_val = []
    for mode in range(1,(size+1)):           
        eigen_val.append(2*np.sqrt(A*B)*np.cos(mode*np.pi/(size+1))-(A+B))
        eigen_vec.append([np.exp(0.5*np.log(A/B)*n)*np.sin(mode*np.pi*n/(size+1)) for n in range(1,size+1)])
        eigen_vec[-1] /= np.linalg.norm(eigen_vec[-1])
    eigen_val = np.array(eigen_val)
    if not np.all(eigen_val<0):
        print("ERROR in eval_transilt:\n  positive eigen values")
    eigen_vec = np.transpose(eigen_vec)
    pass_mat = np.linalg.inv(eigen_vec)
    transit_times = []
    if size%2 == 1:
        met_to_perturb = [np.int(size/2)]
    else:
        met_to_perturb = [np.int(size/2-1),np.int(size/2)]
    for i in met_to_perturb:
        ## Projection of the a perturbation on metabolite "i" of amplitude "1"
        ## in the base of eiven vectors. Equivalent to P.[...,1,...]            
        coeff = pass_mat[:,i]
        temp_ev_mat = eigen_vec*coeff
        def func(t):
            D = (temp_ev_mat*np.exp(eigen_val*t)).sum(1)
            return(np.linalg.norm(D,norm_type))
        transit_times.append(quad(func,0,np.Inf)[0])
    return(np.mean(transit_times)) 
######################################################################


######################################################################
##                       norm_time_series                           ##
######################################################################
def norm_time_series(size , param , nb_points , tmax):
    A = param[0]
    B = param[1]
    eigen_vec = []
    eigen_val = []
    for mode in range(1,(size+1)):           
        eigen_val.append(2*np.sqrt(A*B)*np.cos(mode*np.pi/(size+1))-(A+B))
        eigen_vec.append([np.exp(0.5*np.log(A/B)*n)*np.sin(mode*np.pi*n/(size+1)) for n in range(1,size+1)])
        eigen_vec[-1] /= np.linalg.norm(eigen_vec[-1])
    eigen_val = np.array(eigen_val)
    if not np.all(eigen_val<0):
        print("ERROR in eval_transilt:\n  positive eigen values")
    eigen_vec = np.transpose(eigen_vec)
    pass_mat = np.linalg.inv(eigen_vec)
    transit_times = []
    vec_TS_L1 = []
    vec_TS_L2 = []
    T = np.arange(0,tmax,tmax/nb_points)
    if size%2 == 1:
        met_to_perturb = [np.int(size/2)]
    else:
        met_to_perturb = [np.int(size/2-1),np.int(size/2)]
    for i in met_to_perturb:
        ## Projection of the a perturbation on metabolite "i" of amplitude "1"
        ## in the base of eiven vectors. Equivalent to P.[...,1,...]            
        coeff = pass_mat[:,i]
        temp_ev_mat = eigen_vec*coeff
        def func(t,norm_type):
            D = (temp_ev_mat*np.exp(eigen_val*t)).sum(1)
            return(np.linalg.norm(D,norm_type))
        vec_TS_L1.append(np.array([func(t,1) for t in T]))
        vec_TS_L2.append(np.array([func(t,2) for t in T]))
    vec_TS_L1 = np.array(vec_TS_L1).mean(0)
    vec_TS_L2 = np.array(vec_TS_L2).mean(0)    
    return(T,vec_TS_L1,vec_TS_L2) 
######################################################################


######################################################################
##                      Compute transit time                        ##
######################################################################
def compute_lifetime(jacobian,norm_type=1,tau=10,met_to_perturb=[],method= "inverseJacobian"):
    size = len(jacobian)
    if met_to_perturb==[]:
        met_to_perturb = np.arange(1,size+1,1)
    transit_times = []
    ## Compute T by inverting the jacobian
    if method == "inverseJacobian":
        jacobian = jacobian[1:len(jacobian)-1,1:len(jacobian)-1]        
        inverse = np.linalg.inv(jacobian)       
        Ts = [np.linalg.norm(inverse[:,i],1) for i in range(len(inverse))]
        return(np.mean(Ts))
    elif method=="integration":
        dt =  tau/1000
        mat_dt = linalg.expm(jacobian*dt)
    
        for i in met_to_perturb:
            dconc = np.zeros(size+2)
            dconc[i] = 1
            norm_num = [0]
            norm = [np.linalg.norm(dconc,norm_type)]
            tt = [0]
            while True:
                tt.append(tt[-1] + dt)
                dconc = np.dot(mat_dt,dconc)            
                norm.append(np.linalg.norm(dconc,norm_type))
                norm_num.append(0.5*(norm[-1]+norm[-2])*(tt[-1]-tt[-2]))
                if (norm[-1]/norm[0])< 1e-9:
                    break
            transit_times.append(np.sum(norm_num)/norm[0])
        return(np.max(transit_times))
######################################################################


######################################################################
##                           Compute tau                            ##
######################################################################
def compute_tau(jacobian):
    eigen = np.linalg.eig(jacobian[1:-1,1:-1])    
    eigen_values = np.real(eigen[0])
    eigen_vec = eigen[1][:,np.argsort(eigen_values)]
    eigen_values = np.sort(eigen_values)        
    eigen_vec = eigen_vec[:,eigen_values!=0]
    eigen_values = eigen_values[eigen_values!=0]
    tau_ref = np.max(-1/eigen_values)
    return((tau_ref,eigen_vec[:,-1]))
######################################################################


######################################################################
##                     run_homogenous_model_MA                      ##
######################################################################
## Function that helps running a task in a homogenous network where 
## the rate laws are described by mass action.
## 
## Takes 4 arguments:
## - par: vector [a,b], the forward and backard rates
## - flux: flux value
## - size:  Number of non boundary metabolites
## - first_conc: concentration of the first boundary
## - task: Type of task to do
## Returns a tuple with of the with all the task that have been
#" asked.
def print_matrix(mat):
    for row in mat:
        print("       [ %s ]"%"  ".join(["%1.2f"%k for k in row]))

def run_homogenous_model(conc_st,param_vec,typeofmodel="MA",task_list=[""]):        
    a = param_vec[0]
    b = param_vec[1]    
    size = len(conc_st)-2
    A,B = (0,0)
    S=0
    typeoflaw = ""
    
    if typeofmodel == "MA-conc":
        A = a
        B = b
        S = 1
        typeoflaw = "MA"
    elif typeofmodel == "MM-conc":
        KS = param_vec[2]
        KP = param_vec[3]
        F = (a*conc_st[0]/KS - b*conc_st[1]/KP) / (1 + conc_st[0]/KS + conc_st[1]/KP)
        S = (1 + conc_st[0]/KS + conc_st[1]/KP)
        A = (a-F)/(S*KS)
        B = (b+F)/(S*KP)
        typeoflaw = "MM"
    elif typeofmodel == "MM-tracer":
        KS = param_vec[2]
        KP = param_vec[3]
        F = (a*conc_st[0]/KS - b*conc_st[1]/KP) / (1 + conc_st[0]/KS + conc_st[1]/KP)
        S = (1 + conc_st[0]/KS + conc_st[1]/KP)
        A = a/(S*KS)
        B = b/(S*KP)
    else:
        print("ERROR in function run_homogenous_model:\n   The model %s does not exist."%(typeofmodel))
        
        ## Markers behave like mass action
        typeoflaw = "MA"
    return_vec = ()
    for task in task_list:
        ## param_vec
        if task == "param_vec":
            return_vec = return_vec + (param_vec,)
        ## jacobian
        elif task=="jacobian":
            if typeofmodel == "MA-conc":
                return_vec +=  (compute_jacobian(size,param_vec),)
            elif typeofmodel == "MM-conc":
                return_vec +=  (compute_jacobian_MM(conc_st,param_vec),)
            elif typeofmodel == "MM-tracer":
                return_vec +=  (compute_jacobian(len(conc_st),[A,B]*(size+1)),)
        ## relax_time
        elif task == "relax_time":
            return_vec +=  (1/(A+B - 2*np.sqrt(A*B)*np.cos(np.pi/(size+1))),)
        ## transit_time
        elif task == "transit_time_L1":
            tau = 1/(A+B - 2*np.sqrt(A*B)*np.cos(np.pi/(size+1)))
            T = compute_lifetime(compute_jacobian(size,[A,B]*(size+1)),tau=tau,norm_type=1)
            
            return_vec += (T,)
        elif task == "transit_time_L2":
            return_vec += (eval_transit(size,[A,B],2),)
        else:
            print("ERROR in function run_homogenous_model:\n   The task %s does not exist."%(task))
    if len(return_vec) is 1:
        return(return_vec[0])
    else:
        return(return_vec)
######################################################################


######################################################################
##                        find_tau_iterative                        ##
######################################################################
def find_tau_iterative(M,rtol=1E-3,max_iter=0,x0=[]):
    convergence = False
    L = len(M)    
    if type(M) is not np.ndarray:
        M = np.array(M)
    M = np.linalg.inv(M[1:-1,1:-1])
    if x0 == []:        
        x0 = np.apply_along_axis(lambda vec:np.sum(np.abs(vec)),1,M)
        x0[x0!=0] = np.random.lognormal(1,1, len(x0[x0!=0]))          
    else:
        x0 = x0[1:-1]
    x0 = x0/np.linalg.norm(x0)
    dx = np.ones(L)
    evs = np.ones(L)
    i = 0
    while True:
        i += 1
        dx = np.dot(M,x0) 
        negative_positions = np.squeeze(np.argwhere( np.abs(x0)>1E-10*np.max(np.abs(x0)) ))
        ev = dx[negative_positions]/x0[negative_positions]
        if np.std(ev)<rtol*np.abs(np.mean(ev)):
            convergence = True
            break
        elif max_iter!=0 and i>max_iter:
            break
        x0 = np.array(dx)
        x0 = x0/np.linalg.norm(x0)        
    if np.std(ev)>rtol*np.abs(np.mean(ev)):        
        print("\tThe algorithm did not converge -  max_iter=%d"%max_iter)
    return(convergence,np.mean(-ev),np.append(0,np.append(x0,0)))
######################################################################


######################################################################
##                       confidence_interval                        ##
######################################################################
def confidence_interval(vec,confidence=0.95):
    vec = np.sort(vec)
    L = len(vec)
    center = np.int(L/2)
    index = np.int((L-1)/2)
    remainder = (L-1)%2
    med = (1-remainder)*vec[index] + remainder*vec[index+1]
    index_under = np.int((L-1)*(1-confidence)/2)
    r_under = (L-1)*(1-confidence)/2 - index_under        
    interval_under = (1-r_under)*vec[index_under] + r_under*vec[index_under+1]

    index_over = np.int((L-1)*(1+confidence)/2)
    r_over = (L-1)*(1+confidence)/2 - index_over        
    interval_over = (1-r_over)*vec[index_over] + r_over*vec[index_over+1]

    #import pdb;pdb.set_trace()
    return(med,med-interval_under,interval_over-med)
    
######################################################################

