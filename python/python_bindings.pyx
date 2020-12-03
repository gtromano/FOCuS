# distutils: language = c++
from python_bindings cimport Interval, Quadratic, Info, FOCuS_step_sim, FOCuS_step;

from math import inf
import numpy as np
import matplotlib.pyplot as plt
import time
from IPython.display import clear_output

cdef read_an_interval(Interval Inter):
    return round(Inter.l, 2), round(Inter.u, 2)

cdef read_a_quadratic(Quadratic q):
    a, b, c = round(q.a, 2), round(q.b, 2), round(q.c, 2)
    intervals_list = []
    
    while not q.ints.empty():
        intervals_list.append(read_an_interval(q.ints.front()))
        q.ints.pop_front()
    
    return a, b, c, intervals_list

def evaluate_quadratic(X, a=-1, b=2, c=1):        
    return np.maximum(a*X**2 + b*X + c, 0)

def scan(X, threshold=25, plot=False, pre_change_mean=0):
    #input X, a 1d numpy array of data, significance threshold, pre-change mean (can be None)
    #output the start of change and point at which it was detected at given threshold
    #also optionally plots quadratics in a Jupyter notebook
    
    changepoint = None
    info = PyInfo()
    S_max = np.empty(X.shape)
    stopping_time=len(X)
    
    for t in range(len(X)):
        
        if pre_change_mean is not None:
            info.focus_step_sim(X[t]-pre_change_mean)
        else:
            info.focus_step(X[t])
        
        if plot:
            clear_output(wait=True)
            fig = info.plot_quadratics(threshold=threshold, pre_change_mean=pre_change_mean)
            plt.show()
        
        S_max[t] = info.global_max()
        
        if info.global_max()>threshold:
            changepoint = t+info.time_offset()
            stopping_time = t
            break
    
    return(changepoint, stopping_time, S_max[0:stopping_time+1])

cdef class PyInterval:    
    cdef Interval cinterval
    
    def __cinit__(self, l=-inf, u=inf):
        self.cinterval.l = l
        self.cinterval.u = u
    
    def __repr__(self):
        return f'Interval: [{self.cinterval.l}, {self.cinterval.u}]'
    
    def contains(self, point):
        return self.cinterval.l < point < self.cinterval.u
    
cdef class PyQuadratic:    
    cdef Quadratic cquadratic
    
    def __cinit__(self, a=0, b=0, c=0, interval_list=[(-inf, inf)]):
        self.cquadratic.a = a
        self.cquadratic.b = b
        self.cquadratic.c = c
        
        for interval in interval_list:
            i = PyInterval(interval[0], interval[1])
            self.cquadratic.ints.push_back(i.cinterval)
        
    def __repr__(self):        
        a, b, c, intervals_list = read_a_quadratic(self.cquadratic)        
        return f'Quadratic: {a}x^2+{b}x+{c} defined on intervals ' + repr(intervals_list)
        
cdef class PyInfo:    
    cdef Info cinfo
    
    def __cinit__(self):
        self.cinfo.Q0 = PyQuadratic().cquadratic
        self.cinfo.global_max = 0
        self.cinfo.time_offset = 0
        
    def focus_step_sim(self, double new_point):
        self.cinfo = FOCuS_step_sim(self.cinfo, new_point)
    
    def focus_step(self, double new_point):
        self.cinfo = FOCuS_step(self.cinfo, new_point)
    
    def global_max(self):
        return self.cinfo.global_max
    
    def time_offset(self):
        return int(self.cinfo.time_offset)
    
    def __repr__(self):
        
        out_message = 'FOCuS Info object. \n'
        
        #Q0
        i = 0
        a, b, c, intervals_list = read_a_quadratic(self.cinfo.Q0) 
        out_message += f'Quadratic Q{i}: {a}x^2+{b}x+{c} \n'
        
        #Q1 list of quadratics
        cdef Info temp_info = self.cinfo
        while not temp_info.Q1.empty():
            i += 1
            a, b, c, intervals_list = read_a_quadratic(temp_info.Q1.front())
            out_message += f'Quadratic Q{i}: {a}x^2+{b}x+{c} defined on intervals ' + repr(intervals_list) + '\n'
            temp_info.Q1.pop_front()                        
        
        out_message += f'Global max: {round(self.cinfo.global_max, 2)}'
        
        return out_message
    
    def plot_quadratics(self, pre_change_mean=0, threshold=None):
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("FOCuS Step")
        ax.set_xlabel("$\mu$")
        ax.set_ylabel("$S_{t}(\mu)$", rotation=0)
        
        if pre_change_mean is None:
            a, b, c, intervals_list = read_a_quadratic(self.cinfo.Q0)
            offset = -b/(2*a) #turning point of null quadratic (most likely estimate)
            pre_change_mean=0 #needed to ignore for plotting
        else:
            offset = pre_change_mean
        
        Z = np.linspace(-5, 5, 100) + offset #the x-axis for the plot

        cdef Info temp_info = self.cinfo
        i = 0
        while not temp_info.Q1.empty():
            i += 1
            a, b, c, intervals_list = read_a_quadratic(temp_info.Q1.front())
            ax.plot(Z, evaluate_quadratic(Z - pre_change_mean, a, b, c), color='C0', label=f'Q{i}')
            temp_info.Q1.pop_front() 
        
        if threshold is not None:
            ax.axhline(threshold, color='C1', label="T")

        return fig