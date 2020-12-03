from libcpp.list cimport list as clist


cdef extern from "../src/FOCuS.cpp":
    cdef Info FOCuS_step_sim(Info info, double new_point)
    cdef Info FOCuS_step(Info info, double new_point)
    

cdef extern from "../src/quadratic.cpp":
    pass #Cython needs to read this file

cdef extern from "../src/quadratic.h":
    cdef cppclass Interval:
        double l, u
        
    cdef cppclass Quadratic:
        double a, b, c
        clist[Interval] ints
    
cdef extern from "../src/FOCuS.h":
    
    cdef cppclass Info:
        Quadratic Q0
        clist[Quadratic] Q1
        double global_max
        double time_offset
    
