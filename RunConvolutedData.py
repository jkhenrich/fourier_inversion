# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 09:50:08 2015

@author: jkh
"""

from __future__ import division

import pylab
import numpy as np
from numpy import exp, pi, sqrt
import DataCreation as dc
import FourierTransform as ft


def run_the_transform (e_end, ne, sig):
    number_of_runs = 100
    e = np.linspace(-e_end, e_end, ne)
    t  = np.linspace((-pi*ne)/(4*e_end), (pi*ne)/(4*e_end), ne)
    point_means = np.zeros_like(e, dtype = 'complex')
    point_stds = np.zeros_like(e, dtype = 'complex')
    for i in range(number_of_runs):
        convoi = dc.add_noise(dc.convolve(e,sig))
        iftconvoi = ft.inversetransform(e, convoi, t)
        print iftconvoi
        iftconvoi = abs(iftconvoi)
        point_means +=  iftconvoi 
        point_stds += iftconvoi**2
    point_stds = (1 / (number_of_runs - 1))*(point_stds +(1/number_of_runs)*(point_means**2))
    point_means = (1 / number_of_runs)*point_means
    return t, point_means, point_stds
      
def plot_transform_run(returnedvectors):
    t, point_means, point_stds = returnedvectors
    pylab.errorbar(t, point_means, yerr=point_stds)
    
pylab.clf()        
plot_transform_run(run_the_transform(1.5, 2251, 0.01))
pylab.show()        