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
    number_of_runs = 40
    e = np.linspace(-e_end, e_end, ne)
    t  = np.linspace((-pi*ne)/(4*e_end), (pi*ne)/(4*e_end), ne)
    pretransformedvalues = np.ndarray(shape = (number_of_runs, ne), dtype=float)
    alliftvalues = np.ndarray(shape = (number_of_runs, ne), dtype=float)
    point_means = np.zeros_like(e, dtype = 'complex')
    point_stds = np.zeros_like(e, dtype = 'complex')
    for i in range(number_of_runs):
        convoi = dc.add_noise(dc.convolve(e,sig))
        pretransformedvalues[i] = convoi
        #pylab.plot(e, convoi, '.-')
        iftconvoi = ft.inversetransform(e, convoi, t)
        iftconvoi = abs(iftconvoi)
        alliftvalues[i] = iftconvoi
        pylab.plot(t, iftconvoi, '.')
#        pylab.plot(t, iftconvoi / iftconvoi[(ne-1)/2], '.') #plottingline
        point_means +=  iftconvoi 
        point_stds += iftconvoi**2
    #point_stds = sqrt((1 / (number_of_runs - 1))*(point_stds +(1/number_of_runs)*(point_means**2)))
    point_stds = (1 / (number_of_runs - 1))*(point_stds +(1/number_of_runs)*(point_means**2))
    point_means = (1 / number_of_runs)*point_means
    return number_of_runs, alliftvalues, point_means#t, point_means, point_stds#, allvalues

def find_correlation2rows(number_of_runs, data, point_means):
    firstdesrow = data[2]
    seconddesrow = data[3]
    firstmean = point_means[2]
    secondmean = point_means[3]
    #To calculate the covariance
    firstdesrowmmean = firstdesrow - firstmean
    secdesrowmmean = seconddesrow - secondmean
    tosum = firstdesrowmmean * secdesrowmmean
    covariance = np.sum(tosum) / (number_of_runs-1)
    #To calculate the standarddeviations
    firstrowsqr = firstdesrow * firstdesrow
    secondrowsqr = seconddesrow * seconddesrow
    stdfirstrow = sqrt(np.sum(firstrowsqr)/(number_of_runs-1))
    stdsecondrow = sqrt(np.sum(secondrowsqr)/(number_of_runs-1))
    correlation = covariance / (stdfirstrow*stdsecondrow)
    return correlation
      
def plot_transform_run(returnedvectors):
    t, point_means, point_stds = returnedvectors
    number_of_runs = 100
#    pylab.plot(t, point_stds, '.')
#    for i in np.arange(number_of_runs):
#        pylab.plot(t, allvalues[i], '.')
    pylab.errorbar(t, point_means, yerr=point_stds)
    
#pylab.clf()
#run_the_transform(1.5, 2251, 0.01)
#plot_transform_run(run_the_transform(1.5, 2251, 0.01))
#pylab.show()
#dc.export_data_csv(run_the_transform(1.5, 2251, 0.01), 'originaldistvalues1')
number_of_runs, data, point_means = run_the_transform(1.5, 2251, 0.01)    
print find_correlation2rows(number_of_runs, data, point_means)