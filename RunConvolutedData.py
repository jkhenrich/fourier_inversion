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
        point_stds += np.multiply(iftconvoi, iftconvoi)
    #point_stds = sqrt((1 / (number_of_runs))*(point_stds + (1/number_of_runs)*(point_means**2)))
    #point_stds = sqrt((1 / (number_of_runs - 1))*(point_stds +(1/number_of_runs)*(point_means**2)))
    point_means = (1 / number_of_runs)*point_means
    point_stds = [find_stds(alliftvalues.transpose()[i]-point_means[i], number_of_runs) for i in np.arange(ne)]
    point_stds = np.array(point_stds, dtype='complex')    
    print len(point_stds)
    return number_of_runs, alliftvalues, point_means, t, point_stds#, allvalues

def find_covariance(firstdesrowmmean, secdesrowmmean):
    tosum = firstdesrowmmean * secdesrowmmean
    covariance = np.sum(tosum) / (number_of_runs-1)
    return covariance

def find_stds(rowmmean, number_of_runs):
    firstrowsqr = np.multiply(rowmmean, rowmmean)
    stdfirstrow = sqrt(np.sum(firstrowsqr)/(number_of_runs-1))
    return stdfirstrow


def find_correlation2rows(i, number_of_runs, data, point_means):
    firstdesrow = data[i]
    seconddesrow = data[i+1]
    firstmean = point_means[i]
    secondmean = point_means[i+1]
    #To calculate the covariance
    firstdesrowmmean = firstdesrow - firstmean
    secdesrowmmean = seconddesrow - secondmean
    covariance = find_covariance(firstdesrowmmean, secdesrowmmean)
    #To calculate the standarddeviations
    stdfirstrow = find_stds(firstdesrowmmean, number_of_runs)
    stdsecondrow = find_stds(secdesrowmmean, number_of_runs)
    correlation = covariance / (stdfirstrow*stdsecondrow)
    return correlation
      
def plot_transform_run(number_of_runs, t, point_means, point_stds):
#    pylab.plot(t, point_stds, '.')
#    for i in np.arange(number_of_runs):
#        pylab.plot(t, allvalues[i], '.')
    pylab.errorbar(t, point_means, yerr=point_stds)
    
    
def demo_correlation_btw_neighbors(number_of_runs, data, point_means):
    correlations = np.empty_like(np.arange(len(point_means)-1))
    #for i in np.arange(len(point_means)-1):  
    correlations= [find_correlation2rows(i, number_of_runs, data.transpose(), point_means) for i in np.arange(len(point_means)-1)]
    print correlations

def demo_standard_deviations(number_of_runs, data, point_means):
    standarddeviations = np.empty_like(point_means)
    standarddeviations = [find_stds((data.transpose()[i] - point_means[i]), number_of_runs) for i in np.arange(len(point_means))]
    return np.array(standarddeviations)

def demo_fakedata():
    number_of_runs, data, point_means, t, point_stds = run_the_transform(1.5, 2251, 0.01)
    #demo_correlation_btw_neighbors(number_of_runs, data, point_means)
    stdscalc = demo_standard_deviations(number_of_runs, data, point_means)
    print stdscalc
    print len(stdscalc)
    pylab.clf()
    plot_transform_run(number_of_runs, t, point_means, point_stds)
    pylab.show()
    dc.export_data_csv(data, 'originaldistvalues3')
    dc.export_data_csv(np.vstack([np.float_(point_stds), np.float_(stdscalc)]), 'stdscalc2ways_data3')
    
def demo_realdata():
    importeddata = ft.upload_data()

demo_fakedata()
   