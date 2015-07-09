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
    alliftvalues = np.ndarray(shape = (number_of_runs, ne), dtype=complex)
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
    return number_of_runs, alliftvalues, point_means, t, point_stds#, allvalues

def find_covariance(firstdesrowmmean, secdesrowmmean, number_of_runs):
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
    
def findIQT(t, E, sampledataset, resolutionset):
    Im = ft.inversetransform(E, sampledataset, t)
    RQT = ft.inversetransform(E, resolutionset, t)
    IQT = np.divide(abs(Im), abs(RQT))
#    IQT = np.divide(np.real(Im), np.real(RQT)) #To Do: Make sure real is supposed to be in here 
    return IQT

def plotiftrealdata(e, data, resolution, transformeddata, resolutiont, IQT, t, nt):
    pylab.clf()
    pylab.subplot(1,3,1)
    pylab.plot(e, data, '-o', label='Original Data')
    pylab.plot(e, resolution, '-o', label='Original Resolution')
    pylab.legend()
    pylab.subplot(1,3,2)
    pylab.plot(t, abs(transformeddata), '-o',
                   label='Inverse Fourier Transform')
    pylab.plot(t, abs(resolutiont), '-o', label = 'Resolution')
    pylab.xlabel('Time')
    pylab.ylabel('Intensity (counts/s)')
    pylab.legend()
    pylab.subplot(1,3,3)
    pylab.plot(t, IQT/IQT[len(t)//2], '-o', label = 'IQT')
    pylab.legend()

        
    
def demo_correlation_btw_neighbors(number_of_runs, data, point_means):
    correlations = np.empty_like(np.arange(len(point_means)-1))
    #for i in np.arange(len(point_means)-1):  
    correlations= [find_correlation2rows(i, number_of_runs, data.transpose(), point_means) for i in np.arange(len(point_means)-1)]
    print correlations

def demo_standard_deviations(number_of_runs, data, point_means):
    standarddeviations = np.empty_like(point_means)
    standarddeviations = [find_stds((data.transpose()[i] - point_means[i]), number_of_runs) for i in np.arange(len(point_means))]
    return np.array(standarddeviations)

def demo_stdtest():
    number_of_runs, data, point_means, t, point_stds = run_the_transform(1.5, 2251, 0.01)
    #demo_correlation_btw_neighbors(number_of_runs, data, point_means)
    stdscalc = demo_standard_deviations(number_of_runs, data, point_means)
    print stdscalc
    print point_stds
    print len(stdscalc)
    pylab.clf()
    plot_transform_run(number_of_runs, t, point_means, point_stds)
    pylab.show()
    dc.export_data_csv(data, 'originaldistvalues3')
    dc.export_data_csv(np.vstack([np.float_(point_stds), np.float_(stdscalc)]), 'stdscalc2ways_data3')

def demo_createddata():
    number_of_runs, data, point_means, t, point_stds = run_the_transform(1, 1501, 0.05)
    print len(point_stds)
    assert np.all(abs(np.imag(data)) < 1e-12), 'Cant cast data to float, maxcomplex: %g'%max(abs(np.imag(data)))
    assert np.all(abs(np.imag(point_means)) < 1e-12), 'Cant cast means to float, maxcomplex: %g'%max(abs(np.imag(point_means)))
    assert np.all(abs(np.imag(point_stds)) < 1e-12), 'Cant cast stds to float, maxcomplex: %g'%max(abs(np.imag(point_stds)))
    dc.export_data_csv(np.real(data).transpose(), 'origvalues201578_10')
    dc.export_data_csv(np.vstack([t, np.real(point_means), np.real(point_stds)]).transpose(), 'meanstd201578_10')
    
def graphingofoutput():
    datafile1 = r'D:\Users\jkh\Documents\Python Scripts\fourier_inversion\DataFiles\Created data-differing energy and resolution runs\meanstd201578_06.csv'
    data1st = np.loadtxt(datafile1, delimiter=',', skiprows = 2)
    datafile2 = r'D:\Users\jkh\Documents\Python Scripts\fourier_inversion\DataFiles\Created data-differing energy and resolution runs\meanstd201578_08.csv'
    datasecond = np.loadtxt(datafile2, delimiter=',', skiprows = 2)
    tdata1 = data1st[:, 0]
    tdata2 = datasecond[:, 0]
    meandata1 = data1st[:, 1]
    meandata2 = datasecond[:, 1]
    stdsdata1 = data1st[:, 2]
    stdsdata2 = datasecond[:, 2]
    pylab.clf()
    pylab.plot(tdata1, meandata1, '-o', label='means of data1')
    pylab.plot(tdata2, meandata2, '-o', label='means of data2')
    pylab.xlabel('Time')
    pylab.ylabel('Intensity (counts/s)')
    pylab.legend()

def ft_calc_time(energy):
    """THIS METHOD INVERTS ENERGY <-> TIME FOR FFT

;h = 4.135667516(91)Ã—10âˆ’15 eV.s
;  = 4.135667516(91)Ã—10âˆ’12 meV.s
;  = 4.135667516(91)Ã—10âˆ’3 meV.ns
;  = 4.135667516(91) meV.ps

; #PC CONVERT INPUT ENERGY TO FREQUENCY (ps^-1)
energy = etemp/4.136

; #PC CREATE TIME ARRAY
nenergy = n_elements(energy)
de = energy[1]-energy[0]   ;SET EVEN dE SPACING
n21 = nenergy/2 + 1
; The array of subscripts:
t = indgen(nenergy)
; Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):

; #PC PROVIDE APPROPRIATE SHIFT TO MATCH THE ONE NEEDED FOR FFT
if (nenergy mod 2) eq 0 then val = 2 else val = 1
t[n21] = n21 - nenergy + findgen(n21-val)
; Compute T0 frequency:
t = shift((t/(nenergy*de)),-n21)
return,t"""

    energy = energy/4.136
    nenergy = len(energy)
    de = energy[1]-energy[0]
    return (np.arange(-nenergy/2., nenergy/2.)+1./nenergy)/(nenergy*de)
    _ = """
    n21 = nenergy/2 + 1
    t = np.arange(nenergy)
    #Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):
    #PC PROVIDE APPROPRIATE SHIFT TO MATCH THE ONE NEEDED FOR FFT
    val = 2 if (nenergy % 2) == 0 else 1
    t[n21] = n21 - nenergy + np.arange(n21-val)
    # Compute T0 frequency:
    t = np.fft.fftshift((t/(nenergy*de)))
    print t
    return t
    """

def demo_realdata():
#    datafile = str(raw_input('Enter the name of data file to load:'))
    datafile = r'D:\Users\jkh\Documents\HFBS Data\Reduced Data\20150628_04_pc220Kdyn.txt'
#    resfile = str(raw_input('Enter the name of data file to load:'))
    resfile = r'D:\Users\jkh\Documents\HFBS Data\Reduced Data\20150629_01_pc10Kdyn.txt'
    importedsampledata = np.loadtxt(datafile, dtype='float')
    e = (ft.choose_data_from_file(importedsampledata, 0)) / (10**3)
    I1 = ft.choose_data_from_file(importedsampledata, 5)
    nt = len(e)
    t = np.linspace(-2*(pi * nt) / ((4*abs(e[0]))), 2*(pi * nt) / ((4*abs(e[0]))), nt)
    I1transformed = ft.inversetransform(e, I1, t)
    resolution = ft.choose_data_from_file(np.loadtxt(resfile, dtype='float'), 5)
    restransformed = ft.inversetransform(e, resolution, t)
    IQT1 = findIQT(t, e, I1, resolution)
    t = ft_calc_time(e)
    print t.shape, IQT1.shape
    print t
    plotiftrealdata(e, I1, resolution, I1transformed, restransformed, IQT1, t, nt)
    print I1
#    print I1transformed
#    print restransformed
    print IQT1
    print e.size
    print I1.size
#    importedresdata = np.loadtxt(str(raw_input('Enter the name of data file to load:')),
#                                 dtype='float', delimiter=',')
#    IQT = findIQT()

#demo_createddata()
#demo_realdata()
#testofoutput()
#dc.export_data_csv(np.vstack([np.linspace((-pi*751)/(4*.5), (pi*751)/(4*.5), 751), np.linspace((-pi*1501)/(4*1), (pi*1501)/(4*1),1501), np.linspace((-pi*2251)/(4*1.5), (pi*2251)/(4*1.5), 2251), np.linspace((-pi*5253)/(4*3.5), (pi*5253)/(4*3.5), 5253)]).transpose(), 't201579_01')
graphingofoutput()
print 'orig'
