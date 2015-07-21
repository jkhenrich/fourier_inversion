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
import wsolve as ws


def run_the_transform_artdata (e_end, ne, sig, number_of_runs):
    number_of_runs = 10
    e = np.linspace(-e_end, e_end, ne)
    t  = np.linspace((-pi*ne)/(4*e_end), (pi*ne)/(4*e_end), ne)
    pretransformedvalues = np.ndarray(shape = (number_of_runs, ne), dtype=float)
    alliftvalues = np.ndarray(shape = (number_of_runs, ne), dtype=complex)
    point_means = np.zeros_like(e, dtype = 'complex')
    point_stds = np.zeros_like(e, dtype = 'complex')
    for i in range(number_of_runs):
        convoi = dc.add_noise(dc.convolve(e,sig))
#        convoi = addlinearbackground(0.55, 0.05, convoi, e)
        pretransformedvalues[i] = convoi
#        aproxx = np.append(e[:int(ne/3)], e[(2*int(ne/3)):])
#        aproxy = np.append(convoi[:int(ne/3)], convoi[(2*int(ne/3)):])
#        polymodel = ws.wpolyfit(aproxx, aproxy, degree=1)
#        valuestofft = convoi - polymodel(e)
#        iftconvoi = ft.inversetransform(e, valuestofft, t)
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
    return number_of_runs, alliftvalues, point_means, e, t, point_stds, pretransformedvalues#, allvalues
    
def run_the_transform_redata (e, t, pretransformedvalues, pretransformedres, number_of_runs): #make sure pretranssformed values is in columns
    ne = len(e)
    columnindices = np.arange(number_of_runs)
    indicesfordata = np.where(columnindices %2 == 1)[0]    
    alliftvalues = np.ndarray(shape = (len(indicesfordata), ne), dtype=complex)
    allresvalues = np.empty_like(alliftvalues)
    wobackgroundvalues = np.empty_like(alliftvalues)
    backgroundfunct = []
    portion = 0.1
    for i, xi in enumerate(indicesfordata):
        convoi = ft.choose_data_from_file(pretransformedvalues, xi)
        convoires = ft.choose_data_from_file(pretransformedres, xi)
        '''
        aproxx = np.append(e[:int(ne*portion)], e[-int(ne*portion):])
        aproxy = np.append(convoi[:int(ne*portion)], convoi[-int(ne*portion):])
        polymodel = ws.wpolyfit(aproxx, aproxy, degree=1)
        valuestofft = abs(convoi - polymodel(e))
        '''
        valuestofft = convoi
        wobackgroundvalues[i] = valuestofft
        iftconvoi = ft.inversetransform(e, valuestofft, t)
#        iftconvoi = abs(iftconvoi)
        alliftvalues[i] = iftconvoi
        '''
        aproxyres = np.append(convoires[:int(ne*portion)], convoires[-int(ne*portion):])
        polymodelres = ws.wpolyfit(aproxx, aproxyres, degree=1)
        restofft = abs(convoires - polymodelres(e))
        '''
        restofft = convoires
        iftconvoires = ft.inversetransform(e, restofft, t)
#        iftconvoires = abs(iftconvoires)
        allresvalues[i] = iftconvoires
#        backgroundfunct.append([polymodel, polymodelres])
    return indicesfordata, alliftvalues, allresvalues, backgroundfunct, wobackgroundvalues

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
    
def importrealdataerrors(allrows, rows, datafile):    
    if allrows == 0:
        x = len(rows)
        y = len(datafile[0,:])
        errors = [datafile[xi + 1] for i, xi in enumerate(rows)]
        print errors
        errors = np.array(errors)
    if allrows == 1:
        columnindices = np.arange(len(datafile[0,:]))
        rows = np.where(columnindices %2 == 0)[0]
        x = len(rows)
        y = len(datafile[0, :])
        errors = [datafile[xi + 1] for i, xi in enumerate(rows)]
        errors = np.array(errors)
    
def findIQT(t, E, sampledataset, resolutionset):
    Im = ft.inversetransform(E, sampledataset, t)
    RQT = ft.inversetransform(E, resolutionset, t)
    IQT = np.divide(abs(Im), abs(RQT))
#    IQT = np.divide(np.real(Im), np.real(RQT)) #To Do: Make sure real is supposed to be in here 
    return IQT

#ORIGINAL Function
def plotiftrealdata(e, data, resolution, transformeddata, resolutiont, IQT, t, nt):
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
    pylab.plot(t, IQT, '-o', label = 'IQT')
#    pylab.plot(t, IQT/IQT[len(t)//2], '-o', label = 'IQT')
    pylab.legend()

    
def NOTplotiftrealdata(e, data, resolution, transformeddata, resolutiont, wobackground, IQT, t, nt):
    e = abs(e)
    pylab.subplot(1,3,1)
    pylab.plot(e, data, '-o', label='Original Data')
    pylab.plot(e, resolution, '-o', label='Original Resolution')
    print len(e), len(wobackground)
    pylab.plot(e, wobackground, '-o', label = 'wobackground')
    pylab.legend()
    pylab.subplot(1,3,2)
    pylab.plot(t, abs(transformeddata), '-o',
                   label='Inverse Fourier Transform')
    pylab.plot(t, abs(resolutiont), '-o', label = 'Resolution')
    pylab.xlabel('Time')
    pylab.ylabel('Intensity (counts/s)')
    pylab.legend()
    pylab.subplot(1,3,3)
    pylab.plot(t, IQT, '-o', label = 'IQT')
#    pylab.plot(t, IQT/IQT[len(t)//2], '-o', label = 'IQT')
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
    number_of_runs, data, point_means, e, t, point_stds, pretransformedvalues  = run_the_transform_artdata(1.5, 2251, 0.01, 100)
    #demo_correlation_btw_neighbors(number_of_runs, data, point_means)
    stdscalc = demo_standard_deviations(number_of_runs, data, point_means)
    print 'pretransformedvalues'
    print pretransformedvalues
    pretransformedvaluesmeans = np.sum(pretransformedvalues, axis=0) / number_of_runs
    print pretransformedvaluesmeans
    print stdscalc
    print point_stds
    print len(stdscalc)
    pylab.clf()
    pylab.subplot(2,1,1)
    pylab.plot(e, pretransformedvaluesmeans, label='Original Data')
    pylab.subplot(2,1,2)
    plot_transform_run(number_of_runs, t, point_means, point_stds)
    pylab.show()
#    dc.export_data_csv(data, 'originaldistvalues3')
#    dc.export_data_csv(np.vstack([np.float_(point_stds), np.float_(stdscalc)]), 'stdscalc2ways_data3')

def demo_createddata():
    number_of_runs, data, point_means, e, t, point_stds, pretransformedvalues = run_the_transform_artdata(1, 1501, 0.01, 100)
    print len(point_stds)
    assert np.all(abs(np.imag(data)) < 1e-12), 'Cant cast data to float, maxcomplex: %g'%max(abs(np.imag(data)))
    assert np.all(abs(np.imag(point_means)) < 1e-12), 'Cant cast means to float, maxcomplex: %g'%max(abs(np.imag(point_means)))
    assert np.all(abs(np.imag(point_stds)) < 1e-12), 'Cant cast stds to float, maxcomplex: %g'%max(abs(np.imag(point_stds)))
    dc.export_data_csv(np.real(data).transpose(), 'origvalues2015714_02')
    dc.export_data_csv(np.vstack([t, np.real(point_means), np.real(point_stds)]).transpose(), 'meanstd2015714_02')

def createresolutiondata(e, t, sig):
    gaussian = dc.create_gaussian_norm(e, sig)
    gaussianwbackground = addlinearbackground(0.55, 0.05, gaussian, e)
    gaussianwnoise = dc.add_noise(gaussianwbackground)
    ftgaussianwnoise = ft.inversetransform(e, gaussianwnoise, t)
    return ftgaussianwnoise
    
def addlinearbackground(yintercept, slope, data, e): #make sure data us a 1d horizontal array
    datawbackground = data + (slope*e + yintercept)
    return datawbackground
    
    
def graphingofoutput():#Make sure to know the sigmas of the two data sets
    datafile1 = r'D:\Users\jkh\Documents\Python Scripts\fourier_inversion\DataFiles\Created data-differing energy and resolution runs\meanstd2015710_04.csv'
    data1st = np.loadtxt(datafile1, delimiter=',')
    edata1 = np.linspace(-4.5, 4.5, 2251)
    sigmadata1 =0.01
    datafile2 = r'D:\Users\jkh\Documents\Python Scripts\fourier_inversion\DataFiles\Created data-differing energy and resolution runs\meanstd2015710_04.csv'
    datasecond = np.loadtxt(datafile2, delimiter=',')
    edata2 =np.linspace(-4.5, 4.5, 6753)
    sigmadata2 = 0.01
    tdata1 = data1st[:, 0]
    print len(tdata1)
    tdata2 = datasecond[:, 0]
    print len(tdata2)
    meandata1 = data1st[:, 1]
    meandata2 = datasecond[:, 1]
    stdsdata1 = data1st[:, 2]
    stdsdata2 = datasecond[:, 2]
    gauswbackground = 0
    pylab.clf()
    if gauswbackground == 0:
        pylab.subplot(2,1,1)
        pylab.plot(tdata1, meandata1, '-o',
                   label='Inverse Fourier Transform')
        pylab.plot(tdata2, meandata2, '-o',
                   label='Inverse Fourier Transform')             
#        pylab.plot(tdata1, abs(dc.add_noise(ft.fourier_transform_gaussian(tdata1, sigmadata1))), '-o', label = 'Resolution1')
#        pylab.plot(tdata1, abs(dc.add_noise(ft.fourier_transform_gaussian(tdata2, sigmadata2))), '-o', label = 'Resolution2')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
#        pylab.yscale('log')
        pylab.legend()
#        pylab.plot(tdata1, meandata1, '-o', label='means of data1')
 #       pylab.plot(tdata2, meandata2, '-o', label='means of data2')
        pylab.subplot(2, 1,2)
        pylab.ylim(0,5)
        pylab.plot(tdata1, np.divide(abs(meandata1), abs(ft.fourier_transform_gaussian(tdata1, sigmadata1*2))), '-o', label='means of data1')
        pylab.plot(tdata2, np.divide(abs(meandata2), abs(ft.fourier_transform_gaussian(tdata2, sigmadata2))), '-o', label='means of data2')
#        pylab.plot(tdata1, abs(meandata1) / abs(ft.fourier_transform_gaussian(tdata1, sigmadata1)), '-o', label='means of data1')
#        pylab.plot(tdata2, abs(meandata2) / abs(ft.fourier_transform_gaussian(tdata2, sigmadata2)), '-o', label='means of data2')    
#        pylab.yscale('log')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.legend()
        
    if gauswbackground == 1:
        pylab.subplot(2,1,1)
        pylab.plot(tdata1, meandata1, '-o', label='means of data1')
        pylab.plot(tdata2, meandata2, '-o', label='means of data2')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.legend()
        pylab.subplot(2, 1,2)
        pylab.plot(tdata1, np.divide(abs(meandata1), createresolutiondata(edata1, tdata1, sigmadata1)), '-o', label='means of data1')
        pylab.plot(tdata2, np.divide(abs(meandata2), abs(ft.fourier_transform_gaussian(tdata2, sigmadata2))), '-o', label='means of data2')
        pylab.ylim(0,5)    
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
    """
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

def demo_realdata_nocorr4background():
#    datafile1 = str(raw_input('Enter the name of data file to load:'))
    datafile1 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l6_T240K.txt'
#    datafile2 = str(raw_input('Enter the name of data file to load:'))
    datafile2 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l9_T240K.txt'
#    resfile1 = str(raw_input('Enter the name of data file to load:'))
    resfile1 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l6_T10K.txt'
#    resfile2 = str(raw_input('Enter the name of data file to load:'))
    resfile2 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l9_T10K.txt'
    importedsampledata1 = np.loadtxt(datafile1, dtype='float')
    e1 = (ft.choose_data_from_file(importedsampledata1, 0)) #/ (10**3)
    I1 = ft.choose_data_from_file(importedsampledata1, 5)
    nt1 = len(e1)
    t1 = np.linspace(-2*(pi * nt1) / ((4*abs(e1[0]))), 2*(pi * nt1) / ((4*abs(e1[0]))), nt1)
    t1 = ft_calc_time(e1)
    print len(t1)
    I1transformed = ft.inversetransform(e1, I1, t1)
    resolution1 = ft.choose_data_from_file(np.loadtxt(resfile1, dtype='float'), 5)
    restransformed = ft.inversetransform(e1, resolution1, t1)
    print len(restransformed)
    IQT1 = findIQT(t1, e1, I1, resolution1)
    
    importedsampledata2 = np.loadtxt(datafile2, dtype='float')
    e2 = (ft.choose_data_from_file(importedsampledata2, 0)) #/ (10**3)
    print len(e2)
    I2 = ft.choose_data_from_file(importedsampledata2, 5)
    nt2 = len(e2)
    t2 = np.linspace(-2*(pi * nt2) / ((4*abs(e2[0]))), 2*(pi * nt2) / ((4*abs(e2[0]))), nt2)
    I2transformed = ft.inversetransform(e2, I2, t2)
    resolution2 = ft.choose_data_from_file(np.loadtxt(resfile2, dtype='float'), 5)
    print len(resolution2)
    restransformed2 = ft.inversetransform(e2, resolution2, t2)
    IQT2 = findIQT(t2, e2, I2, resolution2)
    t2 = ft_calc_time(e2)
    print t1.shape, IQT1.shape
    print t1
    pylab.clf()
    plotiftrealdata(e1, I1, resolution1, I1transformed, restransformed, IQT1, t1, nt1)
    plotiftrealdata(e2, I2, resolution2, I2transformed, restransformed2, IQT2, t2, nt2)
    print I1
#    print I1transformed
#    print restransformed
    print IQT1
    print e1.size
    print I1.size
#    importedresdata = np.loadtxt(str(raw_input('Enter the name of data file to load:')),
#                                 dtype='float', delimiter=',')
#    IQT = findIQT()

def demo_realdata_corr4background():
#    datafile1 = str(raw_input('Enter the name of data file to load:'))
    datafile1 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l6_T240K.txt'
#    datafile2 = str(raw_input('Enter the name of data file to load:'))
    datafile2 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l9_T240K.txt'
#    resfile1 = str(raw_input('Enter the name of data file to load:'))
    resfile1 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l6_T10K.txt'
#    resfile2 = str(raw_input('Enter the name of data file to load:'))
    resfile2 = r'D:\Users\jkh\Documents\reduced_data_fromAntonio\fukuoka_l9_T10K.txt'
    importedsampledata1 = np.loadtxt(datafile1, dtype='float')
    e1 = (ft.choose_data_from_file(importedsampledata1, 0)) #/ (10**3) for better viewing
    print e1
    nt1 = len(e1)
    t1 = (2*pi*e1) / 4.1657
    print t1
    t1 = np.linspace(-2*(pi * nt1) / ((4*abs(t1[0]))), 2*(pi * nt1) / ((4*abs(t1[0]))), nt1)
    print t1
    t1_1 = ft_calc_time(e1)
    resolution1 = np.loadtxt(resfile1, dtype='float')
    indicesfordata, alliftvalues, allresvalues, backgroundfunct, wobackground = run_the_transform_redata (e1, t1, importedsampledata1, resolution1, len(importedsampledata1[0]))
    indicesfordata1, alliftvalues1, allresvalues1, backgroundfunct1, wobackground1 = run_the_transform_redata (e1, t1_1, importedsampledata1, resolution1, len(importedsampledata1[0]))
    allIQT = [findIQT(t1, e1, importedsampledata1.transpose()[xi], resolution1.transpose()[xi]) for i,xi in enumerate(indicesfordata)]
    #t1 = ft_calc_time(e1)
    importedsampledata2 = np.loadtxt(datafile2, dtype='float')
    e2 = (ft.choose_data_from_file(importedsampledata2, 0)) / (10**3)
    nt2 = len(e2)
    t2 = np.linspace(-2*(pi * nt2) / ((4*abs(e2[0]))), 2*(pi * nt2) / ((4*abs(e2[0]))), nt2)
    resolution2 = np.loadtxt(resfile2, dtype='float')
    indicesfordata2, alliftvalues2, allresvalues2, backgroundfunct2, wobackground2 = run_the_transform_redata (e2, t2, importedsampledata2, resolution2, len(importedsampledata2[0]))
#    t2 = ft_calc_time(e2)
    allIQT2 = [findIQT(t2, e2, importedsampledata2.transpose()[xi], resolution2.transpose()[xi]) for i,xi in enumerate(indicesfordata2)]
    pylab.clf()
    NOTplotiftrealdata(e1, importedsampledata1.transpose()[5], resolution1.transpose()[5], alliftvalues[2], allresvalues[2], wobackground[2], abs(alliftvalues[2])/abs(allresvalues[2]), t1, nt1)
    NOTplotiftrealdata(e1, importedsampledata1.transpose()[5], resolution1.transpose()[5], alliftvalues1[2], allresvalues1[2], wobackground1[2], abs(alliftvalues[2])/abs(allresvalues[2]), t1_1, nt1)
#    NOTplotiftrealdata(e2, importedsampledata2.transpose()[5], resolution2.transpose()[5], alliftvalues2[2], allresvalues2[2], wobackground2[2], abs(alliftvalues2[2])/(allresvalues2[2]), t2, nt2)


#demo_createddata()
demo_realdata_corr4background()
#demo_stdtest()
#testofoutput()
#dc.export_data_csv(np.vstack([np.linspace((-pi*751)/(4*.5), (pi*751)/(4*.5), 751), np.linspace((-pi*1501)/(4*1), (pi*1501)/(4*1),1501), np.linspace((-pi*2251)/(4*1.5), (pi*2251)/(4*1.5), 2251), np.linspace((-pi*5253)/(4*3.5), (pi*5253)/(4*3.5), 5253)]).transpose(), 't201579_01')
#graphingofoutput()
print 'orig'
