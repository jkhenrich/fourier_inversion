# -*- coding: utf-8 -*-
"""
Created on Fri Jun 05 11:45:39 2015

@author: jkh
"""


from __future__ import division

import pylab
import numpy as np
from numpy import exp, pi


def upload_data():
    """
    upload_data takes data from previously exported csv file, and uploads it
    to be used in future code.
    """
    newdata = np.loadtxt('Trial4.csv', dtype='float', delimiter=',')
#    newdata = np.loadtxt(str(raw_input('Enter the name of data file to load:')),
#                         dtype='float', delimiter=',')
    print newdata.shape
    return newdata


def choose_data_to_transform(data): #data must be a e.size x 4
    """
    choose_data_to_transform takes a multi-dimensional array and reports the
    nth column of the array for future use.
    """
    extracteddata = data
    extracteddata = extracteddata.T
    return extracteddata[0, :]
#    return extracteddata[int(raw_input('Column of Desired Data: ')),:]


def fourier_transform_gaussian(t, sig):
    """
    fourier_transform_gaussian creates an array containing points on a Normal
    Gaussian curve with standard deviation, sigma, and mu = 0 represents the
    Fourier Transform from the same sigma from a Normal Gaussian of
    the form exp(-e**2 / (2 * sig**2)) / (sig * (sqrt(2 * pi))).
    """
    FRQe = exp(-(t**2) * (sig**2) / 2)
    print FRQe
    return FRQe #returns an array of points on the Gaussian curve


def fourier_transform_lorentzian(t):
    """
    fourier_transform_lorentzian creates an array containing points on an
    expontential curve with the same tau as the Lorentzian whose Fourier
    transformation it represents.
    """
    tau = 50
    lft = t.copy()
    lft = exp(-(abs(lft)) / (2*tau))
    print lft
    return lft


def transform_convolution(t, sig):
    """
    tranform_convolution creates an array containing multiplication of the
    points of the fourier transformed gaussian distribution, and fourier
    transformed lorentzian distribution.
    This should be equivalent to the fourier tranform of the convoluted
    gaussian and lorentzian distributions.
    """
    convotrans = fourier_transform_gaussian(t, sig) * fourier_transform_lorentzian(t)
    print convotrans
    return convotrans



def inversetransform(e, data, t):
    """
    inversetransform takes a discrete data set and computes its inverse fourier
    transform over a given time interval. This is accomplished using the
    trapoziodal approximation of the integral.
    """
    initialarraysize = data.size
    indices = np.arange(initialarraysize)
    iftdata = np.empty_like(data, dtype='complex')
    print iftdata.dtype
    for index in indices:
        transformedvalue = np.trapz(data* exp(1j*e*t[index]), x=e)
        #np.put(iftdata, index, transformedvalue)
        iftdata[index] = transformedvalue
    print iftdata
    return iftdata


def plotift(data, iftdata):
    rowused = 0 #Make sure this matchs the integer in return extracteddata[0, :]
                #from choose_data_to_transform(data)
    n = len(data)
    E = np.linspace(-1.5, 1.5, n)
    t = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    sig = 0.01
    pylab.clf()
    pylab.plot(E, data, '-o', label='Original Data')
    if rowused == 0:
        pylab.plot(t, fourier_transform_lorentzian(t), label='Lorentzian Transformed Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Lorentzian signal')
        pylab.legend()
    elif rowused == 1:
        pylab.plot(t, fourier_transform_gaussian(t, sig), '-o', label='Gaussian Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Gaussian signal')
        pylab.legend()
    elif rowused == 2:
        pylab.plot(t, transform_convolution(t, sig), label='Convolution Transformed Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Gaussian signal')
        pylab.legend()
    elif rowused == 3:
        pylab.plot(t, transform_convolution(t, sig), label='Convolution Transformed Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform of Convoluted Data with Noise')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Gaussian signal')
        pylab.legend()
    assert rowused < 4 and rowused >= 0


def fouriertransform(t, data, e):
    """
    fouriertransform takes a discrete data set and computes its fourier
    transform over a given time interval. This is accomplished using the
    trapoziodal approximation of the integral.
    """
    initialarraysize = data.size
    indices = np.arange(initialarraysize)
    ftdata = data.copy()
    for index in indices:
        transformedvalue = np.trapz(data* exp(-1*1j*t*e[index]), x=t)
        np.put(ftdata, index, transformedvalue)
    print ftdata
    return ftdata

def plotfouriertransform(data, fouriertransform):
    n = len(data)
    E = np.linspace(-1.5, 1.5, n)
    pylab.clf()
    pylab.plot(E, data, label='Original Data')
    pylab.plot(E, fouriertransform, '-o',
               label='Original data back from Inverse Fourier Transform')
    pylab.legend()


data_to_use = choose_data_to_transform(upload_data())
n = len(data_to_use)
E = np.linspace(-1.5, 1.5, n)
T = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
result = inversetransform(E, data_to_use, T)
        #fouriertransform(inversetransform(data))
#print "norm",np.linalg.norm(data-result)
plotift(data_to_use, result)
fouriertransform(T, result, E)
#ftreverse = fouriertransform(np.linspace(-(1/3)*pi*n/2 , (1/3)*pi*n/2, n),
#                             result, np.linspace (-1.5, 1.5, n))
#plotfouriertransform(data, ftreverse)
