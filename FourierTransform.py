#!/usr/bin/env python

# -*- coding: utf-8 -*-
#This program is public domain
"""
Created on Fri Jun 05 11:45:39 2015

@author: jkh
"""

from __future__ import division

import pylab
import numpy as np
from numpy import exp, pi, sqrt


def upload_data():
    """
    upload_data takes data from previously exported csv file, and uploads it
    to be used in future code.
    """
#    newdata = np.loadtxt('DataFiles\Trial4.csv', dtype='float', delimiter=',')
    newdata = np.loadtxt(str(raw_input('Enter the name of data file to load:')),
                         dtype='float', delimiter=',')
    print newdata.shape
    return newdata


def choose_data_from_file(data, desiredcolumn): #data must be a e.size x 4
    """
    choose_data_from_file takes a multi-dimensional array and reports the
    nth column of the array for future use.
    """
    extracteddata = data
    extracteddata = extracteddata.T
    return extracteddata[desiredcolumn, :]
#    return extracteddata[int(raw_input('Column of Desired Data: ')),:]


def fourier_transform_gaussian(t, sig):
    """
    fourier_transform_gaussian creates an array containing points on a Normal
    Gaussian curve with standard deviation, sigma, and mu = 0 represents the
    Fourier Transform from the same sigma from a Normal Gaussian of
    the form exp(-e**2 / (2 * sig**2)) / (sig * (sqrt(2 * pi))).
    """
    FRQe = (1/(2*pi))*exp(-(t**2) * (sig**2) / 2)
    print FRQe
    return FRQe #returns an array of points on the Gaussian curve


def fourier_transform_lorentzian(t):#Does not result in I(Q,t), results in (2*pi)*I(Q,t)
    """
    fourier_transform_lorentzian creates an array containing points on an
    expontential curve with the same tau as the Lorentzian whose Fourier
    transformation it represents.
    """
    tau = 50
    lft = t.copy()
    lft = (1/(2*pi))*exp(-(abs(lft)) / (2*tau))
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
    convotrans = fourier_transform_gaussian(t, sig) * (2*pi)*fourier_transform_lorentzian(t)
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
#    print iftdata.dtype
    for index in indices:
        transformedvalue = (1/(2*pi))*np.trapz(data* exp(1j*e*t[index]), x=e)
        #np.put(iftdata, index, transformedvalue)
        iftdata[index] = transformedvalue
#    print iftdata
    return iftdata


def plotift(data, iftdata):
    rowused = 2 #Make sure this matchs the integer in return extracteddata[0, :]
                #from choose_data_from_file(data)
    n = len(data)
    E = np.linspace(-1.5, 1.5, n)
    t = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    sig = 0.01
    pylab.clf()
    pylab.plot(E, data, '-o', label='Original Data')
    if rowused == 0:
        pylab.plot(t, fourier_transform_lorentzian(t),
                   label='Lorentzian Transformed Exact')
        pylab.plot(t, iftdata, '-o',
                   label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Lorentzian signal')
        pylab.legend()
    elif rowused == 1:
        pylab.plot(t, fourier_transform_gaussian(t, sig), '-o',
                   label='Gaussian Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Gaussian signal')
        pylab.legend()
    elif rowused == 2:
        pylab.plot(t, transform_convolution(t, sig),
                   label='Convolution Transformed Exact')
        pylab.plot(t, iftdata, '-o', label='Inverse Fourier Transform')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Convoluted signal')
        pylab.legend()
    elif rowused == 3:
        pylab.plot(t, transform_convolution(t, sig),
                   label='Convolution Transformed Exact')
        pylab.plot(t, iftdata, '-o',
                   label='Inverse Fourier Transform of Convoluted Data with Noise')
        pylab.xlabel('Time')
        pylab.ylabel('Intensity (counts/s)')
        pylab.title('Simulated intensity for a Convoluted signal with noise')
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
    #print indices
    ftdata = np.empty_like(data, dtype='complex')
    for index in indices:
        transformedvalue = np.trapz(data* exp(-1*1j*t*e[index]), x=t)
        ftdata[index] = transformedvalue
#    print ftdata
    return ftdata


def plotfouriertransform(data, fouriertransform):
    n = len(data)
    E = np.linspace(-1.5, 1.5, n)
    pylab.clf()
    pylab.plot(E, data, label='Original Data')
    pylab.plot(E, fouriertransform, '-o',
               label='Original data back from Inverse Fourier Transform')
    pylab.legend()


def demo1():
    data_to_use = choose_data_from_file(upload_data(), 0)
    n = len(data_to_use)
    E = np.linspace(-1.5, 1.5, n)
    T = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    result = inversetransform(E, data_to_use, T)
    ftreverse = fouriertransform(T, result, E)
    plotfouriertransform(data_to_use, ftreverse)
    pylab.show()


def demo2():
    data_to_use = choose_data_from_file(upload_data(), 0)
    n = len(data_to_use)
    E = np.linspace(-1.5, 1.5, n)
    T = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    result = fouriertransform(T, data_to_use, E)
    plotift(data_to_use, result)
    pylab.show()


def demo3():
    data_to_use = choose_data_from_file(upload_data(), 2)
    n = len(data_to_use)
    E = np.linspace(-1.5, 1.5, n)
    T = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    result = inversetransform(E, data_to_use, T)
    print result
    plotift(data_to_use, result)
    pylab.show()


def demo4():
    data_to_use = choose_data_from_file(upload_data(), 3)
    n = len(data_to_use)
    E = np.linspace(-1.5, 1.5, n)
    T = np.linspace(-(1/3)*pi*n/2, (1/3)*pi*n/2, n)
    result = inversetransform(E, data_to_use, T)
    plotift(data_to_use, result)
    pylab.show()


def show_testing(E, t, ynum, Y, yt):
    """
    shows the graphs for test_Gaussian and test_Lorentzian
    """
    if 1:
        pylab.clf()
        pylab.subplot(211)
        pylab.plot(E, abs(ynum), '.', label='data')
        pylab.plot(E, abs(yt), '.-', label='ift(ft(data))')
        pylab.legend()
        pylab.subplot(212)
        pylab.plot(t, abs(Y), '.', label='ft(data)')
        pylab.legend()
        pylab.show()


def test_Gaussian():
    """
    Checks that the fourier transform's value at 0 matches the area under the
    original curve, taking the fourier transform of the inverse transform
    matches the original data, and the inverse transform of the data
    corresponds to the therotical definition of its inverse transform for the
    Gaussian distribution.
    """
    #nE, nt, width = 81, 51, 5
    nE, nt, width = 2251, 2251, 150
    mu, sigma = 0, 0.01
    E = np.linspace(-width*sigma+mu, width*sigma+mu, nE)
    tLoren = np.linspace(-(pi * nt) / ((4*width*sigma+mu)), (pi * nt) / ((4*width*sigma+mu)), nt)
    y1 = exp(-0.5*(E-mu)**2/sigma**2)/sqrt(2*pi*sigma**2)
#    print 'y1:', y1

    # F(0) = \int f(x) e^{-i x 0} dx = \int f(x) e^0 dx = \int f(x)
    ft = fouriertransform(E, y1, tLoren)
    ft_zero = ft[(nt-1)/2] #Index must be the center of the matrix
    y_area = np.trapz(x=E, y=y1)
    assert abs(y_area - ft_zero) < 1e-8, "ft(0): %g, Sy: %g"%(ft_zero, y_area)

    # Fourier transform should be invertible, in that the inverse transform
    # of the forward transform should match the original function.
    Y = inversetransform(E, y1, tLoren)
    yt = fouriertransform(tLoren, Y, E)
    show_testing(E, tLoren, y1, Y, yt)
    assert np.all(abs(y1-yt) < 1e-12), "max err: %g"%max(abs(y1-yt))

    # Fourier transform of a Gaussian using unitary angular frequency
    # is a Gaussian of 1/sigma, scaled so that the peak is 1.  Any shift
    # in the center corresponds to a phase shift of e^{-i mu t}
    Y1theory = exp(-1j*mu*tLoren) * (1/(2*pi)) * exp(-(tLoren**2) * (sigma**2) / 2)
    assert np.all(abs(Y-Y1theory) < 1e-12), "max err: %g"%max(abs(Y-Y1theory))


def test_Lorentzian():
    """
    Checks that the fourier transform's value at 0 matches the area under the
    original curve, taking the fourier transform of the inverse transform
    matches the original data, and the inverse transform of the data
    corresponds to the therotical definition of its inverse transform for the
    Lorentzian distribution.
    """
    nE, nt, width = 2251, 2251, 150
    mu, sigma = 0, 0.01
    E = np.linspace(-width*sigma+mu, width*sigma+mu, nE)
    tLoren = np.linspace(-(pi * nt) / (4*(width*sigma+mu)), (pi * nt) / (4*(width*sigma+mu)), nt)
#    print tLoren
    tau = 50
    gamma = 1 / tau
    y0 = (1/pi)*(gamma/2) / (E**2 + (gamma/2)**2)
    print y0

    #F(0) = \int f(x) e^{-i x 0} dx = \int f(x) e^0 dx = \int f(x)
    ft = fouriertransform(E, y0, tLoren)
    ft_zero = ft[(nt-1)/2] #Index must be the center of the matrix
    y_area = np.trapz(x=E, y=y0)
    assert abs(y_area - ft_zero) < 0.001, "ft(0): %g, Sy: %g"%(ft_zero, y_area)

    # Fourier transform should be invertible, in that the inverse transform
    # of the forward transform should match the original function.
    Y = inversetransform(E, y0, tLoren)
    yt = fouriertransform(tLoren, Y, E)
    show_testing(E, tLoren, y0, Y, yt)
    assert np.all(abs(y0-yt) < 0.001), "max err: %g"%max(abs(y0-yt))


    # Fourier transform of a Lorentzian is an expontential of the form
    # (1/(2*pi)e^(t/(2*tau)).  Any shift in the center corresponds to a phase
    # shift of e^{-i mu t}.
    Y0theory = (1/(2*pi))*exp(-(abs(tLoren)) / (2*tau))
    assert np.all(abs(Y-Y0theory) < 0.001), "max err: %g"%max(abs(Y-Y0theory))


def test_convolution():
    """
    Checks that the fourier transform's value at 0 matches the area under the
    original curve and the inverse transform of the data corresponds to the
    therotical definition of its inverse transform for the
    convoluted curve.
    """
    nE, nt, width = 2251, 2251, 150
    mu, sigma = 0, 0.01
    E = np.linspace(-width*sigma+mu, width*sigma+mu, nE)
    tLoren = np.linspace(-(pi * nt) / (4*(width*sigma+mu)), (pi * nt) / (4*(width*sigma+mu)), nt)
    tau = 50
    gamma = 1 / tau
    y0 = (1/pi)*(gamma/2) / (E**2 + (gamma/2)**2)
    y1 = exp(-0.5*(E-mu)**2/sigma**2)/sqrt(2*pi*sigma**2)
    y2 = np.convolve(y0, y1, 'same')*(E[1]-E[0])
#    print y0, y1, y2
    #F(0) = \int f(x) e^{-i x 0} dx = \int f(x) e^0 dx = \int f(x)
    ft = fouriertransform(E, y2, tLoren)
    ft_zero = ft[(nt-1)/2] #Index must be the center of the matrix
    y_area = np.trapz(x=E, y=y2)
    assert abs(y_area - ft_zero) < 1e-8, "ft(0): %g, Sy: %g"%(ft_zero, y_area)

    # Fourier transform should be invertible, in that the inverse transform
    # of the forward transform should match the original function.
    Y = inversetransform(E, y2, tLoren)
    print 'Y: ', Y
    ygaus = inversetransform(E, y1, tLoren)
    yloren = inversetransform(E, y0, tLoren)
    yconvot = ygaus * yloren *(2*pi)
    if 1:
        pylab.clf()
        pylab.plot(tLoren, abs(yconvot), '.', label='Convoluted divided by Gaussian')
        pylab.plot(tLoren, abs(ygaus), '.', label='Guassian curve')
        pylab.plot(tLoren, abs(yloren), '.', label='Fourier Transformed Laurentian')
        pylab.plot(tLoren, abs(Y), '.', label='ft(data)')
        pylab.legend()
        pylab.show()
    assert np.all(abs(Y - yconvot) < 0.001), "max err: %g"%max(abs(Y - yconvot))


if __name__ == "__main__":
    #demo1()
    #demo2()
    #demo3()
    demo4()
    #test_Gaussian()
    #test_Lorentzian()
    #test_convolution()

