#!/usr/bin/env python

# -*- coding: utf-8 -*-
#This program is public domain
"""
Created on Mon Jun 01 15:39:51 2015

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
    newdata = np.loadtxt('DataFiles\Mo72Cr30_HT4K_dyn', dtype='float', delimiter=',')
#    newdata = np.loadtxt(str(raw_input('Enter the name of data file to load:')),
#                         dtype='float', delimiter=',')
    print newdata.shape
    return newdata


def choose_data_from_file(data, desiredcolumn): #data must be a e.size x 4
    """
    choose_data_to_transform takes a multi-dimensional array and reports the
    nth column of the array for future use.
    """
    extracteddata = data
    extracteddata = extracteddata.T
    return extracteddata[desiredcolumn, :]
#    return extracteddata[int(raw_input('Column of Desired Data: ')),:]

def create_lorentzian(e):
    """
    create_lorentzian creates an array containing points on a Lorentzian curve
    with parameter gamma from points within the energy window requires an input
    e which is an symetric array of energy points.
    """
    tau = 50
    gamma = 1/tau
    sqe = (1/pi)*(gamma/2) / (e**2 + (gamma/2)**2) #Standard equation of
                                                   #a Lorentzian function
    return sqe #returns an array of points on the Lorentzian curve


def create_gaussian_norm(e, sig):
    """
    create_gaussian_norm creates an array containing points on a Normal
    Gaussian curve with standard deviation, sigma, and mu = 0
    (i.e. symmetic about 0) from points within the energy window requires an
    input e which is an symetric array of energy points, and sig which is the
    standard deviation of the Gaussian.
    """
    rqe = exp(-e**2 / (2*sig**2)) / (sig*sqrt(2*pi))
    return rqe #returns an array of points on the Gaussian curve


def convolve(e, sig):
    """
    convolve does a discrete convolution of points on the Lorentzian and
    Gaussian curves, requires an input e which is an symetric array of energy
    points, and sig which is the standard deviation of the Gaussian.
    """
    l1 = create_lorentzian(e)
    g1 = create_gaussian_norm(e, sig)
    deltaE = e[1] - e[0] #the change in e for each step along the curves
    convo = np.convolve(l1, g1, 'same')
    convo = convo*deltaE
    return convo


def add_noise(orig):
    """
    add_noise adds random noise into an array of data and requires an array of
    convoluted values as an input.
    """
    levelofnoise = 0.3
    result = [xi+np.random.randn()*levelofnoise*xi for xi in orig]
        #At each index, a new value, consisiting of the original value plus
        #random noise, i.e. (the level of noise)*(a random number from
        #the normal distribution*the original value)
    return np.array(result)



def show_data(e, s, lsig, gsig, convosig, rn_convosig):
    """
    Plots shows the Lorentzian, Gaussian, Convoluted values, and convoluted
    values with random noise against the energy points used as an input.
    It requires an input e which is an symetric array of energy points, sig
    which is the standard deviation of the Gaussian, and the Lorentzian,
    Gaussian, and Convoluted values, and convoluted values with random noise.
    """
    pylab.clf()
    pylab.plot(e, lsig, '-o', label='Lorentzian')
    pylab.plot(e, gsig, '-o', label='Gaussian')
    pylab.plot(e, convosig, '-o', label='Convolution %g'%s)
#    pylab.plot(e, rn_convosig, '-o', label='RandomNoiseConvolution %g'%s)
    pylab.legend()
    pylab.show()


def print_values(e, lsig, gsig, convosig, rn_convosig):
    deltaE = (e[1]-e[0])
    print 'DeltaE = %g'%deltaE
    print 'Lorentzian ='
    print lsig
    print 'Gaussian ='
    print gsig
    print 'Convoluted Data ='
    print convosig
    print 'Convoluted Data with noise ='
    print rn_convosig
    print 'Size of Lorentzian = %g'%lsig.size
    print 'Size of Gaussian = %g'%gsig.size
    print 'Size of Convoluted Data = %g'%convosig.size
    print 'Size of Convoluted Data = %g'%rn_convosig.size
    print 'Area under Lorentzian = %g'%np.sum(lsig*deltaE)
    print 'Area under Gaussian = %g' %np.sum(gsig*deltaE)
    print 'Area under Convoluted Function = %g'%np.sum(convosig*deltaE)
    print 'Area under Convoluted Function = %g'%np.sum(rn_convosig*deltaE)



def export_data_csv(data, filename):
    """
    export_data requires an array of data as an input and then exports that
    data to a text file.
    """
    np.savetxt(filename+'.csv', data, delimiter=',')
        #use .csv to be able to read with excel


def create_data_to_export(e, s, basename):
    """
    create_data_to_export calls the above functions to report the Lorentzian
    values, Gaussian values, Convoluted values, and convoluted values with
    random noise.
    It requires an input e which is an symetric array of energy points, and an
    array of sigmas which are the standard deviations of the Gaussians.
    """
    for i,sig in enumerate(s):
        lsig = create_lorentzian(e)
        gsig = create_gaussian_norm(e, sig)
        convosig = convolve(e, sig)
        rn_convosig = add_noise(convosig)
        show_data(e, sig, lsig, gsig, convosig, rn_convosig)
        print_values(e, lsig, gsig, convosig, rn_convosig)
        toexport = np.vstack((lsig, gsig, convosig, rn_convosig))
        toexport = toexport.transpose()
        #export_data_csv(toexport, basename+str(i)+'sigma-%g'%sig)



def count_values(e, sig):
    """
    count_values counts the number of the convolted values for a specific
    energy array and sigma above a certain value.
    It is used to calculate the number of values wihin the Gaussian.
    """
    origarray = convolve(e, sig)
    arraywcondition = origarray[np.where(origarray >= 3.2)]
        #One needs to change float value in origarray >= for each sigma to
        #measure the number of points within the Gaussian.
    print arraywcondition.size


def demo1():
    create_data_to_export(np.linspace(-1.5, 1.5, 2251),
                          [0.003, 0.0038, 0.005, 0.035, 0.05],
                          raw_input("Choose a general filename: "))

def demo2():
    create_data_to_export(np.linspace(-1.5, 1.5, 2251), [0.01],
                          raw_input("Choose a general filename: "))

if __name__ == "__main__":
    demo2()


