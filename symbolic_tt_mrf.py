#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:07:51 2022

@author: abouabdallah
"""

import sys
import os


import numpy as np

#import matplotlib.pyplot as plt


def build_symbolic_ttfactors(n, edges):

    m = edges.shape[0]
    facs = [[] for i in range(m)]
    for l in range(m):
        cfacs = ['empty' for i in range(n)]
        i = int(edges[l, 0])
        j = int(edges[l, 1])
        # print(l, list_of_edges[l], i, j)
        for k in range(i):
            cfacs[k] = '1'
        for k in range(j+1, n):
            cfacs[k] = '1'
        cfacs[i] = 'U'
        for k in range(i+1, j):
            cfacs[k] = 'I'
        cfacs[j] = 'V'
        facs[l] = cfacs
    return facs




def loadprod(load) :
    prod =0
    for k in range (len(load)) :
        prod= prod+np.log(load[k])
    return (prod)


def getFactors(n):
    # Create an empty list for factors
    factors=[];

    # Loop over all factors
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    # Return the list of factors
    return factors
