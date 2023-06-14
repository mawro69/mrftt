#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:26:06 2022

@author: abouabdallah
"""



import numpy as np

# fontions de bases

def listoflist(d):
    """ initialise a liste of liste of shape d

   arguments :
       @ d  int : shape of litse

   return :
       @ listeflists : lempty list of list"""

    listoflists = []
    a_list = [1, 5]
    for k in range(d):
        listoflists.append(a_list)
    return(listoflists)


# ----------------------------------------------------------------------


def multi(vec1, vec2):
    """ Compute the elementwise product of two vectors
    arguments :
       @ vec1, vec2 : two vecors

   return :
       @ vec3: vectors elementwise product of vec1 and vec2
    """
    n = len(vec1)
    vec3 = np.ones(n)
    for k in range(n):
        vec3[k] = vec1[k]*vec2[k]
    return vec3

# ======================================================================
    
def get_m(n):
    """ get the number of links m. 
    arguments :
        n : int number of individuals
    return :
        m : int : number of links 

    """
    return int(n*(n-1)/2)


# ======================================================================


def getFactors(m):
    """ get the number the dividors of m. 
    arguments :
        n : int number of individuals
    return :
        m : int : number of links 

    """
    facs = []     # Create an empty list for factors

    # Loop over all factors
    for i in range(1, m + 1):
        if m % i == 0:
            facs.append(i)

    return facs     # Return the list of factors


# ======================================================================


def getFactors_clean(n):
    lis = getFactors(n)
    f = len(lis)
    lis.remove(lis[f-1])
    lis.remove(lis[0])
    return lis
