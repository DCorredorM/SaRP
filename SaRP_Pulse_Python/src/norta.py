#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 12:58:58 2018
Updateed on Wed Jun 19 13:14:19 2019
@author: Daniel Duque
@email: d.duque25@gmail.com
This module implements Normal-to-Anything (NORTA) algorithm
to generate correlated random vectors. The original paper is by
Cario and Nelson (2007). 
"""

import numpy as np
from scipy.linalg import cholesky
import scipy.stats as stats
#from fitter import Fitter

rnd = np.random.RandomState(0)  # Random stream
Z = stats.norm(0, 1)  # Standar normal
GAUSS_SAMPLING = True  # Parameter for montecarlo integration
OUTPUT = 1  # Output flag. 0: no output, 1: process output
MC_SAMPLES = 10000000  # Number of samples to compute the integral


def reset_stream():
    global rnd
    rnd = np.random.RandomState(0)  # Random stream


def find_rho_z(i, j, F_invs, CovX, EX):
    """
    Computes the correlation of the multivariate normal used in the
    generation of random variables. The correlation is found by means of
    a binary search where for each value, the targe covariance between i and j
    is computed as an integral via Monte Carlo. The sampling procedure leverages
    in the bivarite normal distribution embeded in the covariance term.
    
    Args:
        i,j (int): pair of indices for which the correlation is being computed
        F_invs (list of func): list of functions. Each function is the inverse of 
            the marginal distribution of each random variable.
        CovX (ndarray): Covariance matrix of the input
        EX (ndarray): Mean vector of the input
    """
    if OUTPUT == 1:
        print("Computing rhoZ(%i,%i)" % (i, j))
    cor_dem = np.sqrt(CovX[i, i] * CovX[j, j])
    rho_z = CovX[i, j] / cor_dem
    rho_u = 1 if CovX[i, j] > 0 else 0
    rho_d = -1 if CovX[i, j] < 0 else 0
    F_i_inv = F_invs[i]
    F_j_inv = F_invs[j]
    EXi, EXj = EX[i], EX[j]
    while np.abs(rho_u - rho_d) > 1e-4:
        covZ = np.array([[1, rho_z], [rho_z, 1]])
        f = conv_exp(covZ, F_i_inv, F_j_inv, gaussian_sampling=GAUSS_SAMPLING)
        EXiXj = montecarlo_integration(
            f, n=MC_SAMPLES, c=covZ, m=np.zeros(2), gaussian_sampling=GAUSS_SAMPLING
        )
        CXiXj = EXiXj - EXi * EXj
        if OUTPUT == 1:
            print(
                "  rhoZ=%10.4e, C(i,j)=%10.4e, Cov=%10.4e" % (rho_z, CXiXj, CovX[i, j])
            )
        if np.abs(CXiXj - CovX[i, j]) / cor_dem < 1e-4:
            #
            return rho_z
        else:
            if CXiXj > CovX[i, j]:
                rho_u = rho_z
                rho_z = 0.5 * (rho_z + rho_d)
            else:  # rhoC_ij <= rho_ij
                rho_d = rho_z
                rho_z = 0.5 * (rho_z + rho_u)

    return rho_z


def montecarlo_integration(f, m=None, c=None, n=1000000, gaussian_sampling=False):
    """
    Computes the integral for the particular function in NORTA.
    WARNING: This method is not general for other functions as it is.
    """
    if gaussian_sampling:
        assert type(m) != type(None), "Mean and Cov are required for gaussian sampling"
        z_trial = rnd.multivariate_normal(m, c, n)
        integral = np.sum(f(z_trial[:, 0], z_trial[:, 1]))
        return integral / n
    else:
        return montecarlo_integration_uniform(f, n)


def montecarlo_integration_uniform(f, n=1000):
    """
    Basic integration function using uniform sampling. The cube size
    is determined based on the fact that almost all the mass in a bivariate
    standar normal distribution is within -5,5 x -5,5. 
    """
    cube_size = 5
    z1_trials = rnd.uniform(-cube_size, cube_size, n)
    z2_trials = rnd.uniform(-cube_size, cube_size, n)
    V = 2 * cube_size * 2 * cube_size
    integral = np.sum(f(z1_trials, z2_trials))
    return V * integral / n


def PolyArea(x, y):
    """
    Nice function to compute the area enclosed by a sequence of points.
    Not used in this module, but potencialy usefull for other montecarlo 
    integration functions. 
    """
    # https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def conv_exp(covZ, F_i_inv, F_j_inv, gaussian_sampling=True):
    """
    Integrant in NORTA. 
    Args:
        covZ (ndarray): Covariance of the bivariate normal distribution
        F_i_inv (func): Inverse function of the marginal for variable i (j)
        gaussian_sampling (bool): True if the function needs to be modified due
            to the sampling mechanism in the Montecarlo Integration method
    """
    if gaussian_sampling:

        def f(z1, z2):
            return F_i_inv(Z.cdf(z1)) * F_j_inv(
                Z.cdf(z2)
            )  # remove bi_x pdf since montecarlo is sampling according to  bi_z

        return f
    else:
        bi_z = stats.multivariate_normal(cov=covZ)

        def f(z1, z2):
            z1z2 = np.vstack((z1, z2)).transpose()
            return F_i_inv(Z.cdf(z1)) * F_j_inv(Z.cdf(z2)) * bi_z.pdf(z1z2)

        return f


def build_empirical_inverse_cdf(X):
    """
    Builds an inverse CDF function given a sorted vector of values defining a 
    marginal distribution.
    
    Args:
        X:Sorted vector of observations
    """
    n = len(X)

    def f(prob):
        """
        Args:
            prob (ndarray): vector with probablities to compute the inverse
        """
        # assert 0<=prob<=1, 'Argument of inverse function is a probability >=0 and <= 1.'
        return X[np.minimum((n * np.array(prob)).astype(int), n - 1)]

    return f


def fit_NORTA(n, d,data=None,CovX=None,EX=None, F_invs=None, output_flag=0):
    """
    Computes covariance matrix for NORTA algorith. 
    NORTA object assumes emprical marginal distributions.
    
    Args:
        data (ndarray): a n x d array with the data
        n (int): number of observations for each random variable
        d (int): dimanention of the random vector.
        F_invs (list of func): optional parameter to specify the marginal 
            distributions. Default is None, which constructs the marginals from
            the data.
        CovX: (ndarray): optional parameter with a d x d array with the desired covariance matrix
        EX: (ndarray): optional parameter with a 1 x d array with the desired Expected values of each random variable
    Return:
        NORTA_GEN (NORTA): an object that stores the necessary information to 
            generate NORTA random vectors. 
    """
    reset_stream()
    global OUTPUT
    OUTPUT = output_flag

    if type(data)!=type(None):
        assert len(data) == n, "Data needs to be a d x n matrix"
        assert len(data[0]) == d, "Data needs to bo a d x n matrix"
        CovX = np.cov(data, rowvar=False)
    
    if OUTPUT == 1:
        print("Starting NORTA fitting")
        print("Finding %i correlation terms" % (int(d * (d - 1) / 2)))
    
    C = None  # matrix for NORTA
    
    lambda_param = 0.01

    

    VarX = np.diag(np.diag(CovX))

    procedure_done = False
    while procedure_done == False:
        D = np.eye(d)
        if EX==None:
            EX = np.mean(data, axis=0)
            print('esto es exp: ',EX )
        if type(F_invs) != list:
            F_invs = [
                build_empirical_inverse_cdf(np.sort(data[:, i])) for i in range(d)
            ]

        for i in range(d):
            for j in range(i + 1, d):
                D[i, j] = find_rho_z(i, j, F_invs, CovX, EX)
                D[j, i] = D[i, j]
        try:
            C = cholesky(D, lower=True)
            procedure_done = True
        except:
            CovX = (1 - lambda_param) * CovX + lambda_param * VarX
            print("Cholesky factorization failed, starting over")

    NORTA_GEN = NORTA(F_invs, C)
    return NORTA_GEN


def fit_NORTA_CovX(CovX,EX, F_invs,data=None, output_flag=0):
    """
    Computes covariance matrix for NORTA algorith. 
    NORTA object assumes emprical marginal distributions.
    
    Args:
        data (ndarray): a n x d array with the data
        n (int): number of observations for each random variable
        d (int): dimanention of the random vector.
        F_invs (list of func): optional parameter to specify the marginal 
            distributions. Default is None, which constructs the marginals from
            the data.
        CovX: (ndarray): optional parameter with a d x d array with the desired covariance matrix
        EX: (ndarray): optional parameter with a 1 x d array with the desired Expected values of each random variable
    Return:
        NORTA_GEN (NORTA): an object that stores the necessary information to 
            generate NORTA random vectors. 
    """
    reset_stream()
    global OUTPUT
    OUTPUT = output_flag
    d=len(EX)
    
    if type(data)!=type(None):
        assert len(data) == n, "Data needs to be a d x n matrix"
        assert len(data[0]) == d, "Data needs to bo a d x n matrix"
        CovX = np.cov(data, rowvar=False)
    
    if OUTPUT == 1:
        print("Starting NORTA fitting")
        print("Finding %i correlation terms" % (int(d * (d - 1) / 2)))
    
    C = None  # matrix for NORTA
    
    lambda_param = 0.01

    

    VarX = np.diag(np.diag(CovX))

    procedure_done = False
    while procedure_done == False:
        D = np.eye(d)
        if EX==None:
            EX = np.mean(data, axis=0)
            print('esto es exp: ',EX )
        if type(F_invs) != list:
            F_invs = [
                build_empirical_inverse_cdf(np.sort(data[:, i])) for i in range(d)
            ]

        for i in range(d):
            for j in range(i + 1, d):
                D[i, j] = find_rho_z(i, j, F_invs, CovX, EX)
                D[j, i] = D[i, j]
        try:
            C = cholesky(D, lower=True)
            procedure_done = True
        except:
            CovX = (1 - lambda_param) * CovX + lambda_param * VarX
            print("Cholesky factorization failed, starting over")

    NORTA_GEN = NORTA(F_invs, C)
    return NORTA_GEN



class NORTA:
    """
    Class to create a Normal-to-Anything model
    Attributes:
        F (list of func): Marginal inverse CDFs
        C (ndarry): numpy array with the Cholesky factorization
    """

    def __init__(self, Finv, C):
        reset_stream()
        assert len(Finv) == len(C), "Dimension of the marginals and C dont match."
        self.F_inv = Finv
        self.C = C
        reset_stream()

    def gen(self, n=1):
        """
        Generates an array of vectors of where each component follow the
        marginal distribution and the realization are correlated by means of
        the covariance matrix CovZ computed in the fitting process.
        
        Args:
            n (int): number of samples to generate
        """
        d = len(self.F_inv)
        w = rnd.normal(size=(d, n))
        z = self.C.dot(w)

        '''
        print([i.mean() for i in z])

        import matplotlib.pyplot as plt
        p=[Z.cdf(z[i]) for i in range(d)]
        for d,i in enumerate(p):
            plt.hist(i)
            plt.show()
        '''

        X = np.array([self.F_inv[i](Z.cdf(z[i])) for i in range(d)]).transpose()
        return X

