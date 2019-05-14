# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:29:41 2016

@author: bakerh
"""
import numpy as np


class PCA1(object):
    def __init__(self, a, standardize=True):
        """
        compute the SVD of a and store data for PCA.  Use project to
        project the data onto a reduced set of dimensions

        Inputs:

          *a*: a numobservations x numdims array
          *standardize*: True if input data are to be standardized. If False,
          only centering will be carried out.

        Attrs:

          *a* a centered unit sigma version of input a

          *numrows*, *numcols*: the dimensions of a

          *mu*: a numdims array of means of a. This is the vector that points
          to the origin of PCA space.

          *sigma*: a numdims array of standard deviation of a

          *fracs*: the proportion of variance of each of the principal
          components

          *s*: the actual eigenvalues of the decomposition

          *Wt*: the weight vector for projecting a numdims point or array into
          PCA space

          *Y*: a projected into PCA space


        The factor loadings are in the Wt factor, i.e., the factor
        loadings for the 1st principal component are given by Wt[0].
        This row is also the 1st eigenvector.

        """
        n, m = a.shape
        if n > m:
            raise RuntimeError('we assume data in a is organized with '
                               'numrows>numcols')

        self.numrows, self.numcols = n, m
        self.mu = a.mean(axis=0)
        self.sigma = a.std(axis=0)
        self.standardize = standardize

        a = self.center(a)

        self.a = a

        Vh, s, U = np.linalg.svd(a)

        # Note: .H indicates the conjugate transposed / Hermitian.

        # The SVD is commonly written as a = U s V.H.
        # If U is a unitary matrix, it means that it satisfies U.H = inv(U).

        # The rows of Vh are the eigenvectors of a.H a.
        # The columns of U are the eigenvectors of a a.H.
        # For row i in Vh and column i in U, the corresponding eigenvalue is
        # s[i]**2.

        self.Wt = Vh
        Vh1 = Vh.copy()
        for i in range(84):
            Vh1[:, i] = (Vh1[:, i] - np.mean(Vh1[:, i])) / np.std(Vh1[:, i])
        # save the transposed coordinates
        Y = np.dot(a.T, Vh1)
        self.Y = Y
        self.Wt1 = Vh1
        # save the eigenvalues
        self.s = s**2

        # and now the contribution of the individual components
        vars = self.s/float(len(s))
        self.fracs = vars/vars.sum()

    def project(self, x, minfrac=0.):
        '''
        project x onto the principle axes, dropping any axes where fraction
        of variance<minfrac
        '''
        x = np.asarray(x)

        ndims = len(x.shape)

        if (x.shape[-1] != self.numcols):
            raise ValueError('Expected an array with dims[-1]==%d' %
                             self.numcols)

        Y = np.dot(self.Wt, self.center(x).T).T
        mask = self.fracs >= minfrac
        if ndims == 2:
            Yreduced = Y[:, mask]
        else:
            Yreduced = Y[mask]
        return Yreduced

    def center(self, x):
        '''
        center and optionally standardize the data using the mean and sigma
        from training set a
        '''
        if self.standardize:
            return (x - self.mu)/self.sigma
        else:
            return (x - self.mu)

    @staticmethod
    def _get_colinear():
        c0 = np.array([
            0.19294738,  0.6202667,   0.45962655,  0.07608613,  0.135818,
            0.83580842,  0.07218851,  0.48318321,  0.84472463,  0.18348462,
            0.81585306,  0.96923926,  0.12835919,  0.35075355,  0.15807861,
            0.837437,    0.10824303,  0.1723387,   0.43926494,  0.83705486])

        c1 = np.array([
            -1.17705601, -0.513883,   -0.26614584,  0.88067144,  1.00474954,
            -1.1616545,   0.0266109,   0.38227157,  1.80489433,  0.21472396,
            -1.41920399, -2.08158544, -0.10559009,  1.68999268,  0.34847107,
            -0.4685737,   1.23980423, -0.14638744, -0.35907697,  0.22442616])

        c2 = c0 + 2*c1
        c3 = -3*c0 + 4*c1
        a = np.array([c3, c0, c1, c2]).T
        return a
