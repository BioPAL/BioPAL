import logging
import numpy as np
from   scipy.stats import chi2

def get_rho_omega(j,p,n):
    # j is the number of covariance matrices being considered
    # p is the dimension of the covariance matrix (.ie 2 for dual pol, 3 for full pol)
    # n is the number of looks
    rho = 1 - (2*p**2 - 1)*(1 + 1/(j*(j - 1)))/(6*p*n)
    omega2 =-p**2*(1 - 1/rho)**2/4 + p**2*(p**2 - 1)*(1 + (2*j - 1)/(j**2*(j - 1)**2))/(24*n**2*rho**2)
    return rho,omega2


def funR_wishart_SU(Y,X,j,n,p):
    
    # Here X is the next time point and Y is the sum of the previous X's after
    # the last change was detected. j is the total number of data being
    # considered ie. those which have contributed to Y plus one.
    
    lnR = n*(p*(j*np.log(j) - (j-1)*np.log(j-1)) + (j-1)*np.log(np.linalg.det(Y)) + np.log(np.linalg.det(X)) - j*np.log(np.linalg.det(Y+X)))
    lnR = np.real(lnR)
    return lnR


def alg_wishart_SU(X,Y,j,p,n,alpha,verbose):

    change = 0
    rho, omega2 = get_rho_omega(j+1,p,n)
    testStat    = -2*rho*funR_wishart_SU(Y,X,j+1,n,p)
    
    if testStat>0:
        probR = 1 - (chi2.cdf(testStat,p**2) + omega2*(chi2.cdf(testStat,p**2+4) - chi2.cdf(testStat,p**2)))
        if verbose:
            logging.info('Marginal Hypothesis %0.4f {} \n'.format( probR) )

        prob = 1-probR;
        if probR<alpha:
            # Change detected
            Y = X
            j = 1
            change = 1
        else:
            Y = Y + X
            j = j + 1
    else:
        prob = 0

    return Y, j, change, prob, testStat

