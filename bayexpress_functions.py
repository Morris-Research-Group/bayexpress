
# Please read, follow and use our LICENSE and CITATION.

import scipy.special as sc
import pandas as pd
import numpy as np


# Calculating log_10 Bayes factors for differential gene expression (BF_21)

# Optional hyper-parameters (here we use 1 for a flat prior)
u_1 = 1
u_2 = 1

def get_BF_21(N_1, # Total number of reads for genes in experiment 1 
              n_1, # Number of reads for a gene in experiment 1
              N_2, # Total number of reads for genes in experiment 2
              n_2 # Number of reads for a gene in experiment 2
              ):
    
    # betaln returns the natural log of the beta function
    # hence we after convert natural log to log_10
    return (sc.betaln( u_1 + n_1, u_2 + N_1 - n_1) + sc.betaln( u_1 + n_2, u_2 + N_2 - n_2) - sc.betaln( u_1 + n_1 + n_2, u_2 + N_1 - n_1 + N_2 - n_2)) / np.log(10) 


# ratio of expression 
# calculating log_2 fold change

def get_FC(N_1, # Total number of reads for genes in experiment 1 
              n_1, # Number of reads for a gene in experiment 1
              N_2, # Total number of reads for genes in experiment 2
              n_2 # Number of reads for a gene in experiment 2
              ):
    
    rate_1 = (u_1 + n_1) / (u_2 + N_1 - n_1)
    rate_2 = (u_1 + n_2) / (u_2 + N_2 - n_2)

    return np.log2(rate_2 / rate_1)

# calculating q (following Laplace's rule of succession)
def get_q(n, # Number of reads mapping to a gene
          N # Total number of reads in an experiment
          ):

    
    return (n+1)/(N+2)


# calculating log_10 Bayes factors (BF_k1) for consistency 
def get_BF_k1(data):
    # this range is irrelevant if we want to do all replicates
    k = len(data.columns)

    evidence2 = np.full(len(data), 0)

    # iterating over j until k
    for col in data.columns[1:k]: 
        n_j = data[col]
        # print(n_j, 'n_j')
        N_j = sum(data[col])
        # print(N_j, 'N_j')
        evidence2 = evidence2 + sc.betaln(u_1 + n_j, u_2 + N_j - n_j)

    N = sum(data.iloc[:,1:k].sum(axis=0, numeric_only=True))
    n_i = data.iloc[:,1:k].sum(axis=1, numeric_only=True)

    # betaln returns the natural log of the beta function

    evidence1 = sc.betaln( u_1 + n_i, u_2 + N - n_i)

    # convert natural log to log_10
    return (evidence2 - evidence1) / np.log(10) 
