# Using bayexpress

_Bayexpress_ is a framework designed for Differential Gene Expression Analysis of processed RNA-Seq data. It compares sequencing read counts for genes between different experiments.


## Dependencies

All basic python functions to use our Bayesian framework for differential gene expression analysis can be found in [bayexpress_functions.py](bayexpress_functions.py). Their dependencies are the numpy, pandas and scipy python packages.

All basic R functions to use our Bayesian framework for differential gene expression analysis can be found in [bayexpress_functions.R](bayexpress_functions.R). 

$~$


# Minimal working example

See also the jupyter notebook [minimal_working_example.ipynb](minimal_working_example.ipynb).


## Calculating Bayes factors for differential gene expression (BF_21)

### Passing Single Values

    BF_gene = get_BF(N_1, n_1, N_2, n_2)
    print(BF_gene)

### Example test frame

Read counts for 5 genes in two different experiments (cound be conditions, developmental timepoints, ...) with 3 replicates each. 

| genes | exp 1 rep 1 | exp 1 rep 2 | exp 1 rep 3 |  exp 2 rep 1 | exp 2 rep 2 | exp 2 rep 3 |
| :----------- | :------: |  :------: |  :------: | :------: | :------: |------------: |
| gene1 | 10 | 12 | 8 | 100 | 120 | 80 |
| gene2 | 100 | 120 | 80 | 100 | 120 | 80 |
| gene3 | 1000 | 1200 | 1100 | 1200 | 1300 | 1400 |
| gene4 | 10000 | 12000 | 11000 | 12000 | 13000 | 14000 |
| gene5 | 1000 | 1000 | 1000 | 1000 | 1000 |

        
### Passing a data frame or array

    in_data = pd.DataFrame({'genes': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
                            'exp 1 rep 1': [10, 100, 1000, 10000, 1000],
                            'exp 1 rep 2': [12, 120, 1200, 12000, 1000],
                            'exp 1 rep 3': [8, 80, 1100, 11000, 1000],
                            'exp 2 rep 1': [100, 100, 1200, 12000, 1000],
                            'exp 2 rep 2': [120, 120, 1300, 13000, 1000],
                            'exp 2 rep 3': [80, 80, 1400, 14000, 1000]
                            })

    print(in_data)

    out_data = pd.DataFrame({'genes': in_data.genes})

    # n_1 is the sum of all read counts in experiment 1 (replicates 1-3)
    n_1 = in_data.iloc[:,1:4].sum(axis=1)

    # n_2 is the sum of all read counts in experiment 2 (replicates 1-3)
    n_2 = in_data.iloc[:,4:].sum(axis=1)

    # N_1 and N_2 are the total number of reads in experiments 1 and 2, respectively
    N_1 = in_data.iloc[:,1:4].sum(axis=1).sum()
    N_2 = in_data.iloc[:,4:].sum(axis=1).sum()

    out_data['BF'] = get_BF(N_1, n_1, N_2, n_2)
    out_data['FC'] = get_FC(N_1, n_1, N_2, n_2)

    print(out_data)

## Calculating Bayes factors for consistency of replicates (BF_IC)

### Passing a data frame

    in_data = pd.DataFrame({'genes': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
                            'exp 1 rep 1': [10, 100, 1000, 10000, 1000],
                            'exp 1 rep 2': [12, 120, 1200, 12000, 1000],
                            'exp 1 rep 3': [8, 80, 1100, 11000, 1000],
                            'exp 2 rep 1': [100, 100, 1200, 12000, 1000],
                            'exp 2 rep 2': [120, 120, 1300, 13000, 1000],
                            'exp 2 rep 3': [80, 80, 1400, 14000, 1000]
                            })

    print(in_data)

    out_data = pd.DataFrame({'genes': in_data.genes})
    
    # input data is a data frame containing all read counts for all genes (rows) and replicates (colums)

    out_data['BF_IC_1'] = get_BF_IC(in_data.iloc[:,1:4])
    out_data['BF_IC_2'] = get_BF_IC(in_data.iloc[:,4:])

    print(out_data)

