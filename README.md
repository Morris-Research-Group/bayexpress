# Using bayexpress

Bayexpress is a framework designed for Differential Gene Expression Analysis of processed RNA-Seq data. It compares sequencing read counts for genes between different experiments.

## Dependencies

All python functions to use our Bayesian framework for differential gene expression analysis can be found in [bayexpress_functions.py](bayexpress_functions.py). Their dependencies are the numpy and scipy python packages.

All R functions to use our Bayesian framework for differential gene expression analysis can be found in [bayexpress_functions.R](bayexpress_functions.R). Example usage is shown in [SIM.Rmd](SIM.Rmd) which recreates two of the figures from the [manuscript](https://www.biorxiv.org/content/10.1101/2025.01.20.633909v1.full.pdf).

## Minimal working example

If you are interested in calculating Bayes factors using either python or R you can follow our [minimal working example](minimal_working_example.md).

All code is written in Python 3.10.6 (except running _DESeq2_, _edgeR_ in R), the conda environment export can be found in [environment.yml](environment.yml) containing all python and R package necessary.

$~$

# A very friendly introduction to Bayesian statistics

Should you be new to Bayesian concepts and you want to learn about it in an easier read than the papers mentioned below you can read my [friendly introduction to Bayesian statistics](A-very-friendly-Introduction-to-Bayesian-statistics.md).

$~$

-----------------------------

$~$

> [!NOTE]
> There are three manuscripts in preparation related to the work in this repository, here is a description what to find where from which paper.

1. [Hoerbst et al. 2024: Closed Form Solution for the 2-Sample Problem in Differential Gene Expression Analysis](#Closed-Form-Solution-for-the-2-Sample-Problem-in-Differential-Gene-Expression-Analysis)
2. [Hoerbst et al. 2025: A Bayesian framework for ranking genes based on their statistical evidence for differential expression](#A-Bayesian-framework-for-ranking-genes-based-on-their-statistical-evidence-for-differential-expression)
3. [Hoerbst et al. 2025: What is a differentially expressed gene?](#What-is-a-differentially-expressed-gene?)

$~$

-----------------------------

$~$

All code and data used for the bioRxiv preprint
https://www.biorxiv.org/content/10.1101/2025.01.20.633909v1

# A Bayesian framework for ranking genes based on their statistical evidence for differential expression

by 

Franziska Hoerbst, Gurpinder Singh Sidhu, Thelonious Omori, Melissa Tomkins, and Richard J. Morris

$~$


### Reproducing the Figures in the paper

### Figure 2
In Figure 2 we explore the analytical Bayesian framework using synthetic data. All plots from the paper (and more) can be found in [SIM.ipynb](SIM.ipynb).

### Figure 3
Venn diagrams of identified differentially expressed genes between WT and mutant (Snf2) yeast using _DESeq2_, _edgeR_ and our new Bayesian framework _bayexpress_. The Venn diagrams can be produced with [package_comparison_RALL.ipynb](package_comparison_RALL.ipynb). RALL stands for Replicates ALL, meaning all 44/42 clean yeast replicates have been used. The data is imported from [DGE_results](DGE_results) which in turn has been carried out via [do_DGE.ipynb](do_DGE.ipynb).


### Figure 4 and 5
Plot q-plots gene by gene and find the examples from the paper in [example_genes.ipynb](example_genes.ipynb). For more detailed analyses we also created the stalk(genes) function to look at one gene ('LSR1') or a list of genes (e.g. ['LSR1', 'YLL039C', 'YJL098W', 'YLL026W']) to get read counts in WT and Snf2 mutant, q-value visualisation, and more [stalking_genes.ipynb](stalking_genes.ipynb).

### Figure 6 and 7
In Figure 6 and 7 we explore the analytical Bayesian framework using synthetic data a bit further to Figure 2. All plots from the paper (and more) can be found in [SIM.ipynb](SIM.ipynb).

### Figure 8
Can be reproduced via [package_comparison_RALL.ipynb](package_comparison_RALL.ipynb).

### Figure 9
Can be reproduced via [investigating_length_bias.ipynb](investigating_length_bias.ipynb).

### Table 1
Can be reproduced via [package_comparison_RALL.ipynb](package_comparison_RALL.ipynb).

$~$

-----------------------------

$~$

All code and data used for the bioRxiv preprint
https://www.biorxiv.org/content/10.1101/2025.01.31.635902v1.article-metrics

# What is a differentially expressed gene?

by 

Franziska Hoerbst, Gurpinder Singh Sidhu, Melissa Tomkins, and Richard J. Morris

$~$


### Figures 1 (and 4-12)
We found the representation in these Figures so useful, we automized it -- it can be found in [stalking_genes.ipynb](stalking_genes.ipynb). Use the stalk(genes) function to look at one gene ('LSR1') or a list of genes (e.g. ['LSR1', 'YLL039C', 'YJL098W', 'YLL026W']) to get read counts in WT and Snf2 mutant, q-value visualisation, and more. 

In [explore_clean_yeast.ipynb](explore_clean_yeast.ipynb) we identified the example genes we discuss in the manuscript. 

In [example_genes_WIADEG.ipynb](example_genes_WIADEG.ipynb) all of these examples can be found. 

In [example_genes_WIADEG_SNF2.ipynb](example_genes_WIADEG_SNF2.ipynb) we investigate the SNF2 targets discussed in the Appendix. 



### Figure 2 (and 13-16)
For the control experiments between wild-type and wild-type, DGE analysis has only been done in _bayexpress_. The analysis has been done in [do_DGE.ipynb](do_DGE.ipynb), we recycled bootstrapping data created for the comparison experiment. The figures have been created in [CONTROL_R3_BF1.ipynb](CONTROL_R3_BF1.ipynb), [CONTROL_R3_BF10.ipynb](CONTROL_R3_BF10.ipynb), and [CONTROL_R3_BF100.ipynb](CONTROL_R3_BF100.ipynb) for 3 replicates and [CONTROL_R10_BF1.ipynb](CONTROL_R10_BF1.ipynb), [CONTROL_R10_BF10.ipynb](CONTROL_R10_BF10.ipynb), and [CONTROL_R10_BF100.ipynb](CONTROL_R10_BF100.ipynb) for 10 replicates. Figure 16 can be reproduced in [CONTROL_R3_BF1_CIG.ipynb](CONTROL_R3_BF1_CIG.ipynb) and [CONTROL_R10_BF1_CIG.ipynb](CONTROL_R10_BF1_CIG.ipynb).

### Figure 3
Calculating Bayes factors for consistency (BF_k1) and performing bootstrapping experiments to identify consistently inconsistent genes (list of genes marked* as CIG = Consistently Inconsistent Gene) has been carried out using [explore_clean_yeast_consistency.ipynb](explore_clean_yeast_consistency.ipynb). We explored these sets of genes in later parts of the notebook. 

$~$

-----------------------------

$~$

All code and data used for the 2024 arXiv preprint (arXiv:2406.19989 [stat.ME]) https://doi.org/10.48550/arXiv.2406.19989
# Closed Form Solution for the 2-Sample Problem in Differential Gene Expression Analysis

by 

Franziska Hoerbst, Gurpinder Singh Sidhu, Melissa Tomkins, and Richard J. Morris

$~$


### Figure 3
In Figure 3 we explore the analytical Bayesian framework using synthetic data. All plots from the paper (and more) can be found in [SIM.ipynb](SIM.ipynb).

$~$

-----------------------------

$~$


# File by file

[do_DGE.ipynb](do_DGE.ipynb) is a notebook doing all Differential Gene Expression analysis discussed in the paper. For running DESeq2 and edgeR we are calling R scripts (e.g. 'Do_DGE-DESeq2.R') to carry out the analysis. This is where the parameters for the packages are set. 

For bootstrapping experiments we created 100 shuffled data sets consisting of 3, 6, 12, and 20 replicates of the pool of 44/42 which we used for package comparison. This was done in [comparison_data.ipynb](comparison_data.ipynb).


