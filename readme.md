# A501-3DAC-AlphaFold

## Requirements

The required dependencies for the analysis procedures are PyMOL, matplotlib, numpy, pandas and seaborn, you can install the requirements by
```python
# if you have conda environment
conda create -n A501_3DAC python=3.8
conda install pymol matplotlib numpy pandas seaborn biopython
```
You might also need jupyter notebook to run notebook files.

## 1 ParaFold

We used a parallelized version of AlphaFold which called ParaFold for high-throughput structure prediction. The installation and usage of ParaFold is available from [ParaFold](https://github.com/Zuricho/ParallelFold).

Here we provide scripts for submitting multiple structure prediction jobs using ParaFold.

## 2 Structure analysis

We used our python scripts to do structural analysis (hydrogen bond, salt bridge, disulfur bind, secondary structure, SASA, etc.)

## RMSD calculation

We wrote scripts to calculate RMSD between all predicted structures from 7 enzymes in 24 species. The python scripts is located in ./3_rmsd_calculation

## Structure sequence alignment

We conduct structure and sequence alignment for all predicted structures from 7 enzymes in 24 species. 

## Proteome-wide analysis

The data and analysis code is in ./5_proteome_wide_analysis
