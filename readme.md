# A501-3DAC-AlphaFold

## Requirements

The required dependencies for the analysis procedures are PyMOL, matplotlib, numpy, pandas and seaborn, you can install the requirements by
```python
# if you have conda environment
conda create -n A501_3DAC python=3.8
conda install pymol matplotlib numpy pandas seaborn biopython scipy
```
You might also need jupyter notebook to run notebook files.

## 1 ParaFold

We used a parallelized version of AlphaFold which called ParaFold for high-throughput structure prediction. The installation and usage of ParaFold is available from [ParaFold](https://github.com/Zuricho/ParallelFold). Here we provide scripts for submitting multiple structure prediction jobs using ParaFold.

## 2 Result analysis

We used our python scripts (in jupyter notebook format) to do structural analysis (hydrogen bond, salt bridge, disulfur bind, secondary structure, SASA, etc.). In this repo, we provided a demo to run the notebook using 20 example AlphaFold structure.

Detail guide can be found in notebook.

## 3 Structure sequence alignment

We conduct structure and sequence alignment for all predicted structures from 7 enzymes in 24 species. 

First, you should download the prediction results from https://zenodo.org/record/6387901#.Yu15YHYzaUk

Then, the jupyter notebook `3_structure_sequence_alignment/result_analysis.ipynb` will help you calculate the sequence identity and structure RMSD, results saved in `./3_structure_sequence_alignment/result` folder

Finally `3_structure_sequence_alignment/figure.ipynb` is used for plot figures in our context.

## 4 Proteome-wide analysis

The proteome-wide structural feature data in csv format and analysis code is in `./4_proteome_wide_analysis`
