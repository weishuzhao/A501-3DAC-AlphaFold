# ParaFold Scripts

AlphaFold: https://github.com/deepmind/alphafold

ParaFold: https://github.com/Zuricho/ParallelFold

These scripts can submit numerous feature and inference jobs of ParaFold. 

Feature jobs calculate the multiple sequence alignment and search the template, result in a `feature.pkl` file for following AlphaFold neural network prediction. This step is done by 8-core CPUs.

AlphaFold inference jobs take `feature.pkl` file as input, output predicted structure. This step is done by 1*NVIDIA A100 for each job in our project.



Requirement: slurm, ParaFold

## Usage

`batch_feature.sh` will help you automatically modify `template_feature.slurm` to generate a new `sub_feature.slurm`, and then submit this job to slurm.

After finishing all feature process, run `batch_alphafold.sh` will automatically modify `template_alphafold.slurm` to generate a new `sub_alphafold.slurm`, and then submit this job to slurm





