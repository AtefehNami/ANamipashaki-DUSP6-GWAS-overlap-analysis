# DUSP6 GWAS Overlap Analysis

This repository contains the MATLAB code and input data used for GWAS overlap permutation analysis in the DUSP6 study.

## How to run

1. Download and extract this repository.
2. Open MATLAB and set the Current Folder to the main directory containing `startup.m`, `code`, and `data`.
3. Run the following commands:

```matlab
startup
run_DUSPoverlap_permutations
```

The script automatically creates an `output` folder for the generated results.

The analysis is configured to run 100,000 permutations. Because no fixed random seed is specified, permutation p-values may vary slightly between runs.

The original folder structure must be maintained because the scripts use relative file paths.
