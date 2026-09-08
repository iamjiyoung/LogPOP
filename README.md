# Log-polynomial optimization

This folder contains the code used for the numerical experiments in the paper.

## Requirements

- MATLAB R2025b
- MATLAB Optimization Toolbox
- MOSEK 10.0
- GloptiPoly 3
- YALMIP
- SeDuMi
- Python with PySCIPOpt 6.2.1 and SCIP 10.0.2 for the SCIP comparison

## MATLAB setup

Start MATLAB, change to this folder, and add the required packages and local
extraction functions to the MATLAB path. Replace the package paths below with
their locations on your computer.

```matlab
addpath(genpath('path/to/gloptipoly3'))
addpath(genpath('path/to/YALMIP'))
addpath(genpath('path/to/SeDuMi'))
addpath('path/to/MOSEK/10.0/toolbox/r2017a')
addpath(fullfile(pwd,'functions'))
```

Each experiment can then be reproduced by opening the corresponding script and
clicking **Run**, or by entering its file name without `.m` in the Command
Window.

## Files

- `Ex_5_3_woLME.m`, `Ex_5_3_wLME.m`: Example 5.3, standard and LME relaxations
- `Ex_6_1.m`: Example 6.1
- `Ex_6_2_woLME.m`, `Ex_6_2_wLME.m`: Example 6.2, standard and LME relaxations
- `Ex_6_3_1.m`, `Ex_6_3_2.m`: Example 6.3(i) and (ii)
- `Ex_6_4.m`, `Ex_6_5.m`: Examples 6.4 and 6.5
- `Ex_6_6_1.m`--`Ex_6_6_9.m`: the nine paternity instances in Example 6.6
- `Ex_6_7_1.m`--`Ex_6_7_3.m`: the first three LCM instances in Example 6.7
- `Ex_6_7_4.m`: the non-saturated LCM instance and 100-start EM comparison
- `Table_5.m`: the scaling experiment reported in Table 5
- `Table_6.m`: the weight experiment reported in Table 6
- `Table_7_SCIP.py`: the comparison reported in Table 7
- `Section_5_rank2_LME.m`: the rank-two LME illustration in Section 5
- `functions/`: extraction and Table 5 helper functions; these files are not run directly

The active coefficient vector in each `Ex_6_2` script reproduces Example
6.2(i). To reproduce the five instances in Example 6.2(ii), uncomment one of
the five alternative vectors in each script and leave only that assignment
active.

## SCIP comparison

Run all four comparison cases from a terminal with

```text
python Table_7_SCIP.py --example all --time-limit 600 --relative-gap 1e-6
```

Wall-clock times, and the SCIP upper bound at the time limit, can vary with the
computing environment.
