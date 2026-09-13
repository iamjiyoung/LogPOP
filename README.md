# LogPOP

This repository contains code for reproducing the numerical experiments in
[Log-Polynomial Optimization](https://arxiv.org/abs/2601.02797) by Jiyoung Choi,
Jiawang Nie, Xindong Tang, and Suhan Zhong.

The paper studies optimization problems whose objectives are weighted sums of
logarithms of polynomial functions and develops moment relaxations for solving
them.

## Requirements

- MATLAB R2025b
- MOSEK 10.0
- GloptiPoly 3
- YALMIP
- SeDuMi

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

- `Ex_5_3_standard.m`, `Ex_5_3_LME.m`: Example 5.3, standard and LME relaxations
- `Ex_6_1.m`: Example 6.1
- `Ex_6_2_standard.m`, `Ex_6_2_LME.m`: Example 6.2, standard and LME relaxations
- `Ex_6_3_1.m`, `Ex_6_3_2.m`: Example 6.3(i) and (ii)
- `Ex_6_4.m`, `Ex_6_5.m`: Examples 6.4 and 6.5
- `Ex_6_6_1.m`--`Ex_6_6_9.m`: the nine paternity instances in Example 6.6
- `Ex_6_7_1.m`--`Ex_6_7_3.m`: the three LCM instances in Example 6.7
- `Ex_6_7_1_data.m`--`Ex_6_7_3_data.m`: data loaded by the corresponding LCM scripts
- `functions/`: extraction functions used by the example scripts; these files are not run directly

The active coefficient vector in each `Ex_6_2` script reproduces Example
6.2(i). The five alternative vectors, in their listed order, correspond to
instances #1--#5 in Example 6.2(ii). Uncomment one vector in each script and
leave only that assignment active. In `Ex_6_2_standard.m`, use relaxation
orders 4, 4, 5, 6, and 5 for instances #1--#5, respectively. In
`Ex_6_2_LME.m`, use relaxation order 3 for all five instances.

The extracted atoms, and hence gaps evaluated at those atoms, may vary
slightly because they are obtained from a numerical eigendecomposition. The
reported relaxation values, flat-truncation orders, and numerical ranks are
the primary reproducible outputs.
