# Log-Polynomial Optimization

This repository contains the numerical codes for
[Log-Polynomial Optimization](https://arxiv.org/abs/2601.02797) by
Jiyoung Choi, Jiawang Nie, Xindong Tang, and Suhan Zhong.

The paper studies optimization problems of the form

$$
\max_{x\in K}\ \sum_{i=1}^m a_i\log p_i(x),
$$

where $a_i>0$, the functions $p_i$ are polynomials, and $K$ is a
semialgebraic set.

The MATLAB codes require
[GloptiPoly](https://homepages.laas.fr/henrion/software/gloptipoly3/),
[YALMIP](https://yalmip.github.io/),
[SeDuMi](https://github.com/sqlp/sedumi), and
[MOSEK](https://www.mosek.com/). The SQP calculations also require the
MATLAB Optimization Toolbox.

## MATLAB

Set the MATLAB Current Folder to this directory. Open an `Ex_*.m` or
`Table_*.m` file and click **Run**. Each entry script automatically calls
`functions/setup_paths.m` and prints the quantities reported in the paper
at the end of the MATLAB Command Window.

- `Ex_4.m`--`Ex_12_iv.m` reproduce the numerical examples.
- `Ex_12_iv_screening.m` reproduces the candidate-dataset screening
  used to select the stress instance in Example 12(iv).
- `Table_7.m` reproduces the scaling study.
- `Table_10.m` and `Table_10_sage.py` reproduce the coefficient study.
- `Ex_4_to_12_SCIP.py` reproduces the SCIP comparisons in Tables 8--9.
- `Ex_4_to_12_SAGE.py` reproduces the SAGE applicability results in Table 11.
- `Ex_4_11_12_EM.m` reproduces the EM comparisons in Table 12.
- `functions` contains supporting codes and should not be run directly.
  The optimizer-extraction routines are in `functions/extraction`.

The default package locations are specified in
`functions/setup_paths.m`. Modify them if the packages are installed
elsewhere on your computer.

## SAGE

The SAGE calculations for Tables 10 and 11 use Python with
[sageopt](https://pypi.org/project/sageopt/), NumPy, and MOSEK 10.0.
On Windows, the script detects MOSEK in its default installation directory,
`C:\Program Files\Mosek\10.0`. If MOSEK 10.0 is installed elsewhere, set
the `MOSEK10_HOME` environment variable to its installation directory.
From a terminal in this directory, run

```bash
python Table_10_sage.py
```

The script examines all six coefficient profiles, performs the two
applicable SAGE calculations, and prints a final summary table. The
`small_integer` calculation can take approximately 10--15 minutes.
To run only one applicable profile, use

```bash
python Table_10_sage.py --profile unit
python Table_10_sage.py --profile small_integer
```

The codes print their results to the MATLAB Command Window or Python
terminal and do not create separate result files.
Objective values, gaps, and ranks should agree with the paper to the
reported precision; runtimes can vary with the computer and system load.
Each reported MATLAB moment-relaxation runtime is measured on the second
of two consecutive solves of the same instance, formulation, and order.

The SAGE applicability comparison in Table 11 and the seven launched
paternity cases are reproduced by

```bash
python Ex_4_to_12_SAGE.py
```

Use `--case classify` to print only the applicability classification or
`--case 1` through `--case 7` to solve one paternity instance.

## SCIP

The general-purpose global-solver comparisons in Tables 8 and 9 use
[PySCIPOpt](https://pypi.org/project/PySCIPOpt/), which includes SCIP.
Install it with

```bash
python -m pip install pyscipopt
```

Run all Examples 4--12 with

```bash
python Ex_4_to_12_SCIP.py
```

Use `--example 10` for one example, `--example easy` or `--example hard`
for the predefined groups, or `--time-limit 600` to change the
per-instance time limit.
The scripts print the evaluated feasible value, SCIP upper bound,
optimality gap, runtime, and node count.

## Citation

If you use these codes, please cite
[Log-Polynomial Optimization](https://doi.org/10.48550/arXiv.2601.02797).
Citation metadata are provided in `CITATION.cff`.
