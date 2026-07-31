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
- `Table_5.m` and `Table_6.m` reproduce the scaling and coefficient studies.
- `functions` contains supporting codes and should not be run directly.
  The optimizer-extraction routines are in `functions/extraction`.

The default package locations are specified in
`functions/setup_paths.m`. Modify them if the packages are installed
elsewhere on your computer.

## SAGE

The SAGE calculations for Table 6 use Python with
[sageopt](https://pypi.org/project/sageopt/), NumPy, and MOSEK 10.0.
On Windows, the script detects MOSEK in its default installation directory,
`C:\Program Files\Mosek\10.0`. If MOSEK 10.0 is installed elsewhere, set
the `MOSEK10_HOME` environment variable to its installation directory.
From a terminal in this directory, run

```bash
python Table_6_sage.py
```

The script examines all six coefficient profiles, performs the two
applicable SAGE calculations, and prints a final summary table. The
`small_integer` calculation can take approximately 10--15 minutes.
To run only one applicable profile, use

```bash
python Table_6_sage.py --profile unit
python Table_6_sage.py --profile small_integer
```

The codes print their results to the MATLAB Command Window or Python
terminal and do not create separate result files.

## Citation

If you use these codes, please cite
[Log-Polynomial Optimization](https://doi.org/10.48550/arXiv.2601.02797).
Citation metadata are provided in `CITATION.cff`.
