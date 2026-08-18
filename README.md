# pedigree_matrix

`pedigree_matrix` estimates parameter uncertainty for life-cycle inventory
data using the ecoinvent pedigree matrix. It implements the lognormal pedigree
calculation and has no runtime dependencies.

## Installation

Install the package from a local clone:

```console
python -m pip install .
```

For an editable development installation, use `python -m pip install -e .`.

## Usage

A pedigree matrix contains five quality scores in this order: reliability,
completeness, temporal correlation, geographical correlation, and further
technological correlation. Scores must be integers from 1 (best quality) to 5
(worst quality).

```python
from pedigree_matrix import PedigreeMatrix

matrix = PedigreeMatrix(version="expert")
matrix.from_numbers(1, 2, 3, 4, 5, 2)

# Uncertainty factors corresponding to the six quality scores
factors = matrix.get_values()

# Basic uncertainty and pedigree factors are expressed as GSD squared.
log_sigma = matrix.calculate(basic_uncertainty=1.05)
gsd = matrix.calculate(basic_uncertainty=1.05, output="gsd")
gsd_squared = matrix.calculate(
    basic_uncertainty=1.05, output="gsd_squared"
)
```

The package includes two factor tables:

- `version="expert"` (or the legacy `version=1`) selects the original
  expert-judgment factors and accepts an optional sixth sample-size score. If
  omitted, sample size defaults to 1.
- `version="empirical"` (or the legacy `version=2`) selects the later
  empirically based factors. It uses five indicators; sample size is excluded
  because it is expected to be represented in basic uncertainty.

The scores can instead be read from an ecoSpold-style string:

```python
matrix = PedigreeMatrix(version="empirical")
matrix.from_string("(1, 2, 3, 4, 5)")
result = matrix.calculate(basic_uncertainty=1.05, output="gsd_squared")
```

Missing string values written as `-`, `na`, or an empty field default to a
score of 1. `calculate()` returns `ln(GSD)` by default. Select `output="gsd"`
for GSD or `output="gsd_squared"` for the uncertainty factor convention used
in the ecoinvent pedigree literature.

## Scope and assumptions

The calculation combines independent basic and additional uncertainties for a
lognormal distribution. The pedigree approach is a semi-quantitative fallback
for cases where adequate statistical data are unavailable; direct statistical
estimation should be preferred when raw observations exist.

Muller et al. also derive coefficient-of-variation formulas for normal,
uniform, triangular, beta PERT, and gamma distributions. Those extensions are
not currently implemented here.

## References

- Muller, S. et al. (2014), *The application of the pedigree approach to the
  distributions foreseen in ecoinvent v3*,
  <https://doi.org/10.1007/s11367-014-0759-5>.
- Ciroth, A. et al. (2016), *Empirically based uncertainty factors for the
  pedigree matrix in ecoinvent*,
  <https://doi.org/10.1007/s11367-013-0670-5>.

## License

BSD-3-Clause. See [LICENSE](LICENSE).
