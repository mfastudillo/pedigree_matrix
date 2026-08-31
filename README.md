# pedigree_matrix

`pedigree_matrix` estimates parameter uncertainty for life-cycle inventory
(LCI) data using pedigree-matrix scores. The package currently implements the
legacy ecoinvent lognormal calculation and has no runtime dependencies.

## Scientific status

The legacy expert-factor calculation is implemented consistently with
Muller et al. (2014), Equation 2, when the basic uncertainty and all additional
factors are supplied as squared geometric standard deviations (GSD²). The test
suite reproduces the paper's rounded total GSD example of 1.690.

The `empirical` factor table is experimental and should not currently be used
for quantitative results. Its values reproduce the GSD ratios in Table 6 of
Ciroth et al. (2013), but the current calculation treats them as GSD²
contributors. This understates the resulting uncertainty. In addition, the
paper reports no value for score 5 for reliability, completeness, or
geographical correlation; the current table silently repeats score 4 for these
cells. Both issues are scheduled for correction.

This package does not yet represent the complete approach described by the
papers in `literature`. In particular, it does not provide native variance
inputs, basic-uncertainty lookup tables, sector-specific factors, or the
extensions to distributions other than lognormal.

## Installation

Install the package from a local clone:

```console
python -m pip install .
```

For an editable development installation, use `python -m pip install -e .`.

## Current usage

Scores are supplied in this order:

1. reliability;
2. completeness;
3. temporal correlation;
4. geographical correlation;
5. further technological correlation;
6. sample size, only for the legacy ecoinvent v2 expert table.

Each score must be an integer from 1 (best quality) to 5 (worst quality).

```python
from pedigree_matrix import PedigreeMatrix

# The scientifically validated path in the current release is the legacy
# expert table combined with the lognormal calculation.
matrix = PedigreeMatrix(version="expert")
matrix.from_numbers(1, 2, 3, 4, 5)

# In this API, basic uncertainty and expert pedigree factors are GSD².
log_sigma = matrix.calculate(basic_uncertainty=1.05)
gsd = matrix.calculate(basic_uncertainty=1.05, output="gsd")
gsd_squared = matrix.calculate(
    basic_uncertainty=1.05,
    output="gsd_squared",
)
```

`calculate()` returns `ln(GSD)` by default. The other output choices are
`"gsd"` and `"gsd_squared"`.

### Parameterizing a lognormal distribution

For a positive deterministic exchange amount `amount`, treat `amount` as the
median of the lognormal distribution. The two parameters of the underlying
normal distribution are then:

```text
mu    = ln(amount)
sigma = ln(GSD) = matrix.calculate(...)
```

Do not pass the GSD or GSD² output as `sigma`; common Python samplers expect
the standard deviation in natural-log space. For example, the Python standard
library can be used without adding a package dependency:

```python
import math
import random

from pedigree_matrix import PedigreeMatrix

amount = 2.5  # Deterministic amount, interpreted as the median

matrix = PedigreeMatrix(version="expert")
matrix.from_numbers(1, 2, 3, 4, 5)

# basic_uncertainty uses the legacy GSD² convention.
log_sigma = matrix.calculate(basic_uncertainty=1.05)
log_mu = math.log(amount)

sample = random.lognormvariate(log_mu, log_sigma)
```

The equivalent parameterization in NumPy is
`rng.lognormal(mean=log_mu, sigma=log_sigma)`. In SciPy it is
`scipy.stats.lognorm(s=log_sigma, scale=amount)` with the default `loc=0`.
NumPy and SciPy are optional and are not installed by this package.

A two-parameter lognormal distribution cannot represent zero or cross zero.
For a negative exchange, apply the distribution to `abs(amount)` and restore
the fixed negative sign to each sample. A zero amount requires a different
uncertainty model or should remain deterministic.

The arithmetic mean of the resulting positive distribution is
`amount * exp(log_sigma**2 / 2)`; it is larger than the deterministic median
whenever `log_sigma` is nonzero.

#### Using `stats_arrays` and Brightway

Brightway uses [`stats_arrays`](https://github.com/brightway-lca/stats_arrays)
to represent and sample uncertainty distributions. It can be installed from
PyPI as an optional dependency:

```console
python -m pip install stats-arrays
```

For `stats_arrays.LognormalUncertainty`, `loc` is the mean of the underlying
normal distribution and `scale` is its standard deviation. Therefore, this
package maps directly to the expected fields:

```python
import math

from pedigree_matrix import PedigreeMatrix
from stats_arrays import (
    LognormalUncertainty,
    MCRandomNumberGenerator,
    UncertaintyBase,
)

amount = -2.5  # Deterministic amount, interpreted as the median

matrix = PedigreeMatrix(version="expert")
matrix.from_numbers(1, 2, 3, 4, 5)
log_sigma = matrix.calculate(basic_uncertainty=1.05)

# A Brightway exchange dictionary uses "uncertainty type" with a space.
brightway_uncertainty = {
    "amount": amount,
    "uncertainty type": LognormalUncertainty.id,
    "loc": math.log(abs(amount)),
    "scale": log_sigma,
    "negative": amount < 0,
}

# A direct stats_arrays parameter dictionary uses uncertainty_type.
params = UncertaintyBase.from_dicts(
    {
        "uncertainty_type": LognormalUncertainty.id,
        "loc": math.log(abs(amount)),
        "scale": log_sigma,
        "negative": amount < 0,
    }
)

generator = MCRandomNumberGenerator(params, seed=42)
sample = next(generator)[0]
```

The `negative` flag restores a negative sign after sampling the positive
lognormal magnitude. In Brightway, it describes the sign of the stored
exchange amount; it is separate from any sign change applied when an exchange
is inserted into a technosphere matrix.

`stats_arrays` requires a strictly positive lognormal `scale`. If
`log_sigma == 0`, represent the amount as deterministic uncertainty instead of
creating a lognormal parameter row. See the
[`stats_arrays` parameter-array documentation](https://stats-arrays.readthedocs.io/en/latest/)
for heterogeneous arrays, bounds, and additional distribution types.

### Available table aliases

- `version="expert"` or `version=1` selects the legacy ecoinvent v2
  expert-judgment factors. Five scores are accepted; an optional sixth
  sample-size score is supported for compatibility and defaults to 1 when
  omitted.
- `version="empirical"` or `version=2` selects the tentative generic factors
  from Ciroth et al. (2013). This option has the known convention and
  missing-cell problems described above and is not recommended for
  calculations.

These aliases describe historical factor sources; `expert` should not be
interpreted as a claim of compatibility with every ecoinvent v3 release.

### Reading score strings

Scores can be read from a tuple-like string:

```python
matrix = PedigreeMatrix(version="expert")
matrix.from_string("(1, 2, 3, 4, 5)")
```

The current parser converts `-`, `na`, and empty fields to score 1. This is a
legacy behavior and is potentially unsafe: unknown data quality should not be
silently interpreted as ideal data quality. Validate or replace missing scores
before calculating. A future release will require an explicit missing-value
policy.

## Scope and assumptions

The implemented calculation assumes that basic uncertainty and each pedigree
contribution are independent and lognormally distributed. It combines their
log-space variances while preserving the deterministic value.

For a lognormal exchange, the deterministic value is the geometric mean and
median, not the arithmetic mean. Increasing pedigree uncertainty preserves
that median, but increases the arithmetic mean by increasing the log-space
variance. Consequently, a deterministic result based on median exchange
amounts need not equal the mean of Monte Carlo results. This relationship is
exact for an individual positive lognormal variable; it is not a guarantee
about the direction of the difference for a complete product system, where
signed exchanges, aggregation, and matrix calculations also matter.

The pedigree approach is a semi-quantitative fallback for cases where adequate
statistical information is unavailable. Direct estimation from representative
raw observations should be preferred when possible.

Important parts of the published approach that are not implemented include:

- native variance-of-log-transformed-data input and output;
- basic uncertainty factors selected by exchange and process type;
- the updated generic and sector-specific factors from Muller et al. (2016);
- normal, uniform, triangular, beta PERT, and gamma calculations based on
  coefficient of variation, as derived by Muller et al. (2014);
- the Bayesian procedure for updating uncertainty factors;
- qualitative matrix-cell descriptions for deriving scores from source data.

Muller et al. (2017) also shows that the chosen probability distribution can
affect the median, dispersion, and skewness of LCA results. Analyses that
compare product systems should use a consistent uncertainty-distribution
policy across those systems.

See [NEXT_STEPS.md](NEXT_STEPS.md) for the proposed correction and extension
roadmap.

## References

- Ciroth, A., Muller, S., Weidema, B., and Lesage, P. (2013), *Empirically
  based uncertainty factors for the pedigree matrix in ecoinvent*.
  <https://doi.org/10.1007/s11367-013-0670-5>.
- Muller, S., Lesage, P., Ciroth, A., Mutel, C., Weidema, B. P., and Samson,
  R. (2014), *The application of the pedigree approach to the distributions
  foreseen in ecoinvent v3*.
  <https://doi.org/10.1007/s11367-014-0759-5>.
- Muller, S., Lesage, P., and Samson, R. (2016), *Giving a scientific basis for
  uncertainty factors used in global life cycle inventory databases: an
  algorithm to update factors using new information*.
  <https://doi.org/10.1007/s11367-016-1098-5>.
- Muller, S., Mutel, C., Lesage, P., and Samson, R. (2017), *Effects of
  distribution choice on the modeling of life cycle inventory uncertainty: an
  assessment on the ecoinvent v2.2 database*.
  <https://doi.org/10.1111/jiec.12574>.
- Mutel, C. (2013; updated 2018), *Why does the ecoinvent database love the
  lognormal distribution?* <https://chris.mutel.org/ecoinvent-lognormal.html>.
  This is a supplementary conceptual explanation, not the provenance for any
  factor table or a source of current ecoinvent exchange statistics.

Copies of these papers are included in the `literature` directory for project
development and review. The Mutel article is an external web reference.

## License

BSD-3-Clause. See [LICENSE](LICENSE).
