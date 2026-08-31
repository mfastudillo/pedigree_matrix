# Scientific and development roadmap

This roadmap follows a direct comparison of the package with the four papers
in `literature`. It separates correctness fixes from optional extensions so
that the existing API is not expanded on top of ambiguous factor conventions.

## Current assessment

The core legacy lognormal formula is correct when all inputs are GSD². In the
notation of Muller et al. (2014), Equation 2, the implementation calculates:

```text
total GSD² = exp(sqrt(ln(basic GSD²)² + sum(ln(additional GSD²)²)))
ln(GSD)     = 0.5 * sqrt(ln(basic GSD²)²
                         + sum(ln(additional GSD²)²))
```

The existing expert-table example reproduces the rounded GSD of 1.690 in
Muller et al. (2014). The highest-priority problem is therefore not this
formula but the fact that factor tables with different meanings are passed
through it as though they shared one convention.

## Priority 0: restore scientific correctness

### 1. Give every factor table explicit provenance and units

Replace anonymous `version_1` and `version_2` dictionaries with structured
table definitions containing at least:

- stable identifier;
- publication, year, and table number;
- supported indicators;
- factor representation: GSD, GSD², or log-space variance;
- unavailable cells;
- whether sample size is part of the table;
- whether factors are generic or sector-specific.

Suggested stable identifiers are:

```text
ecoinvent_v2_expert
ciroth_2013_empirical
muller_2016_generic
muller_2016_<sector>
```

The aliases `1`, `2`, `expert`, and `empirical` can remain temporarily for
backward compatibility, but should emit a deprecation warning where their
meaning is ambiguous.

### 2. Correct the Ciroth 2013 empirical convention

Ciroth et al. (2013), Table 6, reports GSD ratios. The current implementation
stores these values but combines them as GSD². Make the conversion explicit:

```text
additional GSD² = published GSD ratio²
```

Alternatively, convert every table to log-space variance at load time:

```text
variance = ln(GSD)²
```

The internal variance representation is preferable because contributions can
then be added directly and table conventions remain isolated at the input
boundary.

Add a regression test demonstrating the difference. With the five current
score-5 entries and neutral basic uncertainty, the existing implementation
returns total GSD approximately 1.813, whereas treating the published entries
as GSD gives approximately 3.288. The final expected value must also reflect
the unavailable-cell policy described below.

### 3. Do not invent unavailable empirical cells

Ciroth et al. (2013) reports no score-5 factor for reliability, completeness,
or geographical correlation. Remove the current silent repetition of score 4.

The default behavior should raise a clear error when an unavailable cell is
selected. If imputation is desired, require an explicit policy and return or
record which values were imputed. Possible policies include:

```text
unavailable="error"       # recommended default
unavailable="previous"    # reproduce the historical package behavior
unavailable="expert"      # use a named fallback table after conversion
```

### 4. Replace the missing-score behavior

The parser currently maps `-`, `na`, and empty fields to score 1, which means
perfect quality and no additional uncertainty. This reverses the conservative
meaning of unknown/default information in the published matrix.

Parse missing entries as `None` first, then require an explicit policy:

```text
missing="error"   # recommended default
missing="worst"   # map to score 5
missing="best"    # legacy behavior, explicit only
```

Also replace `[na\.]{1,4}` with exact accepted tokens, reject malformed values,
and decide whether matrices embedded in longer ecoSpold comments should be
found. Enable a real multiple-matrix test after that decision.

### 5. Make the five- versus six-indicator model real

Use version-specific indicator collections. A five-indicator table should
store and return five values; it should not append a neutral sample-size entry.
Keep the sixth indicator only in the explicitly legacy ecoinvent v2 table.

## Priority 1: provide a representation-safe calculation API

### 1. Use log-space variance internally

Add an API that accepts the representation used by ecoinvent v3 and by Muller
et al. (2016): variance of the logarithm.

```python
result = matrix.calculate_from_variance(
    basic_variance=0.0006,
    output="variance",
)
```

Supported outputs should include:

```text
variance     = total log-space variance
log_sigma    = sqrt(variance)
gsd          = exp(sqrt(variance))
gsd_squared  = exp(2 * sqrt(variance))
```

The existing `calculate(basic_uncertainty=...)` method should remain as a GSD²
compatibility wrapper. Do not silently reinterpret its argument.

Validate all numerical inputs as real, finite, and within the domain required
by their representation.

Document that the deterministic value of a lognormal exchange is its geometric
mean and median. Adding pedigree variance preserves that value but increases
the arithmetic mean. Do not generalize the resulting single-variable
relationship into a guaranteed direction of change for complete product-system
or LCIA Monte Carlo results.

### 2. Separate scores, factors, and results

Avoid using `factors` to mean pedigree scores. Prefer explicit objects or
properties such as:

```text
scores
additional_factors
additional_variances
total_variance
```

Expose conversions independently so callers can audit each contribution.

### 3. Test every scientific conversion

Add tests for:

- every score-to-factor cell in each table;
- GSD, GSD², and variance conversions;
- total variance as the sum of independent contributions;
- equivalence of the variance-native API and the legacy GSD² wrapper;
- the Muller et al. (2014) 1.690 example;
- missing and unavailable cells;
- five- and six-indicator table behavior;
- `NaN`, infinity, booleans, malformed strings, and invalid scores.

## Priority 2: complete the published pedigree approach

### 1. Implement basic uncertainty factors

Basic uncertainty is part of the total uncertainty, not an optional pedigree
indicator. Add structured lookup tables for exchange and process types, with
clear provenance.

Muller et al. (2016), Table 4, contains updated posterior basic factors for
some combinations. Where no updated value is available, distinguish an absent
value from a deliberate fallback to a legacy factor. Avoid a global neutral
default of 1.0 unless the caller requests it explicitly.

### 2. Add the Muller 2016 additional factors

Represent both the updated generic factors and the available sector-specific
posterior factors from Muller et al. (2016), Table 3. Sector selection should
be explicit, and the documented fallback order should be:

1. sector-specific factor when available;
2. updated generic factor;
3. a named legacy factor only when explicitly requested.

Preserve the distinction between prior, likelihood, and posterior values. The
runtime calculation should normally use posterior factors, while the other
values may be useful for provenance and future updates.

### 3. Implement qualitative matrix definitions

Store the published description for every indicator and score. This would let
the package represent the pedigree assessment itself instead of only accepting
already-decided numbers. Automated scoring should not be attempted until the
descriptions and ambiguous cases have a documented interpretation.

### 4. Support distributions other than lognormal

Implement the coefficient-of-variation method and parameter transformations
from Muller et al. (2014) for:

- normal;
- uniform;
- triangular;
- beta PERT;
- gamma with a zero location parameter.

The distribution-specific API must preserve the appropriate deterministic
value and explicitly document the assumptions used to preserve distribution
shape. Include the paper's cautions: possible negative values for unbounded
symmetric distributions, bounded-distribution approximations, and the poorer
agreement reported for gamma.

Binomial uncertainty should be identified as not modified by additional
pedigree uncertainty in that framework.

## Priority 3: updating factors and studying distribution choice

### Bayesian factor updates

Muller et al. (2016) describes how prior factors can be updated when new data
become available. Implementing this is a separate statistical feature, not a
requirement for calculating from published tables. If added, it should retain:

- prior and likelihood distribution assumptions;
- posterior distribution parameters, not only their means;
- sample and sector provenance;
- the completeness-indicator assumptions;
- warnings for non-monotonic or data-sparse results.

### Distribution-choice validation

Muller et al. (2017) demonstrates that a default distribution can affect the
median, standard deviation, and skewness of LCIA results even when product
comparisons are relatively stable. Add documentation and integration examples
that use one consistent distribution policy across compared product systems.
Do not present one distribution family as universally correct.

## Documentation and release plan

Before a correctness release:

1. Correct the empirical convention and unavailable cells.
2. Make missing-value behavior explicit.
3. Add provenance and representation metadata to all tables.
4. Add variance-native calculations and regression tests.
5. Document the backward-compatibility behavior and deprecations.

Before claiming compatibility with a particular ecoinvent release, compare
the implemented tables and conventions with that release's official data
quality guidelines. Historical paper compatibility and current ecoinvent
compatibility should be stated separately.

## Acceptance criteria

The package can claim a table or paper is implemented only when:

- every supported value is traceable to a publication and table;
- unavailable values remain distinguishable from fallback values;
- factor units are encoded and tested;
- numerical examples reproduce the source within its reported precision;
- the supported distributions and independence assumptions are explicit;
- missing data cannot silently reduce uncertainty;
- documentation distinguishes implemented behavior from planned behavior.

## Literature reviewed

- Ciroth et al. (2013), *Empirically based uncertainty factors for the
  pedigree matrix in ecoinvent*.
- Muller et al. (2014), *The application of the pedigree approach to the
  distributions foreseen in ecoinvent v3*.
- Muller et al. (2016), *Giving a scientific basis for uncertainty factors
  used in global life cycle inventory databases: an algorithm to update
  factors using new information*.
- Muller et al. (2017), *Effects of distribution choice on the modeling of
  life cycle inventory uncertainty: an assessment on the ecoinvent v2.2
  database*.

Supplementary conceptual source:

- Mutel (2013; updated 2018), *Why does the ecoinvent database love the
  lognormal distribution?* This source explains the lognormal mechanics and
  median-preserving pedigree modifier, but is not used as provenance for factor
  values, current database statistics, or non-lognormal formulas.
