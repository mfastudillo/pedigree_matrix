from .data import version_1, version_2
from .from_string import find_pedigree_matrix
from pprint import pformat
import math
import numbers
import warnings


class PedigreeMatrix(object):
    version_aliases = {
        1: 1,
        2: 2,
        "expert": 1,
        "empirical": 2,
    }

    labels = (
        "reliability",
        "completeness",
        "temporal correlation",
        "geographical correlation",
        "further technological correlation",
        "sample size",
    )

    def __init__(self, version=1):
        if isinstance(version, bool):
            raise ValueError("Version must be 1, 2, 'expert', or 'empirical'")
        try:
            self.version = self.version_aliases[version]
        except (KeyError, TypeError):
            raise ValueError(
                "Version must be 1, 2, 'expert', or 'empirical'"
            ) from None
        self.factors = {}

    def from_numbers(self, *args):
        """Set pedigree scores in matrix order.

        Scores must be integers from 1 (best quality) to 5 (worst quality).
        The sixth, sample-size score belongs only to the expert-judgment table.
        """
        if len(args) not in (5, 6):
            raise ValueError("Must provide either 5 or 6 scores")

        if self.version == 2 and len(args) == 6:
            warnings.warn(
                "Sample size is not part of the empirical five-indicator "
                "matrix; the sixth score is ignored",
                UserWarning,
                stacklevel=2,
            )
            args = args[:5] + (1,)
        if len(args) == 5:
            args = args + (1,)

        for label, score in zip(self.labels, args):
            if (
                isinstance(score, bool)
                or not isinstance(score, numbers.Integral)
                or not 1 <= score <= 5
            ):
                raise ValueError(
                    "Score for '{}' must be an integer from 1 to 5".format(label)
                )

        factors = {}
        for index, factor in enumerate(args):
            factors[self.labels[index]] = factor
        self.factors = factors

    def from_string(self, string):
        factors = find_pedigree_matrix(string)
        if not factors:
            raise ValueError("Can't find Pedigree Matrix factors")
        if self.version == 2:
            # ecoSpold strings can contain the legacy sample-size field. It is
            # deliberately neutral in the empirical five-indicator matrix.
            factors = factors[:5]
        self.from_numbers(*factors)

    def calculate(
        self,
        basic_uncertainty=1.0,
        as_geometric_sigma=None,
        *,
        output="log_sigma",
    ):
        """Combine basic and pedigree uncertainty for a lognormal model.

        ``basic_uncertainty`` is the squared geometric standard deviation
        (GSD squared), as are the factors in the pedigree tables.

        ``output`` can be ``"log_sigma"`` (the default), ``"gsd"``, or
        ``"gsd_squared"``. The legacy ``as_geometric_sigma`` argument is
        retained for compatibility but is deprecated because ``True`` has
        historically returned GSD squared rather than GSD.
        """
        if (
            isinstance(basic_uncertainty, bool)
            or not isinstance(basic_uncertainty, numbers.Real)
            or basic_uncertainty < 1
        ):
            raise ValueError("basic_uncertainty must be a GSD-squared value >= 1")

        if as_geometric_sigma is not None:
            if not isinstance(as_geometric_sigma, bool):
                raise ValueError("as_geometric_sigma must be a boolean")
            if output != "log_sigma":
                raise ValueError(
                    "Use either as_geometric_sigma or output, not both"
                )
            warnings.warn(
                "as_geometric_sigma is deprecated; use output='gsd_squared' "
                "for the historical True behavior",
                DeprecationWarning,
                stacklevel=2,
            )
            output = "gsd_squared" if as_geometric_sigma else "log_sigma"

        if output not in ("log_sigma", "gsd", "gsd_squared"):
            raise ValueError("output must be 'log_sigma', 'gsd', or 'gsd_squared'")

        values = [basic_uncertainty] + self.get_values()
        log_sigma = math.sqrt(sum(math.log(x) ** 2 for x in values)) / 2
        if output == "gsd":
            return math.exp(log_sigma)
        if output == "gsd_squared":
            return math.exp(2 * log_sigma)
        return log_sigma

    def get_values(self):
        if not self.factors:
            raise ValueError("Must provide Pedigree Matrix factors")
        data = version_1 if self.version == 1 else version_2
        return [data[key][index - 1] for key, index in self.factors.items()]

    def __repr__(self):
        if not self.factors:
            return u"Empty Pedigree Matrix"
        else:
            return pformat(self.factors)
