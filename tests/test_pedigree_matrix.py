import math
import unittest
import warnings

from pedigree_matrix import PedigreeMatrix


class PedigreeMatrixTest(unittest.TestCase):
    def test_named_version_aliases(self):
        self.assertEqual(PedigreeMatrix("expert").version, 1)
        self.assertEqual(PedigreeMatrix("empirical").version, 2)

    def test_rejects_invalid_version(self):
        for version in ("unknown", True):
            with self.subTest(version=version), self.assertRaises(ValueError):
                PedigreeMatrix(version)

    def test_rejects_invalid_scores(self):
        for score in (0, 6, 1.5, True):
            with self.subTest(score=score), self.assertRaises(ValueError):
                matrix = PedigreeMatrix()
                matrix.from_numbers(score, 1, 1, 1, 1)

    def test_empirical_matrix_ignores_sample_size(self):
        matrix = PedigreeMatrix("empirical")
        with self.assertWarns(UserWarning):
            matrix.from_numbers(1, 2, 3, 4, 5, 5)
        self.assertEqual(matrix.get_values()[-1], 1.0)

    def test_article_lognormal_example(self):
        matrix = PedigreeMatrix("expert")
        matrix.from_numbers(5, 5, 5, 5, 5)

        # Table 5 reports a total GSD of 1.690 (rounded).
        total_gsd = matrix.calculate(1.279**2, output="gsd")
        self.assertAlmostEqual(total_gsd, 1.690, places=2)

    def test_output_conventions(self):
        matrix = PedigreeMatrix()
        matrix.from_numbers(1, 2, 3, 4, 5)

        log_sigma = matrix.calculate(1.05)
        gsd = matrix.calculate(1.05, output="gsd")
        gsd_squared = matrix.calculate(1.05, output="gsd_squared")

        self.assertAlmostEqual(gsd, math.exp(log_sigma))
        self.assertAlmostEqual(gsd_squared, gsd**2)

    def test_legacy_output_argument(self):
        matrix = PedigreeMatrix()
        matrix.from_numbers(1, 1, 1, 1, 1)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            legacy = matrix.calculate(1.05, as_geometric_sigma=True)
        explicit = matrix.calculate(1.05, output="gsd_squared")

        self.assertEqual(legacy, explicit)
        self.assertTrue(
            any(item.category is DeprecationWarning for item in caught)
        )

    def test_rejects_invalid_calculation_options(self):
        matrix = PedigreeMatrix()
        matrix.from_numbers(1, 1, 1, 1, 1)

        for basic_uncertainty in (0, 0.5, -1, "1.2", True):
            with self.subTest(value=basic_uncertainty), self.assertRaises(
                ValueError
            ):
                matrix.calculate(basic_uncertainty)

        with self.assertRaises(ValueError):
            matrix.calculate(output="variance")


if __name__ == "__main__":
    unittest.main()
