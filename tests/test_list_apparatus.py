"""Unit-testing for computing the list apparatus."""

import unittest
import os
from pyfroger.list_apparatus import ListApparatus


DATASETS = os.path.join(os.path.dirname(__file__), "test_datasets")


class TestListApparatus(unittest.TestCase):
    """Tests that building the list apparatus works properly.
    """

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "normal_tradition.csv"),
            ref_witness="MS_A"
        )

    def test_variant_group(self):
        """Tests that obtaining the variant group works properly."""
        # Check that the variant group at locus_2 is properly computed
        first_error = self.apparatus.error_groups()[1]
        self.assertEqual(
            first_error.reading,
            "vaillant"
        )
        self.assertEqual(
            first_error.witnesses,
            {"MS_C", "MS_D"}
        )

    def test_compute_frequency(self):
        """Tests that computing the frequency behaves as expected."""
        # Check the frequency of the variant group "MS_D, MS_C
        self.assertEqual(
            self.apparatus.group_by_frequency()[frozenset({"MS_D", "MS_C"})]
            ["frequency"],
            3)


if __name__ == "__main__":
    unittest.main()
