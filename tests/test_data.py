"""Unit tests for the ListApparatus class."""
import os
import unittest
from pyfroger.backup.data import ListApparatus


DATASETS = os.path.join(os.path.dirname(__file__), "test_datasets")


class TestFromCSVRawAlignment(unittest.TestCase):
    """Tests for from_csv with a raw alignment matrix."""

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "raw_alignment.csv")
        )

    def test_witnesses_extracted(self):
        self.assertEqual(
            self.apparatus.witnesses,
            {"MS_A", "MS_B", "MS_C", "MS_D"},
        )

    def test_number_of_groups(self):
        # locus_1: 2 readings, locus_2: 2 readings, locus_3: 2 readings
        self.assertEqual(len(self.apparatus.groups), 6)

    def test_grouping_at_locus_1(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_1"]
        readings = {g.reading: g.witnesses for g in groups}
        self.assertEqual(readings["le cheval blanc"],
                         frozenset({"MS_A", "MS_B"}))
        self.assertEqual(readings["le cheval noir"],
                         frozenset({"MS_C", "MS_D"}))

    def test_grouping_at_locus_2(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_2"]
        readings = {g.reading: g.witnesses for g in groups}
        self.assertEqual(readings["il marche"], frozenset({"MS_A", "MS_D"}))
        self.assertEqual(readings["il court"], frozenset({"MS_B", "MS_C"}))

    def test_grouping_at_locus_3(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_3"]
        readings = {g.reading: g.witnesses for g in groups}
        self.assertEqual(readings["vers la mer"],
                         frozenset({"MS_A", "MS_B", "MS_C"}))
        self.assertEqual(readings["vers la montagne"], frozenset({"MS_D"}))

    def test_all_groups_have_locus(self):
        for g in self.apparatus.groups:
            self.assertIsInstance(g.locus, str)
            self.assertTrue(len(g.locus) > 0)

    def test_all_groups_have_witnesses(self):
        for g in self.apparatus.groups:
            self.assertIsInstance(g.witnesses, frozenset)
            self.assertTrue(len(g.witnesses) > 0)

    def test_default_genealogical(self):
        for g in self.apparatus.groups:
            self.assertTrue(g.genealogical)

    def test_default_reading_type(self):
        for g in self.apparatus.groups:
            self.assertEqual(g.reading_type, "")


class TestFromCSVCodedAlignment(unittest.TestCase):
    """Tests for from_csv with a pre-coded alignment matrix."""

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "coded_alignment.csv")
        )

    def test_witnesses_extracted(self):
        self.assertEqual(
            self.apparatus.witnesses,
            {"MS_A", "MS_B", "MS_C", "MS_D"},
        )

    def test_grouping_at_locus_1(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_1"]
        readings = {g.reading: g.witnesses for g in groups}
        self.assertEqual(readings["a"], frozenset({"MS_A", "MS_B"}))
        self.assertEqual(readings["b"], frozenset({"MS_C", "MS_D"}))

    def test_grouping_at_locus_2(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_2"]
        readings = {g.reading: g.witnesses for g in groups}
        self.assertEqual(readings["a"], frozenset({"MS_A", "MS_D"}))
        self.assertEqual(readings["b"], frozenset({"MS_B", "MS_C"}))


class TestFromCSVWithLacunae(unittest.TestCase):
    """Tests for from_csv with missing readings."""

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "with_lacunae.csv")
        )

    def test_lacuna_excluded_from_group(self):
        """A witness with an empty cell should not appear in any group
        at that locus."""
        groups = [g for g in self.apparatus.groups if g.locus == "locus_1"]
        all_witnesses = set()
        for g in groups:
            all_witnesses |= g.witnesses
        self.assertNotIn("MS_C", all_witnesses)

    def test_lacuna_at_locus_2(self):
        groups = [g for g in self.apparatus.groups if g.locus == "locus_2"]
        all_witnesses = set()
        for g in groups:
            all_witnesses |= g.witnesses
        self.assertNotIn("MS_A", all_witnesses)

    def test_no_lacuna_locus_complete(self):
        """A locus with no lacunae should contain all witnesses."""
        groups = [g for g in self.apparatus.groups if g.locus == "locus_3"]
        all_witnesses = set()
        for g in groups:
            all_witnesses |= g.witnesses
        self.assertEqual(all_witnesses, {"MS_A", "MS_B", "MS_C"})


class TestFromCSVNoVariation(unittest.TestCase):
    """Tests when all witnesses agree at every locus."""

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "no_variation.csv")
        )

    def test_one_group_per_locus(self):
        loci = {g.locus for g in self.apparatus.groups}
        for locus in loci:
            groups = [g for g in self.apparatus.groups if g.locus == locus]
            self.assertEqual(
                len(groups), 1,
                f"Expected 1 group at {locus}, got {len(groups)}",
            )

    def test_all_witnesses_in_group(self):
        for g in self.apparatus.groups:
            self.assertEqual(
                g.witnesses,
                frozenset({"MS_A", "MS_B", "MS_C"}),
            )


class TestGroupManipulation(unittest.TestCase):
    """Test that grouping by group of manuscripts behaves
    as expected."""
    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "coded_alignment.csv")
        )

    def test_group_frequency(self):
        self.assertEqual(
            self.apparatus.group_by_frequency(),
            {frozenset({'MS_A', 'MS_B'}):
             {'frequency': 1, 'readings': 'a'},
             frozenset({'MS_D', 'MS_C'}):
                {'frequency': 1, 'readings': 'b'},
                frozenset({'MS_A', 'MS_D'}):
             {'frequency': 1, 'readings': 'a'},
                frozenset({'MS_C', 'MS_B'}):
                    {'frequency': 1, 'readings': 'b'},
                frozenset({'MS_A', 'MS_C', 'MS_B'}):
                    {'frequency': 1, 'readings': 'a'},
                frozenset({'MS_D'}): {'frequency': 1, 'readings': 'b'}})


if __name__ == "__main__":
    unittest.main()
