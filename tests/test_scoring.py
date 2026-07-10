"""Tests that computing the similarities between stemmata
for scoring behaves as expected."""

import unittest
import dendropy
from pyfroger.evaluate_reconstructed_trees import _all_labels, infer_witnesses, _normalize, \
    clade_based_metrics, _build_graph, average_sign_distance, \
    internal_ancestor_identification, evaluate


class TestScoringMethods(unittest.TestCase):
    """Tests that scoring the methods behave as
    expected."""

    def setUp(self) -> None:
        """Setup the unittests by defining required
        examples for evaluating the methods."""
        # Tree that will be used as testcase
        self.newick_1 = "((w7,w8)w6,(w4,w5)w3,w2)w1;"
        self.newick_2 = "((w4,w5)w6,(w7,w8)w3,w2)w1;"
        self.newick_3 = "((w4,w5)w6,(w7,w8)w3,w2)α;"
        self.newick_4 = "((w4,w5),w6)β,(w7,w8)w3,w2)α;"
        self.tree_1 = dendropy.Tree.get(
            data=self.newick_1,
            schema="newick",
            suppress_internal_node_taxa=False)
        self.tree_2 = dendropy.Tree.get(
            data=self.newick_2,
            schema="newick",
            suppress_internal_node_taxa=False)
        self.tree_3 = dendropy.Tree.get(
            data=self.newick_3,
            schema="newick",
            suppress_internal_node_taxa=False)

        return super().setUp()

    def test_all_labels(self):
        """Tests that all the labels are properly returned."""
        expected_labels = set([f"w{ix}" for ix in range(1, 9)])
        self.assertEqual(
            _all_labels(self.tree_1),
            expected_labels
        )

    def test_infer_witnesses(self):
        """Tests that infering witnesses behaves properly,
        in spite of having an hyperarchetype denoted with
        a greek letter."""
        expected_witnesses = set([f"w{ix}" for ix in range(2, 9)])
        self.assertEqual(
            infer_witnesses(self.tree_3),
            expected_witnesses
        )

    # ---- CLADE DETECTION PERFORMANCE

    def test_normalize_all_extant_witnesses(self):
        """Tests that normalization is done properly when
        all witnesses are extant."""
        list_witnesses = set([f"w{ix}" for ix in range(1, 9)])
        _normalize(self.tree_1,
                   list_witnesses)
        self.assertEqual(
            self.tree_1.as_string("newick").strip(),
            "((((w7,w8),w6:0.0),((w4,w5),w3:0.0),w2),w1:0.0);")

    def test_normalize_lost_witnesses(self):
        """Tests that normalization is done properly when
        some witnesses have disappeared."""
        list_witnesses = set([f"w{ix}" for ix in range(2, 9)])
        _normalize(self.tree_3,
                   list_witnesses)
        self.assertEqual(self.tree_3.as_string("newick").strip("\n"),
                         "(((w4,w5),w6:0.0),((w7,w8),w3:0.0),w2);")

    def test_clade_based_metrics(self):
        """Tests that clade based metrics are properly computed by
        comparing tree 1 and tree 2.
        Tree 1 (true) has (non-trivial) bipartitions:
            - [w4, w5]
            - [w3, w4, w5]
            - [w7, w8]
            - [w6, w7, w8]

        Tree 2 (recon) has (non-trivial) bipartitions:
            - [w4, w5]
            - [w6, w4, w5]
            - [w7, w8]
            - [w3, w7, w8]

        Shared bipartitions (2):
            - [w4, w5]
            - [w7, w8]

        Discordant bipartitions (4):
            - [w3, w4, w5]  (T1 only)   -> w3 heads {w4, w5}
            - [w6, w7, w8]  (T1 only)   -> w6 heads {w7, w8}
            - [w6, w4, w5]  (T2 only)   -> w6 heads {w4, w5}
            - [w3, w7, w8]  (T2 only)   -> w3 heads {w7, w8}

        -> RF = 4 // Normalization by 8
        """
        clade_metrics = clade_based_metrics(
            self.newick_1,
            self.newick_2
        )
        self.assertDictEqual(
            {'rf_distance': 4,
             'normalised_rf': 0.5,
             'similarity': 0.5,
             'clade_precision': 0.5,
             'clade_recall': 0.5,
             'clade_f1': 0.5,
             'true_clades_recovered': '2/4',
             'spurious_clades': 2,
             'n_taxa_compared': 8,
             'witnesses_missing_from_recon': [],
             'witnesses_extra_in_recon': [],
             'rooted': False},
            clade_metrics
        )

    # --- ROOS SIGN DISTANCE

    def test_build_graph(self):
        """Tests that a newick string is properly transformed into a graph."""
        graph, wit = _build_graph(self.newick_1,
                                  witnesses=set([f"w{ix}" for ix
                                                 in range(1, 9)]))

        # Check that witness correspondance is properly stored
        self.assertDictEqual(
            {'w1': 0,
             'w6': 1,
             'w7': 2,
             'w8': 3,
             'w3': 4,
             'w4': 5,
             'w5': 6,
             'w2': 7},
            wit
        )

    def test_average_sign_distance(self):
        """Tests that computing the average sign distance
        behaves as expected."""
        computed_average_sign = average_sign_distance(
            self.newick_1,
            self.newick_2,
            witnesses=set([f"w{ix}" for ix
                           in range(1, 9)])
        )
        self.assertEqual(
            computed_average_sign["n_witnesses_scored"],
            len(set([f"w{ix}" for ix
                     in range(1, 9)])),
        )
        self.assertEqual(
            computed_average_sign["skeleton_sign_distance"],
            69.048
        )

    def test_internal_ancestor_identification(self):
        """Tests the evaluation of internal ancestor identification."""
        # Tree that adds wrongly an archetype
        newick_wrong_archetype = "(((w7,w8),w6)β,(w4,w5)w3,w2)w1;"
        self.assertEqual(internal_ancestor_identification(
            self.newick_1,
            newick_wrong_archetype,
            witnesses=set([f"w{ix}" for ix
                           in range(1, 9)])),
                         {'internal_ancestor_recall': 0.6667,
                          'internal_ancestor_precision': 1.0,
                          'internal_ancestor_f1': 0.8,
                          'correctly_identified': '2/3',
                          'missed_ancestors': ['w6'],
                          'spurious_ancestors': []})

    def test_evaluate_all_metrics(self):
        """Tests that computing all metrics behaves as expected."""
        self.assertEqual(
            evaluate(self.newick_1, self.newick_2),
            {'robinson_foulds': {'rf_distance': 4,
                                 'normalised_rf': 0.5,
                                 'similarity': 0.5,
                                 'clade_precision': 0.5,
                                 'clade_recall': 0.5,
                                 'clade_f1': 0.5,
                                 'true_clades_recovered': '2/4',
                                 'spurious_clades': 2,
                                 'n_taxa_compared': 8,
                                 'witnesses_missing_from_recon': [],
                                 'witnesses_extra_in_recon': [],
                                 'rooted': False},
             'roos_sign_distance': {'skeleton_sign_distance': 69.048,
                                    'n_witnesses_scored': 8},
                'internal_ancestors': {'internal_ancestor_recall': 1.0,
                                       'internal_ancestor_precision': 1.0,
                                       'internal_ancestor_f1': 1.0,
                                       'correctly_identified': '3/3',
                                       'missed_ancestors': [],
                                       'spurious_ancestors': []}})


if __name__ == "__main__":
    unittest.main()
