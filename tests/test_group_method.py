"""Tests that the group method for building the unoriented stemma behaves
as expected."""
import unittest
import os
from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

DATASETS = os.path.join(os.path.dirname(__file__), "test_datasets")


class TestGroupMethod(unittest.TestCase):
    """Tests that the group method behaves as expected."""

    @classmethod
    def setUpClass(cls):
        cls.apparatus = ListApparatus.from_csv(
            os.path.join(DATASETS, "normal_tradition.csv"),
            ref_witness="MS_A"
        )
        cls.apparatus_contaminated = ListApparatus.from_csv(
            os.path.join(DATASETS, "contaminated_tradition.csv"),
            ref_witness="MS_A"
        )
        cls.stemma = UnorientedStemma(cls.apparatus)
        # TODO: perform the unit tests in the case of contaminated things
        cls.stemma_contaminated = UnorientedStemma(cls.apparatus_contaminated)

    def test_compute_covering_relations(self):
        """Tests that iteratively computing the covering relations
        between the different groups behaves as expected."""
        # All different variant groups from the apparatus
        groups = [{"MS_A", "MS_B", "MS_C", "MS_D"},
                  {"MS_B"},
                  {"MS_C", "MS_D"},
                  {"MS_D"},
                  {"MS_C"}]
        # Put all items as frozenset for the MS group to become
        # hashable
        groups = frozenset([frozenset(g) for g in groups])
        expected_edges = [(frozenset({'MS_D', 'MS_C'}), frozenset({'MS_D'})),
                          (frozenset({'MS_D', 'MS_C'}), frozenset({'MS_C'})),
                          (frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}),
                           frozenset({'MS_B'})),
                          (frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}),
                           frozenset({'MS_D', 'MS_C'}))]
        expected_children_of = {frozenset({'MS_D', 'MS_C'}):
                                [frozenset({'MS_C'}), frozenset({'MS_D'})],
                                frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}):
                                    [frozenset({'MS_B'}), frozenset({'MS_D', 'MS_C'})]}
        expected_contamination_overlap = []
        # No contamination in the normal tradition csv
        edges, children_of, contamination_overlap = \
            self.stemma._compute_covering_relations(groups)
        self.assertEqual(set(expected_edges), set(edges))
        self.assertEqual(
            {k: set(v) for k, v in expected_children_of.items()},
            {k: set(v) for k, v in children_of.items()}
        )
        self.assertEqual(
            set(expected_contamination_overlap), set(contamination_overlap)
        )

    def test_compute_edges(self):
        """Tests that the edges are properly computed in the case
        of a normal tradition."""
        # Compute the edges of the stemma
        self.stemma.compute_edges()
        # Expected data
        expected_edges = [(frozenset({'MS_D', 'MS_C'}), frozenset({'MS_D'})),
                          (frozenset({'MS_D', 'MS_C'}), frozenset({'MS_C'})),
                          (frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}),
                           frozenset({'MS_B'})),
                          (frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}),
                           frozenset({'MS_D', 'MS_C'}))]
        expected_children_of = {frozenset({'MS_D', 'MS_C'}):
                                [frozenset({'MS_C'}), frozenset({'MS_D'})],
                                frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}):
                                    [frozenset({'MS_B'}),
                                     frozenset({'MS_D', 'MS_C'})]}
        expected_contamination_overlap = []
        # Check the proper assignment of the different computed values
        self.assertEqual(set(expected_edges),
                         set(self.stemma.edges))
        self.assertEqual(
            {k: set(v) for k, v in expected_children_of.items()},
            {k: set(v) for k, v in self.stemma.children_of.items()}
        )
        self.assertEqual(
            set(expected_contamination_overlap),
            set(self.stemma.contamination_overlap)
        )

    def test_compute_covering_relations_contaminated(self):
        """Tests that iteratively computing the covering relations
        between the different groups behaves as expected in
        the case of a contaminated tradition."""
        pass

    def test_build_tree_structure(self):
        """Tests that building the tree structure behaves as expected."""
        children_of = {frozenset({'MS_D', 'MS_C'}):
                       [frozenset({'MS_C'}), frozenset({'MS_D'})],
                       frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'}):
                       [frozenset({'MS_B'}),
                        frozenset({'MS_D', 'MS_C'})]}
        root_witness = "MS_A"
        root_group = frozenset({'MS_B', 'MS_A', 'MS_D', 'MS_C'})
        # All witnesses
        self.assertEqual(self.stemma._build_tree_structure(
            children_of,
            root_group,
            root_witness
        ),
            {'MS_A':
             {'MS_B': {},
              'α': {'MS_C': {}, 'MS_D': {}}}})

    def test_build_stemma(self):
        """Tests that the stemma is correctly built."""
        self.assertDictEqual(self.stemma.build_tree(
            root_witness="MS_A"
        ),
            {'MS_A':
             {'MS_B': {},
              'α': {'MS_C': {}, 'MS_D': {}}}})

    def test_stemma_to_newick(self):
        """Tests that putting the output stemma into a Newick string
        format behaves as expected."""
        self.assertEqual(self.stemma.to_newick(root_witness="MS_A"),
                         "(MS_B,(MS_C,MS_D)α)MS_A;")

    def test_render_tree(self):
        """Tests that the ASCII-art rendering of a nested-dict tree is correct."""
        tree = {'MS_A': {'MS_B': {}, 'α': {'MS_C': {}, 'MS_D': {}}}}
        expected = (
            "MS_A\n"
            "├── MS_B\n"
            "└── α\n"
            "    ├── MS_C\n"
            "    └── MS_D"
        )
        self.assertEqual(UnorientedStemma._render_tree(tree), expected)

    def test_display_tree(self):
        """Tests that display_tree returns the expected ASCII-art string."""
        expected = (
            "MS_A\n"
            "├── MS_B\n"
            "└── α\n"
            "    ├── MS_C\n"
            "    └── MS_D"
        )
        self.assertEqual(self.stemma.display_tree(root_witness="MS_A"), expected)

    def test_build_tree_structure_two_hypotheticals(self):
        """Tests that two Greek-letter hypothetical ancestors are introduced
        when two sub-groups each have no single identifiable ancestor witness."""
        children_of = {
            frozenset({'MS_A', 'MS_B', 'MS_C', 'MS_D', 'MS_E', 'MS_F'}): [
                frozenset({'MS_B', 'MS_C'}),
                frozenset({'MS_D', 'MS_E', 'MS_F'}),
            ],
            frozenset({'MS_B', 'MS_C'}): [
                frozenset({'MS_B'}),
                frozenset({'MS_C'}),
            ],
            frozenset({'MS_D', 'MS_E', 'MS_F'}): [
                frozenset({'MS_D'}),
                frozenset({'MS_E'}),
                frozenset({'MS_F'}),
            ],
        }
        root_group = frozenset({'MS_A', 'MS_B', 'MS_C', 'MS_D', 'MS_E', 'MS_F'})
        self.assertEqual(
            UnorientedStemma._build_tree_structure(children_of, root_group, 'MS_A'),
            {'MS_A': {
                'α': {'MS_B': {}, 'MS_C': {}},
                'β': {'MS_D': {}, 'MS_E': {}, 'MS_F': {}},
            }}
        )

    def test_newick_static(self):
        """Tests the _to_newick_string static method directly on a tree
        that contains two hypothetical ancestors."""
        tree = {
            'MS_A': {
                'α': {'MS_B': {}, 'MS_C': {}},
                'β': {'MS_D': {}, 'MS_E': {}, 'MS_F': {}},
            }
        }
        self.assertEqual(
            UnorientedStemma._to_newick_string(tree),
            "((MS_B,MS_C)α,(MS_D,MS_E,MS_F)β)MS_A;"
        )


if __name__ == "__main__":
    unittest.main()
