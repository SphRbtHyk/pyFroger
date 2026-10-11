# pyFroger

An implementation of Froger's methods for stemmatological reconstruction of manuscript transmission.

PyFroger rebuilds the stemma of a manuscript tradition from a collation, following the
*méthode des groupements* of Dom Jean Froger (*La critique des textes et son automatisation*,
1968): witnesses sharing the same copying errors are grouped, the groups are ordered by
inclusion, and the genealogy is read off that order. Overlapping groups betray a
contamination and are kept as reticulations.

## Documentation

- **[The method of groupings — user guide](docs/methode_des_groupements.html)**
  (in French): the principles of the method, how to prepare a collation, and a step-by-step
  example of tree reconstruction with `list_apparatus.py` and `group_method.py`.
- **[Benchmark report](paper/benchmarking/analysis/results/metrics_report.html)**:
  Robinson–Folds distance and reconstructed-skeleton percentage of pyFroger against
  iqtree, fastme, mrbayes and parsimony, across copy-error rates.

## Quick start

```python
from pyfroger.list_apparatus import ListApparatus
from pyfroger.group_method import UnorientedStemma

apparatus = ListApparatus.from_csv("collation.csv", ref_witness="ref")
stemma = UnorientedStemma(apparatus)
stemma.build_tree()
stemma.display_tree()          # rendered in the console
stemma.to_dot()                # Graphviz output
```

The collation is a CSV file: one row per variant locus, one column per witness, plus the
reference witness (`ref` by default). See the user guide for the full format and options.

## Repository layout

```
src/pyfroger/
├── list_apparatus.py    # collation → apparatus "en liste" (variant/error groups)
├── group_method.py      # Froger's group method → unoriented stemma
├── simulator.py         # synthetic tradition generator (benchmarking)
├── evaluate_reconstructed_trees.py
└── evaluate_contamination_detection.py
paper/benchmarking/      # benchmark data and analysis
docs/                    # documentation
tests/                   # test suite
```
