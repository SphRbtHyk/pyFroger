"""Run parsimony and distnce based methods using BioPoython."""
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

from Bio import AlignIO

aln = AlignIO.read('../../datasets/artificial_scribes/notre_besoin.nexus', 'nexus')

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')
dm = calculator.get_distance(aln)

print(dm)

tree = constructor.build_tree(aln)

print(tree.format("newick"))

# print(tree)

# scorer = ParsimonyScorer()
# searcher = NNITreeSearcher(scorer)
# constructor = ParsimonyTreeConstructor(searcher)
# pars_tree = constructor.build_tree(aln)

# Phylo.draw_ascii(pars_tree.as_phyloxml())