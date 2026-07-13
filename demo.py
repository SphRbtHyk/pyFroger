from io import StringIO
import matplotlib.pyplot as plt
from Bio import Phylo

newick = "(((((A,B)α,(K,L)β)γ,M)δ,Ac)ε,Ab,(Ad,Z)ζ,(Ae,(S,T)η,W)θ,Ba,((Bb,Bd)ι,Be)κ,(C,(Cd,E)λ)μ,Ca,Cb,(Cc,Ce)ν,Cf,(Da,J)ξ,F,G,(H,X)ο,N,P,R,V)O;"
# Parse the tree
tree = Phylo.read(StringIO(newick), "newick")

# Draw it
fig, ax = plt.subplots(figsize=(12, 8))
Phylo.draw(tree, axes=ax, do_show=False)
plt.title("Newick Tree")
plt.show()