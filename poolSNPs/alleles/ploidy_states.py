import numpy as np
from cyvcf2 import VCF, Variant

"""
Drafts for implemeting GL parsing in cyvcf2.

1) Recursive tree for enumerating all the genotype states (non-ordered) 
given the ploidy and the ref/alt alleles at a marker. test with diploid SNP
"""


def recurse_genotypes(node: set, alleles: list, depth: int, ploidy: int) -> set:
    childs = []
    if depth < ploidy:
        for parent in node:
            for a in alleles:
                parent = sorted(parent)
                child = sorted(parent + [a])
                childs.append(tuple(child))
        return recurse_genotypes(childs, alleles, depth + 1, ploidy)
    else:
        return set(node)



ref_alt = ['R', 'A']
states = recurse_genotypes([''], ref_alt, 0, 2)
print(states)
