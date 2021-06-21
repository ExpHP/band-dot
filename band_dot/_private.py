from collections import defaultdict
from functools import reduce
import numpy as np

# returns [(aset, bset)] where aset and bset are a set of indices from
#    each file that have nonzero dot products
def subspaces(sparse_dots):
    sparse_sets = {k:set(v) for (k,v) in sparse_dots.items()}
    unique_b_sets = group_into_non_overlapping_sets(sparse_sets.values())

    def pairs():
        for bbs in unique_b_sets:
            aas = [a for a in sparse_sets if sparse_sets[a] & bbs]
            yield (sorted(aas), sorted(bbs))

    return pairs()

def group_into_non_overlapping_sets(sets):
    sets = [frozenset(s) for s in sets]

    @thru(set)
    def step(sets):
        for s in sets:
            yield reduce(frozenset.union, [x for x in sets if x & s], frozenset())

    return frozenset(fixpoint(sets, step))

# returns [abset] where abset is a set of contiguous indices representing
#  a block diagonal block in the dot product matrix
#
# Compared to 'subspaces', this is less informative and produces larger subspaces.  It exists as a
# crutch for code that is not yet equipped to deal with the added complexity of non-contiguous subspaces.
def block_slices(sparse_abs):
    pairs = subspaces(sparse_abs)
    sets = [  frozenset(range(min(aas), max(aas)+1))
            | frozenset(range(min(bbs), max(bbs)+1))
            for (aas, bbs) in pairs]
    sets = group_into_non_overlapping_sets(sets)
    sets = sorted(sets, key=min)
    return [slice(min(x), max(x)+1) for x in sets]

def thru(outer, *a1, **k1):
    return lambda f: lambda *a2, **k2: outer(f(*a2, **k2), *a1, **k1)

def group_keys_by_matching_value(d):
    inv = defaultdict(set)
    for (k, v) in d.items():
        inv[v].add(k)
    return [(ks, v) for v in inv]

def find_permutation(sparse_dots):
    import scipy.sparse as sparse

    n = len(sparse_dots)
    sparse_abs = {s: {t: abs(v) for (t, v) in vs.items()} for (s, vs) in sparse_dots.items()}

    def coo_matrix_from_triples(triples):
        i, j, x = zip(*ijvpairs())
        return sparse.coo_matrix((list(x), (list(i), list(j))))

    def ijvpairs():
        for i in range(n):
            for j in range(n):
                v = sparse_abs[i].get(j, 0)
                if v:
                    yield (i, j, v)
    mat = coo_matrix_from_triples(ijvpairs()).tocsc()

    perm = np.arange(len(sparse_abs))

    # permute to optimize an objective function within each block subspace
    for slc in block_slices(sparse_abs):
        block = mat[slc, slc]
        # surprise surprise scipy we actually want elementwise multiplication
        block = np.array(block.todense().tolist())

        size = block.shape[0]
        weights = (np.arange(size)[:, None] - np.arange(size)[None, :])**2

        objective = lambda block: (block * weights).sum()

        # a function that tries swapping each pair of indices once, greedily accepting any
        # changes that reduce the objective function
        def step(data):
            block, perm = data
            block, perm = block.arr, perm.arr
            bestBlock, bestPerm, bestValue = block, perm, objective(block)

            for i in range(size):
                for j in range(size):
                    if i == j: continue
                    block, perm = bestBlock.copy(), bestPerm.copy()
                    block[:,[i,j]] = block[:,[j,i]]
                    perm[[i,j]] = perm[[j,i]]

                    if (bestValue is None) or objective(block) < bestValue:
                        bestBlock = block.copy()
                        bestValue = objective(block)
                        bestPerm  = perm.copy()

            return CompareAsBool(bestBlock), CompareAsBool(bestPerm)

        # now we just have to repeatedly apply that function until it stops producing any changes
        _, out = fixpoint((CompareAsBool(block), CompareAsBool(perm[slc])), step)
        perm[slc] = out.arr

    return perm

# ffs numpy
class CompareAsBool:
    def __init__(self, arr): self.arr = arr
    def __eq__(self, other): return (self.arr == other.arr).all()

def fixpoint(x, f):
    while True:
        next = f(x)
        if next == x:
            return x
        x = next
