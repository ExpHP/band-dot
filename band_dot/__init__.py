
# Token required to use the private constructor.
_PAIR_OF_EIGENSYSTEMS__CAN_BREAK_MY_CODE = object()

from collections import defaultdict
import numpy as np
import warnings

from . import _private

def _get_angy_about_private_constructor(klass):
    return TypeError(
        'The constructor of {} is private. Please construct using one of the classmethods instead.'
        .format(klass.__name__)
    )

class PairOfEigensystems:
    '''
    Public python API for band_dot.
    
    This is a type you can construct from two sets of eigenvectors, with methods that provide information
    targeted at the problem of finding matching eigenvectors between the two.
    '''
    def __init__(self, sparse_dots, _private_constructor=None):
        """ Do not use. Construct using a class method instead. """
        if _private_constructor is not _PAIR_OF_EIGENSYSTEMS__CAN_BREAK_MY_CODE:
            raise _get_angy_about_private_constructor(type(self))
        self._sparse_dots = sparse_dots

    @classmethod
    def from_eigenvectors(cls, vectors1, vectors2, threshold=0):
        '''
        Construct from two sets of eigenvectors.
        
        ``vectors1`` and ``vectors2`` are 2D arrays (real or complex) where each row is a normalized
        eigenvector. (this is the transpose of the conventional definition of the eigenvector matrix!).
        Any dot products whose square norm is less than ``threshold`` will be considered zero.
        '''
        vectors1 = np.array(vectors1, copy=False)
        vectors2 = np.array(vectors2, copy=False)

        if vectors1.ndim != 2 or vectors2.ndim != 2:
            raise TypeError('Expected two 2D arrays')

        return cls.from_dense_dots(vectors1.conj() @ vectors2.T, threshold=threshold)

    @classmethod
    def from_dense_dots(cls, dense_dots, threshold=0):
        '''
        Construct from a 2D array of inner products.

        ``dense_dots`` is a 2D array (real or complex) where ``dense_dots[i, j]`` is the inner product
        between the ``i``th vector in the first set, and the ``j``th vector in the second set.

        Any elements whose square norm is less than ``threshold`` will be truncated to zero.
        '''
        dense_dots = np.array(dense_dots, copy=True)
        if dense_dots.ndim != 2:
            raise TypeError('Expected a 2D arrays')

        sparse_dots = defaultdict(dict)
        dense_dots[np.abs(dense_dots)**2 < threshold] = 0  # note: preserves NaN
        for (ai, bi), dot in np.ndenumerate(dense_dots):
            if dot != 0:
                sparse_dots[ai][bi] = abs(dot)
        return cls(sparse_dots, _private_constructor=_PAIR_OF_EIGENSYSTEMS__CAN_BREAK_MY_CODE)

    def nonzero_terms(self):
        '''
        Iterate over nonzero dot products between the vectors.

        Produces tuples ``((i, j), dot)``, where ``i`` is a vectors1 index, ``j`` is a vectors2 index,
        and ``dot`` is a complex value with magnitude >= 0.
        '''
        for (i, inner) in self._sparse_dots.items():
            for (j, dot) in inner.items():
                yield ((i, j), dot)

    def get_matrix_element(self, i, j):
        '''
        Get the inner product between vector ``i`` in ``vectors1`` and vector ``j`` in ``vectors2``.
        (which may have been truncated to zero based on ``threshold``)
        '''
        return self._sparse_dots[i].get(j, 0)

    def subspaces(self):
        '''
        Partitions each set of eigenvectors into subspaces that overlap between the two sets,
        allowing one to identify where mixing occurs.

        Iterates over tuples ``(idxs1, idxs2)``, where ``idxs1`` is a sorted list of indices into
        ``vectors1``, and ``idxs2`` is a sorted list of indices into ``vectors2``.

        The order of the pairs returned is unspecified.

        ``idxs1`` is the union of all indices into ``vectors1`` that have at least one nonzero dot product
        with any of the indices in ``idxs2``, and vice versa. Each index into ``vectors1`` will only appear
        in at most a single ``idxs1``, and likewise for ``vectors2/idxs2``.  ``idxs1`` and ``idxs2`` will
        usually tend to be the same size in any given pair (due to orthonormality conditions), and in general
        will not contain contiguous ranges of integers.

        This method only makes sense if you used a ``threshold > 0``.  (otherwise, it will trivially return
        a single subspace of all eigenvectors in each set)
        '''
        yield from _private.subspaces(self._sparse_dots)

    def permutation(self):
        ''' deprecated '''
        warnings.warn("deprecated; use permutation_2_to_1", DeprecationWarning)
        return self.permutation_2_to_1()

    def inverse_permutation(self):
        ''' deprecated '''
        warnings.warn("deprecated; use permutation_1_to_2", DeprecationWarning)
        return self.permutation_1_to_2()

    def permutation_2_to_1(self):
        '''
        Compute an integer 1D array that permutes vectors2 to resemble vectors1.

        That is, ``vectors2[perm]`` will have its rows ordered to match the rows of ``vector1`` as
        closely as possible.
        '''
        return _private.find_permutation(self._sparse_dots)

    def permutation_1_to_2(self):
        '''
        Compute an integer 1D array that permutes vectors2 to resemble vectors1.

        That is, ``vectors2[perm]`` will have its rows ordered to match the rows of ``vector1`` as
        closely as possible.
        '''
        return _private.find_permutation(self._sparse_dots)
