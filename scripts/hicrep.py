from typing import Union
from contextlib import suppress
import numpy as np
import scipy.sparse as sp

def upperDiagCsr(m: sp.coo_matrix, nDiags: int):
    """Convert an input sp.coo_matrix into a sp.csr_matrix where each row in the
    the output corresponds to one diagonal of the upper triangle of the input.
    Args:
        m (sp.coo_matrix): input matrix
        nDiags (int): output diagonals with index in the range [1, nDiags)
        as rows of the output matrix
    Returns: `sp.csr_matrix` whose rows are the diagonals of the input
    """
    row = m.col - m.row
    idx = np.where((row > 0) & (row < nDiags))
    idxRowp1 = row[idx]
    # the diagonal index becomes the row index
    idxRow = idxRowp1 - 1
    # offset in the original diagonal becomes the column index
    idxCol = m.col[idx] - idxRowp1
    ans = sp.csr_matrix((m.data[idx], (idxRow, idxCol)),
                        shape=(nDiags - 1, m.shape[1]), dtype=m.dtype)
    ans.eliminate_zeros()
    return ans

def varVstran(n: Union[int, np.ndarray]):
    """
    Calculate the variance of variance-stabilizing transformed
    (or `vstran()` in the original R implementation) data. The `vstran()` turns
    the input data into ranks, whose variance is only a function of the input
    size:
        ```
        var(1/n, 2/n, ..., n/n) = (1 - 1/(n^2))/12
        ```
    or with Bessel's correction:
        ```
        var(1/n, 2/n, ..., n/n, ddof=1) = (1 + 1.0/n)/12
        ```
    See section "Variance stabilized weights" in reference for more detail:
    https://genome.cshlp.org/content/early/2017/10/06/gr.220640.117
    Args:
        n (Union(int, np.ndarray)): size of the input data
    Returns: `Union(int, np.ndarray)` variance of the ranked input data with Bessel's
    correction
    """
    with suppress(ZeroDivisionError), np.errstate(divide='ignore', invalid='ignore'):
        return np.where(n < 2, np.nan, (1 + 1.0 / n) / 12.0)


def sccByDiag(m1: sp.coo_matrix, m2: sp.coo_matrix, nDiags: int):
    """Compute diagonal-wise hicrep SCC score for the two input matrices up to
    nDiags diagonals
    Args:
        m1 (sp.coo_matrix): input contact matrix 1
        m2 (sp.coo_matrix): input contact matrix 2
        nDiags (int): compute SCC scores for diagonals whose index is in the
        range of [1, nDiags)
    Returns: `float` hicrep SCC scores
    """
    # convert each diagonal to one row of a csr_matrix in order to compute
    # diagonal-wise correlation between m1 and m2
    m1D = upperDiagCsr(m1, nDiags)
    m2D = upperDiagCsr(m2, nDiags)
    nSamplesD = (m1D + m2D).getnnz(axis=1)
    rowSumM1D = m1D.sum(axis=1).A1
    rowSumM2D = m2D.sum(axis=1).A1
    # ignore zero-division warnings because the corresponding elements in the
    # output don't contribute to the SCC scores
    with np.errstate(divide='ignore', invalid='ignore'):
        cov = m1D.multiply(m2D).sum(axis=1).A1 - rowSumM1D * rowSumM2D / nSamplesD
        rhoD = cov / np.sqrt(
            (m1D.power(2).sum(axis=1).A1 - np.square(rowSumM1D) / nSamplesD ) *
            (m2D.power(2).sum(axis=1).A1 - np.square(rowSumM2D) / nSamplesD ))
        wsD = nSamplesD * varVstran(nSamplesD)
        # Convert NaN and Inf resulting from div by 0 to zeros.
        # posinf and neginf added to fix behavior seen in 4DN datasets
        # 4DNFIOQLTI9G and DNFIH7MQHOR at 5kb where inf would be reported
        # as an SCC score
        wsNan2Zero = np.nan_to_num(wsD, copy=True, posinf=0.0, neginf=0.0)
        rhoNan2Zero = np.nan_to_num(rhoD, copy=True, posinf=0.0, neginf=0.0)

    return rhoNan2Zero @ wsNan2Zero / wsNan2Zero.sum()

