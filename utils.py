from random import randint
from sage.all import *

def gf2n_encode_matrices(F):
    assert F.characteristic() == 2

    n = F.degree()
    mod = F.modulus()

    entry = [[0]*n for _ in range(n)]
    for i in range(n-1):
        entry[i+1][i] = 1
    for i in range(n):
        entry[i][n-1] = mod[i]

    A = matrix(GF(2), entry)
    return [A**i for i in range(n)]

def gf2vec_to_gf2n(vec, encode_mats):
    return sum(b * mat for b, mat in zip(vec, encode_mats))

def gf2n_to_gf2vec(mat):
    return mat.column(0)

def random_nonzero_gf2n_elem(F):
    assert F.characteristic() == 2

    n = F.degree()
    return F.fetch_int(randint(1, 2**n-1))

def random_gf2n_elem(F):
    assert F.characteristic() == 2

    n = F.degree()
    return F.fetch_int(randint(0, 2**n-1))

def random_rain_lin(F):
    """Build matrix from a random linearized polynomial L(X) such that
    1) L^{-1}(X) exists
    2) coefficients of L(X) and L^{-1} are all nonzero
    3) degree of L(X) and L^{-1}(X) are n-1
    """
    assert F.characteristic() == 2
    n = F.degree()

    # construct appropriate linearized polynomial
    while True:
        linpoly_coeffs = [random_nonzero_gf2n_elem(F) for i in range(n)]

        # construct dickson matrix
        entry = [[linpoly_coeffs[j-i]**(2**i) for j in range(n)] for i in range(n)]
        matrix_dickson = matrix(F, entry)

        # check inverse polynomial
        inv_coeffs = matrix_dickson.adjugate().row(0) # skip multiplying 1/det
        if matrix_dickson.is_invertible() and F.zero() not in inv_coeffs:
            break

    # Precompute basis and polynomial evaluations
    basis = [F.fetch_int(2)**i for i in range(n)]
    dual_basis = F.dual_basis()
    poly_evals = [sum([linpoly_coeffs[j] * (basis[i]**(2**j)) for j in range(n)])
                  for i in range(n)]

    entry = [[(dual_basis[i] * poly_evals[j]).trace() for j in range(n)]
             for i in range(n)]

    return matrix(GF(2), entry)

