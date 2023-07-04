from find_quad_eqs import FindEqs
from sage.all import *
from utils import *

def aim_generate_instance(n, exponents):
    F = GF(2**n, 'z')

    l = len(exponents)-1

    # Generate linear layers
    mats = [random_rain_lin(F) for _ in range(l)]

    while True:
        pt = random_nonzero_gf2n_elem(F)

        state = vector(F.zero())
        for (e, mat) in zip(exponents, mats):
            state += mat * vector(pt**e)

        state = F.fetch_int(ZZ(list(state), base=2))
        if state == 0:
            continue

        ct = state**(exponents[-1]) + pt
        return mats, pt, ct, F

def aim_generate_mq_system(n, exponents, mats, ct, F):
    F_encoding_mats = gf2n_encode_matrices(F)
    l = len(exponents)-1
    nvars = n * (l+1)

    B = BooleanPolynomialRing(nvars, "b")
    b = B.gens()

    # y_i = x^(e_i)
    x = vector(b[:n])
    ys = [vector(b[n*i:n*i+n]) for i in range(1,l+1)]

    X = gf2vec_to_gf2n(x, F_encoding_mats)
    Ys = [gf2vec_to_gf2n(y, F_encoding_mats) for y in ys]

    # v = u^(e_star)
    u = sum(mat * y for mat, y in zip(mats, ys))
    v = x + vector(ct)

    U = gf2vec_to_gf2n(u, F_encoding_mats)
    V = gf2vec_to_gf2n(v, F_encoding_mats)

    sbox_inout_vars = [list(x) + list(y) for y in ys] + [list(u) + list(v)]
    Feqs = FindEqs(n)

    eqs_basic = []
    for e, Y in zip(exponents, Ys):
        EQ = X**(e+1) + X*Y
        eqs_basic.extend(gf2n_to_gf2vec(EQ))
    EQ = U**(exponents[-1]+1) + U*V
    eqs_basic.extend(gf2n_to_gf2vec(EQ))

    eqs_full = []
    for inout, e in zip(sbox_inout_vars, exponents):
        eqs_ = Feqs.find_quad(e)
        for eq_ in eqs_:
            eq = sum(prod(inout[ind] for ind in mono_.iterindex()) for mono_ in eq_)
            eqs_full.append(eq)

    return eqs_basic, eqs_full

def rain_generate_instance(n, num_round):
    F = GF(2**n, 'z')

    # Generate linear layers
    mats = [random_rain_lin(F) for _ in range(num_round-1)]
    rcs = [random_nonzero_gf2n_elem(F) for _ in range(num_round)]

    while True:
        pt = random_gf2n_elem(F)
        key = random_gf2n_elem(F)

        state = pt
        for rc, mat in zip(rcs, mats):
            state = state + key + rc
            if state == 0:
                break

            state = state**(-1)
            state = F.fetch_int(ZZ(list(mat * vector(state)), base=2))

        # Last round
        else:
            state =  pt + key + rcs[-1]
            if state == 0:
                continue

            state = state**(-1)
            ct = state + key
            return mats, rcs, key, pt, ct, F

def rain_generate_mq_system(n, mats, rcs, pt, ct, F):
    F_encoding_mats = gf2n_encode_matrices(F)
    I = F_encoding_mats[0]
    num_round = len(rcs)
    nvars = n * num_round

    B = BooleanPolynomialRing(nvars, "b")
    b = B.gens()

    k = vector(b[-n:])
    xs = [k + vector(rcs[i]) for i in range(num_round)]
    ys = [vector(b[n*i:n*i+n]) for i in range(num_round-1)] + [k + vector(ct)]

    xs[0] += vector(pt)
    for r in range(num_round-1):
        xs[r+1] += mats[r] * ys[r]

    Xs = [gf2vec_to_gf2n(x, F_encoding_mats) for x in xs]
    Ys = [gf2vec_to_gf2n(y, F_encoding_mats) for y in ys]

    eqs_basic = []
    eqs_full = []
    for X, Y in zip(Xs, Ys):
        EQ = X*Y + I
        eqs_basic.extend(gf2n_to_gf2vec(EQ))
        eqs_full.extend(gf2n_to_gf2vec(EQ))
        for Z in (X, Y, X**3, Y**3):
            EQ_ = EQ * Z
            eqs_full.extend(gf2n_to_gf2vec(EQ_))

    return eqs_basic, eqs_full

def em_generate_instance(n, e):
    F = GF(2**n, 'z')

    # Generate linear layers
    mats = [random_rain_lin(F) for _ in range(2)]

    while True:
        pt = random_gf2n_elem(F)
        key = random_gf2n_elem(F)

        # Perm = Lin * Exp * Lin
        state = mats[0] * vector(key + pt)
        state = F.fetch_int(ZZ(list(state), base=2))
        if state == 0:
            continue
        state = state**e
        state = mats[1] * vector(state)

        state = F.fetch_int(ZZ(list(state), base=2))
        ct = state + key

        return mats, key, pt, ct, F

def em_generate_mq_system(n, e, mats, pt, ct, F):
    F_encoding_mats = gf2n_encode_matrices(F)
    I = F_encoding_mats[0]
    nvars = n

    B = BooleanPolynomialRing(nvars, "b")
    b = B.gens()

    k = vector(b)
    x = mats[0] * (k + vector(pt))
    y = mats[1]**(-1) * (k + vector(ct))

    X = gf2vec_to_gf2n(x, F_encoding_mats)
    Y = gf2vec_to_gf2n(y, F_encoding_mats)

    inout = list(x) + list(y)

    EQ = X**(e+1) + X*Y
    eqs_basic = list(gf2n_to_gf2vec(EQ))

    eqs_full = []
    if e == -1:
        eqs_full.extend(eqs_basic)
        for Z in (X, Y, X**3, Y**3):
            EQ_ = EQ * Z
            eqs_full.extend(gf2n_to_gf2vec(EQ_))

    else:
        Feqs = FindEqs(n)
        eqs_ = Feqs.find_quad(e)
        for eq_ in eqs_:
            eq = sum(prod(inout[ind] for ind in mono_.iterindex()) for mono_ in eq_)
            eqs_full.append(eq)

    return eqs_basic, eqs_full
