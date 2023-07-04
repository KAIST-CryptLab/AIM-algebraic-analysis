from itertools import combinations, compress
from sage.all import *
import subprocess
import os
import re

def run_xl(equations):
    """Perform eXtended Linearization on boolean quadratic system

    For each iteration, the system is updated as following
        1) extend the equations by multiplying degree-1 monomials on them
           (skip for the initial iteration)
        2) construct and echlonize the coefficient matrix
        3) remove useless equations
        4) rebuild the system

    Args:
        equations: a list of multivariate boolean quadratic polynomials

    Yields:
        deg: degree of the system
        rank: rank of the coefficient matrix
        num_eqs: number of equations
        num_mono: number of appearing monomials (except '1')
    """
    eqs = equations[:]
    ring = equations[0].ring()
    variables = ring.gens()

    # prepare monomials
    monomials = [ring.one()]
    monomials.extend(variables)
    monomials.extend(map(prod, combinations(variables,2)))
    mono_to_ind = {mono: ind for ind, mono in enumerate(monomials)}
    mono_exist = [0] * len(monomials)

    # construct coefficient matrix
    coeff_mat = matrix(GF(2), len(eqs), len(monomials), sparse=False)
    for i, eq in enumerate(eqs):
        for mono in eq:
            ind = mono_to_ind[mono]
            coeff_mat[i,ind] = 1
            mono_exist[ind] = 1

    # echelonize the coefficient matrix and remove useless equations
    coeff_mat.echelonize()
    pivots = coeff_mat.pivots()
    rank = len(pivots)

    yield 2, rank, len(eqs), sum(mono_exist)-mono_exist[0]

    for deg in range(3, len(variables)+1):
        # Rebuild the system
        eqs = [sum(compress(monomials, coeff_mat[row])) for row in range(rank)]

        # extend the equations
        eqs_ = eqs[:]
        eqs.extend(map(mul, cartesian_product_iterator((variables, eqs_))))

        # prepare monomials
        nmono = len(monomials)
        for ind, mono in enumerate(map(prod, combinations(variables, deg))):
            monomials.append(mono)
            mono_to_ind[mono] = nmono+ind
        mono_exist = [0] * len(monomials)

        # construct coefficient matrix
        coeff_mat = matrix(GF(2), len(eqs), len(monomials), sparse=False)
        for i, eq in enumerate(eqs):
            for mono in eq:
                ind = mono_to_ind[mono]
                coeff_mat[i,ind] = 1
                mono_exist[ind] = 1

        # echelonize the coefficient matrix and remove useless equations
        coeff_mat.echelonize()
        pivots = coeff_mat.pivots()
        rank = len(pivots)

        yield deg, rank, len(eqs), sum(mono_exist)-mono_exist[0]

def run_gb_with_magma(equations, fname):
    """Compute groebner basis (gb) of boolean quadratic system using Magma

    Steps are following:
        1) generate magma script and write to './MagmaScripts/[fname]'
        2) run magma with generated script

    Args:
        equations: a list of multivariate boolean quadratic polynomials
        fname: file name to save magma script

    Returns:
        deg: solving degree (maximal degree reached during gb computation)
        steps: gb computation steps
    """
    ring = equations[0].ring()
    variables = ring.gens()
    nvars = len(variables)

    eq_str = ',\n'.join(str(eq) for eq in equations)
    for i, var in enumerate(variables):
        eq_str = re.sub(r'\b%s\b' % str(var), f'b[{i+1}]', eq_str)

    script = f'''
F<[b]> := BooleanPolynomialRing({nvars}, "grevlex");
B := [{eq_str}];
I := Ideal(B);

gb, degs := GroebnerBasis(I);
printf "%o, %o", Max(degs), #degs;
quit;
'''

    path = f'MagmaScripts/{fname}.m'
    out = _run_magma(script, path)

    data = out.split('\n')[4]
    deg, steps = data.split(', ')

    return int(deg), int(steps)

def bench_gb_with_magma(equations, fname, numWarm=1, numRun=3, numOutIter=4):
    """Benchmark groebner basis (gb) computation with magma

    Steps are following:
        1) generate magma script and write to './MagmaScripts/bench/[fname]'
        2) run magma with generated script

    Args:
        equations: a list of multivariate boolean quadratic polynomials
        fname: file name to save magma script
        numWarm: number of warming iterations
        numRun: number of benchmarking iterations

    Returns:
        deg: solving degree (maximal degree reached during gb computation)
        steps: gb computation steps
        cpi: average clock cycles
    """
    ring = equations[0].ring()
    variables = ring.gens()
    nvars = len(variables)

    eq_str = ',\n'.join(str(eq) for eq in equations)
    for i, var in enumerate(variables):
        eq_str = re.sub(r'\b%s\b' % str(var), f'b[{i+1}]', eq_str)

    script = f'''
F<[b]> := BooleanPolynomialRing({nvars}, "grevlex");
B := [{eq_str}];

Is := [];
Js := [];
for i in [1..{numWarm}] do
    Append(~Is, Ideal(B));
end for;
for i in [1..{numRun}] do
    Append(~Js, Ideal(B));
end for;

for i in [1..{numWarm}] do
    gb, degs := GroebnerBasis(Is[i]);
end for;

cycleStart := ClockCycles();
for i in [1..{numRun}] do
    gb, degs := GroebnerBasis(Js[i]);
end for;
cycleEnd := ClockCycles();
printf "%o, %o, %.2o", Max(degs), #degs, Real((cycleEnd - cycleStart)/{3});
quit;
'''
    path = f'MagmaScripts/bench/{fname}.m'

    times = [0] * numOutIter
    for i in range(numOutIter):
        out = _run_magma(script, path)

        data = out.split('\n')[4]
        deg, steps, time = data.split(', ')
        times[i] = float(time)

    cut = numOutIter // 4
    times.sort()
    time = times[cut:-cut] # remove extreme data

    avg_time = sum(times) / len(times)
    return int(deg), int(steps), avg_time

def _run_magma(script, path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(path, 'w') as f:
        f.write(script)

    proc = subprocess.run(['magma', path], capture_output=True, text=True)
    return proc.stdout

