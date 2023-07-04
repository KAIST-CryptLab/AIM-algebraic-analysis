from algebraic_analysis import *
from primitives import *
import sys
import os

aim_params = {
        2: # l = 2
        ((7, (3,5,5)),
        (11, (3,9,5)),
        (12, (3,9,5)),
        (13, (3,11,5)),
        (14, (3,11,5)),
        (15, (3,11,5)),
        (16, (3,11,5)),
        (17, (3,11,5)),
        (18, (3,11,5)),
        (19, (3,9,5)),
        (20, (3,9,5))),
        3: # l = 3
        ((11, (3,9,5,5)),
        (12, (3,9,7,5)),
        (13, (3,11,9,5)),
        (14, (3,11,9,5)),
        (15, (3,11,7,5)),
        (16, (3,11,7,5)),
        (17, (3,11,7,5)),
        (18, (3,11,7,5)),
        (19, (3,9,7,5)),
        (20, (3,9,7,5)))}

def em_params(nu):
    assert nu in (2,3,5)

    # NGG
    if nu == 2:
        for n in range(8, 33, 2):
            m = n//2
            e = 2**(m+1) + 2**(m-1) - 1
            yield (n,e)

    # Mersenne
    if nu == 3:
        for n in range(5, 37):
            for k in range(3,2**n,2):
                Feqs = FindEqs(n)
                eqs_ = Feqs.find_quad(2**k-1)
                if len(eqs_) == 3*n:
                    yield (n, 2**k-1)
                    break

    # Inverse
    if nu == 5:
        for n in range(5, 37):
            yield (n,-1)

    ''' # Gold
    if nu == 1:
        for n in range(9, 33, 2):
            s = (n-1)//2
            if is_even(s):
                e = 2**s + 2**(s//2) - 1
            else:
                e = 2**s + 2**((3*s-1)//2) - 1

            yield (n,e)
    '''

def analysis_XL(n, eqs, primitive, param, eqsgen_type, attack_type):
    dirname = f"Data/XL/{primitive}{param}/"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    path = dirname + f"{eqsgen_type}_{n}"

    with open(path, 'w') as f:
        txt = "n, deg: #eq, #mono, rank"
        f.write(txt + '\n')
        print(txt)

        for deg, rank, neqs, nmono in run_xl(eqs):
            txt = f"{n}, {deg}: {neqs}, {nmono}, {rank}"
            f.write(txt + '\n')
            f.flush()
            print(txt)

            if nmono - rank < n:
                break

def analysis_GB(n, eqs, primitive, param, eqsgen_type, attack_type):
    dirname = f"Data/GB/"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    magma_fname = f"{primitive}{param}_{eqsgen_type}_{n}"
    path = dirname + f"{primitive}{param}_{eqsgen_type}"
    if not os.path.exists(path):
        with open(path, 'w') as f:
            f.write("n: deg, #step")
            print("n: deg, #step")

    with open(path, 'a') as f:
        deg, step = run_gb_with_magma(eqs, magma_fname)
        f.write(f"{n}: {deg}, {step}")
        print(f"{n}: {deg}, {step}")

def analysis_benchGB(n, eqs, primitive, param, eqsgen_type, attack_type):
    dirname = f"Data/benchGB/"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    magma_fname = f"{primitive}{param}_{eqsgen_type}_{n}"
    path = dirname + f"{primitive}{param}_{eqsgen_type}"
    if not os.path.exists(path):
        with open(path, 'w') as f:
            f.write("n: deg, #step, time")
            print("n: deg, #step, time")

    with open(path, 'a') as f:
        deg, step, time = bench_gb_with_magma(eqs, magma_fname)
        f.write(f"{n}: {deg}, {step}, {time:.2f}")
        print(f"{n}: {deg}, {step}, {time:.2f}")

def gen_eqs_aim(ell, eqsgen_type):
    for n, exps_ in aim_params[ell]:
        exps = tuple(2**e - 1 for e in exps_)
        mats, pt, ct, F = aim_generate_instance(n, exps)
        eqs_basic, eqs_full = aim_generate_mq_system(n, exps, mats, ct, F)

        if eqsgen_type == "basic":
            eqs = eqs_basic
        else:
            eqs = eqs_full

        yield n, eqs

def gen_eqs_rain(rnds, eqsgen_type):
    for n in range(5, 31):
        mats, rcs, key, pt, ct, F = rain_generate_instance(n, rnds)
        eqs_basic, eqs_full = rain_generate_mq_system(n, mats, rcs, pt, ct, F)

        if eqsgen_type == "basic":
            eqs = eqs_basic
        else:
            eqs = eqs_full

        yield n, eqs

def gen_eqs_em(nu, eqsgen_type):
    for n, exp in em_params(nu):
        mats, key, pt, ct, F = em_generate_instance(n, exp)
        eqs_basic, eqs_full = em_generate_mq_system(n, exp, mats, pt, ct, F)

        if eqsgen_type == "basic":
            eqs = eqs_basic
        else:
            eqs = eqs_full

        yield n, eqs

def main():
    usage = '''
Usage: python3 main.py [primitive] [param] [eqsgen_type] [attack_type]

primitive: aim / rain / em
param: ell(=2or3) of aim / rnds(=3or4) of rain / nu(=2or3or5) of em
eqsgen_type: basic / full
attack_type: XL / GB / benchGB
'''
    if len(sys.argv) < 5:
        print(usage)
        exit()

    primitive = sys.argv[1]
    param = sys.argv[2]
    eqsgen_type = sys.argv[3]
    attack_type = sys.argv[4]

    primitive_params = {
            "aim": ("2", "3"),
            "rain": ("3", "4"),
            "em": ("2", "3", "5")
    }

    if (primitive not in ("aim", "rain", "em")
        or param not in primitive_params[primitive]
        or eqsgen_type not in ("basic", "full")
        or attack_type not in ("XL", "GB", "benchGB")):
        print(usage)
        exit()

    gen_eqs_map = {
            "aim": gen_eqs_aim,
            "rain": gen_eqs_rain,
            "em": gen_eqs_em
    }
    analyzer_map = {
            "XL": analysis_XL,
            "GB": analysis_GB,
            "benchGB": analysis_benchGB
    }

    analyzer = analyzer_map[attack_type]
    gen_eqs = gen_eqs_map[primitive]

    for n, eqs in gen_eqs(int(param), eqsgen_type):
        analyzer(n, eqs, primitive, param, eqsgen_type, attack_type)

main()
