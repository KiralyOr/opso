import sympy
import time
import pickle
import csv
from sympy import Matrix, Integer
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from collections import defaultdict
from functools import partial

def gen_seq(c, m):
    if m == 0:
        return [[]]
    else:
        arr1 = [[x] + subseq for subseq in gen_seq(c, m - 1) for x in range(1, c)]
        arr2 = [[c] + subseq for subseq in gen_seq(c + 1, m - 1)]
        return arr1 + arr2

def to_pairs(s):
    return list(zip(s[::2], s[1::2]))

def gen_pairs(m):
    return [to_pairs(s) for s in gen_seq(1, 2 * m)]

def dfa_transition_matrix(m, S):
    U = {u for u, _ in S}
    Q = list(U.union({0, float('inf')}))
    Sigma = set(range(1, m + 1))
    state_index = {state: idx for idx, state in enumerate(Q)}
    num_states = len(Q)
    transition_matrix = [[Integer(0)] * num_states for _ in range(num_states)]

    for u in Q:
        for v in Sigma:
            if u == float('inf') or (u, v) in S:
                transition_matrix[state_index[u]][state_index[float('inf')]] += 1
            elif (u, v) not in S and v in U:
                transition_matrix[state_index[u]][state_index[v]] += 1
            else:
                transition_matrix[state_index[u]][state_index[0]] += 1

    return Matrix(transition_matrix)

def sympy_perm(num, k):
    return sympy.factorial(num) // sympy.factorial(num - k)

def calculate_single_cached(args, V, L):
    try:
        S_frozen, k, count = args
        S = list(S_frozen)
        A = dfa_transition_matrix(V, S)
        path_num = (A ** L)[0, -1]
        return count * path_num * sympy_perm(V, k)
    except Exception as e:
        print(f"[ERROR] Worker failed on input {args}: {e}")
        return Integer(0)

def calculate_total_sum_cached(m, V, L, all_pairs):
    group_counts = defaultdict(int)
    k_max_for_group = {}

    for s in all_pairs:
        S_key = frozenset(s)
        group_counts[S_key] += 1
        k = max(v for pair in s for v in pair)
        k_max_for_group[S_key] = max(k_max_for_group.get(S_key, 0), k)

    inputs = [(S_key, k_max_for_group[S_key], group_counts[S_key]) for S_key in group_counts]

    with ProcessPoolExecutor() as executor:
        worker = partial(calculate_single_cached, V=V, L=L)
        results = list(tqdm(executor.map(worker, inputs), total=len(inputs)))

    return sum(results, Integer(0))

def get_e(m, total_sum, precision=50):
    V = sympy.Integer(2)**10
    L = sympy.Integer(2)**21

    factor = sympy.Integer(m * 2)
    denominator = V**(L + factor)
    result = 1 - (total_sum / denominator).evalf(precision)
    e_x = (result * V ** (factor)).evalf(precision)
    return e_x, result

if __name__ == "__main__":
    V = Integer(2**10)
    L = Integer(2**21)
    log_file = "summary_log.csv"

    with open(log_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["m", "pair_count", "unique_S", "runtime_sec", "pkl_file", "result", "e_x"])

        for m in range(1, 5):
            print(f"\nCalculating total_sum for m = {m}")
            all_pairs = gen_pairs(m)
            pair_count = len(all_pairs)

            start_time = time.time()
            total_sum = calculate_total_sum_cached(m, V, L, all_pairs)
            e_x, result = get_e(m, total_sum, precision=50)
            duration = round(time.time() - start_time, 2)
            filename = f"total_sum_m{m}.pkl"
            with open(filename, "wb") as f:
                pickle.dump(total_sum, f)

            unique_s = len(set(frozenset(s) for s in all_pairs))
            print(f"Saved {filename} (pairs: {pair_count}, time: {duration}s)")
            print(f"result (1-p): {result.evalf(10)}")
            print(f"e_x: {e_x.evalf(10)}")

            writer.writerow([m, pair_count, unique_s, duration, filename, result.evalf(50), e_x.evalf(50)])