import sympy
import time
import pickle
import csv
from sympy import Matrix, Integer
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

def gen_seq(c, n):
    if n == 0:
        return [[]]
    else:
        arr1 = [[x] + subseq for subseq in gen_seq(c, n - 1) for x in range(1, c)]
        arr2 = [[c] + subseq for subseq in gen_seq(c + 1, n - 1)]
        return arr1 + arr2

def to_pairs(s):
    return list(zip(s[::2], s[1::2]))

def gen_pairs(n):
    return [to_pairs(s) for s in gen_seq(1, 2 * n)]

def dfa_transition_matrix(n, S):
    U = {u for u, _ in S}
    Q = list(U.union({0, float('inf')}))
    Sigma = set(range(1, n + 1))
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

def sympy_perm(n, k):
    return sympy.factorial(n) // sympy.factorial(n - k)

def calculate_single(args):
    try:
        s, V, L = args
        s = list(s)
        A = dfa_transition_matrix(V, s)
        path_num = (A ** L)[0, -1]
        k = max(v for pair in s for v in pair)
        return path_num * sympy_perm(V, k)
    except Exception as e:
        print(f"[ERROR] Worker failed on input {args}: {e}")
        return Integer(0)

def calculate_total_sum(n, V, L, all_pairs):
    inputs = [(s, V, L) for s in all_pairs]

    with ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(calculate_single, inputs), total=len(inputs)))

    return sum(results, Integer(0))

if __name__ == "__main__":
    V = Integer(2**10)
    L = Integer(2**21)
    log_file = "summary_log.csv"

    with open(log_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["n", "pair_count", "runtime_sec", "pkl_file"])

        for n in range(1,5):
            print(f"\nCalculating total_sum for n = {n}")
            all_pairs = gen_pairs(n)
            pair_count = len(all_pairs)

            start_time = time.time()
            total_sum = calculate_total_sum(n, V, L, all_pairs)
            duration = round(time.time() - start_time, 2)

            filename = f"total_sum_n{n}.pkl"
            with open(filename, "wb") as f:
                pickle.dump(total_sum, f)

            print(f"Saved {filename} (pairs: {pair_count}, time: {duration}s)")
            writer.writerow([n, pair_count, duration, filename])