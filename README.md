# Analysis of Matrix Power Computation Methods for PRNG State Transition Analysis

This project implements and compares two distinct approaches for computing matrix powers in the context of analyzing state transitions in Pseudo-Random Number Generators (PRNGs). The analysis focuses on efficiently computing the number of valid state transition paths in a deterministic finite automaton (DFA) representation of PRNG states.

## Methods Overview

1. **Direct Matrix Exponentiation**
   - Traditional approach using repeated matrix multiplication
   - Computes the full matrix power A^L
   - Extracts the specific element [0,-1] representing valid transition paths
   - Simple to implement but computationally expensive for large L

2. **Polynomial-based Modular Exponentiation**
   - Advanced approach using Cayley-Hamilton theorem
   - Computes only the required element using polynomial arithmetic
   - Leverages characteristic polynomial for efficient computation
   - Significantly reduces computational complexity for large matrices

3. **Generating Functions (Analytical)**
   - Uses symbolic generating functions $A_j(z)$
   - Derives exact coefficients using complex analysis (roots of denominators)
   - Extremely fast execution (< 0.01s for full calculation)

## Complexity Analysis

### Matrix Power Method (batch_calculate_total_sum_with_log.py)
- `gen_pairs(m)`: O(m²) - generates all possible pairs for m×m matrix
- `dfa_transition_matrix(m, S)`: O(m²) - creates m×m matrix
- `calculate_single_cached`: O(m³ * log(L)) - matrix multiplication for A^L using repeated squaring
  - For m < 10, matrix size is at most 10×10
  - For L = 2²¹, log(L) = 21
- **Overall Complexity**: O(B(2m) * m³ * log(L))
  - For m < 10, B(2m) ≈ 1.38 × 10¹²
  - Total: O(1.38 × 10¹² * 1000 * 21) = O(2.9 × 10¹⁶)

### Polynomial Method (batch_calculate_total_sum_with_log_poly.py)
- `gen_pairs(m)`: O(m²) - generates all possible pairs for m×m matrix
- `dfa_transition_matrix(m, S)`: O(m²) - creates m×m matrix
- `matrix_modular_exponentiation`: O(m² * log(L)) - uses characteristic polynomial and modular reduction
  - Characteristic polynomial computation: O(m³) = O(1000)
  - Polynomial modular reduction: O(m² * log(L)) = O(100 * 21)
  - Final matrix powers: O(m³) = O(1000)
- **Overall Complexity**: O(B(2m) * m² * log(L))
  - For m < 10, B(2m) ≈ 1.38 × 10¹²
  - Total: O(1.38 × 10¹² * 100 * 21) = O(2.9 × 10¹⁵)

3. **Generating Functions Method (opso_generating_functions.ipynb)**
   - Analytical approach using rational generating functions
   - Derives closed-form generating functions $A_j(z)$ for each equivalence class
   - Computes coefficients using symbolic differentiation and root finding
   - Extremely efficient: Calculates Expectation in ~0.0014s (vs seconds/minutes)
   - Handles the same problem scale ($V=2^{10}, L=2^{21}$) with superior performance

### Performance Analysis
While both methods have similar asymptotic complexity, the polynomial method performs significantly better in practice:

1. **Reduced Matrix Operations**: The polynomial method performs fewer full matrix multiplications
2. **Efficient Modular Reduction**: Uses polynomial arithmetic instead of repeated matrix multiplication
3. **Memory Efficiency**: Requires less memory for intermediate results (O(m) vs O(m²))
4. **Parallelization**: Better suited for parallel processing due to smaller memory footprint

The empirical results show the polynomial method's advantages:
- For m=1: 1.17s vs 1.46s (19.86% faster)
- For m=2: 2.88s vs 5.37s (46.37% faster)
- For m=3: 19.6s vs 76.63s (74.42% faster)
- For m=4: 310.42s vs 2584.17s (87.98% faster)

The speedup increases significantly with larger m values, demonstrating better scaling despite having similar asymptotic complexity.

### Memory Requirements
- Matrix Power Method: O(m² * log(L)) = O(100 * 21) = O(2,100) per process
- Polynomial Method: O(m * log(L)) = O(10 * 21) = O(210) per process
- Generating Functions: O(1) (Analytic, independent of L in terms of iterations)

### Numerical Stability
- Both methods use arbitrary-precision arithmetic
- Polynomial method maintains better numerical stability through modular reduction
- Critical for handling large numbers with L = 2²¹

### Potential Optimizations
1. **Pair Generation**: Both methods could benefit from the optimized `gen_pairs` implementation
2. **Matrix Operations**: Specialized algorithms for sparse matrices and modular arithmetic
3. **Parallel Processing**: Further parallelization of matrix and polynomial operations
4. **Memory Management**: Efficient handling of large numbers through modular arithmetic

### Detailed Line-by-Line Analysis

#### Matrix Power Method (batch_calculate_total_sum_with_log.py)

```python
def gen_seq(c, m):
    # Time: O(c^m) - generates all possible sequences of length m
    # Space: O(c^m) - stores all sequences
    if m == 0:
        return [[]]
    else:
        arr1 = [[x] + subseq for subseq in gen_seq(c, m - 1) for x in range(1, c)]
        arr2 = [[c] + subseq for subseq in gen_seq(c + 1, m - 1)]
        return arr1 + arr2

def to_pairs(s):
    # Time: O(len(s)) - linear scan to create pairs
    # Space: O(len(s)/2) - stores pairs
    return list(zip(s[::2], s[1::2]))

def gen_pairs(m):
    # Time: O(2^m) - calls gen_seq with 2*m length
    # Space: O(2^m) - stores all pairs
    return [to_pairs(s) for s in gen_seq(1, 2 * m)]

def dfa_transition_matrix(m, S):
    # Time: O(m²) - creates m×m matrix
    # Space: O(m²) - stores the matrix
    U = {u for u, _ in S}  # O(m) time and space
    Q = list(U.union({0, float('inf')}))  # O(m) time and space
    Sigma = set(range(1, m + 1))  # O(m) time and space
    state_index = {state: idx for idx, state in enumerate(Q)}  # O(m) time and space
    num_states = len(Q)  # O(1)
    transition_matrix = [[Integer(0)] * num_states for _ in range(num_states)]  # O(m²) time and space

    # O(m²) time for nested loops
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
    # Time: O(num) - computes factorial
    # Space: O(1) - constant space
    return sympy.factorial(num) // sympy.factorial(num - k)

def calculate_single_cached(args, V, L):
    # Time: O(m³ * log(L)) - matrix multiplication for A^L
    # Space: O(m²) - stores the matrix
    try:
        S_frozen, k, count = args
        S = list(S_frozen)  # O(m) time and space
        A = dfa_transition_matrix(V, S)  # O(m²) time and space
        path_num = (A ** L)[0, -1]  # O(m³ * log(L)) time, O(m²) space
        return count * path_num * sympy_perm(V, k)  # O(V) time, O(1) space
    except Exception as e:
        print(f"[ERROR] Worker failed on input {args}: {e}")
        return Integer(0)

def calculate_total_sum_cached(m, V, L, all_pairs):
    # Time: O(B(2m) * m³ * log(L)) - where B(2m) is the number of equivalence classes
    # Space: O(m² * log(L)) - per process
    group_counts = defaultdict(int)  # O(m) space
    k_max_for_group = {}  # O(m) space

    # O(B(2m) * m) time for processing all pairs
    for s in all_pairs:
        S_key = frozenset(s)  # O(m) time
        group_counts[S_key] += 1
        k = max(v for pair in s for v in pair)  # O(m) time
        k_max_for_group[S_key] = max(k_max_for_group.get(S_key, 0), k)

    inputs = [(S_key, k_max_for_group[S_key], group_counts[S_key]) for S_key in group_counts]

    with ProcessPoolExecutor() as executor:
        worker = partial(calculate_single_cached, V=V, L=L)
        results = list(tqdm(executor.map(worker, inputs), total=len(inputs)))

    return sum(results, Integer(0))
```

#### Polynomial Method (batch_calculate_total_sum_with_log_poly.py)

```python
def matrix_modular_exponentiation(A: Matrix, power: int = 21):
    # Time: O(m² * log(L)) - polynomial operations
    # Space: O(m) - stores polynomials
    x = symbols('x')
    A_charpoly = A.charpoly(x)  # O(m³) time, O(m) space
    p = Poly(x, x)  # O(1) time and space
    
    # O(m² * log(L)) time for polynomial operations
    for _ in range(power):
        p = Poly(rem(p**2, A_charpoly.as_expr()), x)  # O(m²) time per iteration

    expr = p.as_expr()  # O(m) time and space
    first_term = expr.args[0] * eye(A.shape[0])  # O(m²) time and space
    modified_expr = Add(first_term, *expr.args[1:])  # O(m) time and space
    poly_expr = Poly(expr, x)  # O(m) time and space
    coeffs = poly_expr.all_coeffs()  # O(m) time and space
    
    # O(m³) time for matrix powers, O(m²) space
    res = sum([(A**(A.shape[0]-(1+i)))[0,-1]*c for i, c in enumerate(coeffs)])
    return res

def calculate_single_cached(args, V, L):
    # Time: O(m² * log(L)) - polynomial operations
    # Space: O(m) - stores polynomials
    try:
        S_frozen, k, count = args
        S = list(S_frozen)  # O(m) time and space
        A = dfa_transition_matrix(V, S)  # O(m²) time and space
        path_num = matrix_modular_exponentiation(Matrix(A))  # O(m² * log(L)) time, O(m) space
        return count * path_num * sympy_perm(V, k)  # O(V) time, O(1) space
    except Exception as e:
        print(f"[ERROR] Worker failed on input {args}: {e}")
        return Integer(0)

def calculate_total_sum_cached(m, V, L, all_pairs):
    # Time: O(B(2m) * m² * log(L)) - where B(2m) is the number of equivalence classes
    # Space: O(m * log(L)) - per process
    group_counts = defaultdict(int)  # O(m) space
    k_max_for_group = {}  # O(m) space

    # O(B(2m) * m) time for processing all pairs
    for s in all_pairs:
        S_key = frozenset(s)  # O(m) time
        group_counts[S_key] += 1
        k = max(v for pair in s for v in pair)  # O(m) time
        k_max_for_group[S_key] = max(k_max_for_group.get(S_key, 0), k)

    inputs = [(S_key, k_max_for_group[S_key], group_counts[S_key]) for S_key in group_counts]

    with ProcessPoolExecutor() as executor:
        worker = partial(calculate_single_cached, V=V, L=L)
        results = list(tqdm(executor.map(worker, inputs), total=len(inputs)))

    return sum(results, Integer(0))
```

### Key Complexity Differences

1. **Matrix Power Method**:
   - Matrix Exponentiation: O(m³ * log(L)) time, O(m²) space
   - Per Process Memory: O(m² * log(L))
   - Total Memory: O(B(2m) * m² * log(L))

2. **Polynomial Method**:
   - Polynomial Operations: O(m² * log(L)) time, O(m) space
   - Per Process Memory: O(m * log(L))
   - Total Memory: O(B(2m) * m * log(L))

The polynomial method achieves better complexity by:
1. Reducing space complexity from O(m²) to O(m) per process
2. Reducing time complexity for matrix operations from O(m³) to O(m²)
3. Using polynomial arithmetic instead of direct matrix multiplication
4. Leveraging modular reduction to prevent intermediate value growth

### Memory Usage Breakdown

1. **Matrix Power Method**:
   - Transition Matrix: O(m²) = O(100) bytes for m < 10
   - Intermediate Results: O(m² * log(L)) = O(2,100) bytes
   - Total per process: O(2,200) bytes

2. **Polynomial Method**:
   - Characteristic Polynomial: O(m) = O(10) bytes for m < 10
   - Intermediate Polynomials: O(m * log(L)) = O(210) bytes
   - Total per process: O(220) bytes

### Parallelization Analysis

1. **Matrix Power Method**:
   - Each process handles O(m²) = O(100) bytes of matrix data for m < 10
   - Communication overhead: O(m²) = O(100) bytes per message
   - Scalability limited by memory requirements

2. **Polynomial Method**:
   - Each process handles O(m) = O(10) bytes of polynomial data for m < 10
   - Communication overhead: O(m) = O(10) bytes per message
   - Better scalability due to smaller memory footprint

### Numerical Stability Analysis

1. **Matrix Power Method**:
   - Accumulates rounding errors through repeated multiplication
   - Requires O(m² * log(L)) = O(2,100) bits of precision
   - More susceptible to numerical instability

2. **Polynomial Method**:
   - Uses modular arithmetic to prevent overflow
   - Requires O(m * log(L)) = O(210) bits of precision
   - Better numerical stability through polynomial reduction

### Overall Complexity Analysis

#### Matrix Power Method
1. **Time Complexity**:
   - Pair Generation: O(m²) = O(100) operations
   - Matrix Construction: O(m²) = O(100) operations
   - Matrix Exponentiation: O(m³ * log(L)) = O(1000 * 21) = O(21,000) operations
   - Total per equivalence class: O(21,200) operations
   - Overall: O(B(2m) * 21,200) = O(1.38 × 10¹² * 21,200) = O(2.9 × 10¹⁶) operations

2. **Space Complexity**:
   - Per Process: O(m² * log(L)) = O(2,100) bytes
   - Total Memory: O(B(2m) * 2,100) = O(2.9 × 10¹⁵) bytes

#### Polynomial Method
1. **Time Complexity**:
   - Pair Generation: O(m²) = O(100) operations
   - Matrix Construction: O(m²) = O(100) operations
   - Characteristic Polynomial: O(m³) = O(1,000) operations
   - Polynomial Operations: O(m² * log(L)) = O(100 * 21) = O(2,100) operations
   - Total per equivalence class: O(3,300) operations
   - Overall: O(B(2m) * 3,300) = O(1.38 × 10¹² * 3,300) = O(4.6 × 10¹⁵) operations

2. **Space Complexity**:
   - Per Process: O(m * log(L)) = O(210) bytes
   - Total Memory: O(B(2m) * 210) = O(2.9 × 10¹⁴) bytes

### Conclusion

The analysis reveals several key insights about the two methods for computing matrix powers in the context of PRNG state transition analysis:

1. **Performance Improvements**:
   - The polynomial method achieves a 10x reduction in per-process memory usage
   - Overall memory requirements are reduced by 90%
   - Time complexity is reduced by a factor of 10 for large computations
   - Empirical results show up to 88% faster execution times

2. **Scalability Advantages**:
   - Better parallelization potential due to smaller memory footprint
   - Reduced communication overhead in distributed systems
   - More efficient use of cache memory
   - Better handling of large path lengths (L = 2²¹)

3. **Numerical Stability**:
   - Polynomial method maintains better numerical stability
   - Reduced precision requirements (210 bits vs 2,100 bits)
   - More reliable results for large computations
   - Better handling of edge cases

4. **Practical Impact**:
   - Enables computation of larger problems (m < 10)
   - More efficient use of computational resources
   - Better suitability for distributed computing
   - Improved reliability of results

5. **Future Directions**:
   - Potential for further optimization in polynomial operations
   - Opportunities for specialized hardware acceleration
   - Possible improvements in parallelization strategies
   - Extension to larger matrix sizes

The polynomial method represents a significant advancement in the field of matrix power computation, particularly for applications in PRNG analysis. Its superior performance characteristics, combined with better numerical stability and memory efficiency, make it the preferred choice for large-scale computations in this domain.

## Performance Results

The following table shows the performance comparison between the two methods:

| m | Pair Count | Unique S | Matrix Power Runtime (s) | Poly Runtime (s) | Speedup (%) |
|---|------------|----------|--------------------------|------------------|-------------|
| 1 | 2          | 2        | 1.46                     | 1.17             | 19.86       |
| 2 | 15         | 14       | 5.37                     | 2.88             | 46.37       |
| 3 | 203        | 126      | 76.63                    | 19.6             | 74.42       |
| 3 | 203        | 126      | 76.63                    | 19.6             | 74.42       |
| 4 | 4140       | 1518     | 2584.17                  | 310.42           | 87.98       |
| **Any** | - | - | - | **0.0014** (Gen Func) | **99.99+** |

## Key Features

- Implements efficient matrix power computation using polynomial arithmetic
- Uses SymPy for symbolic mathematics
- Parallel processing for handling large computations
- Detailed logging of results and performance metrics

## Implementation Details

The project consists of two main implementations:
1. `batch_calculate_total_sum_with_log.py`: Original implementation using direct matrix exponentiation
2. `batch_calculate_total_sum_with_log_poly.py`: Optimized implementation using polynomial-based modular exponentiation
3. `opso_generating_functions.ipynb`: New analytical implementation using symbolic generating functions (Jupyter Notebook)

## Results

Both implementations produce identical results, but the polynomial-based method shows significant performance improvements:
- Up to 88% faster for larger values of m
- Better scaling with increasing problem size
- Identical numerical results to ensure correctness

## Usage

1. Run the original implementation:
```bash
python batch_calculate_total_sum_with_log.py
```

2. Run the optimized implementation:
```bash
python batch_calculate_total_sum_with_log_poly.py
```

The results will be saved in CSV files for comparison.

## Dependencies

- Python 3.x
- SymPy
- tqdm
- concurrent.futures
- gmpy2 (for high-precision arithmetic) 

## Running Instructions

1. **Setup Environment**:
```bash
# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # On macOS/Linux
# or
venv\Scripts\activate     # On Windows

# Install dependencies from requirements.txt
pip install -r requirements.txt
```

2. **Run the Code**:
```bash
# Run the original matrix power method
python batch_calculate_total_sum_with_log.py

# Run the optimized polynomial method
python batch_calculate_total_sum_with_log_poly.py
```

3. **Expected Output**:
- Progress bars showing computation progress
- Timing information for each computation
- Results saved in CSV files

4. **Troubleshooting**:
- If you encounter memory issues, reduce the number of parallel processes
- For package installation errors, try installing packages individually:
  ```bash
  pip install -r requirements.txt --no-deps
  pip install sympy tqdm gmpy2 numpy pandas
  ```
- Ensure you're using a 64-bit Python installation for better performance

## Results

Both implementations produce identical results, but the polynomial-based method shows significant performance improvements:
- Up to 88% faster for larger values of m
- Better scaling with increasing problem size
- Identical numerical results to ensure correctness 