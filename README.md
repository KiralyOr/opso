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

## Complexity Analysis

### Matrix Power Method (`batch_calculate_total_sum_with_log.py`)
- `gen_pairs(n)`: O(2ⁿ) - generates all possible pairs
- `dfa_transition_matrix(n, S)`: O(n²) - creates n×n matrix
- `calculate_single_cached`: O(n³) - matrix multiplication for A^L
- Overall Complexity: O(2ⁿ * n³) - dominated by matrix power computation

### Polynomial Method (`batch_calculate_total_sum_with_log_poly.py`)
- `gen_pairs(n)`: O(2ⁿ) - same as matrix power method
- `dfa_transition_matrix(n, S)`: O(n²) - same as matrix power method
- `matrix_modular_exponentiation`: O(n²) - polynomial operations instead of matrix multiplication
- Overall Complexity: O(2ⁿ * n²) - improved by replacing O(n³) matrix operations with O(n²) polynomial operations

The polynomial method achieves better asymptotic complexity by:
- Replacing matrix multiplication (O(n³)) with polynomial arithmetic (O(n²))
- Maintaining the same exponential growth in pair generation (O(2ⁿ))
- Reducing the overall complexity from O(2ⁿ * n³) to O(2ⁿ * n²)

## Performance Results

The following table shows the performance comparison between the two methods:

| n | Pair Count | Unique S | Matrix Power Runtime (s) | Poly Runtime (s) | Speedup (%) |
|---|------------|----------|--------------------------|------------------|-------------|
| 1 | 2          | 2        | 1.46                     | 1.17             | 19.86       |
| 2 | 15         | 14       | 5.37                     | 2.88             | 46.37       |
| 3 | 203        | 126      | 76.63                    | 19.6             | 74.42       |
| 4 | 4140       | 1518     | 2584.17                  | 310.42           | 87.98       |

## Key Features

- Implements efficient matrix power computation using polynomial arithmetic
- Uses SymPy for symbolic mathematics
- Parallel processing for handling large computations
- Detailed logging of results and performance metrics

## Implementation Details

The project consists of two main implementations:
1. `batch_calculate_total_sum_with_log.py`: Original implementation using direct matrix exponentiation
2. `batch_calculate_total_sum_with_log_poly.py`: Optimized implementation using polynomial-based modular exponentiation

## Results

Both implementations produce identical results, but the polynomial-based method shows significant performance improvements:
- Up to 88% faster for larger values of n
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