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


- **Performance Optimization**: Develop and compare efficient methods for computing matrix powers in large state spaces
- **Correctness Verification**: Ensure both methods produce identical results while optimizing performance
- **Scalability Analysis**: Study how different methods scale with increasing problem size
- **Practical Application**: Apply optimized methods to PRNG state transition analysis

## Performance Results

The following table shows the performance comparison between the two methods:

| n | Pair Count | Unique S | Original Runtime (s) | Poly Runtime (s) | Speedup (%) |
|---|------------|----------|----------------------|------------------|-------------|
| 1 | 2          | 2        | 1.46                 | 1.17             | 19.86       |
| 2 | 15         | 14       | 5.37                 | 2.88             | 46.37       |
| 3 | 203        | 126      | 76.63                | 19.6             | 74.42       |
| 4 | 4140       | 1518     | 2584.17              | 310.42           | 87.98       |

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