import sympy as sp
import time
import pickle

def calculate_second_moment():
    print("Starting Second Moment Calculation...")
    t0 = time.time()
    
    z = sp.Symbol("z")
    alpha = sp.Rational(1024)
    t = 2
    
    ms = []
    ms.append(sp.Matrix([[1, 0], [0, 1]]))
    ms.append(sp.Matrix([[1 + z, 0], [0, 1]]))
    ms.append(sp.Matrix([[1, z], [0, 1]]))
    ms.append(sp.Matrix([[1, 0], [z, 1 + z]]))
    ms.append(sp.Matrix([[1 + z, 0], [0, 1 + z]]))
    ms.append(sp.Matrix([[1, z], [z, 1]]))

    def R_func(M):
        return M[0, 0] + M[1, 1] - M[0, 1] - M[1, 0]
    
    def Q_func(M):
        return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    
    def A_jk(M):
        z_alpha = z / alpha
        Q_eval = Q_func(M).subs(z, z_alpha)
        R_eval = R_func(M).subs(z, z_alpha)
        numerator = Q_eval
        denominator = (1 - z) * Q_eval + z_alpha**t * R_eval
        return numerator / denominator

    def nth_coefficient_from_R(R_expr, n_val=2**21):
        z = sp.symbols('z')
        R_expr = sp.simplify(R_expr)
        
        P_expr, Q_expr = sp.fraction(R_expr)
        P = sp.poly(P_expr, z)
        Q = sp.poly(Q_expr, z)

        # Solve Q(z) = 0
        z_roots = sp.solve(Q, z)
        
        # rho_k = 1 / z_k
        rhos = [1 / root for root in z_roots]
        Q_prime = Q.diff()

        a_values = []
        for rho in rhos:
            z_inv = 1 / rho
            numerator = -rho * P.as_expr().subs(z, z_inv)
            denominator = Q_prime.as_expr().subs(z, z_inv)
            a_k = sp.simplify(numerator / denominator)
            a_values.append(a_k)

        # Final result: sum(a_k * rho_k^n)
        # We perform the sum symbolically then evalf closely
        result = sum(a * rho**n_val for a, rho in zip(a_values, rhos))
        
        # KEY IMPROVEMENT: 
        # 1. Use higher precision (100 digits)
        # 2. Explicitly take the real part sp.re() to remove imaginary noise
        return sp.re(result.evalf(100))

    # Generate the 6 generating functions
    A_jks = [A_jk(M) for M in ms]

    print("Calculating coefficients (this may take a moment)...")
    A_jks_coefficient = [nth_coefficient_from_R(A_jk) for A_jk in A_jks]
    
    # Ns definition
    sigma = 2**10
    n1 = sigma * (sigma - 1)**2 * (sigma - 2)
    n2 = 2 * sigma * (sigma - 1) * (sigma - 2)
    n3 = 2 * sigma * (sigma - 1) * (sigma - 2)
    n4 = 4 * sigma * (sigma - 1)
    n5 = sigma * (sigma - 1)
    n6 = sigma * (sigma - 1)
    Ns = [n1, n2, n3, n4, n5, n6]
    
    weighted_score = [coeff * n for coeff, n in zip(A_jks_coefficient, Ns)]
    sum_weighted_score = sum(weighted_score)
    
    t1 = time.time()
    duration = t1 - t0
    
    print(f"Sum of Weighted Scores: {sum_weighted_score}")
    print(f"Time for Second Moment Sum: {duration:.4f} seconds")
    
    results = {
        "sum_weighted_score": sum_weighted_score,
        "time_seconds": duration
    }
    
    with open("moment2.pkl", "wb") as f:
        pickle.dump(results, f)
    print("Saved results to moment2.pkl")

if __name__ == "__main__":
    calculate_second_moment()
