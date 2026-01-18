import sympy as sp
import time
import pickle

def calculate_first_moment():
    print("Starting First Moment Calculation...")
    t0 = time.time()
    
    z = sp.Symbol("z")
    alpha = sp.Rational(1024)
    t = 2
    
    # Define A_j generator
    def A_j(P):
        z_alpha = z / alpha
        num = P.subs(z, z_alpha)
        denom = z_alpha**t + (1 - z) * num
        return num / denom

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

    P_distinct = sp.Poly(1, z)
    P_repeat = sp.Poly(1 + z, z)
    Aj_distinct = A_j(P_distinct.as_expr())
    Aj_repeat = A_j(P_repeat.as_expr())

    print("Calculating coefficients...")
    repeat_coeff = nth_coefficient_from_R(Aj_repeat)
    distinct_coeff = nth_coefficient_from_R(Aj_distinct)
    
    E_x = (alpha**2 - alpha) * distinct_coeff + alpha * repeat_coeff
    
    t1 = time.time()
    duration = t1 - t0
    print(f"Expectation (E_x): {E_x}")
    print(f"Time for Expectation: {duration:.4f} seconds")
    
    results = {
        "E_x": E_x,
        "time_seconds": duration
    }
    
    with open("moment1.pkl", "wb") as f:
        pickle.dump(results, f)
    print("Saved results to moment1.pkl")

if __name__ == "__main__":
    calculate_first_moment()
