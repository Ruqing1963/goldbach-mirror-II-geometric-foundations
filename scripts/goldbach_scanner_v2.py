"""
Project Mirror v2: Goldbach Comet Geometric Scanner
====================================================
Corrected version. Changes from v1:

BUG FIX: v1 used the ORDERED Hardy-Littlewood formula (2*C2*...)
but counted UNORDERED pairs (p ≤ N). This produced ratios ≈ 0.56
instead of ≈ 1.0. v2 counts ORDERED pairs (both (p,q) and (q,p))
to match the standard HL prediction.

TERMINOLOGY FIX: "Sato-Tate Vol" renamed to "HL Prediction" since
the Hardy-Littlewood singular series is a number-theoretic object,
not a Sato-Tate volume.
"""
import math
import time


def sieve_of_eratosthenes(limit):
    """Return a bytearray where is_prime[i] = 1 iff i is prime."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return is_prime


def wormhole_factor(n, primes):
    """
    Compute the static conduit factor: ∏_{p|n, p>2} (p-1)/(p-2).
    This is the Hardy-Littlewood singular series correction for 2N = n.
    """
    factor = 1.0
    temp = n
    limit = int(math.isqrt(n))
    for p in primes:
        if p > limit:
            break
        if p > 2 and temp % p == 0:
            factor *= (p - 1) / (p - 2)
            while temp % p == 0:
                temp //= p
    if temp > 2:  # remaining large prime factor
        factor *= (temp - 1) / (temp - 2)
    return factor


def classify_orbit(n):
    """Classify 2N by its mod-6 structure (determines comet band)."""
    if n % 6 == 0:
        return "Resonant (mod 6)"   # upper band: 3 | N
    elif n % 3 == 0:
        return "Stable (mod 3)"
    else:
        return "Generic"            # lower band: 3 ∤ N


def run_project_mirror():
    print("=" * 70)
    print("  Project Mirror v2: Goldbach Comet Geometric Scanner")
    print("  Static Conduit Verification (HL formula corrected)")
    print("=" * 70)

    # --- Configuration ---
    START_2N = 10_000_000
    RANGE_SIZE = 2000
    END_2N = START_2N + RANGE_SIZE

    # 1. Sieve
    print(f"[*] Sieve up to {END_2N}...")
    t0 = time.time()
    is_prime = sieve_of_eratosthenes(END_2N)
    primes = [i for i, x in enumerate(is_prime) if x]
    print(f"[*] Found {len(primes)} primes. ({time.time()-t0:.2f}s)")
    print("-" * 70)

    # 2. Twin prime constant
    C2 = 0.6601618158468

    # 3. Scan
    print(f"[*] Scanning 2N from {START_2N} to {END_2N}...")
    print(f"{'2N':<12} | {'G(2N)':<10} | {'HL Pred':<12} | {'Ratio':<8} | {'Wormhole':<8} | {'Orbit Type'}")
    print("-" * 70)

    results = []

    for n in range(START_2N, END_2N + 2, 2):
        # A. Count ORDERED Goldbach representations
        #    G(2N) = #{ordered (p,q): p+q=2N, p,q prime}
        #    We count p from 2 to n-2 (both orderings counted)
        count = 0
        for p in primes:
            if p >= n:
                break
            if p > 1 and n - p > 1 and is_prime[n - p]:
                count += 1
        # count now includes both (p, n-p) and (n-p, p) when p ≠ n-p

        # B. Hardy-Littlewood prediction (ORDERED count)
        #    G(2N) ~ 2 * C2 * ∏_{p|N,p>2} (p-1)/(p-2) * 2N/(ln 2N)^2
        log_n = math.log(n)
        integral = n / (log_n ** 2)
        wf = wormhole_factor(n, primes)
        predicted = 2 * C2 * wf * integral

        ratio = count / predicted if predicted > 0 else 0
        orbit = classify_orbit(n)

        # Print selected rows
        if (n - START_2N) < 20 or (n - START_2N) % 200 == 0:
            print(f"{n:<12} | {count:<10} | {predicted:<12.1f} | {ratio:<8.4f} | {wf:<8.4f} | {orbit}")

        results.append((n, count, predicted, wf, orbit, ratio))

    print("-" * 70)

    # 4. Statistics by orbit type
    print("\n[*] Statistics by orbit type:")
    for otype in ["Resonant (mod 6)", "Stable (mod 3)", "Generic"]:
        subset = [r for r in results if r[4] == otype]
        if subset:
            avg_ratio = sum(r[5] for r in subset) / len(subset)
            avg_wf = sum(r[3] for r in subset) / len(subset)
            print(f"    {otype:<20}: n={len(subset):<5} avg_ratio={avg_ratio:.6f}  avg_wormhole={avg_wf:.4f}")

    avg_all = sum(r[5] for r in results) / len(results)
    print(f"\n[*] Overall average ratio: {avg_all:.6f}")
    print("[*] If HL formula is correct, ratio should be close to 1.0")
    print("[*] Deviations from 1.0 reflect the integral approximation")
    print("    (Li_2(n) vs n/(ln n)^2), which converges slowly.")

    # 5. Discriminant verification for a small example
    print("\n" + "=" * 70)
    print("[*] Discriminant verification (Theorem 2.1)")
    print("=" * 70)
    N_test, p_test = 15, 7  # 2N=30, p=7, 2N-p=23 (both prime!)
    roots = [0, p_test, -p_test, 2*N_test - p_test, -(2*N_test - p_test)]
    disc = 1
    for i in range(len(roots)):
        for j in range(i+1, len(roots)):
            disc *= (roots[i] - roots[j])**2
    
    # Predicted: 2^12 * p^6 * (2N-p)^6 * (N-p)^4 * N^4
    p, q = p_test, 2*N_test - p_test
    Nm = N_test
    predicted_disc = (2**12) * (p**6) * (q**6) * ((Nm-p)**4) * (Nm**4)
    
    print(f"  2N=30, p=7, 2N-p=23")
    print(f"  Roots: {roots}")
    print(f"  Computed Δ = {disc}")
    print(f"  Formula  Δ = 2^12 · 7^6 · 23^6 · 8^4 · 15^4 = {predicted_disc}")
    print(f"  Match: {'✓' if disc == predicted_disc else '✗ — MISMATCH'}")
    if disc != predicted_disc:
        print(f"  Ratio: {disc / predicted_disc}")


if __name__ == "__main__":
    run_project_mirror()
