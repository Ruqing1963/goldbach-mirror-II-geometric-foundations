"""
verify_discriminant.py — Independent Verification of Goldbach-Frey Discriminant
================================================================================
Verifies Theorem 3.2 and Theorem 5.1 of:
  R. Chen, "The Goldbach Mirror: Conductor Rigidity and the Static Conduit in GSp(4)"

Theorem 3.2:  Δ(f) = 2^12 · p^6 · (2N-p)^6 · (N-p)^4 · N^4
Theorem 5.1(iii):  ord_r(Δ) = 4·ord_r(N) when r ∤ p(2N-p)(N-p)
"""


def compute_discriminant(roots):
    """Compute Δ = ∏_{i<j} (e_i - e_j)^2 from a list of roots."""
    disc = 1
    n = len(roots)
    for i in range(n):
        for j in range(i + 1, n):
            disc *= (roots[i] - roots[j]) ** 2
    return disc


def ord_p(n, p):
    """Compute the p-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v


def verify_theorem_3_2():
    """Verify Δ = 2^12 · p^6 · (2N-p)^6 · (N-p)^4 · N^4."""
    print("=" * 60)
    print("Theorem 3.2: Discriminant Formula")
    print("Δ = 2^12 · p^6 · (2N-p)^6 · (N-p)^4 · N^4")
    print("=" * 60)

    test_cases = [
        (15, 7),    # 2N=30, p=7, 2N-p=23
        (10, 3),    # 2N=20, p=3, 2N-p=17
        (25, 7),    # 2N=50, p=7, 2N-p=43
        (50, 13),   # 2N=100, p=13, 2N-p=87
        (100, 29),  # 2N=200, p=29, 2N-p=171
        (49, 3),    # 2N=98, p=3 (N=49=7², tests square factors)
    ]

    all_pass = True
    for N, p in test_cases:
        q = 2 * N - p
        roots = [0, p, -p, q, -q]
        disc = compute_discriminant(roots)
        formula = (2**12) * (p**6) * (q**6) * ((N - p)**4) * (N**4)
        match = disc == formula
        if not match:
            all_pass = False
        print(f"  2N={2*N:>4}, p={p:>3}, 2N-p={q:>3}: "
              f"{'✓' if match else '✗ MISMATCH'}")

    print(f"\n  Result: {'ALL PASSED ✓' if all_pass else 'SOME FAILED ✗'}\n")
    return all_pass


def verify_theorem_5_1():
    """Verify ord_r(Δ) = 4·ord_r(N) at static conduit primes."""
    print("=" * 60)
    print("Theorem 5.1(iii): ord_r(Δ) = 4·ord_r(N)")
    print("for r > 2 with r ∤ p(2N-p)(N-p)")
    print("=" * 60)

    test_cases = [
        # (N, p, r, expected_ord_r_N)
        (15, 7, 3, 1),     # N=15, 3|N, ord_3(N)=1
        (15, 7, 5, 1),     # N=15, 5|N, ord_5(N)=1
        (49, 3, 7, 2),     # N=49=7², ord_7(N)=2
        (125, 7, 5, 3),    # N=125=5³, ord_5(N)=3
        (21, 5, 3, 1),     # N=21=3·7, r=3
        (21, 5, 7, 1),     # N=21=3·7, r=7
        (27, 5, 3, 3),     # N=27=3³, ord_3(N)=3
    ]

    all_pass = True
    for N, p, r, expected_v in test_cases:
        q = 2 * N - p
        # Check hypotheses: r ∤ p(2N-p)(N-p)
        if p % r == 0 or q % r == 0 or (N - p) % r == 0:
            print(f"  2N={2*N:>4}, p={p:>3}, r={r}: SKIP (hypothesis violated)")
            continue

        roots = [0, p, -p, q, -q]
        disc = compute_discriminant(roots)
        v_r = ord_p(disc, r)
        v_N = ord_p(N, r)
        expected = 4 * v_N
        match = v_r == expected
        if not match:
            all_pass = False
        print(f"  2N={2*N:>4}, p={p:>3}, r={r}: "
              f"ord_r(N)={v_N}, ord_r(Δ)={v_r}, "
              f"4·ord_r(N)={expected}  "
              f"{'✓' if match else '✗ MISMATCH'}")

    print(f"\n  Result: {'ALL PASSED ✓' if all_pass else 'SOME FAILED ✗'}\n")
    return all_pass


def verify_conduit_uniformity():
    """Verify Theorem 5.1(ii): conductor at r|N is independent of p."""
    print("=" * 60)
    print("Theorem 5.1(ii): Uniform conductor at static conduit")
    print("ord_r(Δ) is the same for all p with r ∤ p(2N-p)(N-p)")
    print("=" * 60)

    N, r = 15, 5  # r=5 divides N=15
    print(f"  N={N}, r={r}, ord_r(N)={ord_p(N, r)}")
    print(f"  Expected ord_r(Δ) = 4·{ord_p(N, r)} = {4*ord_p(N, r)}")
    print()

    results = []
    for p in range(3, 2 * N, 2):
        q = 2 * N - p
        if p <= 1 or q <= 1:
            continue
        if p % r == 0 or q % r == 0 or (N - p) % r == 0:
            continue
        roots = [0, p, -p, q, -q]
        disc = compute_discriminant(roots)
        v = ord_p(disc, r)
        results.append((p, v))
        print(f"    p={p:>3}: ord_{r}(Δ) = {v}")

    vals = set(v for _, v in results)
    uniform = len(vals) == 1
    print(f"\n  Uniform: {'✓' if uniform else '✗'} "
          f"(values: {vals})\n")
    return uniform


if __name__ == "__main__":
    print()
    r1 = verify_theorem_3_2()
    r2 = verify_theorem_5_1()
    r3 = verify_conduit_uniformity()

    print("=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"  Theorem 3.2 (discriminant formula):  "
          f"{'PASS ✓' if r1 else 'FAIL ✗'}")
    print(f"  Theorem 5.1(iii) (ord_r = 4·ord_r(N)): "
          f"{'PASS ✓' if r2 else 'FAIL ✗'}")
    print(f"  Theorem 5.1(ii) (uniformity):          "
          f"{'PASS ✓' if r3 else 'FAIL ✗'}")
    print()
