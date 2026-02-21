"""
export_table_data.py — Export Table 1 data to CSV for the repository
"""
import math, csv

def sieve(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return is_prime

def odd_radical(n):
    if n <= 1: return 1
    r, temp = 1, n
    while temp % 2 == 0: temp //= 2
    d = 3
    while d * d <= temp:
        if temp % d == 0:
            r *= d
            while temp % d == 0: temp //= d
        d += 2
    if temp > 1: r *= temp
    return r

LIMIT = 20000
is_prime = sieve(LIMIT)
N_t = 8192

rows = []
for p in range(3, N_t // 2 + 1):
    q = N_t - p
    if q <= 1: continue
    ro = odd_radical(p) * odd_radical(q)
    cond = ro * ro  # rad_odd(M)=1 for N=2^13
    rho = math.log(cond) / math.log(N_t) if cond > 1 else 0.0
    if is_prime[p] and is_prime[q]:
        typ = "Goldbach"
    elif not is_prime[p] and not is_prime[q]:
        typ = "Composite"
    else:
        typ = "Mixed"
    rows.append([N_t, p, q, typ, ro, cond, round(rho, 6)])

with open('/home/claude/repo/data/table1_N8192.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'p', 'q', 'type', 'rad_odd_pq', 'conductor_proxy', 'rho'])
    w.writerows(rows)

# Also export the stability band scan data
rows2 = []
for N_even in range(100, 10002, 2):
    rhos = []
    for p in range(3, N_even // 2 + 1, 2):
        q = N_even - p
        if is_prime[p] and is_prime[q]:
            ro = odd_radical(p) * odd_radical(q) * odd_radical(N_even//2)**2
            cond = ro * ro
            rhos.append(math.log(cond) / math.log(N_even) if cond > 1 else 0.0)
    if rhos:
        import numpy as np
        rows2.append([N_even, len(rhos), round(min(rhos), 6), round(float(np.mean(rhos)), 6), round(max(rhos), 6)])

with open('/home/claude/repo/data/stability_band_scan.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'num_goldbach_pairs', 'rho_min', 'rho_mean', 'rho_max'])
    w.writerows(rows2)

print(f"Exported {len(rows)} rows to table1_N8192.csv")
print(f"Exported {len(rows2)} rows to stability_band_scan.csv")
