"""
generate_figures_v2.py — Corrected Conductor Computation
=========================================================
Conductor proxy = (rad_odd(p) * rad_odd(q) * rad_odd(N/2)^2)^2
matching the values in the user's draft table.

For N = 2^k, rad_odd(N/2) = 1, so conductor = (rad_odd(p)*rad_odd(q))^2.
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─── Utilities ────────────────────────────────────────────────────────────────

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

def conductor_proxy(N_even, p):
    """Squared odd-radical conductor: N_proxy = (rad_odd(p)*rad_odd(q)*rad_odd(N/2)^2)^2"""
    q = N_even - p
    half_N = N_even // 2
    base = odd_radical(p) * odd_radical(q) * odd_radical(half_N)**2
    return base * base  # squared

def chens_ratio(N_even, p):
    cond = conductor_proxy(N_even, p)
    if cond <= 1 or N_even <= 1: return 0.0
    return math.log(cond) / math.log(N_even)

# ─── Setup ────────────────────────────────────────────────────────────────────

LIMIT = 20000
print("[*] Sieving...")
is_prime = sieve(LIMIT)

# ═══════════════════════════════════════════════════════════════════════════════
# TABLE 1: N = 8192
# ═══════════════════════════════════════════════════════════════════════════════
N_t = 8192
print(f"\n{'='*78}")
print(f"TABLE 1: Conductor Compression Metrics for N = {N_t}")
print(f"{'='*78}")

gb, comp, mixed = [], [], []
for p in range(3, N_t // 2 + 1):
    q = N_t - p
    if q <= 1: continue
    rho = chens_ratio(N_t, p)
    ro = odd_radical(p) * odd_radical(q)
    cond = conductor_proxy(N_t, p)
    row = (p, q, ro, cond, rho)
    if is_prime[p] and is_prime[q]:
        gb.append(row)
    elif not is_prime[p] and not is_prime[q]:
        comp.append(row)
    else:
        mixed.append(row)

gb.sort(key=lambda x: x[4])
comp.sort(key=lambda x: x[4], reverse=True)  # highest-ρ composites first
mixed.sort(key=lambda x: x[4])

print(f"\n{'Type':<12} {'(p, q)':<18} {'rad_odd(pq)':<14} {'Conductor':<16} {'ρ':<8}")
print("-"*78)
for row in gb[:5]:
    p,q,ro,cond,rho = row
    print(f"{'Goldbach':<12} ({p}, {q}){'':<{12-len(f'({p}, {q})')}} {ro:<14} {cond:<16.4e} {rho:<8.4f}")
print("  ...")
for row in gb[-3:]:
    p,q,ro,cond,rho = row
    print(f"{'Goldbach':<12} ({p}, {q}){'':<{12-len(f'({p}, {q})')}} {ro:<14} {cond:<16.4e} {rho:<8.4f}")
print("-"*78)
for row in comp[:5]:
    p,q,ro,cond,rho = row
    print(f"{'Composite':<12} ({p}, {q}){'':<{12-len(f'({p}, {q})')}} {ro:<14} {cond:<16.4e} {rho:<8.4f}")
print("-"*78)
for row in mixed[:3]:
    p,q,ro,cond,rho = row
    which = "p prime" if is_prime[p] else "q prime"
    print(f"{'Mixed':<12} ({p}, {q}){'':<{12-len(f'({p}, {q})')}} {ro:<14} {cond:<16.4e} {rho:<8.4f}")

gb_rhos = [r[4] for r in gb]
comp_rhos = [r[4] for r in comp if r[4] > 0]
mixed_rhos = [r[4] for r in mixed if r[4] > 0]

print(f"\n  Goldbach (n={len(gb)}):   ρ ∈ [{min(gb_rhos):.4f}, {max(gb_rhos):.4f}], mean={np.mean(gb_rhos):.4f}")
if comp_rhos:
    print(f"  Composite (n={len(comp)}): ρ ∈ [{min(comp_rhos):.4f}, {max(comp_rhos):.4f}], mean={np.mean(comp_rhos):.4f}")
if mixed_rhos:
    print(f"  Mixed (n={len(mixed)}):    ρ ∈ [{min(mixed_rhos):.4f}, {max(mixed_rhos):.4f}], mean={np.mean(mixed_rhos):.4f}")

# Key structural observation
print(f"\n  Structural observation:")
print(f"    Goldbach pairs concentrate in band ρ ∈ [{min(gb_rhos):.2f}, {max(gb_rhos):.2f}]")
print(f"    Composite pairs spread over ρ ∈ [{min(comp_rhos):.2f}, {max(comp_rhos):.2f}]")
print(f"    Goldbach band width: {max(gb_rhos)-min(gb_rhos):.4f}")
print(f"    Composite spread:    {max(comp_rhos)-min(comp_rhos):.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# Multi-N scan for the stability band plot
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[*] Scanning ρ across N = 100..10000 for Figure 1...")

N_values = list(range(100, 10002, 2))
avg_rho_list, min_rho_list, max_rho_list, N_plot = [], [], [], []

for N_even in N_values:
    rhos = []
    for p in range(3, N_even // 2 + 1, 2):
        q = N_even - p
        if is_prime[p] and is_prime[q]:
            rhos.append(chens_ratio(N_even, p))
    if rhos:
        avg_rho_list.append(np.mean(rhos))
        min_rho_list.append(min(rhos))
        max_rho_list.append(max(rhos))
        N_plot.append(N_even)

N_arr = np.array(N_plot)
avg_arr = np.array(avg_rho_list)
min_arr = np.array(min_rho_list)
max_arr = np.array(max_rho_list)

# 2^k points
powers_of_2 = [128, 256, 512, 1024, 2048, 4096, 8192]
p2_N, p2_avg, p2_min = [], [], []
for pk in powers_of_2:
    if pk in N_plot:
        idx = N_plot.index(pk)
        p2_N.append(pk)
        p2_avg.append(avg_rho_list[idx])
        p2_min.append(min_rho_list[idx])

print("  2^k points:")
for i,pk in enumerate(p2_N):
    print(f"    N=2^{int(math.log2(pk))}={pk}: avg_ρ={p2_avg[i]:.4f}, min_ρ={p2_min[i]:.4f}")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Static Conduit Stability Band
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[*] Generating Figure 1...")

fig, ax = plt.subplots(figsize=(12, 6))

# Rolling average
w = 40
avg_smooth = np.convolve(avg_arr, np.ones(w)/w, mode='valid')
N_smooth = N_arr[w//2:w//2+len(avg_smooth)]

# Rigidity band envelope
band_lo = np.percentile(avg_arr, 5)
band_hi = np.percentile(avg_arr, 95)
ax.fill_between(N_arr, band_lo, band_hi, color='#3366CC', alpha=0.07,
                label=f'GSp(4) rigidity band $[{band_lo:.2f},\\,{band_hi:.2f}]$')

# Scatter (subsample for clarity)
ax.scatter(N_arr[::4], avg_arr[::4], s=1.5, color='#6699DD', alpha=0.3, zorder=2)

# Smoothed trend
ax.plot(N_smooth, avg_smooth, color='#1a3a6e', linewidth=1.8, alpha=0.85,
        label='Rolling mean $\\langle\\rho\\rangle$ (Goldbach pairs)', zorder=3)

# 2^k markers
ax.scatter(p2_N, p2_avg, s=90, color='#CC3333', marker='v', zorder=5,
           edgecolors='#881111', linewidth=0.8,
           label='$N = 2^k$ anchor points')
for i, pk in enumerate(p2_N):
    ax.annotate(f'$2^{{{int(math.log2(pk))}}}$',
                (pk, p2_avg[i]), textcoords="offset points",
                xytext=(10, -10), fontsize=9, color='#CC3333', fontweight='bold')

med = np.median(avg_arr)
ax.axhline(y=med, color='#444444', linestyle=':', linewidth=0.7, alpha=0.5,
           label=f'Median $\\rho = {med:.3f}$')

ax.set_xlabel('Even number $N$', fontsize=13)
ax.set_ylabel("Chen's Ratio  $\\rho = \\log\\mathcal{N} / \\log N$", fontsize=13)
ax.set_title('Figure 1 — Static Conduit Stability Band in the Siegel Moduli Space', fontsize=14)
ax.legend(loc='lower right', fontsize=10, framealpha=0.92)
ax.grid(True, alpha=0.15)
ax.set_xlim(50, 10300)

plt.tight_layout()
plt.savefig('/home/claude/figure1.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/figure1.png', dpi=200, bbox_inches='tight')
print("    Saved figure1.pdf/.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: The Barrier at N = 8192
# ═══════════════════════════════════════════════════════════════════════════════
print("[*] Generating Figure 2...")

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5),
                                 gridspec_kw={'width_ratios': [3, 1.2]})

# Collect all decompositions
p_gb, rho_gb = [], []
p_comp, rho_comp = [], []
p_mix, rho_mix = [], []

for p in range(3, N_t // 2 + 1):
    q = N_t - p
    if q <= 1: continue
    rho = chens_ratio(N_t, p)
    if rho <= 0: continue
    if is_prime[p] and is_prime[q]:
        p_gb.append(p); rho_gb.append(rho)
    elif not is_prime[p] and not is_prime[q]:
        p_comp.append(p); rho_comp.append(rho)
    else:
        p_mix.append(p); rho_mix.append(rho)

# Left: scatter
ax1.scatter(p_comp, rho_comp, s=2, color='#CC8855', alpha=0.15,
            label=f'Both composite ($n={len(p_comp)}$)', zorder=1)
ax1.scatter(p_mix, rho_mix, s=3, color='#88AA55', alpha=0.2,
            label=f'Mixed ($n={len(p_mix)}$)', zorder=2)
ax1.scatter(p_gb, rho_gb, s=12, color='#2255BB', alpha=0.7,
            label=f'Goldbach pairs ($n={len(p_gb)}$)', zorder=3)

# Mean lines
mean_gb = np.mean(rho_gb)
mean_comp = np.mean(rho_comp)
ax1.axhline(y=mean_gb, color='#2255BB', linestyle='--', linewidth=1.2, alpha=0.6)
ax1.axhline(y=mean_comp, color='#CC8855', linestyle='--', linewidth=1.2, alpha=0.6)

# Annotate gap
mid = (mean_gb + mean_comp) / 2
ax1.annotate('', xy=(N_t//2 - 200, mean_gb), xytext=(N_t//2 - 200, mean_comp),
             arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.5))
ax1.text(N_t//2 - 100, mid, f'$\\Delta\\rho = {abs(mean_gb-mean_comp):.3f}$',
         fontsize=10, va='center', color='#333333')

ax1.set_xlabel('Summand $p$', fontsize=12)
ax1.set_ylabel("$\\rho(8192, p)$", fontsize=12)
ax1.set_title('Decomposition landscape at $N = 2^{13}$', fontsize=13)
ax1.legend(loc='lower right', fontsize=9, framealpha=0.92)
ax1.grid(True, alpha=0.15)

# Right: histogram
all_rho = rho_gb + rho_comp + rho_mix
lo, hi = min(all_rho) - 0.05, max(all_rho) + 0.05
bins = np.linspace(lo, hi, 55)

ax2.hist(rho_gb, bins=bins, orientation='horizontal', color='#2255BB',
         alpha=0.6, label='Goldbach', density=True)
ax2.hist(rho_comp, bins=bins, orientation='horizontal', color='#CC8855',
         alpha=0.35, label='Composite', density=True)
ax2.hist(rho_mix, bins=bins, orientation='horizontal', color='#88AA55',
         alpha=0.25, label='Mixed', density=True)

ax2.set_xlabel('Density', fontsize=11)
ax2.set_ylabel("$\\rho$", fontsize=11)
ax2.set_title('Distribution', fontsize=12)
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.15)

plt.suptitle('Figure 2 — The Conductor Barrier at $N = 2^{13} = 8192$',
             fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('/home/claude/figure2.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/figure2.png', dpi=200, bbox_inches='tight')
print("    Saved figure2.pdf/.png")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary for LaTeX
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*78}")
print("SUMMARY FOR LATEX")
print(f"{'='*78}")
print(f"Global stability band: [{band_lo:.4f}, {band_hi:.4f}]")
print(f"Median ρ: {med:.4f}")
print(f"N=8192 Goldbach: mean={mean_gb:.4f}, band=[{min(rho_gb):.4f}, {max(rho_gb):.4f}]")
print(f"N=8192 Composite: mean={mean_comp:.4f}, band=[{min(rho_comp):.4f}, {max(rho_comp):.4f}]")
print(f"N=8192 gap: {abs(mean_gb-mean_comp):.4f}")
print(f"Ground state: p={gb[0][0]}, q={gb[0][1]}, ρ={gb[0][4]:.4f}")

# Verify user's table values
print(f"\nVerification of user's table:")
for pp in [13, 31]:
    qq = 8192 - pp
    ro = odd_radical(pp)*odd_radical(qq)
    cc = ro*ro  # since rad_odd(N/2)=1 for N=2^13
    rr = math.log(cc)/math.log(8192)
    print(f"  p={pp}, q={qq}: rad_odd={ro}, cond={cc:.4e}, ρ={rr:.4f}")
