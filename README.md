# The Goldbach Mirror II: Geometric Foundations of Conductor Rigidity and the Static Conduit in GSp(4)

**Author:** Ruqing Chen, GUT Geoservice Inc., Montréal, QC, Canada

## Paper

> R. Chen, *The Goldbach Mirror II: Geometric Foundations of Conductor Rigidity and the Static Conduit in GSp(4)*, 2026.

This is the eighth paper in the conductor rigidity series, building on:

| # | Paper | Zenodo |
|---|-------|--------|
| 1 | Conductor Incompressibility for Frey Curves Associated to Prime Gaps | [10.5281/zenodo.18682375](https://zenodo.org/records/18682375) |
| 2 | Density Thresholds for Equidistribution in Prime-Indexed Geometric Families | [10.5281/zenodo.18682721](https://zenodo.org/records/18682721) |
| 3 | Weil Restriction Rigidity and Prime Gaps via Genus 2 Hyperelliptic Jacobians | [10.5281/zenodo.18683194](https://zenodo.org/records/18683194) |
| 4 | On Landau's Fourth Problem: Conductor Rigidity and Sato–Tate Equidistribution | [10.5281/zenodo.18683712](https://zenodo.org/records/18683712) |
| 5 | The 2-2 Coincidence: Conductor Rigidity for Primes in Arithmetic Progressions | [10.5281/zenodo.18684151](https://zenodo.org/records/18684151) |
| 6 | The Genesis of Prime Constellations: Weil Restriction on GSp(8) | [10.5281/zenodo.18684352](https://zenodo.org/records/18684352) |
| 7 | The Goldbach Mirror (companion paper) | [10.5281/zenodo.18684892](https://zenodo.org/records/18684892) |
| 8 | **The Goldbach Mirror II** (this paper) | *forthcoming* |

**Companion repository:** [goldbach-mirror-conductor-rigidity](https://github.com/Ruqing1963/goldbach-mirror-conductor-rigidity) (Paper #7)

## Abstract

Building on the algebraic framework of Paper #7, we introduce **Chen's ratio** ρ(N,p) = log 𝒩(J_{N,p}) / log N as a normalised measure of the conductor of the Goldbach–Frey Jacobian in the Siegel moduli space 𝒜₂. Computational scanning over N ∈ [10², 10⁴] reveals that Goldbach pairs are rigidly confined to a narrow ρ-band—the **static conduit stability band**—while composite decompositions spread over a range roughly three times wider. At N = 2^k, the odd radical of the static conduit vanishes, producing dramatic **conductor dips** that isolate the pure boundary-prime structure. We formalise these observations through a **geometric obstruction analysis** and identify the precise analytic gap (effective GSp(4) Sato–Tate equidistribution) required to convert this framework into a proof of the Goldbach conjecture.

## Key Results

| Result | Section | Status |
|--------|---------|--------|
| Chen's ratio ρ and conductor proxy 𝒩_proxy | §2 | Defined and justified (via Ogg–Saito) |
| Radical rigidity: Goldbach pairs satisfy rad(p) = p | Prop 2.4 | Proved unconditionally |
| Composite spread: rad_odd(p) can collapse arbitrarily | Prop 2.5 | Proved unconditionally |
| Static conduit stability band (Figure 1) | §3 | Computational (N ∈ [100, 10000]) |
| 2^k anchor phenomenon | §4, Prop 4.1 | Proved unconditionally |
| N = 8192 conductor barrier (Figure 2, Table 1) | §4.2 | Computational (76 Goldbach pairs) |
| Geometric obstruction analysis | §5 | Framework established |
| Singular series = aperture of static conduit | Obs 6.1 | Sieve-theoretic interpretation |
| Analytic gap identification | Problem 5.3, Remark 5.5 | Open problem formulated |

## Repository Contents

```
├── README.md
├── LICENSE
├── paper/
│   ├── Goldbach_Mirror_II_Geometric_Foundations.tex    # LaTeX source
│   └── Goldbach_Mirror_II_Geometric_Foundations.pdf    # Compiled paper (9 pages)
├── figures/
│   ├── figure1.pdf            # Static Conduit Stability Band
│   ├── figure1.png            # (PNG version for preview)
│   ├── figure2.pdf            # The Conductor Barrier at N = 2^13
│   └── figure2.png            # (PNG version for preview)
├── scripts/
│   ├── generate_figures.py    # Generates both figures and Table 1 data
│   ├── goldbach_scanner_v2.py # Hardy-Littlewood verification (from Paper #7)
│   ├── verify_discriminant.py # Discriminant formula verification (from Paper #7)
│   └── export_table_data.py   # Exports CSV data files
└── data/
    ├── table1_N8192.csv          # All 4094 decompositions of N=8192
    └── stability_band_scan.csv   # ρ statistics for N ∈ [100, 10000]
```

## Figures

### Figure 1 — Static Conduit Stability Band

![Figure 1](figures/figure1.png)

The mean Chen's ratio ⟨ρ⟩_N over all Goldbach pairs of each even N, scanned from 100 to 10,000. The shaded region marks the 5th–95th percentile envelope. Red triangles mark N = 2^k (k = 7, …, 13), which fall dramatically below the band due to the vanishing of rad_odd(M).

### Figure 2 — The Conductor Barrier at N = 2¹³ = 8192

![Figure 2](figures/figure2.png)

**Left:** Chen's ratio for all decompositions of N = 8192. Goldbach pairs (blue, n = 76) cluster in a narrow band; composite decompositions (orange, n = 3144) spread across the full range. **Right:** Density histogram confirming the clustering. Mean gap Δρ = 0.558.

## Quick Start

All scripts use only the Python standard library plus NumPy and Matplotlib:

```bash
# Install dependencies
pip install numpy matplotlib

# Generate figures and table data
python scripts/generate_figures.py

# Run the Hardy-Littlewood scanner (from Paper #7)
python scripts/goldbach_scanner_v2.py

# Verify discriminant formula (from Paper #7)
python scripts/verify_discriminant.py

# Export CSV data files
python scripts/export_table_data.py
```

### Compiling the paper

```bash
cd paper
# Copy figures into the paper directory
cp ../figures/figure1.pdf ../figures/figure2.pdf .
pdflatex Goldbach_Mirror_II_Geometric_Foundations.tex
pdflatex Goldbach_Mirror_II_Geometric_Foundations.tex   # twice for cross-refs
```

## Data Format

### table1_N8192.csv

| Column | Description |
|--------|-------------|
| `N` | Even number (8192) |
| `p` | First summand |
| `q` | Second summand (N − p) |
| `type` | `Goldbach` (both prime), `Composite` (both composite), or `Mixed` |
| `rad_odd_pq` | rad_odd(p) × rad_odd(q) |
| `conductor_proxy` | (rad_odd(p) × rad_odd(q))² (since rad_odd(M) = 1 for N = 2^k) |
| `rho` | Chen's ratio: log(conductor_proxy) / log(N) |

### stability_band_scan.csv

| Column | Description |
|--------|-------------|
| `N` | Even number |
| `num_goldbach_pairs` | Count of Goldbach pairs |
| `rho_min` | Minimum ρ over Goldbach pairs |
| `rho_mean` | Mean ρ over Goldbach pairs |
| `rho_max` | Maximum ρ over Goldbach pairs |

## Requirements

- Python 3.6+
- NumPy
- Matplotlib
- LaTeX (for compiling the paper; texlive recommended)

## Citation

```bibtex
@article{Chen2026GMII,
  author  = {Ruqing Chen},
  title   = {The Goldbach Mirror {II}: Geometric Foundations of Conductor
             Rigidity and the Static Conduit in {GSp(4)}},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/Ruqing1963/goldbach-mirror-II-geometric-foundations}
}
```

## License

This work is released under the [MIT License](LICENSE).
