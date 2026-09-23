# Numerical Solution of the Schrödinger Equation — Infinite Square Well

A Python implementation that solves the time-independent Schrödinger equation for a particle
in an infinite square well using the finite-difference method, then verifies the numerical
results against the known analytic solution.

Originally written as a final-year physics project at the University of Botswana.

## Overview

The time-independent Schrödinger equation is:

$$-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi(x) = E\psi(x)$$

For an infinite square well of width `W`, the potential `V(x)` is zero inside the well and
infinite outside it, which forces the boundary conditions ψ(0) = ψ(W) = 0. This project
discretizes that equation on a grid, turns it into a matrix eigenvalue problem, and solves
for the allowed energies and wavefunctions numerically.

## Method

1. **Discretize space** into `N` grid points across the well, excluding the two boundary
   points (where ψ is already known to be zero).
2. **Build the Hamiltonian matrix** using the standard second-order central-difference
   approximation of the second derivative, which produces a tridiagonal matrix.
3. **Diagonalize the Hamiltonian** with `scipy.linalg.eigh` to get all eigenvalues (energy
   levels) and eigenvectors (wavefunctions) in one call.
4. **Normalize** each wavefunction so that the total probability integrates to 1.
5. **Compare** the numerical energy levels against the analytic formula:

$$E_n = \frac{n^2\pi^2\hbar^2}{2mW^2}, \quad n = 1, 2, 3, \dots$$

6. **Visualize** the probability density |ψ(x)|² for the lowest few energy levels, in both
   a 3D surface plot and a 2D line plot.

Atomic units are used throughout (ħ = 1, m = 1), so all energies are in Hartree and all
lengths are in Bohr radii.

## Example output

With `W = 2` and `N = 100`, the first four numerical energy levels match the analytic
values to several decimal places, confirming that the finite-difference approximation
converges correctly. The plots show the familiar sinusoidal probability density patterns
for each quantum number `n`, with `n` nodes for the `n`-th excited state.

*(Add your own screenshots of the 3D and 2D plots here once you run the script — see
"Adding screenshots" below.)*

## Requirements

- Python 3.8+
- numpy
- scipy
- matplotlib

Install everything with:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python schrodinger_square_well.py
```

This will:
- print the first four numerical energy eigenvalues alongside the analytic values
- open a 3D plot of probability density vs. position and energy level
- open a 2D plot comparing the probability densities directly

You can change the well width, number of grid points, or number of energy levels shown by
editing the `W`, `N`, and `num_levels` variables near the top of the script.

## Adding screenshots

`plt.show()` opens an interactive window but doesn't save an image. To capture plots for
this README, add a line like this before each `plt.show()` call:

```python
plt.savefig("figures/probability_density_3d.png", dpi=150)
```

Then reference the saved image in this README with:

```markdown
![3D probability density](figures/probability_density_3d.png)
```

## Project structure

```
.
├── schrodinger_square_well.py   # Main script
├── requirements.txt             # Python dependencies
├── README.md                    # This file
└── LICENSE                      # MIT license
```

## Background

This project was completed as part of a BSc Physics research project at the University of
Botswana, applying numerical linear algebra methods to a foundational problem in quantum
mechanics.

## Author

**Montle Fredah Segomotso**

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
