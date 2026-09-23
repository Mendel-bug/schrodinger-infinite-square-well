'''Proving the Schrodinger equation for an infinite square well using python
By Montle Fredah Segomotso 202107060'''

# Import necessary libraries
import numpy as np                          # For numerical operations
import matplotlib.pyplot as plt              # For plotting
from mpl_toolkits.mplot3d import Axes3D      # For 3D plotting
from scipy.linalg import eigh                # For eigenvalue problems

# ─────────────────────────────────────────────
# Build the Hamiltonian matrix for the infinite square well
# ─────────────────────────────────────────────
def schrodinger_hamiltonian(N, W):
    """
    Constructs the finite-difference Hamiltonian matrix for a particle
    in an infinite square well of width W, discretized on N-1 interior
    grid points (the walls at x=0 and x=W are excluded since psi=0 there).

    Uses the standard second-order central-difference approximation of
    the kinetic energy operator:  -hbar^2/(2M) * d^2/dx^2.
    """
    M = 1.0      # Mass in atomic units
    Hbar = 1.0   # Reduced Planck's constant in atomic units

    dx = W / (N - 1)                       # Correct grid spacing (parentheses fixed)
    D = -(Hbar**2) / (2 * M * dx**2)       # Correct off-diagonal coefficient (grouping fixed)

    H = np.zeros((N - 1, N - 1))

    # Fill the tridiagonal matrix correctly: lower, main, and upper diagonals
    for x in range(N - 1):
        H[x, x] = -2 * D                   # Main diagonal (no energy term needed; eigh finds E directly)
        if x > 0:
            H[x, x - 1] = D                # Lower diagonal
        if x < N - 2:
            H[x, x + 1] = D                # Upper diagonal (this was missing entirely before)

    # Infinite walls are already enforced by excluding the boundary points
    # from the grid (Dirichlet boundary condition psi(0)=psi(W)=0), so no
    # artificial large penalty value is needed.
    return H


# ─────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────
W = 2     # Infinite square well width
N = 100   # Number of grid points
num_levels = 4  # Number of energy levels/wavefunctions to show

# ─────────────────────────────────────────────
# Solve for eigenvalues (energies) and eigenvectors (wavefunctions) directly.
# eigh on a single, correctly-built Hamiltonian returns ALL eigenvalues
# sorted in ascending order in one shot — no bisection search needed.
# ─────────────────────────────────────────────
H = schrodinger_hamiltonian(N, W)
all_eigenvalues, all_eigenvectors = eigh(H)

energy_levels = all_eigenvalues[:num_levels]
print(f'Found energy eigenvalues: {energy_levels}')

# Compare against the analytic solution E_n = n^2 * pi^2 * hbar^2 / (2*M*W^2)
analytic = [ (n**2) * np.pi**2 / (2 * W**2) for n in range(1, num_levels + 1) ]
print(f'Analytic energy eigenvalues: {analytic}')

# ─────────────────────────────────────────────
# Build wavefunctions and probability densities for each level
# ─────────────────────────────────────────────
x_interior = np.linspace(0, W, N + 1)[1:-1]   # Interior grid points (excludes the walls)
X, Y = np.meshgrid(x_interior, energy_levels)  # Meshgrid for 3D plotting
Z = np.zeros_like(X)

for i in range(num_levels):
    wavefunction = all_eigenvectors[:, i].copy()     # i-th eigenvector = i-th energy level's wavefunction

    # Normalize so that the integral of |psi|^2 dx = 1
    dx = W / (N - 1)
    wavefunction /= np.sqrt(np.sum(wavefunction**2) * dx)

    Z[i, :] = wavefunction**2   # Probability density

# ─────────────────────────────────────────────
# 3D Visualization
# ─────────────────────────────────────────────
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='Greens')

for i in range(num_levels):
    ax.plot(x_interior, [energy_levels[i]] * len(x_interior), Z[i, :],
            label=f'n={i+1}', linewidth=2)

ax.set_title('3D Wave Function Probability Density For Infinite Square Well')
ax.set_xlabel('Position (x)')
ax.set_ylabel('Energy Levels')
ax.set_zlabel('Probability Density')
ax.legend()
plt.tight_layout()
plt.show()

# ─────────────────────────────────────────────
# 2D Visualization for clarity
# ─────────────────────────────────────────────
plt.figure(figsize=(10, 6))
for i in range(num_levels):
    plt.plot(x_interior, Z[i, :], label=f'n={i+1}')
plt.title('Wave Function Probability Density for Infinite Square Well')
plt.xlabel('Position (x)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
