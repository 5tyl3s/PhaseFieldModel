import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Constants
# -----------------------------
k  = 1.380649e-23      # Boltzmann (J/K)
h  = 6.62607015e-34    # Planck (J·s)
NA = 6.02214076e23     # Avogadro (1/mol)
R  = 8.314462618       # Gas constant (J/mol·K)

# -----------------------------
# Material / Model Parameters
# (match your config)
# -----------------------------
gamma = 0.0292  # J/m^2

drivingForceSlope = 0.0212 * 1000   # J/mol/
drivingForceIntercept = -61.3952 * 1000  # J/mol

molarMass = 95.95 / 1000   # kg/mol
density = 10280            # kg/m^3
molarVolume = molarMass / density  # m^3/mol

Q_atom = 5e-19  # J per atom (your diffusion activation energy)

# Atomic density
atomicDensity = NA / molarVolume

# -----------------------------
# Temperature range
# -----------------------------
T = np.linspace(2600, 2900, 500)  # K

rates = []

for temp in T:
    # Driving force (J/m^3)
    dGv = (drivingForceSlope * temp + drivingForceIntercept) / molarVolume

    # No nucleation above equilibrium
    if dGv >= 0:
        rates.append(0)
        continue

    # Critical barrier
    dGstar = (16 * np.pi * gamma**3) / (3 * dGv**2)

    # Attempt frequency
    nu = (k * temp) / h

    # Diffusion term (converted to molar form)
    diffTerm = np.exp(-(Q_atom * NA) / (R * temp))

    # Nucleation rate (1/m^3/s)
    rate = atomicDensity * nu * diffTerm * np.exp(-dGstar / (k * temp))

    rates.append(rate)

rates = np.array(rates)

# -----------------------------
# Plot
# -----------------------------
plt.figure()
plt.semilogy(T, rates)
plt.xlabel("Temperature (K)")
plt.ylabel("Nucleation Rate (1/m³·s)")
plt.title("Homogeneous Nucleation Rate vs Temperature")
plt.grid(True)
plt.show()