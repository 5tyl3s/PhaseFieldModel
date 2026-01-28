# Inoculation Phase Field Model - Complete Equation Reference

This document contains all mathematical equations used in the phase field simulation model.

---

## 1. Temperature Evolution

### Temperature Field
The temperature at each node evolves based on an initial vertical gradient and linear cooling:

$$T(x, y, t) = T_{\text{start}} + y \cdot G_t \cdot \Delta x - G_t \cdot C_{\text{rate}} \cdot \Delta t \cdot t$$

where:
- $T_{\text{start}}$ = starting temperature (K)
- $y$ = vertical position (lattice coordinates)
- $G_t$ = temperature gradient (K/m)
- $\Delta x$ = grid spacing (m)
- $C_{\text{rate}}$ = cooling rate (K/s)
- $\Delta t$ = timestep (s)
- $t$ = current timestep number

---

## 2. Phase Field Evolution

### Phase Field Update
The phase field $\phi$ evolves via the Allen-Cahn equation:

$$\frac{\partial \phi}{\partial t} = -M_\phi \frac{\delta F}{\delta \phi}$$

where the driving force is:

$$\frac{\delta F}{\delta \phi} = \frac{dE_{\text{local}}}{d\phi} + \nabla^2 E_{\text{grad}}$$

Discrete form:
$$\phi^{n+1} = \phi^n - \frac{\Delta t}{(\Delta x)^2} \times \text{diffEnergy}_\phi \times 3 \times 10^{-16}$$

### Phase Local Energy Density
The local phase energy includes contributions from undercooling, driving forces, grain interaction, and particle interaction:

$$E_{\text{local}} = (-2 + 2\phi) \cdot u_c \cdot \Delta G_{\text{sol}} + 2\phi(1 - u_c) \Delta G_{\text{liq}}$$

$$+ A_\phi \left[2(\phi - 1) \sum_g \psi_g + 2\phi \left(1 - \sum_g \psi_g\right) + (2\phi - 2) \sum_g \psi_g \right]$$

$$+ \phi \cdot c_p^2 \cdot \gamma_s \cdot f_s + (2\phi - 2) \cdot c_p^2 \cdot \gamma_l \cdot f_l$$

where:
- $u_c$ = undercooling fraction (dimensionless, 0 to 1)
- $\Delta G_{\text{sol}}$ = driving force to solidification (J/m³)
- $\Delta G_{\text{liq}}$ = driving force to liquification (J/m³)
- $A_\phi$ = phase coefficient
- $\psi_g$ = grain phase field
- $c_p$ = particle composition
- $\gamma_s, \gamma_l$ = particle-solid and particle-liquid interfacial energies
- $f_s, f_l$ = surface-to-volume ratio of particles

### Phase Gradient Energy
Using 4-point Laplacian (von Neumann stencil):

$$E_{\text{grad}} = K_\phi \nabla^2 \phi \times 10$$

$$\nabla^2 \phi = \sum_{\text{neighbors}} \phi_{\text{nbr}} - 4 \cdot \phi_{\text{center}}$$

---

## 3. Grain Field Evolution

### Grain Phase Field Update
Each grain field $\psi_g$ evolves independently:

$$\frac{\partial \psi_g}{\partial t} = -M_g \frac{\delta F}{\delta \psi_g}$$

Discrete form:
$$\psi_g^{n+1} = \psi_g^n - \frac{\Delta t}{(\Delta x)^2} \times \text{diffEnergy}_g \times 5 \times 10^{-15}$$

### Grain Local Energy Density
The grain energy includes growth driving force and competition with other grains:

$$E_{\text{grain}} = C_{\text{pre}} \left[\left(\psi_g^3 - \psi_g\right) \cdot u_c + \psi_g \left(\sum_{j \neq g} \psi_j\right)^2 \cdot \kappa_g \right]$$

where:
- $C_{\text{pre}}$ = grain energy pre-coefficient
- $\kappa_g$ = grain interaction width
- $u_c$ = undercooling fraction

### Grain Gradient Energy
The gradient term includes grain boundary energy calculation:

$$E_{\text{grad}} = -K_g \cdot \gamma_{GB}(\mathbf{n}) \cdot \nabla^2 \psi_g / (\Delta x)^2$$

where $\gamma_{GB}$ is the orientation-dependent grain boundary energy computed using spherical harmonics (see Section 6).

---

## 4. Driving Forces

### Undercooling Function
The undercooling is a smooth transition from liquid (above melt) to solid (below melt):

$$u_c(T) = 0.5 \left(1 - \tanh\left(10^6 \left(\frac{T}{T_m} - 1\right)\right)\right)$$

where:
- $T$ = current temperature (K)
- $T_m$ = melt temperature (K)

### Driving Force to Solidification
Linear fit to thermodynamic data:

$$\Delta G_{\text{sol}} = 1000 \times \frac{-0.0212 \cdot T + 61.3952}{\bar{V}}$$

where $\bar{V}$ = molar volume (m³/mol)

### Driving Force to Liquification
$$\Delta G_{\text{liq}} = 1000 \times \frac{0.0349 \cdot T - 101.0704}{\bar{V}}$$

---

## 5. Particle Composition Evolution

### Particle Composition Update
Particles move via Cahn-Hilliard dynamics:

$$\frac{\partial c_p}{\partial t} = \nabla \cdot \left(M(c_p, \phi) \nabla \mu\right)$$

Discrete form:
$$c_p^{n+1} = c_p^n + \Delta t \cdot \text{flux}$$

### Chemical Potential
The variational chemical potential drives particle motion:

$$\mu = \frac{\partial f}{\partial c_p} = 2c_p \phi^2 \gamma_s f_s + 2c_p(1-\phi)^2 \gamma_l f_l$$

where:
- $\gamma_s, \gamma_l$ = particle interfacial energies
- $f_s = 0.66 \times 4\pi r_p^2 / \left(\frac{4}{3}\pi r_p^3\right)$ = surface-to-volume ratio

### Particle Mobility
The mobility is phase-dependent (higher in liquid):

$$M(\phi) = M_{\text{liquid}} (1 - \phi) + M_{\text{solid}} \phi$$

where:
- $M_{\text{liquid}} = 10^{-12}$ m²/J·s (typical)
- $M_{\text{solid}} = 10^{-24}$ m²/J·s (typical)

### Particle Flux and Laplacian
$$\text{flux} = M(\phi) \cdot \nabla^2 \mu / (\Delta x)^2$$

$$\nabla^2 \mu = \sum_{\text{4-neighbors}} \mu_{\text{nbr}} - 4 \mu_{\text{center}}$$

---

## 6. Grain Boundary Energy (Spherical Harmonics)

### Surface Normal Calculation
From the grain phase gradient $\nabla \psi_g$:

$$\mathbf{n} = \frac{(-\partial_x \psi_g, -\partial_y \psi_g, -\partial_z \psi_g)}{|\nabla \psi_g|}$$

### Gradient Calculation (2D)
Horizontal gradient:
$$g_x = \left(\psi_g^{\text{right}} - \psi_g^{\text{left}}\right) \cdot \frac{1}{2\Delta x} + 0.5 \left(\psi_g^{\text{UR}} - \psi_g^{\text{UL}} + \psi_g^{\text{DR}} - \psi_g^{\text{DL}}\right) \cdot \frac{1}{2\Delta x}$$

Vertical gradient:
$$g_y = \left(\psi_g^{\text{down}} - \psi_g^{\text{up}}\right) \cdot \frac{1}{2\Delta x} + 0.5 \left(\psi_g^{\text{DR}} - \psi_g^{\text{UR}} + \psi_g^{\text{DL}} - \psi_g^{\text{UL}}\right) \cdot \frac{1}{2\Delta x}$$

### Euler Angle Rotation
The gradient normal is rotated by the grain orientation (Bunge convention):

$$\mathbf{R} = R_z(\phi_2) \cdot R_x(\Phi) \cdot R_z(\phi_1)$$

where $(\phi_1, \Phi, \phi_2)$ are the Euler angles. Element-by-element:

$$R_{00} = \cos\phi_1 \cos\phi_2 - \sin\phi_1 \sin\phi_2 \cos\Phi$$
$$R_{01} = -\cos\phi_1 \sin\phi_2 - \sin\phi_1 \cos\phi_2 \cos\Phi$$
$$R_{02} = \sin\phi_1 \sin\Phi$$

(and so on for other elements)

### Rotated Normal in Spherical Coordinates
Convert rotated normal to $(r, \theta, \phi)$:

$$\theta = \arccos\left(\frac{z}{r}\right), \quad \phi = \arctan2(y, x)$$

### Grain Boundary Energy via Spherical Harmonics
$$\gamma_{GB}(\theta, \phi) = \sum_{l=0}^{L_{\max}} \sum_{m=-l}^{l} c_{lm} Y_{lm}^{\text{real}}(\theta, \phi)$$

where $Y_{lm}^{\text{real}}$ are real spherical harmonics and $c_{lm}$ are coefficients read from a data file.

### Real Spherical Harmonics
$$Y_{lm}^{\text{real}}(\theta, \phi) = \begin{cases}
N_l \cdot P_l^0(\cos\theta) & \text{if } m = 0 \\
\sqrt{2} N_l P_l^m(\cos\theta) \cos(m\phi) & \text{if } m > 0 \\
\sqrt{2} N_l P_l^{|m|}(\cos\theta) \sin(|m|\phi) & \text{if } m < 0
\end{cases}$$

where:
- $N_l = \sqrt{\frac{2l+1}{4\pi} \frac{(l-|m|)!}{(l+|m|)!}}$ = normalization factor
- $P_l^m(x)$ = associated Legendre polynomial

---

## 7. Nucleation

### Homogeneous Nucleation Probability
Uses classical nucleation theory:

$$P_{\text{hom}} = 1 - \exp\left(-(\Delta x)^3 \cdot I_{\text{hom}} \cdot \Delta t\right)$$

where:
$$I_{\text{hom}} = \frac{k_B T}{h} N_{\text{atoms}} \exp\left(-\frac{Q_d}{k_B T}\right) \exp\left(-\frac{\Delta G^*}{k_B T}\right)$$

### Critical Nucleus Free Energy (3D Spherical)
$$\Delta G^* = \frac{16\pi \gamma^3}{3 (\Delta G_{\text{driving}})^2}$$

where:
- $\gamma = 0.05$ J/m² = surface/interfacial energy
- $\Delta G_{\text{driving}}$ = volumetric driving force

### Nucleation Rate Components
$$N_{\text{atoms}} = \rho_A \times V_{\text{cell}}$$

where $\rho_A$ = atomic number density = $N_A / \bar{V}$

$$\text{DiffusionTerm} = \exp\left(-\frac{Q_d}{R T}\right)$$

where:
- $Q_d$ = diffusion activation energy
- $R = 8.314$ J/(mol·K) = gas constant

Constants:
- $k_B = 1.380649 \times 10^{-23}$ J/K = Boltzmann constant
- $h = 6.62607015 \times 10^{-34}$ J·s = Planck constant
- $N_A = 6.02214076 \times 10^{23}$ mol⁻¹ = Avogadro number

### Heterogeneous Nucleation
Occurs when undercooling exceeds a threshold:

$$u_c < 0 \quad \text{and} \quad P_{\text{het}} = \Delta t \cdot c_p$$

where the threshold is computed from:
$$u_c = \frac{1000}{\bar{V}} \left(-0.0212 T + 61.3952\right) + \Delta T_{\text{het,uc}}$$

---

## 8. Grain Propagation & Activation

### Grain Addition to Pending List
When a grain reaches sufficient phase (> 0.01) in a neighbor, it's marked for addition:

$$\psi_g^{\text{neighbor}} > 0.01 \implies \text{add grain to propagation list}$$

The grain with the **highest phase** from each neighbor is selected to avoid duplicate activations.

### Orientation Storage During Propagation
When grain $g$ is added to the pending list, its orientation is also stored:

$$\mathbf{O}_g^{\text{stored}} = \mathbf{O}_g^{\text{neighbor}}$$

### Grain Activation
At the end of the timestep, pending grains are activated:

$$\psi_g^{\text{new}} \gets 0 \quad \text{(initialized but will evolve)}$$
$$\mathbf{O}_g^{\text{new}} \gets \mathbf{O}_g^{\text{stored}}$$

---

## 9. Numerical Parameters

### Discretization
- Grid spacing: $\Delta x = 1 \text{ μm}$ (typical)
- Timestep: $\Delta t = 10^{-8}$ to $10^{-12}$ s (typical)

### Coefficients (Typical Values)
- Phase gradient coefficient: $K_\phi \approx 1$ (tunable)
- Grain gradient coefficient: $K_g \approx 1$ (tunable)
- Phase coefficient: $A_\phi \approx 0.25$
- Grain coefficient: $A_g \approx 0.125$
- Grain energy pre-coefficient: $C_{\text{pre}} \approx \text{tunable}$

### Clamping & Safety
All computed diffusion energies are clamped to $[-10^{12}, 10^{12}]$ J/m³ to prevent numerical overflow.

---

## 10. Euler Angle Convention

The model uses **Bunge Euler angles** with the convention:

- $\phi_1$ = first rotation about Z-axis (0 to 2π)
- $\Phi$ = rotation about new X-axis (0 to π)
- $\phi_2$ = second rotation about new Z-axis (0 to 2π)

Rotation matrix composition:
$$\mathbf{R} = R_z(\phi_2) \circ R_x(\Phi) \circ R_z(\phi_1)$$

---

## 11. Summary of Field Variables

| Variable | Symbol | Units | Range |
|----------|--------|-------|-------|
| Phase field | $\phi$ | - | [0, 1] |
| Grain phase | $\psi_g$ | - | [0, 1] |
| Temperature | $T$ | K | [min, max] |
| Particle composition | $c_p$ | - | [0, 1] |
| Chemical potential | $\mu$ | J/m³ | ℝ |
| Undercooling | $u_c$ | - | [0, 1] |
| Euler angles | $(\phi_1, \Phi, \phi_2)$ | rad | [0, 2π] or [0, π] |

---

## References

- **Phase field modeling**: Wheeler et al., Phys. Rev. E 47, 1893 (1993)
- **Grain boundary energy**: Read & Shockley, Phys. Rev. 78, 275 (1950)
- **Classical nucleation theory**: Kelton, Solid State Physics 45, 75 (1991)
- **Spherical harmonics**: Varshalovich, Moskalev & Khersonskii, Quantum Theory of Angular Momentum (1988)
