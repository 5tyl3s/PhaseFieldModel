# Inoculation Phase Field Model - Complete Equation Reference

This document contains all mathematical equations used in the current phase field simulation model, including the updated nucleation mechanics with particle redistribution.

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

Discrete form (with explicit tuning coefficient $3 \times 10^{-16}$):
$$\phi^{n+1} = \phi^n - \frac{\Delta t}{(\Delta x)^2} \times \text{diffEnergy}_\phi \times 3 \times 10^{-16}$$

**Note**: The coefficient $3 \times 10^{-16}$ is a tuning parameter to control phase field evolution speed.

### Phase Local Energy Density
The local phase energy includes contributions from undercooling, grain interaction, and particle interaction:

$$E_{\text{local}} = (-2 + 2\phi) \cdot u_c \cdot \Delta G_{\text{sol}} + 2\phi(1 - u_c) \Delta G_{\text{liq}}$$

$$+ A_\phi \left[2(\phi - 1) \sum_g \psi_g + 2\phi \left(1 - \sum_g \psi_g\right) + (2\phi - 2) \sum_g \psi_g \right]$$

$$+ \phi \cdot c_p^2 \cdot \gamma_s \cdot f_s + (2\phi - 2) \cdot c_p^2 \cdot \gamma_l \cdot f_l$$

where:
- $u_c$ = undercooling fraction (dimensionless, 0 to 1)
- $\Delta G_{\text{sol}}$ = driving force to solidification (J/m³)
- $\Delta G_{\text{liq}}$ = driving force to liquification (J/m³)
- $A_\phi$ = phase coefficient (computed from interfacial energy and barrier height)
- $\psi_g$ = grain phase field
- $c_p$ = particle composition
- $\gamma_s, \gamma_l$ = particle-solid and particle-liquid interfacial energies
- $f_s, f_l$ = surface-to-volume ratio of particles

### Phase Gradient Energy
Using 4-point Laplacian (von Neumann stencil):

$$E_{\text{grad}} = K_\phi \nabla^2 \phi \times 10$$

$$\nabla^2 \phi = \sum_{\text{neighbors}} \phi_{\text{nbr}} - 4 \cdot \phi_{\text{center}}$$

**Note**: The coefficient $10$ is a tuning multiplier for gradient energy strength.

---

## 3. Grain Field Evolution

### Grain Phase Field Update
Each grain field $\psi_g$ evolves independently:

$$\frac{\partial \psi_g}{\partial t} = -M_g \frac{\delta F}{\delta \psi_g}$$

Discrete form (with explicit tuning coefficient $5 \times 10^{-15}$):
$$\psi_g^{n+1} = \psi_g^n - \frac{\Delta t}{(\Delta x)^2} \times \text{diffEnergy}_g \times 5 \times 10^{-15}$$

**Note**: The coefficient $5 \times 10^{-15}$ is a tuning parameter to control grain field evolution speed.

### Grain Local Energy Density
The grain energy includes growth driving force and competition with other grains:

$$E_{\text{grain}} = C_{\text{pre}} \left[\left(\psi_g^3 - \psi_g\right) \cdot u_c + \psi_g \left(\sum_{j \neq g} \psi_j\right)^2 \cdot \kappa_g \right]$$

where:
- $C_{\text{pre}}$ = grain energy pre-coefficient (computed from grain interfacial energy)
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

$$\Delta G_{\text{sol}} = 1000 \times \frac{m_{\text{slope}} \cdot T + b_{\text{intercept}}}{\bar{V}}$$

where:
- $m_{\text{slope}} = -0.0212$ K⁻¹ (tuning parameter)
- $b_{\text{intercept}} = 61.3952$ (tuning parameter)
- $\bar{V}$ = molar volume (m³/mol)

### Driving Force to Liquification
$$\Delta G_{\text{liq}} = 1000 \times \frac{-m_{\text{slope}} \cdot T - b_{\text{intercept}}}{\bar{V}}$$

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
- $M_{\text{liquid}} = 1.0 \times 10^{-15}$ m²/J·s (tuning parameter)
- $M_{\text{solid}} = 1.0 \times 10^{-24}$ m²/J·s (tuning parameter)

### Particle Flux and Laplacian
$$\text{flux} = -M(\phi) \cdot \nabla^2 \mu / (\Delta x)^2$$

$$\nabla^2 \mu = \sum_{\text{4-neighbors}} \mu_{\text{nbr}} - 4 \mu_{\text{center}}$$

**Boundary Condition**: No-flux boundary at top and bottom (neighbors with `exists == 0` are skipped in Laplacian calculation).

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

where $(\phi_1, \Phi, \phi_2)$ are the Euler angles.

### Rotated Normal in Spherical Coordinates
Convert rotated normal to $(r, \theta, \phi)$:

$$\theta = \arccos\left(\frac{z}{r}\right), \quad \phi = \arctan2(y, x)$$

### Grain Boundary Energy via Spherical Harmonics
$$\gamma_{GB}(\theta, \phi) = \sum_{l=0}^{L_{\max}} \sum_{m=-l}^{l} c_{lm} Y_{lm}^{\text{real}}(\theta, \phi)$$

where:
- $Y_{lm}^{\text{real}}$ = real spherical harmonics
- $c_{lm}$ = coefficients read from data file
- $L_{\max} = 8$ (maximum harmonic degree)

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

### Homogeneous Nucleation
**Probability**:
$$P_{\text{hom}} = 1 - \exp\left(-(\Delta x)^3 \cdot I_{\text{hom}} \cdot \Delta t\right)$$

where $I_{\text{hom}}$ is computed via classical nucleation theory.

**Particle Redistribution** (NEW): When homogeneous nucleation occurs at node $n$:
1. Set nucleation node: $c_p(n) \gets 0$
2. Calculate total liquid fraction across grid: $L_{\text{total}} = \sum_{\text{all nodes}} (1 - \phi)$
3. Distribute particles proportionally to liquid fraction:
$$c_p(i) \gets c_p(i) + \frac{(1 - \phi(i))}{L_{\text{total}}}$$
4. Clamp all values to [0, 1]

### Heterogeneous Nucleation
**Condition**: Occurs when undercooling exceeds threshold AND particles available:

$$\Delta G_{\text{het}} = 1000 \times \frac{m_{\text{slope}} \cdot T + b_{\text{intercept}}}{\bar{V}} > \Delta T_{\text{het,uc}}$$

AND

$$P_{\text{het}} = c_p \cdot \Delta t \quad (\text{particle consumption rate})$$

**Particle Redistribution** (NEW): When heterogeneous nucleation occurs at node $n$:
1. Set nucleation node: $c_p(n) \gets 1$ (maximum particle absorption)
2. Initialize total loss: $L = 0$
3. Expand radially through neighboring nodes until $L \geq 1$:
   - For each neighbor at distance $d$:
   $$\text{Remove} = (1 - \phi_{\text{neighbor}}) \times c_p(\text{neighbor})$$
   $$L \gets L + \text{Remove}$$
   $$c_p(\text{neighbor}) \gets c_p(\text{neighbor}) - \text{Remove}$$
4. Stop when $L \geq 1$ or all nodes processed
5. Clamp all values to [0, 1]

---

## 8. Grain Propagation & Activation

### Grain Addition to Pending List
When a grain reaches sufficient phase (> 0.01) in a neighbor, it's marked for addition:

$$\psi_g^{\text{neighbor}} > 0.01 \implies \text{add grain to propagation list}$$

The grain with the **highest phase** from each neighbor is selected.

### Orientation Storage During Propagation
When grain $g$ is added to the pending list, its orientation is stored:

$$\mathbf{O}_g^{\text{stored}} = \mathbf{O}_g^{\text{neighbor}}$$

### Grain Activation
At the end of timestep, pending grains are activated and grain boundary energy is computed based on stored orientation.

---

## 9. Numerical Parameters & Coefficients

### Discretization Parameters

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Grid spacing | $\Delta x$ | 1 | μm |
| Timestep | $\Delta t$ | 1×10⁻⁸ | s |

### Material & Thermodynamic Coefficients

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Melt temperature | $T_m$ | 2896 | K | |
| Starting temperature | $T_{\text{start}}$ | 2891 | K | |
| Minimum temperature | $T_{\text{min}}$ | 100 | K | |
| Cooling rate | $C_{\text{rate}}$ | 0.3 | K/s | |
| Temperature gradient | $G_t$ | 100000 | K/m | |
| Specific heat capacity | $c_p$ | 251 | J/(kg·K) | |
| Density | $\rho$ | 10280 | kg/m³ | |
| Molar mass | $M$ | 95.95 | g/mol | |
| Molar volume | $\bar{V}$ | $M / (1000 \rho)$ | m³/mol | Computed |

### Interfacial Energy Coefficients

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Liquid-solid interfacial energy | $\gamma_{\text{LS}}$ | 2.92 | J/m² | |
| Liquid-solid interfacial width | $w_{\text{LS}}$ | 1 | μm | Converted to m in code |
| Grain boundary interfacial energy | $\gamma_{GB}$ | 2.92 | J/m² | |
| Grain boundary interfacial width | $w_{GB}$ | 1 | μm | Converted to m in code |
| Particle-solid interfacial energy | $\gamma_{p,s}$ | 2.5 | J/m² | Tuning parameter |
| Particle-liquid interfacial energy | $\gamma_{p,l}$ | 1.95 | J/m² | Tuning parameter |

### Barrier Height Coefficients (Tuning Parameters)

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Phase field barrier height | $h_\phi$ | 0.25 | Controls phase transition steepness |
| Grain field barrier height | $h_g$ | 0.125 | Controls grain transition steepness |

### Computed Phase Field Coefficients

| Parameter | Symbol | Computed Formula | Value | Units |
|-----------|--------|------------------|-------|-------|
| Phase coefficient | $A_\phi$ | $0.75 \times \gamma_{\text{LS}} / (w_{\text{LS}} \times h_\phi)$ | ≈ 9.76 | J/m² |
| Phase pre-coefficient | $C_{\phi,\text{pre}}$ | $0.75 \times \gamma_{\text{LS}} / (h_\phi \times w_{\text{LS}})$ | ≈ 9.76 | J/m² |
| Phase gradient coefficient | $K_\phi$ | $0.75 \times \gamma_{\text{LS}} \times w_{\text{LS}}$ | ≈ 2.19×10⁻⁶ | J/m |

### Computed Grain Field Coefficients

| Parameter | Symbol | Computed Formula | Value | Units |
|-----------|--------|------------------|-------|-------|
| Grain pre-coefficient | $C_{g,\text{pre}}$ | $0.75 \times \gamma_{GB} / (h_g \times w_{GB})$ | ≈ 19.5 | J/m² |
| Grain gradient coefficient | $K_g$ | $0.5 \times w_{GB}$ | ≈ 5×10⁻⁷ | m |

### Particle Evolution Coefficients (Tuning Parameters)

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Particle mobility (liquid) | $M_{\text{liq}}$ | 1×10⁻¹⁵ | m²/(J·s) | Tuning parameter |
| Particle mobility (solid) | $M_{\text{solid}}$ | 1×10⁻²⁴ | m²/(J·s) | Tuning parameter |
| Initial particle volume fraction | $c_p^0$ | 0.01 | - | Initial condition |
| Particle diameter | $d_p$ | 1.0 | μm | |

### Phase & Grain Field Evolution Tuning Coefficients (Critical!)

| Parameter | Value | Units | Notes |
|-----------|-------|-------|-------|
| Phase field evolution multiplier | 3×10⁻¹⁶ | - | Controls phase field speed; **TUNING PARAMETER** |
| Grain field evolution multiplier | 5×10⁻¹⁵ | - | Controls grain field speed; **TUNING PARAMETER** |
| Phase gradient energy multiplier | 10 | - | Gradient term strength; **TUNING PARAMETER** |

### Driving Force Coefficients (Linear Regression Fit)

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Driving force slope | $m$ | -0.0212 | K⁻¹ | **TUNING PARAMETER** |
| Driving force intercept | $b$ | 61.3952 | - | **TUNING PARAMETER** |
| Diffusion activation energy | $Q_d$ | 1.5 | (arbitrary model units) | **TUNING PARAMETER** |

### Nucleation Coefficients

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Homogeneous nucleation base coefficient | $C_{\text{hom}}$ | 1×10²⁰ | (model units) |
| Heterogeneous nucleation undercooling threshold | $\Delta T_{\text{het,uc}}$ | 0 | (computed from driving force) |

### Spherical Harmonics

| Parameter | Value | Notes |
|-----------|-------|-------|
| Maximum harmonic degree | $L_{\max}$ | 8 |
| Number of coefficients | $(L_{\max}+1)^2$ | 81 |

---

## 10. Summary of Tuning Parameters

The following parameters are **adjustable for model calibration**:

1. **Evolution Speed** (critical for stability):
   - Phase evolution multiplier: $3 \times 10^{-16}$
   - Grain evolution multiplier: $5 \times 10^{-15}$
   
2. **Interfacial & Mobility Properties**:
   - Particle mobility (liquid): $1 \times 10^{-15}$ m²/(J·s)
   - Particle mobility (solid): $1 \times 10^{-24}$ m²/(J·s)
   - Barrier heights: $h_\phi = 0.25$, $h_g = 0.125$

3. **Thermodynamic Driving Force**:
   - Slope: $m = -0.0212$ K⁻¹
   - Intercept: $b = 61.3952$
   - Diffusion activation energy: $Q_d = 1.5$

4. **Particle Interaction**:
   - Particle-solid interfacial energy: $2.5$ J/m²
   - Particle-liquid interfacial energy: $1.95$ J/m²

---

## 11. Field Variables Summary

| Variable | Symbol | Units | Range | Notes |
|----------|--------|-------|-------|-------|
| Phase field | $\phi$ | - | [0, 1] | 0 = liquid, 1 = solid |
| Grain phase | $\psi_g$ | - | [0, 1] | Per grain |
| Temperature | $T$ | K | [min, max] | Spatially/temporally varying |
| Particle composition | $c_p$ | - | [0, 1] | Volume fraction |
| Chemical potential | $\mu$ | J/m³ | ℝ | Variational derivative |
| Undercooling | $u_c$ | - | [0, 1] | 0 = fully liquid, 1 = fully solid |
| Liquid fraction | $(1-\phi)$ | - | [0, 1] | Used in nucleation redistribution |
| Euler angles | $(\phi_1, \Phi, \phi_2)$ | rad | [0, 2π]/[0, π] | Bunge convention |

---

## References

- **Phase field modeling**: Wheeler et al., Phys. Rev. E 47, 1893 (1993)
- **Grain boundary energy**: Read & Shockley, Phys. Rev. 78, 275 (1950)
- **Classical nucleation theory**: Kelton, Solid State Physics 45, 75 (1991)
- **Spherical harmonics**: Varshalovich, Moskalev & Khersonskii, Quantum Theory of Angular Momentum (1988)

