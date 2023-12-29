# Examples of how to use APT

The following examples are provided for the interested user. These examples cover the cases discussed in this reference (https://arxiv.org/abs/2307.04843).

Input decks in which there is an analytic pulse are named with `_apt` and the comparison cases in which the pulse is launched from the simulation boundary and only updated using the internal field solver are named with `_internal`. Figure numbers refer to the above reference. 

Caution: some of these examples can be fairly computationally expensive.

## Benchmarking examples

### Vacuum implementation

Input decks: `vacuum_plane1D`

Reference figure: Figure 2

Precomplier flag: APT_VACUUM

About: used to assess the difference in the total fields between simulations with an analytic pulse and one introduced from the simulation boundary. Simulation consists of either vacuum or a plasma slab with vacuum on either side (reflection case).

Additional details: the provided input deck is configured for the reflection case. The vacuum case can be obtained by commenting out the species blocks for electrons and ions. The only parameter in the input deck necessary to change for the scaling study is `dxi` (the number of cells per laser wavelength). The data points in Fig. 2 are given by $\langle E_A-E_I \rangle$, where $E_A$ and $E_I$ are the total electric field in the laser polarization direction for APT and the case using only the internal field solver, and $\langle E (x) \rangle \equiv (2 \int_{x-\lambda/2}^{x+\lambda/2} \mathrm{d}x' E^2)^{1/2}$ defines the moving cycle-average electric field amplitude. 

### Plasma (analytic current) implementation

Input decks: `plasma_plane1D` and `plasma_ramp1D`

Reference figure: Figure 8

Precomplier flag: APT_PLASMA

About: used to assess the difference in the total fields between simulations with an analytic pulse (including analytic current) and one introduced from the simulation boundary. The APT simulation consists of either a uniform plasma (`plane1D`) or one with a short density ramp down to vacuum (`ramp1D`).

Additional details: in addition to resolution (`dxi`), the relevant parameters to change for the scaling study are `nppc`, (the number of macroparticles per cell), `den_max` (the density, normalized to the critical density), and the small parameters in the phase (`ph0`) and amplitude (`den_a`, ramp case only) that compensate for changes in the initial particle placement with changing resolution. These parameters may need to be adjusted by the user for each case, due to an initialization bug affecting the random particle noise. The analysis is similar to the preceding test.

### Scaling (APT_VACUUM)

Input decks: `vacuum_plane3D`

Reference figure: Figure 4

Precomplier flag: APT_VACUUM

About: used to assess the scaling of the code with a varying number of MPI tasks.

## Examples for applications of APT

### Approximate solutions to Maxwell's equations (paraxial Gaussian beam)

Precomplier flag: APT_VACUUM_GAUSS_2D

#### Example 1: Laser wakefield

Input decks: `gauss2D_lwfa`

Reference figure: Figure 5a

About: tests the validity of the paraxial approximation leading to a Gaussian solution to the wave equation, under conditions relevant to laser wakefield acceleration

#### Example 2: Radiation pressure acceleration

Input decks: `gauss2D_rpa`

Reference figure: Figure 5b/c

About: tests the validity of the paraxial approximation leading to a Gaussian solution to the wave equation, under conditions relevant to radiation pressure acceleration

Additional details: the spot size is controlled by the parameter `w0`.

### Reduced dimensionality simulation

Precomplier flag: APT_VACUUM_GAUSS

Input decks: `compton`

Reference figure: Figure 6

About: demonstrates that a series of 1D simulations using APT can reproduce results of an expensive, fully 3D simulation, in the case of the inverse Compton scattering of high energy electrons from a counter-propagating laser pulse

Additional details: the total spectra were obtained by discrete integration over the results from 1D simulations at a number of different evaluation locations, defined by `y0` and `z0` in the input deck. In this case, the photon emission is sufficiently azimuthally symmetric that purely radial integration can be performed (i.e., performing the radial scan over `y0` > 0 with fixed `z0` = 0).

## Testing stability and accuracy

### Numerical phase velocity (APT_VACUUM)

Precomplier flag: APT_VACUUM

Input decks: `vacuum_gradient1D`

Reference figure: Figure 3

About: tests whether numerical dispersion must be included in analytic pulses for APT to be stable and accurate

Additional details: requires modification of the source code to replace the assignment of $v_\phi$ and $B_0$ with the zeroth order $v_\phi=c$ and $B_0=E_0/c$ (these are found in the subroutine `setup_analytic_pulses`, in `analytic_pulse.F90`).

### Considerations for APT_PLASMA

Input decks: `plasma_heat1D`

Reference figure: Figure 9

About: tests the required accuracy for the analytic pulse and analytic current, including the effect of particle shape on numerical dispersion, and finite duration corrections

Additional details: the various cases in Fig. 9 can be obtained by modifying the source code. The "wrong $\zeta$" case uses `s0 = 0` (in the subroutine `setup_analytic_pulses`, in `analytic_pulse.F90`). The case without the finite duration correction ("no FD") can be reproduced by setting `phdt`, `adt1`, and `adt2` to zero everywhere they appear in `analytic_pulse.F90`.
