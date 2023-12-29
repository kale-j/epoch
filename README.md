# EPOCH with implementation of the analytic pulse technique (APT)

This branch of EPOCH contains a separation of Maxwell's equations into analytic and computational parts, henceforth referred to as the analytic pulse technique (APT). Details of APT and examples of its applications are given in this reference: https://arxiv.org/abs/2307.04843

The following documentation is intended to serve as a supplement to the reference. It is not intended to be sufficient on its own to understand what this code does.

## About APT

The linearity of Maxwell's equations in the current and field terms allows the fields and current to be decomposed into multiple parts. One can use this linearity to treat part of the field/current analytically and part using standard computational electromagnetics approaches (e.g., FDTD). The combination of analytic and computational solutions can be used to extend the capabilities of existing methods for simulating electromagnetics. This can be highly beneficial. In the PIC context, it can allow the user to study new configurations, can be used to bypass certain numerical issues, and can sometimes significantly reduce simulation cost. Several examples of this in practice are given in the reference.

In EPOCH it is convenient to consider either solutions in which only the fields have an analytic part or solutions in which both the fields and the current have an analytic part. These are referred to as APT_VACUUM (and its variants, see precompiler flag section) and APT_PLASMA, respectively.

## Examples

Examples of how to use this code for both benchmarking and physics problems are provided in the `analytic_pulse_usage` directory. This directory also contains information on how to reproduce the results included in the reference.

## Implementation notes

This implementation of APT is based on evaluating the analytic component of the electromagnetic fields and current on the grid. This is significantly faster than evaluation at the particle position when there many particles in the simulation.

The current implementation is geared towards modeling electromagnetic pulses (like laser pulses). The pulses propagate in the +x-direction and are polarized in y.

### Precompiler flags

A few different configurations for the analytic pulse/current are available, using the following precompiler flags:

```
APT_VACUUM

Plane wave implementation of APT without analytic current
```

```
APT_VACUUM_GAUSS

Spatially Gaussian/temporally super-Gaussian pulse without analytic current. Includes support for moving focal location (flying focus)
```


```
APT_VACUUM_GAUSS_2D

Spatially Gaussian/temporally super-Gaussian pulse without analytic current. This version uses 2D focusing (in y). Includes support for moving focal location (flying focus)
```


```
APT_PLASMA

Plane wave implementation of APT with analytic current. At the moment, this only works in 1D (see known issues).
```


### Outputs

New outputs are introduced to represent the total fields (e.g., ex_total) and the particle current minus the analytic current (e.g., jx_diff, in APT_PLASMA only). These variables can be written out following the same approach as is ordinarily used.

NOTE: with APT, the "ordinary" fields (e.g., ex) represents the fields advanced internally by the FDTD field update, i.e., ex is equal to ex_total minus the analytically specified part.

## General disclaimer and known issues

This code serves as an example for how APT can be implemented in a PIC code. An attempt has been made to ensure the code will work under several common use cases (for example, those in the reference), however, it has not been tested in all possible configurations.

You may encounter numerical instabilities or other incorrect behavior if the conditions under which this implementation is valid are violated. Users are advised to perform thorough testing before using this code for research.

### Known issues

Initialization bug: APT_PLASMA currently only works in 1D, and only if load balancing does not occur before the first step. It is advised to use_pre_balance=F and balance_first=F. Attempting to use APT_PLASMA in higher dimensions or with initial load balancing is likely to cause out of bounds memory access.

Particle shape: the phase velocity used in APT_PLASMA depends on the particle shape and is currently only correct for the triangle particle shape. In testing, using the "wrong" particle shape contribution to the phase velocity appears fairly stable under conditions where APT_PLASMA is fairly stable, but can increase the growth rate of the numerical instability when APT_PLASMA is marginally stable or unstable. Use other particle shapes with caution.

## Acknowledgements

In research derived from this work, please cite: https://arxiv.org/abs/2307.04843