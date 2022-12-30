# abaqus-greenStrain
Abaqus script to compute Green strain from Nominal strain

## Overview

Green's strain is sometimes used in large-rotation, small-strain formulations
for such problems as shell buckling. It is not available directly as an Abaqus
field output nor is the displacement gradient which is commonly used in its
calculation.

This script computes and stores Green strain (GE) field output for all odb step
frames which include Nominal strain (NE).

### Reference

1. [Continuum Mechanics](https://www.continuummechanics.org/greenstrain.html)
1. [Abaqus Theory](https://help.3ds.com/2023/english/dssimulia_established/SIMACAETHERefMap/simathe-c-strainmeas.htm)

## Download

Get the [latest release](https://github.com/costerwi/abaqus-greenStrain/releases/latest)

## Execution

From the command line

```
abaqus python greenStrain.py Job-1.odb [Job-2.odb ...]
```

### Optional automatic execution

Copy the `onJobCompletion()` method into a local `abaqus_v6.env` file for
automatic execution whenever an Abaqus job completes. Make sure this
script is in the working directory or the `PYTHONPATH`

## Derivation

Green strain is calculated based on the relationships

$${\bf GE} = {1 \over 2} \left( {\bf F} \cdot {\bf F}^T - {\bf I} \right)$$

$${\bf F} \cdot {\bf F}^T - {\bf I} = {\bf V} \cdot {\bf V}^T - {\bf I} = {\bf V} \cdot {\bf V} - {\bf I} $$

$${\bf V} = {\bf NE} + {\bf I}$$

Where
- ${\bf GE}$ is the Green strain tensor in current configuration
- ${\bf F}$ is the displacement gradient tensor
- ${\bf V}$ is the stretch tensor in current configuration
- ${\bf NE}$ is the Nominal strain tensor in current configuration
