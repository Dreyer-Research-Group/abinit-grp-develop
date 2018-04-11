---
description: To to calculate the temperature dependence of the electronic structure
authors: SP
---

This page gives hints on how to calculate the temperature dependence of the electronic structure with the ABINIT package.

## Introduction

The electronic structure changes with temperature. In most materials, such
changes are mainly driven by the electron-phonon interaction, which is also
present at zero Kelvin, inducing the so-called zero-point motion
renormalization (ZPR) of the eigenvalues. These effects can be computed thanks
to the Allen-Heine-Cardona (AHC) theory [[cite:Allen1976]],
[[cite:Allen1981]], [[cite:Allen1983]], which is based on diagrammatic method
of many-body perturbation theory. An extension to the standard AHC theory also
gives access to the electronic lifetime and decay rates. These physical
properties are available from ABINIT since v7.10.4.

The AHC formalism and the implemented equations can be found in
[[cite:Ponce2014a]]. An extended verification and validation study (also
versus other first-principle codes) of the ABINIT implementation can be found
in [[cite:Ponce2014]]. The AHC implementation can be used with any XC
functional working with the response-function (RF) part of the code, and
requires the use of norm-conserving pseudopotentials. NetCDF support is
mandatory.

The AHC implementation in ABINIT is built on a Sternheimer approach to
efficiently compute the sum over highly energetic bands appearing in the AHC
equations [[cite:Gonze2011]]. Such behavior is controlled by the input
variable [[ieig2rf]].

The **k** -point convergence can be strongly improved by restoring the charge
neutrality through the reading of the Born effective charge and dielectric
tensor (controlled by the input variable [[getddb]]). More information on the
importance of charge neutrality fulfillment can be found in
[[cite:Ponce2015]]. The value of [[elph2_imagden]] sets the imaginary shifts
used to smooth numerical instabilities in the denominator of the sum-over-
states expression.

We have checked that the implementation correctly holds for arbitrarily small
[[elph2_imagden]] parameters, [[cite:Ponce2015]]. The input variable
[[smdelta]] triggers the calculation of the electronic lifetime and the value
of the smearing delta function can be specified through [[esmear]].

A double grid can be used to speed-up the calculations with [[getwfkfine]] or
[[irdwfkfine]]. The variable [[getgam_eig2nkq]] gives the contribution at Γ so
that the Debye-Waller term can be computed. This variable is only relevant for
calculations of AHC using the abinit program only. It is nonetheless
recommended to use the provided python post-processing script
temperature_para.py with its module rf_mods.py in the directory
scripts/post_processing/ to allow for more flexibility. The python scripts
support multi-threading.

The following steps are required to perform an AHC calculation:

* Perform a response function calculation at **q** =Γ with electric field perturbation.
* Perform phonon calculations and produce the EPC for a large set of wavevectors **q** , reading the Born effective charge and dielectric tensor with [[getddb]].
* Gather and compute the impact of the electron-phonon coupling on the electronic eigenenergies using the temperature_para.py python script.

The outputs of the script are provided in text and NetCDF format to allow for
later reading inside ABINIT. This could be used in the future developments of
ABINIT to compute temperature-dependent optical properties for example.

For the temperature dependence of the Fermi energy, see [[topic:ElPhonTransport]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* A lesson has been developed on [[lesson:tdepes|the temperature dependence of the electronic structure]]:.

