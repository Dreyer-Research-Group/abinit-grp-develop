---
description: How to set parameters for a parallel calculation
authors: MT,FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters for a parallel calculation with the ABINIT package.

## Introduction

* For ground-state calculations, the code has been parallelized (MPI-based parallelism) on the k-points, on the spins, on the spinor components, on the bands, and on the FFT grid and plane wave coefficients. For the k-point and spin parallelisations (using MPI), the communication load is generally very small. This allows it to be used on a cluster of workstations. However, the number of nodes that can be used in parallel might be small, and depends strongly on the physics of the problem. A combined FFT / band parallelisation (LOBPCG)is available [[cite:Bottin2008]], and has shown very large speed up (>1000) on powerful computers with a large number of processors and high-speed interconnect. The combination of FFT / band / k point and spin parallelism is also available, and quite efficient for such computers. Available for norm-conserving as well as PAW cases. Automatic determination of the best combination of parallelism levels is available. Use of MPIIO is mandatory for the largest speed ups to be observed. 
* Chebyshev filtering (Chebfi) is a method to solve the linear eigenvalue problem, and can be used as a SCF solver in Abinit. It is implemented in Abinit [[cite:Levitt2015]]. The design goal is for Chebfi to replace LOBPCG as the solver of choice for large-scale computations in Abinit. By performing less orthogonalizations and diagonalizations than LOBPCG, scaling to higher processor counts is possible. A manual to use Chebfi is available [[pdf:howto_chebfi.pdf|here]] 
* For ground-state calculations, with a set of images (e.g. nudged elastic band method, the string method, the path-integral molecular dynamics, the genetic algorithm), MPI-based parallelism is used. The work load for the different images has been distributed, and this parallelism can be combined with the parallelism described in point hereabove, leading to speed-up beyond 5000. 
* For ground-state calculations, GPU can be used. This is based on CUDA+MAGMA. 
  

* For ground-state calculations, the wavelet part of ABINIT (BigDFT) is also very well parallelized : MPI band parallelism, combined with GPUs. 
* For response calculations, the code has been parallelized (MPI-based parallelism) on k-points, spins, bands, as well as on perturbations. For the k-points, spins and bands parallelisation, the communication load is rather small also, and, unlike for the GS calculations, the number of nodes that can be used in parallel will be large, nearly independently of the physics of the problem. Parallelism on perturbations is very similar to the parallelism on images in the ground state case (so, very efficient), although the load balancing problem for perturbations with different number of k points is not adressed at present. Use of MPIIO is mandatory for the largest speed ups to be observed. 
  

* GW calculations are MPI-parallelized over k-points. They are also parallelized over transitions (valence to conduction band pairs), but the two parallelisation cannot be used currently at present. The transition parallelism has been show to allow speed ups as large as 300. 
  

* Ground state, response function, and GW parallel calculations can be done also by using OpenMP parallelism, even combined with MPI parallelism. 



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[lesson:basepar|An introduction on ABINIT in Parallel]] should be read before going to the next lessons about parallelism. One simple example of parallelism in ABINIT will be shown.
* [[lesson:paral_gspw|Parallelism for ground-state calculations, with plane waves]] presents the combined k-point (K), plane-wave (G), band (B), spin/spinor parallelism of ABINIT (so, the "KGB" parallelism), for the computation of total energy, density, and ground state properties 
* [[lesson:paral_moldyn|Parallelism for molecular dynamics calculations]]
* [[lesson:paral_images|Parallelism based on "images", e.g. for the determination of transitions paths (NEB, string method) or PIMD]], that can be activated on top of the "KGB" parallelism for force calculations.
* [[lesson:paral_gswvl|Parallelism for ground-state calculations, with wavelets]] presents the parallelism of ABINIT, when wavelets are used as a basis function instead of planewaves, for the computation of total energy, density, and ground state properties
* [[lesson:paral_dfpt|Parallelism of response-function calculations]] - you need to be familiarized with the calculation of linear-response properties within ABINIT, see the tutorial [[lesson:rf1|Response-Function 1 (RF1)]]
* [[lesson:paral_mbt|Parallelism of Many-Body Perturbation calculations (GW)]] allows to speed up the calculation of accurate electronic structures (quasi-particle band structure, including many-body effects).

