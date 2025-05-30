# Thin-Film Smoothed Particle Hydrodynamics Fluid

[Mengdi Wang](https://wang-mengdi.github.io/), [Yitong Deng](https://yitongdeng.github.io/), Xiangxin Kong, Aditya H. Prasad, [Shiying Xiong](https://shiyingxiong.github.io/), [Bo Zhu](https://faculty.cc.gatech.edu/~bozhu/)

[![webpage](https://img.shields.io/badge/Project-Homepage-green)](https://wang-mengdi.github.io/proj/21-thin-film-sph/)
[![paper](https://img.shields.io/badge/Paper-Public-red)](https://cs.dartmouth.edu/~bozhu/papers/sph_bubble.pdf)
[![code](https://img.shields.io/badge/Source_Code-Github-blue)](https://github.com/wang-mengdi/sph-thin-film)

This repo stores the source code of our SIGGRAPH 2021 paper **Thin-Film Smoothed Particle Hydrodynamics Fluid**.

<figure>
  <img src="./represent3000.jpg" align="left" width="100%" style="margin: 0% 5% 2.5% 0%">
  <figcaption>Several thin-film phenomena as simulated using our proposed method. Counting from left to right: (1) surface flow on a oscillating soap bubble, (2) pinched-off droplets between two circular rims pulling away from one another, (3) a soap bubble bursting into tiny droplets and filaments after being poked from the right, (4) the vibrant, opal-like color pattern caused by Rayleigh-Taylor instability on a large-deforming, wet thin-film surface.</figcaption>
</figure>

## Abstract

We propose a particle-based method to simulate thin-film fluid that jointly facilitates aggressive surface deformation and vigorous tangential flows. We build our dynamics model from the surface tension driven Navier-Stokes equation with the dimensionality reduced using the asymptotic lubrication theory and customize a set of differential operators based on the weakly compressible Smoothed Particle Hydrodynamics (SPH) for evolving pointset surfaces. The key insight is that the compressible nature of SPH, which is unfavorable in its typical usage, is helpful in our application to co-evolve the thickness, calculate the surface tension, and enforce the fluid incompressibility on a thin film. In this way, we are able to two-way couple the surface deformation with the in-plane flows in a physically based manner. We can simulate complex vortical swirls, fingering effects due to Rayleigh-Taylor instability, capillary waves, Newton’s interference fringes, and the Marangoni effect on liberally deforming surfaces by presenting both realistic visual results and numerical validations. The particle-based nature of our system also enables it to conveniently handle topology changes and codimension transitions, allowing us to marry the thin-film simulation with a wide gamut of 3D phenomena, such as pinch-off of unstable catenoids, dripping under gravity, merging of droplets, as well as bubble rupture.

## Usage

1. Install the newest version of [`xmake`](https://xmake.io/#/) build tool.
2. Build the project (you may be prompted to install a few packages):

        $ cd sph-thin-film
        $ xmake -v sph_bubble
3. Run the half bubble example code

        $ xmake r sph_bubble -test 18 -o output

The point positions etc. will be saved as binary files in `./output`. Feel free to [contact me](https://wang-mengdi.github.io/) if you need assistance with visualization.

For those interested in simulating bubbles, we encourage you to follow our subsequent work:

[https://github.com/wang-mengdi/melp-thin-film](https://github.com/wang-mengdi/melp-thin-film)