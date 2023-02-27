---
permalink: /gallery/
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Gallery

This is a display of systems Faunus can be used to simulate.

## Minimal example

This is a minimalistic example of a simple Monte Carlo simulation of charged Lennard-Jones particles in
a high-dielectric continuum solvent.

<a href="https://github.com/mlund/faunus/blob/master/examples/minimal/minimal.yml" class="btn btn--primary">Input</a>

## Liquid Water

![image-right](https://www.dropbox.com/s/2m3499h4n7gp6q6/water-small.png?raw=1){: .align-right}

Here we use the extended simple-point-charge (SPC/E) model for water to simulate the liquid phase in the
_NPT_-ensemble. This also illustrates how to generate a trajectory and perform on-the-fly
analysis. Long range electrostatic interactions are handled using an effective pair-potential.
Ewald summation is also available.

<a href="https://github.com/mlund/faunus/blob/master/examples/water.yml" class="btn btn--primary">Input</a>

## Coarse Grained Lipid Bilayer

![image-right](https://www.dropbox.com/s/7ymnpf9t4w1nt48/membrane-3bead.jpg?raw=1){: .align-right}

Here, a bilayer is assembled from a three-bead lipid model with zero tension using _isochoric_ moves
of the simulation contained dimensions under constant volume. This example also illustrated how to
set up intra-molecular bonds.

<a href="https://github.com/mlund/faunus/blob/master/examples/membrane/membrane.yml" class="btn btn--primary">Input</a>

