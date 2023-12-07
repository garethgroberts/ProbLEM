# ProbLEM-Probabilistic Landscape Evolution Modelling

This directory contains a forward model for evolving longitudinal river profiles that incorporates stochastic erosion.  

It was used to generate Figure 11 (below with caption) in the paper Roberts and Wani (Date TBC), A theory of stochastic fluvial landscape evolution, preprint link: https://eartharxiv.org/repository/view/5578/. 

The ProbLEM algorithm is super simple, most of the effort is in the plotting of results! The code evolves a longitudinal river profile (elevation as a function of distance, $`z(x)`$) via the simple scheme:

```math
z^{x}_{t+1} = 
	\begin{cases}
		z^{x}_{t}	& \text{if} \quad F \le c^*\\
		z^{x}_{t} -\delta_z	& \text{if} \quad F > c^*
	\end{cases}
```
where $`t`$ and $`x`$ indicate time and space indicies, respectively. $\delta_z$ is erosion (units of length), $F$ is a (stochastic) driving force, parameterised in the code as an Ornstein-Uhlenbeck process, and $c^*$ is a critical threshold (force) that must be exceeded for erosion to occur. Here it is assumed to be constant, but it could easily be modified to be stochastic.  The following schematic used to explain the approach is from Roberts and Wani (TBC). 

<p align="center">
<img src="https://github.com/garethgroberts/ProbLEM/assets/11752321/bc55b160-3a1e-4cf4-8cfd-0f8adc5b1a69" width="650" height="300" />
</p>

<br/>
<br/>
<br/>

The following figure shows output from the algorithm, it is Figure 11 in Roberts & Wani (TBC). It is a demonstration of local erosional complexity and emergent simplicity for an evolving theoretical river profile subject to stochastic forcing. 

![alt_text](https://github.com/garethgroberts/ProbLEM/files/13598650/rivevo_ou_staged.pdf)

(a-b) Examples of inserted forces (grey), associated expected values and variance (solid and dashed), critical threshold (red), and resultant analytical probabilities of erosion (thick black; see Figure 10a-b for extended explanation). (c) Example of profile evolution in one simulation. Grey and black curves = profiles every 50 and 200 time steps, respectively. (d) Zoom into panel c. Black outlined bars = river profile at time step 600. (e) Zoom into a separate simulation that has different random driving forces (but same distribution). Black outlined bars = time step 600. Note differences between evolution and resultant profiles in panels d and e. (f) ‘Fuzzy’ black curves = profiles at 200, 400, 600 and 800 time steps (e.g. black curves in panel c) from 10 simulations with different random driving forces (but same distribution). Straight line = starting condition. Circles = expected displacements of erosional front originating at the mouth of the ‘river’ at (1000, 0); colors = time steps (see panel c for scale); calculated variance is smaller than symbol size. Small grey square centred at (190, 250) = position of zoomed-in region shown in panels d, e and g. (g) Calculated positions of longitudinal profiles at time step 600 for the 10 simulations shown in panel f. Colors and line widths simply indicate profiles from different simulations for clarity. Note their variability.

<br/>
Note that analytical calculations of expected displacement values and variance are plotted in panel (f), see Roberts and Wani (TBC) for details. 
<br/>
<br/>
<br/>
---

To run a single simulation of the stochastic model use:

> python ProbLEM_force_OU_evo_single_fig.py

A shell script (bash) is provided for running multiple simulations, which calls the python script referred to above, to run:

> ./loop.sh

it deposits output in directories called run_*. Examples used to create Figure 11 are included in the directories named run_*.

There are matplotlib (python) plotting scripts included in ProbLEM_force_OU_evo_single_fig.py, but I find gmt (https://docs.generic-mapping-tools.org/6.4/index.html) offers a bit more flexibility and hence Figure 11 was produced by running:

> ./plot_results.gmt

Also included in this directory are various gmt colour palette files used for plotting, note that these are generated within the plot_results.gmt script. 

The individual scripts have comments to explain variables and values used. Do feel free to get in touch if anything isn't clear. 

Scripts were genereated and run on mac osx (13.5.2), python = python3.11.6, gmt = gmt6.4.0.  

Gareth Roberts, Dec 2023. 
(https://www.garethgroberts.com/home)https://www.garethgroberts.com/home

