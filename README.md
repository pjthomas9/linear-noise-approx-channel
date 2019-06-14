# linear-noise-approx-channel
Code for linear noise approximation analysis of signal transduction information channel project with Greg Hessler and Andrew Eckford

# Figure Generation Files

Fig 1
Fig 2a
Fig 2b
Fig 3
etc.

# Matlab Utility Files

## White noise input.  

Calculate SE^FULL (spectral efficiency, fully observed case), SE^PART_0
(SE for partially observed, low freq limit), SE^PART_\infty (SE for
partially observed, high freq limit).

MI_3state_chain.m		3-state chain
MI_3state_ring.m		3-state ring (e.g. Channelrhodopsin)
MI_ACh.m			Acetylcholine

These files take p as an input (probability of high intensity input for
IID input signal). In these files, baseline transition rate
(\alpha_{ij}^0), sensitivity (\alpha_{ij}^1), observable (C) are
defined internally.  Edit files to change these parameters.  
They are called by 

SE_sweep_p.m



## Colored noise input.  

MI_3state_chain_colored.m  
MI_3state_colored_sweep_p.m  
MI_ACh_colored.m
MI_ACh_colored_sweep_p.m
SE_3state_chain_plot.m
norms_3state_chain.m
norms_3state_ring.m
norms_ACh.m



# Mathematica Files

To be added...

# To be deleted:

MI_3state_sweep_p.m (superceded by SE_sweep.m)

