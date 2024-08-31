# Source_Wavefield_Simulation

## Overview
This project generates synthetic ground motions for a 1D soil profile using a fast and simplified approach. The method idealizes the soil as an N+1 layered medium with an underlying half-space. Each layer is assumed to be elastic, homogeneous, and isotropic, with constant material properties across the layer, and free of body forces. The source is modeled as a strike-slip fault located in the jth layer.


The P, SV, and SH amplitudes are defined for each layer, and boundary conditions (radiation at infinite depth, displacement continuity, traction continuity, free surface) are used to relate those amplitudes and derive the recurrence relation and the Green's function. Detailed results generated using this method are available in this paper, [Generation of synthetic ground motions for geothermally induced earthquakes](https://drive.google.com/file/d/1WC_JmtXk0oV-6j6XtYoxS9hF3_m4_bzS/view).

## Approach
![image](https://github.com/Aditi-Kumawat/Source_Wavefield_Simulation/assets/72736535/4b977241-f0f9-44fe-8aeb-b3d0c3883af6)

## Visualization
![](src/velocity_animation.gif)


## Code Running Sequence (Located in `/src`)

1. **`sc_save_fSRT_ll.m`**:
   - Generates and saves soil medium and source-depth dependent parameters $f_{SS}$, $f_{SR}$, $f_{RS}$, $f_{RR}$, and $f_{TT}$ (Equation 41) and their derivatives (in Appendix).


2. **`sc_save_I_paras.m`**:
   - Saves $I_r^{(1)}(r,\omega)$, $I_r^{(2)}(r,\omega)$, $I_r^{(3)}(r,\omega)$, $I_\phi^{(1)}(r,\omega)$, $I_\phi^{(2)}(r,\omega)$, $I_z^{(1)}(r,\omega)$, $I_z^{(2)}(r,\omega)$, $I_z^{(3)}(r,\omega)$ (Equations 56-63).



3. **`sc_velCYLcoordi.m`**:
   - Generates ground motion displacement using Equations (46, 48-51) at given epicentral distances and source properties.
 
 
*Equation numbers correspond to the study by Singla and Gupta (2019).

## References
The apporach used for the analysis is described in more detail in the following papers:
- Singla and Gupta (2019): [Surface rotations due to kinematic shear dislocation point source in a multilayered elastic medium](https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/109/1/433/568060/Surface-Rotations-Due-to-Kinematic-Shear?redirectedFrom=fulltext)
- Pei et al. (2008): [Improvements on computation of phase velocities of Rayleigh waves based on the generalized R/T coefficient method](https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/98/1/280/341880/Improvements-on-Computation-of-Phase-Velocities-of?redirectedFrom=fulltext)

By following the steps outlined in the code running sequence, users can generate and analyze synthetic ground motions for various soil profiles (**`fns_inptMatPara.m`**) and source characteristics.
